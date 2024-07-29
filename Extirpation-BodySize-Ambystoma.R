##### READ IN MODERN SALAMANDER DATA ####
library(curl)
library(readr)
library(dplyr)
f <- curl("https://raw.githubusercontent.com/TIMAVID/Extirpation-body-size-Ambystoma-HallsCave/main/Data/Modern_Ambystoma_measurements.csv")
Modern_Ambystoma <- read_csv(f) # READ MODERN MEASUREMENT DATA IN

summary(Modern_Ambystoma)

Tiger_sals <- filter(Modern_Ambystoma,grepl('tigrinum|mavortium|ordinaria|subsalsum',Species)) # MAKING SEPARATE DATAFRAME FOR SPECIES IN THE TIGER SALAMANDER SPECIES COMPLEX
Tiger_sals_noLarval <- filter(Tiger_sals,!grepl("Larva",`Life Stage`)) # MAKING SEPARATE DATAFRAME FOR FOR SPECIES IN THE TIGER SALAMANDER SPECIES COMPLEX WITHOUT SPECIMENS NOTED AS BEING LARVA

###### Bones function TO CREATE LINEAR MODELS FOR PREDICTING SVL------------
varlist <- names(Modern_Ambystoma)[11:13] #list all the different linear measurements taken

bones.lm<-function(variables,data){ # THIS FUNCTION CREATES LINEAR MODELS FOR EACH MEASUREMENT AND ALSO PLOTS THE MEASUREMENTS WITH SVL
  require(ggplot2)
  require(gridExtra)
  figs<-lapply(variables, function(x) {
    ggplot(data = data, 
           aes(log(SVL_P),log(get(x)),)) + geom_point()+ ggtitle(x)+ theme_classic()+ ylab("log(Measurement)") +
      geom_smooth(aes(x = log(SVL_P), y = log(get(x))),alpha = .2, method = "lm", inherit.aes = FALSE )
  })
  do.call(grid.arrange, c(figs, ncol=3, top = deparse(substitute(data))))
  
  models <- lapply(variables, function(x) { #function to perform linear regression on all measurements for each dataset
    lm(substitute(log(SVL_P) ~ log(i), list(i = as.name(x))), data = data)
  })
  names(models) <- variables
  (sum <- (lapply(models, summary)))
  b<-(lapply(sum, function (x)x$coefficients[1]))
  m<-(lapply(sum, function (x)x$coefficients[2]))
  R<-(lapply(sum, function (x)x$adj.r.squared))
  P<-(lapply(sum, function (x)x$coefficients[2,4])) #may need to p.adjust
  MSE <- (lapply(sum, function (x)(mean(x$residuals^2))))
  Con <- (lapply(models, confint))
  C<-(lapply(Con, function (x)x[2,]))
  out<-list(m,b,R,P,MSE,models,C)
  names(out)<-c("slope","Y-intercept","adj.R-squared","P-value","MSE","models","ConfidenceInt")
  out <- do.call(rbind, out)
  out <- t(out)
  out <- as.data.frame(out)
  return(out)
}

Ambystoma_lm <-bones.lm(varlist, Modern_Ambystoma) #liner models of all measurements based on all species

Tiger_sals_lm <-bones.lm(varlist, Tiger_sals) #liner models of all measurements based on tiger salamander species
Tiger_sals_noLarval_lm <-bones.lm(varlist, Tiger_sals_noLarval) #liner models of all measurements based on tiger salamander species w/o larva

##### READ IN AND CREATE Fossil dataframes--------------------
library(tidyr)
library(dplyr)

g <- curl("https://raw.githubusercontent.com/TIMAVID/Extirpation-body-size-Ambystoma-HallsCave/main/Data/Ambystoma_fossil_measurements_total.csv")
fossil_Ambystoma_long <- read_csv(g)

Fossils <- unique(fossil_Ambystoma_long[c("Specimen_Number","Specimen_Letter")])

fossil_Ambystoma <- spread(fossil_Ambystoma_long, Measurement_Type, `Measurement_(mm)`)
fossil_Ambystoma <- fossil_Ambystoma %>% dplyr::relocate(Hum_L, Hum_DW, Fem_L)



## 5 cm intervals # THIS CODE IS ASSIGNING AN AGE TO EACH FOSSIL BASED ON THE DEPTH OF THE FOSSIL IN THE CAVE DEPOSIT
library(data.table)

classes5 <- seq(from = 1, to = 70, by = 1)
upper5 <- seq(from = 5, to = 350, by = 5)
lower5 <- seq(from = 0, to = 345, by = 5)

fossil_Ambystoma5cmBIN = copy(fossil_Ambystoma)
setDT(fossil_Ambystoma5cmBIN)

interval_lookup5 <- data.table(classes5, upper5,lower5)
fossil_Ambystoma5cmBIN<- fossil_Ambystoma5cmBIN[interval_lookup5, Bin:=classes5, on=c("Level_min >= lower5","Level_max <= upper5")]

###### CREATING AGE-DEPTH MODEL ######
library(rbacon)

Bacon('HallsCave', youngest.age = 31, depths=seq(0, 350,
                                                 length=350)) # set the youngest age as 31 ybp (aka 1993 CE)

times5cm_all<-Bacon.hist(seq(from = 2.5, to = 350, by = 5)) # V1 = min, V2 = max, V3 = median, V4 = mean
times5cmborder<-Bacon.hist(seq(from = 0, to = 350, by = 5)) # V1 = min, V2 = max, V3 = median, V4 = mean
times5cmborder[,2] <- c(times5cmborder[,2][-1], NA)



times5cm <- times5cm_all[,4] # AVERAGE AGE FOR EACH 5 CM BIN BASED ON AGE-DEPTH MODEL FROM BACON
setDT(fossil_Ambystoma5cmBIN)
meanage5cm<-(times5cm)
interval_lookup5age <- data.table(meanage5cm, upper5,lower5)

fossil_Ambystoma5cmBIN<- fossil_Ambystoma5cmBIN[interval_lookup5age, Age:=meanage5cm, on=c("Level_min >= lower5","Level_max <= upper5")]

fossil_Ambystoma5cmBIN_long <- fossil_Ambystoma5cmBIN %>%                                   # HERE WE ARE MAKING THE DATA INTO LONG FORMAT WHICH MAKES IT EASIER TO PLOT EVERTHING TOGETHER
  gather(Measurement_type, Measurement, - c(Specimen_Number:Age))


##### Fossil estimates ######

######lm(Phrynosomatids) -predict> fossils------------
estimate_SVL <- function(lm, data, x, y #need to check specific columns [,X:Y] in data match the order in the lm and are the ones to estimate from
) #function to estimate SVL FOR FOSSILS
{
  t<- t(lm)
  t <- as.data.frame(t)
  slopes <-as.vector(t[1,])
  intercepts <- as.vector(t[2,])
  estimates <- mapply(function(m,b,data){{(exp(as.numeric(m)*
                                                 log(data) + as.numeric(b)))}}, 
                      data= data[,x:y], m=slopes, b=intercepts)
  data<-list(estimates, data$Age, data$Specimen_Number, data$Specimen_Letter, data$Family, data$Genus, data$Details)
  names(data)<-c("Estimated_SVL","Age", "Specimen", "Letter", "Family", "Genus", "Details")
  return((data))
}

fossil_Ambystoma_estimates <- estimate_SVL(Ambystoma_lm, fossil_Ambystoma5cmBIN, x=1, y = 3) # ESTAIMATING SVL FOR FOSSIL SPECIMENS BASED ON THE LINEAR MODEL PREVIOUSLY CREATED BASED ON THE MODERN SPECIMENS
fossil_Ambystoma_estimates_tiger <- estimate_SVL(Tiger_sals_noLarval_lm, fossil_Ambystoma5cmBIN, x=1, y = 3) # ESTAIMATING SVL FOR FOSSIL SPECIMENS BASED ON THE LINEAR MODEL PREVIOUSLY CREATED BASED ON TIGER SALAMANDER SPECIMENS

data_formatting <-function(data) # HERE WE ARE CONVERTING THE FOSSILING ELEMENTS FROM A LIST INTO A DATAFRAME AND CLEANING UP THE DATA A LITTLE BIT
{
  require(Hmisc)
  Fossil_estimates<-llist(data)
  Fossil_estimates<- unlist(Fossil_estimates,recursive=FALSE)
  Fossil_estimates<- lapply(Fossil_estimates, data.frame, stringsAsFactors = FALSE)
  Fossil_estimates<-bind_cols(Fossil_estimates, .name_repair = c("universal_quiet"))
  Fossil_estimates<- Fossil_estimates %>% 
    rename("Age" = "X..i.....4",
           "Specimen_Number" = "X..i.....5",
           "Specimen_Letter" = "X..i.....6",
           "Family" = "X..i.....7",
           "Genus" = "X..i.....8",
           "Details" = "X..i.....9")
  # HERE WE ARE MAKING THE DATA INTO LONG FORMAT WHICH MAKES IT EASIER TO PLOT EVERTHING TOGETHER
  Fossil_estimates <- Fossil_estimates %>%  
    gather(Measurement_type, SVL_estimate, - c(Age:Details))
}


fossil_Ambystoma_estimates <- data_formatting(fossil_Ambystoma_estimates) # CONVERT TO DATAFRAME
fossil_Ambystoma_estimates_tiger <- data_formatting(fossil_Ambystoma_estimates_tiger) # CONVERT TO DATAFRAME


####### PLOTTING BODY SIZE PREDICTIONS FOR FOSSILS THROUGH TIME ####### 
library(ggplot2)
library("ggpubr")
library(mdthemes)
library(ggpmisc)


fossil_Ambystoma_estimates %>% # SUMMARY OF BODY SIZE ESTIMATES FOR FOSSILS
  group_by(Measurement_type) %>% drop_na(SVL_estimate) %>% summarise(Mean = mean(SVL_estimate), sd = sd(SVL_estimate), iqr = IQR(SVL_estimate))

fossil_Ambystoma_estimates_tiger %>% # SUMMARY OF BODY SIZE ESTIMATES FOR FOSSILS BASED ON TIGER SALAMANDER MODELS
  group_by(Measurement_type) %>% drop_na(SVL_estimate) %>% summarise(Mean = mean(SVL_estimate), sd = sd(SVL_estimate), iqr = IQR(SVL_estimate))


fossil_Ambystoma_estimates_p <- ggplot(fossil_Ambystoma_estimates, aes(x = Age, y = SVL_estimate)) + # PLOT OF BODY SIZE CHANGES THROUGH TIME FOR DIFFERENT MEASUREMENTS TAKEN ON FOSSILS
  geom_point(aes(color="#FFC857"), size = 4, position=position_jitter(width=30,height=0)) +
  geom_smooth(aes(x = Age, y = SVL_estimate),alpha = .2, method = "lm", inherit.aes = FALSE ) +
  theme_bw() +scale_color_manual(values ="#FFC857", na.value = "#000000") + ylim(70,160) +
  labs(x = "Years BP", y = "SVL estimates", title = "*Ambystoma*") + mdthemes::md_theme_classic(base_size = 16) +
  stat_cor(aes(x = Age, y = SVL_estimate), formula = y ~ x, inherit.aes = FALSE ,label.y = 137)+ theme(panel.spacing.x = unit(10, "mm"), legend.position = "none")+
  stat_regline_equation(aes(x = Age, y = SVL_estimate), formula = y ~ x, inherit.aes = FALSE ,label.y = 135) + facet_wrap("Measurement_type")

Hum_L <- fossil_Ambystoma_estimates %>% 
  dplyr::filter(Measurement_type == "Hum_L") %>% 
  lm(SVL_estimate ~ Age, 
     data = ., 
     na.action = na.omit)
summary(Hum_L) # SUMMARY STATISTICS FOR LINEAR MODEL (BODY_SIZE ~ AGE) BASED ON HUM_L


fossil_Ambystoma_estimates_tiger_p <- ggplot(fossil_Ambystoma_estimates_tiger, aes(x = Age, y = SVL_estimate)) + # PLOT OF BODY SIZE CHANGES THROUGH TIME FOR DIFFERENT MEASUREMENTS TAKEN ON FOSSILS USING TIGER SALAMANDER MODEL
  geom_point(aes(color="#FFC857"), size = 4, position=position_jitter(width=30,height=0)) +
  geom_smooth(aes(x = Age, y = SVL_estimate),alpha = .2, method = "lm", inherit.aes = FALSE ) +
  theme_bw() +scale_color_manual(values ="#FFC857", na.value = "#000000") + ylim(70,160) +
  labs(x = "Years BP", y = "SVL estimates", title = "*Ambystoma*") + mdthemes::md_theme_classic(base_size = 16) +
  stat_cor(aes(x = Age, y = SVL_estimate), formula = y ~ x, inherit.aes = FALSE ,label.y = 137)+ theme(panel.spacing.x = unit(10, "mm"), legend.position = "none")+
  stat_regline_equation(aes(x = Age, y = SVL_estimate), formula = y ~ x, inherit.aes = FALSE ,label.y = 135) + facet_wrap("Measurement_type")


library(cowplot)
plot_grid(fossil_Ambystoma_estimates_p, fossil_Ambystoma_estimates_tiger_p, ncol = 2, labels = "auto", align = "hv", rel_heights = c(1, 1),
          axis = "lr")


##### READ IN PASTCLIM DATA #####
library(pastclim)
# get time series for location
coords<-data.frame(longitude=-99.533333,latitude=c(30.133333)) # coordinates of Hall's Cave
download_dataset(dataset="Beyer2020", bio_variables = c("bio01","bio10", "bio12", "bio16"))
time_series <-location_series(x=coords,#time_bp=animals$time,
                              bio_variables=c("bio01","bio10","bio12", "bio16"),
                              dataset="Beyer2020")

# subset times on study interval
time_series$time_bp <- -(time_series$time_bp) # MAKE TIMES POSITIVE
time_series_sub <- time_series[time_series$time %in% seq(from = 0, to = 20000, by = 1000),] # TAKE DATAPOINTS FROM 0 TO 20,000 YBP BE 1,000 YEARS

degreeC <- expression(degree~C)
bio01<-ggplot()+
  geom_line(data = time_series_sub, aes(time_bp, bio01)) +theme_classic() + labs(
    title = "Annual Mean Temperature", x = "YPB",
    y = degreeC)
bio10 <- ggplot()+
  geom_line(data = time_series_sub, aes(time_bp, bio10)) +theme_classic() + labs(
    title = "Mean temperature of hottest quarter", x = "YPB",
    y = degreeC)
bio12 <- ggplot()+
  geom_line(data = time_series_sub, aes(time_bp, bio12)) +theme_classic() + labs(
    title = "Annual precipitation", x = "YPB",
    y = "mm per year") +xlim(0, 21000)
bio16 <- ggplot()+
  geom_line(data = time_series_sub, aes(time_bp, bio16)) +theme_classic() + labs(
    title = "Precipitation of wettest quarter", x = "YPB",
    y = "mm per quarter")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### MODEL PALEOCLIM DATA  ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library("mgcv")
library("scam")
library("ggplot2")
library("cowplot")
library("tidyr")
#devtools::install_github("gavinsimpson/gratia")
library("gratia")

bio01_mod <- gam(bio01 ~ s(time_bp, k = 15), gamma = .5, data = time_series_sub,
                 correlation = corCAR1(form = ~ time_bp), 
                 method = "REML")
summary(bio01_mod)
par(mfrow = c(2, 2))
gam.check(bio01_mod)


bio12_mod <- gam(log(bio12) ~ s(time_bp, k = 15), gamma = .7, data = time_series_sub,
                 correlation = corCAR1(form = ~ time_bp), 
                 method = "REML")
summary(bio12_mod)
par(mfrow = c(2, 2))
gam.check(bio12_mod)

meanage5cm_df <- data.frame(time_bp= meanage5cm)


bio16_mod <- gam(log(bio16) ~ s(time_bp, k = 15), gamma = .7, data = time_series_sub,
                 correlation = corCAR1(form = ~ time_bp), 
                 method = "REML")
summary(bio16_mod)
par(mfrow = c(2, 2))
gam.check(bio16_mod)

#.............................bio01.............................
fitbio01 <- predict(bio01_mod, meanage5cm_df, se.fit = TRUE)
bio01crit.t <- qt(0.975, df = df.residual(bio01_mod))
bio01_newGCV <- data.frame(Age = meanage5cm,
                           bio01_fit = fitbio01$fit,
                           se.fit = fitbio01$se.fit)
bio01_newGCV <- transform(bio01_newGCV,
                          upper = bio01_fit + (bio01crit.t * se.fit),
                          lower = bio01_fit - (bio01crit.t * se.fit))
ggplot(bio01_newGCV, aes(Age, bio01_fit))+geom_point()+geom_line() +
  geom_point(data = time_series_sub, aes(x = time_bp, y = bio01), color = "blue")

#.............................bio12.............................
fitbio12 <- predict(bio12_mod, meanage5cm_df, se.fit = TRUE)
bio12crit.t <- qt(0.975, df = df.residual(bio12_mod))
bio12_newGCV <- data.frame(Age = meanage5cm,
                           bio12_fit = fitbio12$fit,
                           se.fit = fitbio12$se.fit)
bio12_newGCV <- transform(bio12_newGCV,
                          upper = bio12_fit + (bio12crit.t * se.fit),
                          lower = bio12_fit - (bio12crit.t * se.fit))
ggplot(bio12_newGCV, aes(Age, bio12_fit))+geom_point()+geom_line() +
  geom_point(data = time_series_sub, aes(x = time_bp, y = log(bio12)), color = "blue")

#.............................bio16.............................
fitbio16 <- predict(bio16_mod, meanage5cm_df, se.fit = TRUE)
bio16crit.t <- qt(0.975, df = df.residual(bio16_mod))
bio16_newGCV <- data.frame(Age = meanage5cm,
                           bio16_fit = fitbio16$fit,
                           se.fit = fitbio16$se.fit)
bio16_newGCV <- transform(bio16_newGCV,
                          upper = bio16_fit + (bio16crit.t * se.fit),
                          lower = bio16_fit - (bio16crit.t * se.fit))
ggplot(bio16_newGCV, aes(Age, bio16_fit))+geom_point()+geom_line() +
  geom_point(data = time_series_sub, aes(x = time_bp, y = log(bio16)), color = "blue")

###### COMBINE BODY SIZE AND PALEOCLIM DATA #####

GAM_data <- merge(fossil_Ambystoma_estimates_tiger, bio01_newGCV, by.x = "Age", by.y = "Age", suffixes = c(".sal",".bio01"), no.dups = TRUE)
GAM_data <- merge(GAM_data, bio12_newGCV, by.x = "Age", by.y = "Age", suffixes = c(".bio01",".bio12"), no.dups = TRUE)

GAM_data$Measurement_type <- as.factor(GAM_data$Measurement_type)

GAM_data %>% # count number of measurements
  dplyr::filter(Measurement_type == "Hum_L") %>% drop_na(SVL_estimate) %>% summarise(n())
GAM_data %>% # count number of measurements
  dplyr::filter(Measurement_type == "Hum_DW") %>% drop_na(SVL_estimate) %>% summarise(n())
GAM_data %>% # count number of measurements
  dplyr::filter(Measurement_type == "Fem_L") %>% drop_na(SVL_estimate) %>% summarise(n())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### LINEAR MODELS TESTING RELATIONSHIP BETWEEN BODY SIZE AND TEMP ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hum_L_lm_bio01 <- GAM_data %>% # HUMERUS LENGTH AND MAT
  dplyr::filter(Measurement_type == "Hum_L") %>% 
  lm(SVL_estimate ~ bio01_fit, 
     data = ., 
     na.action = na.omit)
summary(Hum_L_lm_bio01)


Hum_DW_lm_bio01 <- GAM_data %>% # HUMERUS DISTAL WIDTH AND MAT
  dplyr::filter(Measurement_type == "Hum_DW") %>% 
  lm(SVL_estimate ~ bio01_fit, 
     data = ., 
     na.action = na.omit)
base::summary(Hum_DW_lm_bio01)
confint(Hum_DW_lm_bio01)

Fem_L_lm_bio01 <- GAM_data %>% # FEMUR LENGTH AND MAT
  dplyr::filter(Measurement_type == "Fem_L") %>% 
  lm(SVL_estimate ~ bio01_fit, 
     data = ., 
     na.action = na.omit)
base::summary(Fem_L_lm_bio01)


###### PLOTTING LINEAR RELATIONSHIP BETWEEN RECONSTRUCTED BODY SIZE AND TEMPERATURE ######
GAM_data %>% group_by(Measurement_type) %>% # COUNT NUMBER OF MEASURMENTS FOR EACH SIDE
  filter(Details=="1R") %>% drop_na(SVL_estimate)%>% summarise(n())
GAM_data %>% group_by(Measurement_type) %>% # COUNT NUMBER OF MEASURMENTS FOR EACH SIDE
  filter(Details=="1L") %>% drop_na(SVL_estimate) %>% summarise(n())


Body_size_tem_p<-ggplot(GAM_data, aes(x = bio01_fit, y = SVL_estimate)) +
  geom_point(aes(color="#FFC857"), size = 4, position=position_jitter(width=0,height=0)) +
  geom_smooth(aes(x = bio01_fit, y = SVL_estimate),alpha = .3, color = "black", method = "lm", inherit.aes = FALSE ) +
  theme_bw() +scale_color_manual(values ="#FFC857", na.value = "#000000") +
  labs(x = '°C', y = "SVL estimates", title = "*Ambystoma*") + scale_x_continuous(limits = c(12, 19), breaks = c(12.5, 15, 17.5)) +
  mdthemes::md_theme_classic(base_size = 16) + theme(panel.spacing.x = unit(10, "mm"), legend.position = "none")+ facet_wrap("Measurement_type")
# stat_cor(aes(x = bio01_fit, y = SVL_estimate), formula = y ~ x, inherit.aes = FALSE ,label.y = 82, method = "pearson")+ theme(panel.spacing.x = unit(10, "mm"), legend.position = "none")+
# stat_regline_equation(aes(x = bio01_fit, y = SVL_estimate), formula = y ~ x, inherit.aes = FALSE ,label.y = 80) + facet_wrap("Measurement_type")


GAM_data %>% 
  filter(Details=="1L") %>% # PLOT ONLY LEFT SIDE ELEMENTS
  ggplot(., aes(x=bio01_fit, y=SVL_estimate)) +
  geom_point(aes(color="#FFC857"), size = 4, position=position_jitter(width=0,height=0)) +
  geom_smooth(aes(x = bio01_fit, y = SVL_estimate),alpha = .2, method = "lm", inherit.aes = FALSE ) +
  theme_bw() +scale_color_manual(values ="#FFC857", na.value = "#000000") +
  labs(x = degreeC, y = "SVL estimates", title = "*Ambystoma* L elements") + mdthemes::md_theme_classic(base_size = 16) +
  stat_cor(aes(x = bio01_fit, y = SVL_estimate), formula = y ~ x, inherit.aes = FALSE ,label.y = 137)+ theme(panel.spacing.x = unit(10, "mm"), legend.position = "none")+
  stat_regline_equation(aes(x = bio01_fit, y = SVL_estimate), formula = y ~ x, inherit.aes = FALSE ,label.y = 135) + facet_wrap("Measurement_type")


Hum_DW_L_lm_bio01 <- GAM_data %>% # MODEL ONLY INCLUDES LEFT SIDE ELEMENTS
  dplyr::filter(Measurement_type == "Hum_DW") %>% filter(Details=="1L") %>% 
  lm(SVL_estimate ~ bio01_fit, 
     data = ., 
     na.action = na.omit)
base::summary(Hum_DW_L_lm_bio01)
confint(Hum_DW_L_lm_bio01)

GAM_data %>% 
  filter(Details=="1R") %>% 
  ggplot(., aes(x=bio01_fit, y=SVL_estimate)) +
  geom_point(aes(color="#FFC857"), size = 4, position=position_jitter(width=0,height=0)) +
  geom_smooth(aes(x = bio01_fit, y = SVL_estimate),alpha = .2, method = "lm", inherit.aes = FALSE ) +
  theme_bw() +scale_color_manual(values ="#FFC857", na.value = "#000000") +
  labs(x = degreeC, y = "SVL estimates", title = "*Ambystoma* R elements") + mdthemes::md_theme_classic(base_size = 16) +
  stat_cor(aes(x = bio01_fit, y = SVL_estimate), formula = y ~ x, inherit.aes = FALSE ,label.y = 137)+ theme(panel.spacing.x = unit(10, "mm"), legend.position = "none")+
  stat_regline_equation(aes(x = bio01_fit, y = SVL_estimate), formula = y ~ x, inherit.aes = FALSE ,label.y = 135) + facet_wrap("Measurement_type")

Hum_DW_R_lm_bio01 <- GAM_data %>% # MODEL ONLY INCLUDES LEFT SIDE ELEMENTS
  dplyr::filter(Measurement_type == "Hum_DW") %>% filter(Details=="1R") %>% 
  lm(SVL_estimate ~ bio01_fit, 
     data = ., 
     na.action = na.omit)
base::summary(Hum_DW_R_lm_bio01)
confint(Hum_DW_R_lm_bio01)


###### SUBSAMPLING HUMERUS DISTAL WIDTH DATASET
Hum_DW_sub <- GAM_data %>% 
  dplyr::filter(Measurement_type == "Hum_DW") %>% drop_na(SVL_estimate)
Hum_DW_noNA <- GAM_data %>% 
  dplyr::filter(Measurement_type == "Hum_DW") %>% drop_na(SVL_estimate)
set.seed(123)
Hum_DW_sub <- Hum_DW_sub[sample(nrow(Hum_DW_sub), 57), ] # here we are subsampling the humerus_dw data to 57 rows, the same size as the femur length dataset
Hum_DW_sub_lm <- lm(SVL_estimate ~ bio01_fit, 
                    data = Hum_DW_sub, 
                    na.action = na.omit)
base::summary(Hum_DW_sub_lm) # examine whether we can still see significant relationship with subsample
confint(Hum_DW_sub_lm)

library(car)
Hum_DW_sub.boot <- Boot(Hum_DW_sub_lm, R=10000) # bootstrat replicates with the subsample
print(confint(Hum_DW_sub.boot, level=.95, type="norm")) # examine whether we can still see significant relationship with subsample

f <- function () {
  fit <- lm(SVL_estimate ~ bio01_fit, data = Hum_DW_noNA, subset = sample(nrow(Hum_DW_noNA), 57))
  coef(fit)
}
set.seed(0); z <- t(replicate(1000, f()))
par(mfrow = c(1,2))
hist(z[,1], main = "intercept")
hist(z[,2], main = "slope")


Fem_L_noNA <- GAM_data %>% 
  dplyr::filter(Measurement_type == "Fem_L") %>% drop_na(SVL_estimate)
Fem_L_noNA_lm <- lm(SVL_estimate ~ bio01_fit, 
                    data = Fem_L_noNA, 
                    na.action = na.omit)
base::summary(Fem_L_noNA_lm)
confint(Fem_L_noNA_lm)
Fem_L.boot <- Boot(Fem_L_noNA_lm, R=10000)
print(confint(Fem_L.boot, level=.95, type="norm"))


##### BAYESIAN CHANGE POINT #######
library(bcp)

Hum_DW_noNA_ID <- Hum_DW_noNA %>%
  group_by(Age) %>%
  mutate(Age_id = cur_group_id())  

levels(as.factor(Hum_DW_noNA_ID$Age_id))
levels(as.factor(Hum_DW_noNA_ID$Age))

Hum_DW_noNA_bcp.ri <- bcp(as.vector(Hum_DW_noNA_ID$SVL_estimate), id = Hum_DW_noNA_ID$Age_id)
plot(Hum_DW_noNA_bcp.ri)
data_bcp <- data.frame(postmean = Hum_DW_noNA_bcp.ri$posterior.mean, prob = Hum_DW_noNA_bcp.ri$posterior.prob, Age=as.numeric(levels(as.factor(Hum_DW_noNA_ID$Age))))


bchange_plot<-ggplot(data=data_bcp, aes(x=Age, y=X1)) +
  geom_line(alpha = 0.8, colour = "black",linewidth =3)+ 
  geom_point(data = Hum_DW_noNA_ID,
             mapping = aes(x = Age, y = SVL_estimate),
             inherit.aes = FALSE,
             size = 2)+theme_bw() +
  geom_line(aes(x = Age, y = prob/(1/200), alpha = 0.2, colour = "red", linewidth =3),
            inherit.aes = TRUE)+
  scale_y_continuous(
    
    # Features of the first axis
    name = "Body size posterior mean (mm)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*(1/200), name="Posterior probability")
  ) + coord_cartesian(xlim = c(6000,21000)) + theme(legend.position="none") +xlab("Years BP") + scale_x_continuous(breaks =seq(0,100000,2000))

bchange_plot

bio10_P <- ggplot()+
  geom_line(data = time_series_sub, aes(time_bp, bio10)) +theme_classic() + labs(
    title = "Mean temperature of hottest quarter", x = "YPB",
    y = degreeC)+ coord_cartesian(xlim = c(6000,21000)) + theme_bw()+ scale_x_continuous(breaks =seq(0,100000,2000))

bio01_p <- ggplot()+
  geom_line(data = time_series_sub, aes(time_bp, bio01)) +theme_classic() + labs(
    title = "Mean annual temperature", x = "YPB",
    y = degreeC)+ coord_cartesian(xlim = c(6000,21000)) + theme_bw() + scale_x_continuous(breaks =seq(0,100000,2000))

bio12_p <-ggplot()+
  geom_line(data = time_series_sub, aes(time_bp, bio12)) +theme_classic() + labs(
    title = "Anual precipitation", x = "YPB",
    y = degreeC)+ coord_cartesian(xlim = c(6000,21000)) + theme_bw() + scale_x_continuous(breaks =seq(0,100000,2000))

library(cowplot)
plot_grid(bchange_plot, bio01_p, ncol = 1, labels = "auto", align = "hv", rel_heights = c(1.7, 1, 1),
          axis = "lr")


######### GAM BODY SIZE TEMP MODELS #########

M_sal_AR4 <- gamm(formula = SVL_estimate ~
                    s(Measurement_type, bs = "re") +
                    s(Age, bs="tp", by = Measurement_type, k = 15) + 
                    s(bio01_fit, bs="tp", by = Measurement_type, k = 15)
                  ,
                  data = GAM_data, method="ML", 
                  correlation = corARMA(form = ~1 | Age, p = 4),
                  select = TRUE,
                  family= "gaussian")
M_sal_AR3 <- gamm(formula = SVL_estimate ~ 
                    s(Measurement_type, bs = "re") +
                    s(Age, bs="tp", by = Measurement_type, k = 15) + 
                    s(bio01_fit, bs="tp", by = Measurement_type, k = 15)
                  ,
                  data = GAM_data, method="ML", 
                  correlation = corARMA(form = ~1 | Age, p = 3),
                  select = TRUE,
                  family= "gaussian")
M_sal_AR2 <- gamm(formula = SVL_estimate ~ 
                    s(Measurement_type, bs = "re") +
                    s(Age, bs="tp", by = Measurement_type, k = 15) + 
                    s(bio01_fit, bs="tp", by = Measurement_type, k = 15)
                  ,
                  data = GAM_data, method="ML", 
                  correlation = corARMA(form = ~1 | Age, p = 2),
                  select = TRUE,
                  family= "gaussian")
M_sal_AR1 <- gamm(formula = SVL_estimate ~ 
                    s(Measurement_type, bs = "re") +
                    s(Age, bs="tp", by = Measurement_type, k = 15) + 
                    s(bio01_fit, bs="tp", by = Measurement_type, k = 15)
                  ,
                  data = GAM_data, method="ML", 
                  correlation = corARMA(form = ~1 | Age, p = 1),
                  select = TRUE,
                  family= "gaussian")
M_sal_AR0 <- gamm(formula = SVL_estimate ~ 
                    s(Measurement_type, bs = "re") +
                    s(Age, bs="tp", by = Measurement_type, k = 15) + 
                    s(bio01_fit, bs="tp", by = Measurement_type, k = 15)
                  ,
                  data = GAM_data, method="ML",
                  select = TRUE,
                  family= "gaussian")
dev.new()
k.check(M_sal_AR2$gam) # MODEL EVALUATION
par(mfrow = c(2, 2))
gam.check(M_sal_AR2$gam)
qq_plot(M_sal_AR2$gam, method = "simulate") +
  labs(title =NULL, subtitle =NULL)
par(mfrow = c(2, 2))
acf(resid(M_sal_AR4$lme), lag.max = 36, main = "ACF")
pacf(resid(M_sal_AR4$lme), lag.max = 36, main = "pACF")
acf(resid(M_sal_AR2$lme), lag.max = 36, main = "ACF")
pacf(resid(M_sal_AR2$lme), lag.max = 36, main = "pACF")
anova(M_sal_AR0$lme, M_sal_AR1$lme, M_sal_AR2$lme, M_sal_AR4$lme) # MODEL COMPARISON



layout(matrix(1:2, ncol = 2))
res <- resid(M_sal_AR2$lme, type = "normalized") # EXTRACT RESIDUALS
acf(res, lag.max = 36, main = "ACF - AR(2) errors") # MODEL RESIDUAL EVALUATION
pacf(res, lag.max = 36, main = "pACF- AR(2) errors") # MODEL RESIDUAL EVALUATION
layout(1)



summary(M_sal_AR2$gam) #  SUMMARY OF BEST FITTING MODEL
gratia::draw(M_sal_AR2, residuals = TRUE, parametric = FALSE) & theme_bw() # PLOT PARTIAL EFFECTS OF EACH SMOOTH TERM
parital_effects_p <- gratia::draw(M_sal_AR2, select = c("s(bio01_fit)"), ncol = 1, residuals = TRUE, partial_match = TRUE) &
  labs(x = '°C') & theme_bw(base_size = 12) # PLOT PARTIAL EFFECTS OF TEMPRERATURE SMOOTH TERM

dev.new()
plot_grid(Body_size_tem_p, parital_effects_p, ncol = 2, labels = "auto", align = "hv", rel_heights = c(1, 1.4),
          axis = "lr") # PLOT LINEAR MODEL AND GAM MODEL RESULTS

library(tidygam)
M_sal_AR2$gam %>%
  predict_gam(length_out = 20, series = "bio01_fit") %>%
  plot(comparison = "Measurement_type") +
  scale_color_brewer(type = "qual") + scale_fill_brewer(type = "qual") + theme_bw() & theme_bw(base_size = 12)

##### READ IN FOSSIL HERP DATABASE FOR HALL'S CAVE #####
h <- curl("https://raw.githubusercontent.com/TIMAVID/Extirpation-body-size-Ambystoma-HallsCave/main/Data/HallsCave_Lizard_Salamander.csv")
Lizard_SAl_database <- read_csv(h)

###### ASSIGN AGES TO FOSSILS######
head(meanage5cm)

minage5cm<-c(times5cmborder[,1])
minage5cm<-minage5cm[-71]
maxage5cm<-c(times5cmborder[,2])
maxage5cm<-maxage5cm[-71]
interval_lookup5age <- data.table(meanage5cm, upper5,lower5)
interval_lookup5Minage <- data.table(minage5cm, upper5,lower5)
interval_lookup5Maxage <- data.table(maxage5cm, upper5,lower5)

Lizard_SAl5cmBIN = copy(Lizard_SAl_database)
setDT(Lizard_SAl5cmBIN)

Lizard_SAl5cmBIN<- Lizard_SAl5cmBIN[interval_lookup5age, Age:=meanage5cm, on=c("Level_min >= lower5","Level_max <= upper5")]
Lizard_SAl5cmBIN<- Lizard_SAl5cmBIN[interval_lookup5Minage, MinAge:=minage5cm, on=c("Level_min >= lower5","Level_max <= upper5")]
Lizard_SAl5cmBIN<- Lizard_SAl5cmBIN[interval_lookup5Maxage, MaxAge:=maxage5cm, on=c("Level_min >= lower5","Level_max <= upper5")]

######EXAMINE ELEMENT REPRESENTATION ###### 
Lizard_SAl5cmBIN$Element <- ifelse(Lizard_SAl5cmBIN$Element == "pubis+ichium", # combine common elements
                                   "ichium_pubis", Lizard_SAl5cmBIN$Element)
Lizard_SAl5cmBIN$Element <- ifelse(Lizard_SAl5cmBIN$Element == "pubis, ichium", # combine common elements
                                   "ichium_pubis", Lizard_SAl5cmBIN$Element)
Lizard_SAl5cmBIN$Element <- ifelse(Lizard_SAl5cmBIN$Element == "pubis_ichiumm", # combine common elements
                                   "ichium_pubis", Lizard_SAl5cmBIN$Element)
Lizard_SAl5cmBIN$Element <- ifelse(Lizard_SAl5cmBIN$Element == "pubis_ichium", # combine common elements
                                   "ichium_pubis", Lizard_SAl5cmBIN$Element)
Lizard_SAl5cmBIN$Element <- ifelse(Lizard_SAl5cmBIN$Element == "phalanx", # combine common elements
                                   "phalanges", Lizard_SAl5cmBIN$Element)
Lizard_SAl5cmBIN$Element <- ifelse(Lizard_SAl5cmBIN$Element == "metacarples?", # combine common elements
                                   "phalanges", Lizard_SAl5cmBIN$Element)

SalTotal <- Lizard_SAl5cmBIN %>% # TOTAL NUMBER OF FOSSILS
  filter(grepl('salamander', Higher_Classification)) %>% summarise(NISP = sum(NISP))

SAL_element_totals <- Lizard_SAl5cmBIN %>% # NISP OF EACH ELEMENT
  filter(grepl('salamander', Higher_Classification)) %>% 
  group_by(Element) %>% summarise(NISP = sum(NISP))

SAL_element_byage <- Lizard_SAl5cmBIN %>% drop_na(Age) %>% # NISP by AGE
  filter(grepl('salamander', Higher_Classification))%>% 
  group_by(Age) %>% summarise(NISP = sum(NISP))


Element_rep <- Lizard_SAl5cmBIN %>% # NUMBER OF ELEMENTS IN EACH TIME BIN
  filter(grepl('salamander', Higher_Classification)) %>% 
  aggregate(.,   # Applying aggregate to count # of unique elements per bin
            Element ~ Age,
            function(x) length(unique(x)))

Element_rep <- data.frame(yr=Element_rep$Age, #make dataframe
                          Num=(as.vector(as.matrix(Element_rep$Element))),
                          NISP=(as.vector(as.matrix(SAL_element_byage$NISP))))

Element_rep_plot<-ggplot(Element_rep )+ # NUMBER OF ELEMENTS AND NISP IN EACH TIME BIN
  geom_line(aes(yr, Num), colour="black")+
  geom_bar(aes(x=yr, y=NISP),stat="identity", fill="cyan",colour="#006000")+
  scale_x_reverse(breaks =seq(0,100000,2000))+
  scale_y_continuous(breaks =seq(0,500,20))+
  xlab("Age (cal. BP)")+ylab("NISP")+
  coord_flip() +
  theme_classic(base_size = 17)
Element_rep_plot

ggscatter(Element_rep, x = "NISP", y = "Num", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "log(NISP)", ylab = "log(Num)")

######EXAMINE SALAMANDER PERISITENCE THROUGH TIME######
SAl <- Lizard_SAl5cmBIN %>% 
  #filter(!is.na(Details))  %>%
  filter(grepl('salamander', Higher_Classification)) %>% 
  filter(!is.na(Higher_Classification)) %>% 
  group_by(Higher_Classification,Age,MinAge,MaxAge, .drop=FALSE) %>% summarise(NISP = n(), Prescence = any(NISP>0)) %>%  filter(!is.na(Age))

SAl_plot <- ggplot(SAl, aes(Age,sqrt(NISP)))+
  # geom_pointrange(aes(xmin=MinAge, xmax=MaxAge, size=(NISP)))+ scale_size(range = c(.5,2))+
  geom_area() +
  xlim(0, 21000) +
  xlab("Years BP")+ylab("%")+
  #coord_flip()+
  theme_classic2() + theme(legend.position="bottom") +ylab(expression(paste("NISP"^(1/2))))
# theme(axis.title.y=element_blank())
#       # axis.text.y=element_blank(),
#       # axis.ticks.y=element_blank()) +

#####READ IN HALL'S CAVE PALEOPROXIE DATASETS#####
i <- curl("https://raw.githubusercontent.com/TIMAVID/Extirpation-body-size-Ambystoma-HallsCave/main/Data/Mammal_Toomey_NISP.csv")
Mammal_Toomey_NISP <- read_csv(i)

Mammal_Toomey_NISP = copy(Mammal_Toomey_NISP)
setDT(Mammal_Toomey_NISP)

Mammal_Toomey_NISP<- Mammal_Toomey_NISP[interval_lookup5age, Age:=meanage5cm, on=c("Interval == lower5")]

library(colorBlindness)

Mammal_Toomey_NISP_gopher_plot <- Mammal_Toomey_NISP %>% 
  filter(grepl('Geomys_NISP|Thomomys_NISP', Taxon)) %>% 
  ggplot(aes(fill=Taxon, y = NISP, x=Age)) + 
  geom_area(position='fill')+
  # theme_classic(base_size = 15) +
  theme_classic2() +
  xlim(0, 21000) +
  scale_fill_manual(name = "Taxon", values=c(SKY_BLUE, ORANGE)) +
  ylab("Relative Abundace") + xlab("Years BP")+theme(legend.position="bottom")

Mammal_Toomey_NISP_shrew_plot <- Mammal_Toomey_NISP %>% 
  filter(grepl('Notiosorex_NISP|Cryptotis_NISP', Taxon)) %>% 
  ggplot(aes(fill=Taxon, y = NISP, x=Age)) + 
  geom_area(position='fill')+
  # theme_classic(base_size = 15) +
  theme_classic2() +
  xlim(0, 21000) +
  scale_fill_manual(name = "Taxon", values=c(BLUE, VERMILLION)) +
  ylab("Relative Abundace") + xlab("Years BP")+theme(legend.position="bottom")


j <- curl("https://raw.githubusercontent.com/TIMAVID/Extirpation-body-size-Ambystoma-HallsCave/main/Data/Cooke_strontium.csv")
Cooke_strontium <- read_csv(j)
Cooke_strontium$Material <- ifelse(Cooke_strontium$Material == "vole_enamel", # combine common elements
                                   "rodent_enamel", Cooke_strontium$Material)
Cooke_strontium$Material <- ifelse(Cooke_strontium$Material == "pocket_gopher_enamel", # combine common elements
                                   "rodent_enamel", Cooke_strontium$Material)

k <- curl("https://raw.githubusercontent.com/TIMAVID/Extirpation-body-size-Ambystoma-HallsCave/main/Data/Sun_Paleoclim_data_Ti.csv")
Paleoclim_data_Ti <- read_csv(k)
times_Ti <-Bacon.hist(Paleoclim_data_Ti$AvgDepth) # V1 = min, V2 = max, V3 = median, V4 = mean
times_Ti <- data.frame(times_Ti, Depth = Paleoclim_data_Ti$AvgDepth)
Paleoclim_data_Ti <- merge(Paleoclim_data_Ti, times_Ti, by.x = "AvgDepth", by.y = "Depth")
Paleoclim_data_Ti <- Paleoclim_data_Ti %>% 
  rename(
    Age_old = Age,
    Age = X4
  )

j <- curl("https://raw.githubusercontent.com/TIMAVID/Extirpation-body-size-Ambystoma-HallsCave/main/Data/Sun_Paleoclim_data_C_D.csv")
Paleoclim_C_D <- read_csv(j)
times_C_D<-Bacon.hist(seq(from = .5, to = 285.5, by = 1)) # V1 = min, V2 = max, V3 = median, V4 = mean
times_C_D<- data.frame(times_C_D, Depth = seq(from = .5, to = 285.5, by = 1))
Paleoclim_C_D <- merge(Paleoclim_C_D, times_C_D, by.x = "AvgDepth", by.y = "Depth")
Paleoclim_C_D <- Paleoclim_C_D %>% 
  rename(
    Age_old = Age,
    Age = X4
  )


d13c_label <- expression(delta^{13}*C)
dD_label <- expression(delta*D[wax])
Sr_label <- expression(Sr^{87}/Sr^{86})


Ti_plt <- ggplot(Paleoclim_data_Ti, aes(x = Age, y = Ti_pct)) +
  geom_point() + geom_line(color = "red")+ xlim(0, 21000) +
  labs(y = "Ti %", x = "Years Before Present") +theme_classic2()

d13c_plt <- ggplot(Paleoclim_C_D, aes(x = Age, y = d13C)) +
  geom_point() + geom_line(color = "blue")+ xlim(0, 21000)+scale_y_reverse()+
  labs(y = d13c_label, x = "Years Before Present") +theme_classic2()

dD_plt <- ggplot(Paleoclim_C_D, aes(x = Age, y = dD_corr)) +
  geom_point() + geom_line(color = "forestgreen")+ xlim(0, 21000)+ scale_y_reverse()+
  labs(y = dD_label, x = "Years Before Present") +theme_classic2()

dSr_plt <- ggplot(Cooke_strontium, aes(x = YBP, y = Sr87_Sr86, shape = Material, color = Material)) +
  geom_point(size = 2) + scale_color_manual(values=c("#490092", "#004949"))  + stat_smooth(method="lm", se=TRUE,
                                                                                           formula=y ~ poly(x, 2, raw=TRUE), colour="black", alpha = .2)+ 
  labs(y = Sr_label, x = "Years Before Present") +theme_classic2() +theme(legend.position="bottom")+
  geom_hline(yintercept=0.708408, linetype="dashed", color = ORANGE)+annotate("text", x=0, y=0.70835, label="Thin soils")+
  geom_hline(yintercept=0.709968, linetype="dashed", color = SKY_BLUE)+annotate("text", x=0, y=0.710026, label="Thick soils") + xlim(0, 21000)

bio12_p2 <- ggplot()+
  geom_line(data = time_series_sub, aes(time_bp, bio12)) +theme_classic2() + labs(
    title = "Annual precipitation", x = "YPB",
    y = "mm per year") +xlim(0, 21000)

######PLOT ALL PALEOPROXIES WITH SALAMANDER PERSISTENCE######
library(cowplot)
plot_grid(SAl_plot, dSr_plt, Ti_plt,Mammal_Toomey_NISP_gopher_plot, ncol = 1, labels = "auto", align = "hv", rel_heights = c(1, 1.7, 1,1.3),
          axis = "lr")
plot_grid(SAl_plot, bio12_p2, d13c_plt, dD_plt, ncol = 1, labels = "auto", align = "hv", rel_heights = c(.7, 1,1,1),
          axis = "lr")

plot_grid(SAl_plot, dSr_plt,Mammal_Toomey_NISP_gopher_plot,bio12_p2, ncol = 1, labels = "auto", align = "hv", rel_heights = c(1.3, 1.3, 1,1.3),
          axis = "lr")


##### AMBYSTOMA EXIRPATION VARIABLE TESTING #####
library("mgcv")
library("scam")
library("ggplot2")
library("cowplot")
library("tidyr")
# packageurl <- "https://github.com/gavinsimpson/gratia/archive/refs/tags/v0.8.1.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
#devtools::install_github("gavinsimpson/gratia")
library("gratia")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ MODEL PALEOCLIM DATA  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## fit Paleoclim data GAM using gamm() with a CAR(1)

#.............................Carbon.............................----

d13C_mod <- gamm(d13C ~ s(Age, k = 45), data = Paleoclim_C_D,
                 correlation = corCAR1(form = ~ Age), method = "REML")
summary(d13C_mod$gam)
par(mfrow = c(2, 2))
gam.check(d13C_mod$gam)

## MODEL PALEOPROXIE FOR LEVELS WITH SALAMANDERS FROM HALL'S CAVE
meanage5cm_df2 <- NULL
meanage5cm_df2$Age <- meanage5cm_df$time_bp

fitd13C <- predict(d13C_mod$gam, meanage5cm_df2, se.fit = TRUE)
d13Ccrit.t <- qt(0.975, df = df.residual(d13C_mod$gam))
d13C_newGCV <- data.frame(Age = meanage5cm_df2,
                          d13C_fit = fitd13C$fit,
                          se.fit = fitd13C$se.fit)
d13C_newGCV <- transform(d13C_newGCV,
                         upper = d13C_fit + (d13Ccrit.t * se.fit),
                         lower = d13C_fit - (d13Ccrit.t * se.fit))

ggplot(d13C_newGCV, aes(Age, d13C_fit))+geom_point()+geom_line() +
  geom_point(data = Paleoclim_C_D, aes(x = Age, y = d13C), color = "blue")

#.............................dD wax.............................----

dD_mod <- gamm(dD_corr ~ s(Age, k = 40), data = Paleoclim_C_D,
               correlation = corCAR1(form = ~ Age), method = "REML")
summary(dD_mod$gam)
par(mfrow = c(2, 2))
gam.check(dD_mod$gam)

## MODEL PALEOPROXIE FOR LEVELS WITH SALAMANDERS FROM HALL'S CAVE

fitdD <- predict(dD_mod$gam, meanage5cm_df2, se.fit = TRUE)
dDcrit.t <- qt(0.975, df = df.residual(dD_mod$gam))
dD_newGCV <- data.frame(Age = meanage5cm_df2,
                        dD_fit = fitdD$fit,
                        se.fit = fitdD$se.fit)
dD_newGCV <- transform(dD_newGCV,
                       upper = dD_fit + (dDcrit.t * se.fit),
                       lower = dD_fit - (dDcrit.t * se.fit))

ggplot(dD_newGCV, aes(Age, dD_fit))+geom_point()+geom_line() +
  geom_point(data = Paleoclim_C_D, aes(x = Age, y = dD_corr), color = "blue")

#.............................dSr.............................----
Cooke_strontium$Material <- as.factor(Cooke_strontium$Material)

Cooke_strontium_plant <- Cooke_strontium %>% 
  filter(Material == "hackberry_aragonite")

colnames(Cooke_strontium_plant)[4] <- "Age"

dSr_mod <- gam(Sr87_Sr86 ~ s(Age, k = 20), data = Cooke_strontium_plant,
               method = "REML")

summary(dSr_mod)
par(mfrow = c(2, 2))
gam.check(dSr_mod)
library(tidygam)

dSr_mod %>%
  predict_gam(length_out = 20, series = "Age") %>%
  plot(serier = "Age")

## MODEL PALEOPROXIE FOR LEVELS WITH SALAMANDERS FROM HALL'S CAVE

fitdSr <- predict(dSr_mod, meanage5cm_df2, se.fit = TRUE)
dSrcrit.t <- qt(0.975, df = df.residual(dSr_mod))
dSr_newGCV <- data.frame(Age = meanage5cm_df2,
                         dSr_fit = fitdSr$fit,
                         se.fit = fitdSr$se.fit)
dSr_newGCV <- transform(dSr_newGCV,
                        upper = dSr_fit + (dSrcrit.t * se.fit),
                        lower = dSr_fit - (dSrcrit.t * se.fit))

ggplot(dSr_newGCV, aes(Age, dSr_fit))+geom_point()+geom_line() +
  geom_point(data = Cooke_strontium_plant, aes(x = Age, y = Sr87_Sr86), color = "blue")


####### MAKE SALAMANDER PRESCENCE/ABSENCE DATAFRAME ######
SAl_P_A <- Lizard_SAl5cmBIN %>% 
  dplyr::group_by(Age,Higher_Classification, .drop=FALSE) %>% summarise(NISP = n()) %>% 
  ungroup %>%
  complete(Age,Higher_Classification,
           fill = list(NISP = 0)) %>% drop_na() %>% mutate(Prescence = case_when(NISP > 0 ~ 1, NISP == 0 ~ 0)) %>% filter(Higher_Classification == "salamander")

GAM_data_extirpate <- merge(SAl_P_A, d13C_newGCV, by.x = "Age", by.y = "Age", suffixes = c(".sal",".d13C"), no.dups = TRUE) # MERGE SALAMANDER PRESCENCE/ABSENCE WITH PALEOPROXY DATAFRAMES
GAM_data_extirpate <- merge(GAM_data_extirpate, dD_newGCV, by.x = "Age", by.y = "Age", suffixes = c(".d13C",".dD"), no.dups = TRUE)
GAM_data_extirpate <- merge(GAM_data_extirpate, dSr_newGCV, by.x = "Age", by.y = "Age", suffixes = c(".dD",".dSr"), no.dups = TRUE)
GAM_data_extirpate <- merge(GAM_data_extirpate, bio01_newGCV, by.x = "Age", by.y = "Age", suffixes = c(".dSr",".bio01"), no.dups = TRUE)
GAM_data_extirpate <- merge(GAM_data_extirpate, bio12_newGCV, by.x = "Age", by.y = "Age", suffixes = c(".bio01",".bio12"), no.dups = TRUE)
GAM_data_extirpate <- merge(GAM_data_extirpate, bio16_newGCV, by.x = "Age", by.y = "Age", suffixes = c(".bio12",".bio16"), no.dups = TRUE)
# GAM_data_extirpate_long<- tidyr::gather(GAM_data_extirpate, Prescence,  factor_key=TRUE)
GAM_data_extirpate2 <- GAM_data_extirpate

GAM_data_extirpate[,5] <- c(scale(GAM_data_extirpate[,5])) # SCALE PALEOPROXY DATA
GAM_data_extirpate[,9] <- c(scale(GAM_data_extirpate[,9])) # SCALE PALEOPROXY DATA
GAM_data_extirpate[,13] <- c(scale(GAM_data_extirpate[,13])) # SCALE PALEOPROXY DATA
GAM_data_extirpate[,17] <- c(scale(GAM_data_extirpate[,17])) # SCALE PALEOPROXY DATA
GAM_data_extirpate[,21] <- c(scale(GAM_data_extirpate[,21])) # SCALE PALEOPROXY DATA
GAM_data_extirpate[,25] <- c(scale(GAM_data_extirpate[,25])) # SCALE PALEOPROXY DATA


dev.new()
plot(GAM_data_extirpate[,c(3, 4,5,9,13, 17,21,25)]) # show correlation between variables
library("PerformanceAnalytics")
my_data <- GAM_data_extirpate[,c(3, 4,5,9,13, 17,21,25)]
chart.Correlation(my_data, histogram=TRUE, pch=19)

###### PENALIZED LOGISTIC REGRESSION BETWEEN PRESCENCE/ABSENCE AND PALEOPROXIES ######
library(arm)
fit_bayes <- bayesglm(Prescence ~ bio01_fit+d13C_fit*dSr_fit, data=GAM_data_extirpate, family=binomial(link =logit))
display(fit_bayes)
library(logistf)
fit_firth <- logistf(Prescence ~ bio01_fit+d13C_fit*dSr_fit, data=GAM_data_extirpate) # Firth's Bias-Reduced Logistic Regression
summary(fit_firth)
library(brglm)
fit_brglm <- brglm(Prescence ~ bio01_fit+d13C_fit*dSr_fit, data=GAM_data_extirpate, family=binomial(link =logit))
summary(fit_brglm)

exp(fit_brglm$coefficients)
library(ggplot2)
library(sjPlot)
plot_model(fit_brglm, type = "pred", terms = "dSr_fit [all]") + theme_classic()
plot_model(fit_brglm, type = "pred", terms = "d13C_fit [all]")+ theme_classic()
plot_model(fit_brglm, type = "int", terms = "d13C_fit [all]")+ theme_classic()


nd <- expand.grid(d13C_fit = seq(min(GAM_data_extirpate$d13C_fit), max(GAM_data_extirpate$d13C_fit), length = 10),
                  dSr_fit  = seq(min(GAM_data_extirpate$dSr_fit), max(GAM_data_extirpate$dSr_fit), length = 10))

preds <- predict(fit_brglm, newdata = nd, se.fit = TRUE)
preds$fit_prob <- exp(preds$fit)/(1+exp(preds$fit))
nd$pred  <- plogis(preds$fit)

ggplot(nd, aes(d13C_fit, dSr_fit, fill = pred)) +
  geom_tile() +
  scale_fill_viridis_c("Probability of\nAnimal presence") +
  coord_fixed(1) +
  theme_minimal(base_size = 16)


ggplot(GAM_data_extirpate2, aes(x = d13C_fit, y=Prescence)) + geom_point() +
  stat_smooth(method = "glm", method.args = list(family=binomial), se = TRUE) + xlab("d13C") +
  ylab("Probability of presence") +
  ggtitle("Probability of presence of Ambystoma d13C")+theme_classic()

ggplot(GAM_data_extirpate2, aes(x = dSr_fit, y=Prescence)) + geom_point() +
  stat_smooth(method = "glm", method.args = list(family=binomial), se = TRUE) + xlab("dSr") +
  ylab("Probability of presence") +
  ggtitle("Probability of presence of Ambystoma dSr")+theme_classic()

ggplot(GAM_data_extirpate2, aes(x = bio12_fit, y=Prescence)) + geom_point() +
  stat_smooth(method = "glm", method.args = list(family=binomial), se = TRUE) + xlab("bio12") +
  ylab("Probability of presence") +
  ggtitle("Probability of presence of Ambystoma bio12")+theme_classic()

















