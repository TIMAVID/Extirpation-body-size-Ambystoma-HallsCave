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

##### READ IN HERP DATA #####
Halls_herp_database <- read_csv("Data/Halls_herp_database.csv")
levels(as.factor(Halls_herp_database$Higher_Classification))

Lizard_SAl_database <- Halls_herp_database %>% filter(grepl('lizard|salamander', Higher_Classification)) # FILTER FOR ONLY LIZARDS AND SALAMANDERS
Lizard_SAl_database$Level_min <- as.numeric(Lizard_SAl_database$Level_min)
Lizard_SAl_database$Level_max <- as.numeric(Lizard_SAl_database$Level_max)


write.csv(Lizard_SAl_database)


write.csv(Lizard_SAl_database,"Data/HallsCave_Lizard_Salamander.csv", row.names = FALSE)

















