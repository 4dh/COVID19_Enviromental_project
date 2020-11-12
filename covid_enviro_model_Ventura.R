library(ggplot2)
library(mgcv)
library("imputeTS")
#Moritz, Steffen, and Bartz-Beielstein, Thomas. "imputeTS: Time Series Missing Value Imputation in R." R Journal 9.1 (2017). doi: 10.32614/RJ-2017-009.
library(lubridate)
library(data.table)
library(tscount)

rm(list=ls())

#_______________________MSPE function_______________________________
#Mean Squared Prediction Error
#Arguments are 2 vectors for observed, y and predicted observations hat{y} from 
#the test component of the  model. The mspe formula is mean((y_i - \hat{y}_i)^2)
mspe_val <- function(y, y_hat){
  mspe <- 0
  total_df <- data.frame(y, y_hat)
  #total_df$y <- y
  #total_df$y_hat <- y_hat
  total_df$pred_error_sq <- (abs(total_df$y_hat-total_df$y))^2
  mspe <- (sum(total_df$pred_error_sq))/length(total_df$pred_error_sq)
  return(mspe)
}


#_______________________VENTURA___________________________________
Ventura_df = read.csv("C:\\Users\\dheym\\OneDrive\\Documents\\GitHub\\csci708\\src\\Ventura_merged.csv")

#sort and format dates
Ventura_df <- Ventura_df[ order(Ventura_df$Date ),]
Ventura_df$Date  <- ymd(Ventura_df$Date)
str(Ventura_df)

#date range with the most data: 2/15/2020 - 5/31/2020
Ventura_s <- subset(Ventura_df, Ventura_df$Date > "2020-02-15" & Ventura_df$Date < "2020-05-31")
#subset only columns for what I'm interested in: 
#date(2), logaqi(6), boxcoxmeanozone (8),  boxcoxmeanpm2.5 (12),parks% (13), DiffI (16), DiffD (17)
Ventura_s <- Ventura_s[, c(2, 6, 8, 12, 13, 16, 17)]

#________________________Filling in NAs________________________

for (a in 2:ncol(Ventura_s)){
  Ventura_s[a] <- na_kalman(Ventura_s[a], model = "auto.arima", smooth = FALSE)
}

Ventura_t <- transform(Ventura_s, ndate = as.numeric(Date),
                          nyear  = as.numeric(format(Date, '%Y')),
                          nmonth = as.numeric(format(Date, '%m')),
                          day    = as.numeric(format(Date, '%j')))


##############
# MODEL: GLM #
##############


#dataframe for glm with 1 or 2 lags:
Ventura_lags <- subset(Ventura_t)
Ventura_lags$I_diff_t_1 <- Ventura_lags$I_diff
Ventura_lags$I_diff_t_1 <- shift(Ventura_lags$I_diff_t_1, n=1, fill=NA, type="lag")
Ventura_lags$D_diff_t_1 <- Ventura_lags$D_diff
Ventura_lags$D_diff_t_1 <- shift(Ventura_lags$D_diff_t_1, n=1, fill=NA, type="lag")
Ventura_lags$I_diff_t_2 <- Ventura_lags$I_diff
Ventura_lags$I_diff_t_2 <- shift(Ventura_lags$I_diff_t_2, n=2, fill=NA, type="lag")
Ventura_lags$D_diff_t_2 <- Ventura_lags$D_diff
Ventura_lags$D_diff_t_2 <- shift(Ventura_lags$D_diff_t_2, n=2, fill=NA, type="lag")
#get rid of first 2 rows with NAs now 
Ventura_lags <- na.omit(Ventura_lags)

#########train and test sets
V_glm_train_dat <- Ventura_lags[c(1:96),]
V_glm_test_dat <- Ventura_lags[c(97:103),]


#________________________MODEL FITTING______________________________
################        SACRAMENTO       ##################
#response: death_Diff; linear glm model
#glm
mod_V_D_glm <- glm((D_diff ~  D_diff_t_1 + D_diff_t_2 + log_aqi + boxcox_mean_Ozone_ppm + boxcox_mean_PM2_5_µg_m3_LC +
                      parks_percent_change_from_baseline), family = poisson(),
                   data=V_glm_train_dat)

#store 
V_glm_train_dat$resid_D <- resid(mod_V_D_glm)
V_glm_train_dat$fitted_D <- fitted(mod_V_D_glm)

#diagnostic summary and storing std error
summary(mod_V_D_glm)
V_glm_train_dat_SE_D <- summary(mod_V_D_glm)$coefficients[, 2]
plot(mod_V_D_glm)
preds <- predict.glm(mod_V_D_glm, type = "response", interval = "confidence")




#response: infect_Diff; linear glm model
#glm
mod_V_I_glm <- glm((I_diff ~  I_diff_t_1 + I_diff_t_2 + log_aqi + boxcox_mean_Ozone_ppm + boxcox_mean_PM2_5_µg_m3_LC +
                      parks_percent_change_from_baseline), family = poisson(),
                   data=V_glm_train_dat)

#store residual and fitted values
V_glm_train_dat$resid_I<- resid(mod_V_I_glm)
V_glm_train_dat$fitted_I <- fitted(mod_V_I_glm)

#diagnostic summary and storing std error
summary(mod_V_I_glm)
V_glm_train_dat_SE_I <- summary(mod_V_I_glm)$coefficients[, 2]
plot(mod_V_I_glm)


#____________Plotting, Predicting, and Evaluation_________________
#plotting model fit vs observed              
ggplot(V_glm_train_dat , aes(x=Date, y = Values)) + 
  geom_line(aes(y=fitted_I, color="fitted_I")) +
  geom_line(aes(y=I_diff, color="actual_I")) +
  geom_line(aes(y=fitted_D, color="fitted_D")) + 
  geom_line(aes(y=D_diff, color="actual_D"))  +ggtitle("GLM Model Fit- Ventura, CA")


#predict I
newdata_I <- subset(V_glm_test_dat, select = c(I_diff_t_1, I_diff_t_2,
                                               log_aqi, boxcox_mean_Ozone_ppm, 
                                               boxcox_mean_PM2_5_µg_m3_LC,
                                               parks_percent_change_from_baseline))
V_glm_test_dat$I_pred <- predict(mod_V_I_glm, newdata_I, type = "response")


#df is results df1 with confidence interval value columns     
#99% CI
critval <- 1.96
#95% CI 
critval_2 <- 2.576
pred_I <- predict(mod_V_I_glm, type="response", se.fit = TRUE)
df1 <-  subset(V_glm_train_dat, select = c(I_diff, Date))
df1 <- cbind(df1, pred_I = pred_I$fit)
df1 <- cbind(df1, se = pred_I$se.fit) 
df1 <- cbind(df1, ucl = df1$pred_I + critval*df1$se)
df1 <- cbind(df1, lcl = df1$pred_I - critval*df1$se)
df1 <- cbind(df1, ucl_2 = df1$pred_I + critval_2*df1$se)
df1 <- cbind(df1, lcl_2 = df1$pred_I - critval_2*df1$se)

pred_I2 <- predict(mod_V_I_glm, newdata_I, type="response", se.fit = TRUE)
df1a <-  subset(V_glm_test_dat, select = c(I_diff, Date))
df1a <- cbind(df1a, pred_I2 = pred_I2$fit)
df1a <- cbind(df1a, se = pred_I2$se.fit) 
df1a <- cbind(df1a, ucl = df1a$pred_I2 + critval*df1a$se)
df1a <- cbind(df1a, lcl = df1a$pred_I2 - critval*df1a$se)
df1a <- cbind(df1a, ucl_2 = df1a$pred_I2 + critval_2*df1a$se)
df1a <- cbind(df1a, lcl_2 = df1a$pred_I2 - critval_2*df1a$se)

#plotting model fit vs observed   : Infection_diff   

p3 <- ggplot(df1) + 
  #99% CI
  geom_ribbon(data = df1, aes(x = Date, ymin = lcl_2, ymax = ucl_2,  fill = "99% CI"), alpha = .6) +
  #95% CI
  geom_ribbon(data = df1, aes(x = Date, ymin = lcl, ymax = ucl, fill = "95% CI"), alpha = .9)+
  #observed
  geom_point(data=df1, aes(x=Date, y=I_diff, color = "obs_I"), size = 1, shape = 23) +
  #fit
  geom_point(data=df1, aes(x=Date, y=pred_I), size = 1, shape = 23) + geom_line(aes(x = Date, y=pred_I, color="fitted_I")) +
  
  #PREDICTED WEEK
  #99% CI
  geom_ribbon(data = df1a, aes(x = Date, ymin = lcl_2, ymax = ucl_2,  fill = "99% CI"), alpha = .6) +
  #95% CI
  geom_ribbon(data = df1a, aes(x = Date, ymin = lcl, ymax = ucl, fill = "95% CI"), alpha = .9)+
  #observed
  geom_point(data=df1a, aes(x=Date, y=I_diff, color = "obs_I"), size = 1, shape = 23) +
  #predicted
  geom_point(data=df1a, aes(x=Date, y=pred_I2), size = 1, shape = 23) + geom_line(data = df1a, aes(x = Date, y=pred_I2, color="predicted_I")) +
  #vertical line to show prediction region
  geom_vline(xintercept=as.Date("2020-05-24"), colour="orangered1") +
  geom_text(aes(x=as.Date("2020-05-24"), label="prediction", y=50), colour="orangered3", angle=90, vjust = 1.2)



p3 + theme(legend.title=element_blank()) +ggtitle("GLM Infection_diff Model, Ventura, CA")+ theme_bw() + scale_fill_manual(values = c("lightskyblue", "darkolivegreen2")) + scale_color_manual(values = c("mediumvioletred", "royalblue", "orangered3")) 


#predict D
newdata_D <- subset(V_glm_test_dat, select = c(D_diff_t_1, D_diff_t_2,
                                               log_aqi, boxcox_mean_Ozone_ppm, 
                                               boxcox_mean_PM2_5_µg_m3_LC,
                                               parks_percent_change_from_baseline))
V_glm_test_dat$D_pred <- predict(mod_V_D_glm, newdata_D,  type = "response")


#df is results df1 with confidence interval value columns     
#99% CI
critval <- 1.96
#95% CI 
critval_2 <- 2.576
pred_D <- predict(mod_V_D_glm, type="response", se.fit = TRUE)
df2 <-  subset(V_glm_train_dat, select = c(D_diff, Date))
df2 <- cbind(df2, pred_D = pred_D$fit)
df2 <- cbind(df2, se = pred_D$se.fit) 
df2 <- cbind(df2, ucl = df2$pred_D + critval*df2$se)
df2 <- cbind(df2, lcl = df2$pred_D - critval*df2$se)
df2 <- cbind(df2, ucl_2 = df2$pred_D + critval_2*df2$se)
df2 <- cbind(df2, lcl_2 = df2$pred_D - critval_2*df2$se)
pred_D2 <- predict(mod_V_D_glm, newdata_D, type="response", se.fit = TRUE)
df2a <-  subset(V_glm_test_dat, select = c(D_diff, Date))
df2a <- cbind(df2a, pred_D2 = pred_D2$fit)
df2a <- cbind(df2a, se = pred_D2$se.fit) 
df2a <- cbind(df2a, ucl = df2a$pred_D2 + critval*df2a$se)
df2a <- cbind(df2a, lcl = df2a$pred_D2 - critval*df2a$se)
df2a <- cbind(df2a, ucl_2 = df2a$pred_D2 + critval_2*df2a$se)
df2a <- cbind(df2a, lcl_2 = df2a$pred_D2 - critval_2*df2a$se)


#plotting model fit vs observed   : Death_diff   
#add titles
p3 <- ggplot(df2) + 
  #99% CI
  geom_ribbon(data = df2, aes(x = Date, ymin = lcl_2, ymax = ucl_2,  fill = "99% CI"), alpha = .6) +
  #95% CI
  geom_ribbon(data = df2, aes(x = Date, ymin = lcl, ymax = ucl, fill = "95% CI"), alpha = .9)+
  #observed
  geom_point(data=df2, aes(x=Date, y=D_diff, color = "obs_D"), size = 1, shape = 23) +
  #fit
  geom_point(data=df2, aes(x=Date, y=pred_D), size = 1, shape = 23) + geom_line(aes(x = Date, y=pred_D, color="fitted_D")) +
  
  #PREDICTED WEEK
  #99% CI
  geom_ribbon(data = df2a, aes(x = Date, ymin = lcl_2, ymax = ucl_2,  fill = "99% CI"), alpha = .6) +
  #95% CI
  geom_ribbon(data = df2a, aes(x = Date, ymin = lcl, ymax = ucl, fill = "95% CI"), alpha = .9)+
  #observed
  geom_point(data=df2a, aes(x=Date, y=D_diff, color = "obs_D"), size = 1, shape = 23) +
  #predicted
  geom_point(data=df2a, aes(x=Date, y=pred_D2), size = 1, shape = 23) + geom_line(data = df2a, aes(x = Date, y=pred_D2, color="predicted_D")) +
  #vertical line to show prediction region
  geom_vline(xintercept=as.Date("2020-05-24"), colour="orangered1") +
  geom_text(aes(x=as.Date("2020-05-24"), label="prediction", y=1.8), colour="orangered3", angle=90, vjust = 1.2)


p3 + theme(legend.title=element_blank()) +ggtitle("GLM Death_diff Model, Ventura, CA") + theme_bw() + scale_fill_manual(values = c("lightskyblue", "darkolivegreen2")) + scale_color_manual(values = c("mediumvioletred", "royalblue", "orangered3")) 




#simple plot for the test week data
ggplot(V_glm_test_dat , aes(x=Date, y = Values)) + 
  geom_line(aes(y=I_diff, color="I_obs")) + 
  geom_line(aes(y=D_diff, color="D_obs")) +
  geom_line(aes(y=I_pred, color="I_pred")) + 
  geom_line(aes(y=D_pred, color="D_pred"))


##########
## MSPE ##
##########

mspe_val(V_glm_test_dat$I_diff, V_glm_test_dat$I_pred)
#156.8412
mspe_val(V_glm_test_dat$D_diff, V_glm_test_dat$D_pred)
#0.5457606

#______________________TSGLM____________________________
V_tsglm_train_dat <- Ventura_t[c(1:98),]
V_tsglm_test_dat <- Ventura_t[c(99:105),]

V_tsglm_reg_mat <-  data.matrix(subset(V_tsglm_train_dat, select = c(log_aqi, boxcox_mean_Ozone_ppm, 
                                                                     boxcox_mean_PM2_5_µg_m3_LC,
                                                                     parks_percent_change_from_baseline)))

#INFECTION MODEL  
mod_V_I_tsglm <- tsglm(V_tsglm_train_dat$I_diff,
                       model = list(past_obs = 1, past_mean = 1), link = "log", distr = "poisson",
                       xreg = V_tsglm_reg_mat)

summary(mod_V_I_tsglm)
#store 
V_tsglm_train_dat$resid_I <- residuals(mod_V_I_tsglm)
V_tsglm_train_dat$fitted_I <- fitted(mod_V_I_tsglm)
plot(mod_V_I_tsglm)
scoring(mod_V_I_tsglm)
# logarithmic   quadratic   spherical    rankprob      dawseb      normsq     sqerror 
# 3.1211526  -0.2831711  -0.4245136   2.9336820   3.7094821   2.9293893  46.2424099 

#____________Plotting, Predicting, and Evaluation_________________
#plotting model fit vs observed              
ggplot(V_tsglm_train_dat , aes(x=Date, y = Values)) + 
  geom_line(aes(y=fitted_I, color="fitted_I")) +
  geom_line(aes(y=I_diff, color="actual_I")) 
# +geom_line(aes(y=fitted_D, color="fitted_D")) + 
# geom_line(aes(y=D_diff, color="actual_D"))


#predict I
newdata_I2 <- subset(V_tsglm_test_dat, select = c(log_aqi, boxcox_mean_Ozone_ppm, 
                                                  boxcox_mean_PM2_5_µg_m3_LC,
                                                  parks_percent_change_from_baseline))
V_tsglm_test_dat$I_pred <- predict(mod_V_I_tsglm, n.ahead = 7, newxreg = newdata_I2)$pred
df1 <-  subset(V_tsglm_train_dat)
df1a <-  subset(V_tsglm_test_dat)

#df is results df1 with confidence interval value columns     
#99% CI
critval <- 1.96
#95% CI 
critval_2 <- 2.576
# pred_I <- predict(mod_SB_I_tsglm, type="response", se.fit = TRUE)
# df1 <-  subset(SB_tsglm_train_dat, select = c(I_diff, Date))
# df1 <- cbind(df1, pred_I = pred_I$fit)
# df1 <- cbind(df1, se = pred_I$se.fit) 
# df1 <- cbind(df1, ucl = df1$pred_I + critval*df1$se)
# df1 <- cbind(df1, lcl = df1$pred_I - critval*df1$se)
# df1 <- cbind(df1, ucl_2 = df1$pred_I + critval_2*df1$se)
# df1 <- cbind(df1, lcl_2 = df1$pred_I - critval_2*df1$se)
# 
# pred_I2 <- predict(mod_SB_I_glm, newdata_I, type="response", se.fit = TRUE)
# df1a <-  subset(SB_tsglm_test_dat, select = c(I_diff, Date))
# df1a <- cbind(df1a, pred_I2 = pred_I2$fit)
# df1a <- cbind(df1a, se = pred_I2$se.fit) 
# df1a <- cbind(df1a, ucl = df1a$pred_I2 + critval*df1a$se)
# df1a <- cbind(df1a, lcl = df1a$pred_I2 - critval*df1a$se)
# df1a <- cbind(df1a, ucl_2 = df1a$pred_I2 + critval_2*df1a$se)
# df1a <- cbind(df1a, lcl_2 = df1a$pred_I2 - critval_2*df1a$se)

#plotting model fit vs observed   : Infection_diff   


w <- ggplot(df1) + 
  #99% CI
  #  geom_ribbon(data = df1, aes(x = Date, ymin = lcl_2, ymax = ucl_2,  fill = "99% CI"), alpha = .6) +
  #95% CI
  #  geom_ribbon(data = df1, aes(x = Date, ymin = lcl, ymax = ucl, fill = "95% CI"), alpha = .9)+
  #observed
  geom_point(data=df1, aes(x=Date, y=I_diff, color = "obs_I"), size = 1, shape = 23) +
  #fit
  geom_point(data=df1, aes(x=Date, y=fitted_I), size = 1, shape = 23) + geom_line(aes(x = Date, y=fitted_I, color="fitted_I")) +
  
  #PREDICTED WEEK
  #99% CI
  #  geom_ribbon(data = df1a, aes(x = Date, ymin = lcl_2, ymax = ucl_2,  fill = "99% CI"), alpha = .6) +
  #95% CI
  #  geom_ribbon(data = df1a, aes(x = Date, ymin = lcl, ymax = ucl, fill = "95% CI"), alpha = .9)+
  #observed
  geom_point(data=df1a, aes(x=Date, y=I_diff, color = "obs_I"), size = 1, shape = 23) +
  #predicted
  geom_point(data=df1a, aes(x=Date, y=I_pred), size = 1, shape = 23) + geom_line(data = df1a, aes(x = Date, y=I_pred, color="predicted_I")) +
  #vertical line to show prediction region
  geom_vline(xintercept=as.Date("2020-05-24"), colour="orangered1") +
  geom_text(aes(x=as.Date("2020-05-24"), label="prediction", y=50), colour="orangered3", angle=90, vjust = 1.2)



w + theme(legend.title=element_blank()) +ggtitle("TSGLM Infection_diff Model, Ventura, CA")+ theme_bw() + scale_fill_manual(values = c("lightskyblue", "darkolivegreen2")) + scale_color_manual(values = c("mediumvioletred", "royalblue", "orangered3")) 






#DEATH MODEL
mod_V_D_tsglm <- tsglm(V_tsglm_train_dat$D_diff,
                       model = list(past_obs = 1, past_mean = 1), link = "log", distr = "poisson",
                       xreg = V_tsglm_reg_mat)
summary(mod_V_D_tsglm)
#store 
V_tsglm_train_dat$resid_D <- residuals(mod_V_D_tsglm)
V_tsglm_train_dat$fitted_D <- fitted(mod_V_D_tsglm)
plot(mod_V_D_tsglm)
scoring(mod_V_D_tsglm)
# ogarithmic   quadratic   spherical    rankprob      dawseb      normsq     sqerror 
# 0.5923363  -0.6868232  -0.8128928   0.2407830  -1.1094379   0.7542459   0.5012795 

#____________Plotting, Predicting, and Evaluation_________________
#plotting model fit vs observed              
ggplot(V_tsglm_train_dat , aes(x=Date, y = Values)) + 
  geom_line(aes(y=fitted_I, color="fitted_I")) +
  geom_line(aes(y=I_diff, color="actual_I")) +
  geom_line(aes(y=fitted_D, color="fitted_D")) + 
  geom_line(aes(y=D_diff, color="actual_D"))


#predict D
newdata_D2 <- subset(V_tsglm_test_dat, select = c(log_aqi, boxcox_mean_Ozone_ppm, 
                                                  boxcox_mean_PM2_5_µg_m3_LC,
                                                  parks_percent_change_from_baseline))
V_tsglm_test_dat$D_pred <- predict(mod_V_D_tsglm, n.ahead = 7, newxreg = newdata_D2)$pred
df1 <-  subset(V_tsglm_train_dat)
df1a <-  subset(V_tsglm_test_dat)

#df is results df1 with confidence interval value columns     
#99% CI
critval <- 1.96
#95% CI 
critval_2 <- 2.576
# pred_I <- predict(mod_SB_I_tsglm, type="response", se.fit = TRUE)
# df1 <-  subset(SB_tsglm_train_dat, select = c(I_diff, Date))
# df1 <- cbind(df1, pred_I = pred_I$fit)
# df1 <- cbind(df1, se = pred_I$se.fit) 
# df1 <- cbind(df1, ucl = df1$pred_I + critval*df1$se)
# df1 <- cbind(df1, lcl = df1$pred_I - critval*df1$se)
# df1 <- cbind(df1, ucl_2 = df1$pred_I + critval_2*df1$se)
# df1 <- cbind(df1, lcl_2 = df1$pred_I - critval_2*df1$se)
# 
# pred_I2 <- predict(mod_SB_I_glm, newdata_I, type="response", se.fit = TRUE)
# df1a <-  subset(SB_tsglm_test_dat, select = c(I_diff, Date))
# df1a <- cbind(df1a, pred_I2 = pred_I2$fit)
# df1a <- cbind(df1a, se = pred_I2$se.fit) 
# df1a <- cbind(df1a, ucl = df1a$pred_I2 + critval*df1a$se)
# df1a <- cbind(df1a, lcl = df1a$pred_I2 - critval*df1a$se)
# df1a <- cbind(df1a, ucl_2 = df1a$pred_I2 + critval_2*df1a$se)
# df1a <- cbind(df1a, lcl_2 = df1a$pred_I2 - critval_2*df1a$se)

#plotting model fit vs observed   : Infection_diff   


w2 <- ggplot(df1) + 
  #99% CI
  #  geom_ribbon(data = df1, aes(x = Date, ymin = lcl_2, ymax = ucl_2,  fill = "99% CI"), alpha = .6) +
  #95% CI
  #  geom_ribbon(data = df1, aes(x = Date, ymin = lcl, ymax = ucl, fill = "95% CI"), alpha = .9)+
  #observed
  geom_point(data=df1, aes(x=Date, y=D_diff, color = "obs_D"), size = 1, shape = 23) +
  #fit
  geom_point(data=df1, aes(x=Date, y=fitted_D), size = 1, shape = 23) + geom_line(aes(x = Date, y=fitted_D, color="fitted_D")) +
  
  #PREDICTED WEEK
  #99% CI
  #  geom_ribbon(data = df1a, aes(x = Date, ymin = lcl_2, ymax = ucl_2,  fill = "99% CI"), alpha = .6) +
  #95% CI
  #  geom_ribbon(data = df1a, aes(x = Date, ymin = lcl, ymax = ucl, fill = "95% CI"), alpha = .9)+
  #observed
  geom_point(data=df1a, aes(x=Date, y=D_diff, color = "obs_D"), size = 1, shape = 23) +
  #predicted
  geom_point(data=df1a, aes(x=Date, y=D_pred), size = 1, shape = 23) + geom_line(data = df1a, aes(x = Date, y=D_pred, color="predicted_D")) +
  #vertical line to show prediction region
  geom_vline(xintercept=as.Date("2020-05-24"), colour="orangered1") +
  geom_text(aes(x=as.Date("2020-05-24"), label="prediction", y=.9), colour="orangered3", angle=90, vjust = 1.2)



w2 + theme(legend.title=element_blank()) +ggtitle("TSGLM Death_diff Model, Ventura, CA")+ theme_bw() + scale_fill_manual(values = c("lightskyblue", "darkolivegreen2")) + scale_color_manual(values = c("mediumvioletred", "royalblue", "orangered3")) 





##########
## MSPE ##
##########

mspe_val(V_tsglm_test_dat$I_diff, V_tsglm_test_dat$I_pred)
#120.8618
mspe_val(V_tsglm_test_dat$D_diff, V_tsglm_test_dat$D_pred)
#0.5922987





#___________TSGLM y = Death_diff with another covariate, I_diff with lag (TBD in thr function)___________
#This model is similar to the previous TSGLM model, but now we focus only on death_diff as the response
#variable and also introduced a new covariate, I_diff with a number of days as a lag (from 1-15 days).
#I will loop through model generation, altering the column associated with lagged I_diff so that for each
#iteration, the lag value will increase by 1 day (starting from -15). The mspe will be saved to a table,
#and then it will be clear which lag time produces the most effective model.

#new data frame to alter
V_tot_dat<- Ventura_t

#function to generate model for given data sets and lag amount for the infection diff; lags will be from 1-15
#days prior to date, t
model_compare <- function(dat_df){
  #results dataframe
  total_table <- data.frame(matrix(ncol = 2, nrow = 0))
  x_names <- c("lag_days", "MSPE")
  colnames(total_table) <- x_names
  
  for (lag_amt in 1:15){
    #create lag
    #create covariate column for lagged infection_diff
    dat_df$I_diff_t_n <- dat_df$I_diff
    dat_df$I_diff_t_n <- shift(dat_df$I_diff_t_n, n=lag_amt, fill=NA, type="lag")
    
    tmp_train <- dat_df[c(1:98),]
    tmp_test <- dat_df[c(99:105),]
    
    #get rid of first rows with NAs due to the lag created
    tmp_train <- na.omit(tmp_train)
    
    #covariates matrix
    tmp_reg_mat <-  data.matrix(subset(tmp_train, select = c(I_diff_t_n,log_aqi, boxcox_mean_Ozone_ppm, 
                                                             boxcox_mean_PM2_5_µg_m3_LC,
                                                             parks_percent_change_from_baseline)))
    #MODEL  
    tmp_mod_tsglm <- tsglm(tmp_train$D_diff,
                           model = list(past_obs = 1, past_mean = 1), link = "log", distr = "poisson",
                           xreg = tmp_reg_mat)
    
    #predictive/test component
    newdat <- subset(tmp_test, select = c(I_diff_t_n, log_aqi, boxcox_mean_Ozone_ppm, 
                                          boxcox_mean_PM2_5_µg_m3_LC,
                                          parks_percent_change_from_baseline))
    tmp_test$D_pred <- predict(tmp_mod_tsglm, n.ahead = 7, newxreg = newdat)$pred
    tmp_mspe <- mspe_val(tmp_test$D_diff, tmp_test$D_pred)
    
    total_table <- rbind(total_table, c(lag_amt, tmp_mspe))   
    
  }
  
  
  return(total_table)
}

model_compare(V_tot_dat)


#It is evident that the lowest mspe (0.4257723) is that of model with lag of 13 days for infection. Therefore, I will
#use this as a covariate in the model.
#create lag
#create covariate column for lagged infection_diff
V_tot_dat$I_diff_t_n <- V_tot_dat$I_diff
V_tot_dat$I_diff_t_n <- shift(V_tot_dat$I_diff_t_n, n=13, fill=NA, type="lag")

V_tot_train <- V_tot_dat[c(1:98),]
V_tot_test <- V_tot_dat[c(99:105),]

#get rid of first rows with NAs due to the lag created
V_tot_train<- na.omit(V_tot_train)

#covariates matrix
tmp_reg_mat <-  data.matrix(subset(V_tot_train, select = c(I_diff_t_n,log_aqi, boxcox_mean_Ozone_ppm, 
                                                           boxcox_mean_PM2_5_µg_m3_LC,
                                                           parks_percent_change_from_baseline)))
#MODEL  
tmp_mod_tsglm <- tsglm(V_tot_train$D_diff,
                       model = list(past_obs = 1, past_mean = 1), link = "log", distr = "poisson",
                       xreg = tmp_reg_mat)
V_tot_train$fitted_D <- fitted(tmp_mod_tsglm)
#predictive/test component
newdat <- subset(V_tot_test, select = c(I_diff_t_n,log_aqi, boxcox_mean_Ozone_ppm, 
                                        boxcox_mean_PM2_5_µg_m3_LC,
                                        parks_percent_change_from_baseline))
V_tot_test$D_pred <- predict(tmp_mod_tsglm, n.ahead = 7, newxreg = newdat)$pred
tmp_mspe <- mspe_val(V_tot_test$D_diff, V_tot_test$D_pred)
tmp_mspe
#0.4257723


df1 <-  subset(V_tot_train)
df1a <-  subset(V_tot_test)
#plotting model fit vs observed   : Death_diff with infection lagged as covariate   

summary(tmp_mod_tsglm)

w3 <- ggplot(df1) + 
  #99% CI
  #  geom_ribbon(data = df1, aes(x = Date, ymin = lcl_2, ymax = ucl_2,  fill = "99% CI"), alpha = .6) +
  #95% CI
  #  geom_ribbon(data = df1, aes(x = Date, ymin = lcl, ymax = ucl, fill = "95% CI"), alpha = .9)+
  #observed
  geom_point(data=df1, aes(x=Date, y=D_diff, color = "obs_D"), size = 1, shape = 23) +
  #fit
  geom_point(data=df1, aes(x=Date, y=fitted_D), size = 1, shape = 23) + geom_line(aes(x = Date, y=fitted_D, color="fitted_D")) +
  
  #PREDICTED WEEK
  #99% CI
  #  geom_ribbon(data = df1a, aes(x = Date, ymin = lcl_2, ymax = ucl_2,  fill = "99% CI"), alpha = .6) +
  #95% CI
  #  geom_ribbon(data = df1a, aes(x = Date, ymin = lcl, ymax = ucl, fill = "95% CI"), alpha = .9)+
  #observed
  geom_point(data=df1a, aes(x=Date, y=D_diff, color = "obs_D"), size = 1, shape = 23) +
  #predicted
  geom_point(data=df1a, aes(x=Date, y=D_pred), size = 1, shape = 23) + geom_line(data = df1a, aes(x = Date, y=D_pred, color="predicted_D")) +
  #vertical line to show prediction region
  geom_vline(xintercept=as.Date("2020-05-24"), colour="orangered1") +
  geom_text(aes(x=as.Date("2020-05-24"), label="prediction", y=.9), colour="orangered3", angle=90, vjust = 1.2)



w3 + theme(legend.title=element_blank())  +ggtitle("TSGLM Death_diff Model, full Ventura, CA")+ theme_bw() + scale_fill_manual(values = c("lightskyblue", "darkolivegreen2")) + scale_color_manual(values = c("mediumvioletred", "royalblue", "orangered3")) 






