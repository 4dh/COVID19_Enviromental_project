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


#_______________________SACRAMENTO_______________________________
Sacramento_df = read.csv("C:\\Users\\dheym\\OneDrive\\Documents\\GitHub\\csci708\\src\\Sacramento_merged.csv")

#sort and format dates
Sacramento_df <- Sacramento_df[ order(Sacramento_df$Date ),]
Sacramento_df$Date  <- ymd(Sacramento_df$Date)
str(Sacramento_df)

#date range with the most data: 2/15/2020 - 5/31/2020
Sacramento_s <- subset(Sacramento_df, Sacramento_df$Date > "2020-02-15" & Sacramento_df$Date < "2020-05-31")
#subset only columns for what I'm interested in: 
#date(2), boxcoxaqi(6), meanozone (7),  boxcoxmeanpm2.5 (11), parks% (12), DiffI (15), DiffD (16)
Sacramento_s <- Sacramento_s[, c(2, 6, 7, 11, 12, 15, 16)]


#________________________Filling in NAs________________________

for (a in 2:ncol(Sacramento_s)){
  Sacramento_s[a] <- na_kalman(Sacramento_s[a], model = "auto.arima", smooth = FALSE)
}

Sacramento_t <- transform(Sacramento_s, ndate = as.numeric(Date),
                          nyear  = as.numeric(format(Date, '%Y')),
                          nmonth = as.numeric(format(Date, '%m')),
                          day    = as.numeric(format(Date, '%j')))


##############
# MODEL: GLM #
##############


#dataframe for glm with 1 or 2 lags:
Sacramento_lags <- subset(Sacramento_t)
Sacramento_lags$I_diff_t_1 <- Sacramento_lags$I_diff
Sacramento_lags$I_diff_t_1 <- shift(Sacramento_lags$I_diff_t_1, n=1, fill=NA, type="lag")
Sacramento_lags$D_diff_t_1 <- Sacramento_lags$D_diff
Sacramento_lags$D_diff_t_1 <- shift(Sacramento_lags$D_diff_t_1, n=1, fill=NA, type="lag")
Sacramento_lags$I_diff_t_2 <- Sacramento_lags$I_diff
Sacramento_lags$I_diff_t_2 <- shift(Sacramento_lags$I_diff_t_2, n=2, fill=NA, type="lag")
Sacramento_lags$D_diff_t_2 <- Sacramento_lags$D_diff
Sacramento_lags$D_diff_t_2 <- shift(Sacramento_lags$D_diff_t_2, n=2, fill=NA, type="lag")
#get rid of first 2 rows with NAs now 
Sacramento_lags <- na.omit(Sacramento_lags)

#########train and test sets
S_glm_train_dat <- Sacramento_lags[c(1:96),]
S_glm_test_dat <- Sacramento_lags[c(97:103),]


#________________________MODEL FITTING______________________________
################        SACRAMENTO       ##################
#response: death_Diff; linear glm model
#glm
mod_S_D_glm <- glm((D_diff ~  D_diff_t_1 + D_diff_t_2 + boxcox_aqi + mean_Ozone_ppm + boxcox_mean_PM2_5_µg_m3_LC +
                        parks_percent_change_from_baseline), family = poisson(),
                    data=S_glm_train_dat)

#store 
S_glm_train_dat$resid_D <- resid(mod_S_D_glm)
S_glm_train_dat$fitted_D <- fitted(mod_S_D_glm)

#diagnostic summary and storing std error
summary(mod_S_D_glm)
S_glm_train_dat_SE_D <- summary(mod_S_D_glm)$coefficients[, 2]
plot(mod_S_D_glm)
preds <- predict.glm(mod_S_D_glm, type = "response", interval = "confidence")




#response: infect_Diff; linear glm model
#glm
mod_S_I_glm <- glm((I_diff ~  I_diff_t_1 + I_diff_t_2 + boxcox_aqi + mean_Ozone_ppm + boxcox_mean_PM2_5_µg_m3_LC +
                      parks_percent_change_from_baseline), family = poisson(),
                    data=S_glm_train_dat)

#store residual and fitted values
S_glm_train_dat$resid_I<- resid(mod_S_I_glm)
S_glm_train_dat$fitted_I <- fitted(mod_S_I_glm)

#diagnostic summary and storing std error
summary(mod_S_I_glm)
S_glm_train_dat_SE_I <- summary(mod_S_I_glm)$coefficients[, 2]
plot(mod_S_I_glm)


#____________Plotting, Predicting, and Evaluation_________________
#plotting model fit vs observed              
ggplot(S_glm_train_dat , aes(x=Date, y = Values)) + 
  geom_line(aes(y=fitted_I, color="fitted_I")) +
  geom_line(aes(y=I_diff, color="actual_I")) +
  geom_line(aes(y=fitted_D, color="fitted_D")) + 
  geom_line(aes(y=D_diff, color="actual_D"))+ggtitle("GLM Model Fit- Sacramento, CA")


#predict I
newdata_I <- subset(S_glm_test_dat, select = c(I_diff_t_1, I_diff_t_2,
                                                boxcox_aqi, mean_Ozone_ppm, 
                                                boxcox_mean_PM2_5_µg_m3_LC,
                                                parks_percent_change_from_baseline))
S_glm_test_dat$I_pred <- predict(mod_S_I_glm, newdata_I, type = "response")


#df is results df1 with confidence interval value columns     
#99% CI
critval <- 1.96
#95% CI 
critval_2 <- 2.576
pred_I <- predict(mod_S_I_glm, type="response", se.fit = TRUE)
df1 <-  subset(S_glm_train_dat, select = c(I_diff, Date))
df1 <- cbind(df1, pred_I = pred_I$fit)
df1 <- cbind(df1, se = pred_I$se.fit) 
df1 <- cbind(df1, ucl = df1$pred_I + critval*df1$se)
df1 <- cbind(df1, lcl = df1$pred_I - critval*df1$se)
df1 <- cbind(df1, ucl_2 = df1$pred_I + critval_2*df1$se)
df1 <- cbind(df1, lcl_2 = df1$pred_I - critval_2*df1$se)

pred_I2 <- predict(mod_S_I_glm, newdata_I, type="response", se.fit = TRUE)
df1a <-  subset(S_glm_test_dat, select = c(I_diff, Date))
df1a <- cbind(df1a, pred_I2 = pred_I2$fit)
df1a <- cbind(df1a, se = pred_I2$se.fit) 
df1a <- cbind(df1a, ucl = df1a$pred_I2 + critval*df1a$se)
df1a <- cbind(df1a, lcl = df1a$pred_I2 - critval*df1a$se)
df1a <- cbind(df1a, ucl_2 = df1a$pred_I2 + critval_2*df1a$se)
df1a <- cbind(df1a, lcl_2 = df1a$pred_I2 - critval_2*df1a$se)

#plotting model fit vs observed   : Infection_diff   

p2 <- ggplot(df1) + 
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



p2 + theme(legend.title=element_blank()) +ggtitle("GLM Infection_diff Model, Sacramento, CA") + theme_bw() + scale_fill_manual(values = c("lightskyblue", "darkolivegreen2")) + scale_color_manual(values = c("mediumvioletred", "royalblue", "orangered3")) 


#predict D
newdata_D <- subset(S_glm_test_dat, select = c(D_diff_t_1, D_diff_t_2,
                                                boxcox_aqi, mean_Ozone_ppm, 
                                                boxcox_mean_PM2_5_µg_m3_LC,
                                                parks_percent_change_from_baseline))
S_glm_test_dat$D_pred <- predict(mod_S_D_glm, newdata_D,  type = "response")


#df is results df1 with confidence interval value columns     
#99% CI
critval <- 1.96
#95% CI 
critval_2 <- 2.576
pred_D <- predict(mod_S_D_glm, type="response", se.fit = TRUE)
df2 <-  subset(S_glm_train_dat, select = c(D_diff, Date))
df2 <- cbind(df2, pred_D = pred_D$fit)
df2 <- cbind(df2, se = pred_D$se.fit) 
df2 <- cbind(df2, ucl = df2$pred_D + critval*df2$se)
df2 <- cbind(df2, lcl = df2$pred_D - critval*df2$se)
df2 <- cbind(df2, ucl_2 = df2$pred_D + critval_2*df2$se)
df2 <- cbind(df2, lcl_2 = df2$pred_D - critval_2*df2$se)
pred_D2 <- predict(mod_S_D_glm, newdata_D, type="response", se.fit = TRUE)
df2a <-  subset(S_glm_test_dat, select = c(D_diff, Date))
df2a <- cbind(df2a, pred_D2 = pred_D2$fit)
df2a <- cbind(df2a, se = pred_D2$se.fit) 
df2a <- cbind(df2a, ucl = df2a$pred_D2 + critval*df2a$se)
df2a <- cbind(df2a, lcl = df2a$pred_D2 - critval*df2a$se)
df2a <- cbind(df2a, ucl_2 = df2a$pred_D2 + critval_2*df2a$se)
df2a <- cbind(df2a, lcl_2 = df2a$pred_D2 - critval_2*df2a$se)


#plotting model fit vs observed   : Death_diff   
#add titles
p2 <- ggplot(df2) + 
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


p2 + theme(legend.title=element_blank())  +ggtitle("GLM Death_diff Model, Sacramento, CA")+ theme_bw() + scale_fill_manual(values = c("lightskyblue", "darkolivegreen2")) + scale_color_manual(values = c("mediumvioletred", "royalblue", "orangered3")) 




#simple plot for the test week data
ggplot(S_glm_test_dat , aes(x=Date, y = Values)) + 
  geom_line(aes(y=I_diff, color="I_obs")) + 
  geom_line(aes(y=D_diff, color="D_obs")) +
  geom_line(aes(y=I_pred, color="I_pred")) + 
  geom_line(aes(y=D_pred, color="D_pred"))



##########
## MSPE ##
##########

mspe_val(S_glm_test_dat$I_diff, S_glm_test_dat$I_pred)
#64.32438
mspe_val(S_glm_test_dat$D_diff, S_glm_test_dat$D_pred)
#0.2364802

#______________________TSGLM____________________________
S_tsglm_train_dat <- Sacramento_t[c(1:98),]
S_tsglm_test_dat <- Sacramento_t[c(99:105),]

S_tsglm_reg_mat <-  data.matrix(subset(S_tsglm_train_dat, select = c(boxcox_aqi, mean_Ozone_ppm, 
                                                                       boxcox_mean_PM2_5_µg_m3_LC,
                                                                       parks_percent_change_from_baseline)))

#INFECTION MODEL  
mod_S_I_tsglm <- tsglm(S_tsglm_train_dat$I_diff,
                        model = list(past_obs = 1, past_mean = 1), link = "log", distr = "poisson",
                        xreg = S_tsglm_reg_mat)

summary(mod_S_I_tsglm)
#store 
S_tsglm_train_dat$resid_I <- residuals(mod_S_I_tsglm)
S_tsglm_train_dat$fitted_I <- fitted(mod_S_I_tsglm)
plot(mod_S_I_tsglm)
scoring(mod_S_I_tsglm)
# logarithmic   quadratic   spherical    rankprob      dawseb      normsq     sqerror 
# 6.3579425  -0.1631946  -0.3102764   6.8991870   9.6420367   7.9814486 192.8957237 

#____________Plotting, Predicting, and Evaluation_________________
#plotting model fit vs observed              
ggplot(S_tsglm_train_dat , aes(x=Date, y = Values)) + 
  geom_line(aes(y=fitted_I, color="fitted_I")) +
  geom_line(aes(y=I_diff, color="actual_I")) 
# +geom_line(aes(y=fitted_D, color="fitted_D")) + 
# geom_line(aes(y=D_diff, color="actual_D"))


#predict I
newdata_I2 <- subset(S_tsglm_test_dat, select = c(boxcox_aqi, mean_Ozone_ppm, 
                                                   boxcox_mean_PM2_5_µg_m3_LC,
                                                   parks_percent_change_from_baseline))
S_tsglm_test_dat$I_pred <- predict(mod_S_I_tsglm, n.ahead = 7, newxreg = newdata_I2)$pred
df1 <-  subset(S_tsglm_train_dat)
df1a <-  subset(S_tsglm_test_dat)

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



w + theme(legend.title=element_blank()) +ggtitle("TSGLM Infection_diff Model, Sacramento, CA") + theme_bw() + scale_fill_manual(values = c("lightskyblue", "darkolivegreen2")) + scale_color_manual(values = c("mediumvioletred", "royalblue", "orangered3")) 






#DEATH MODEL
mod_S_D_tsglm <- tsglm(S_tsglm_train_dat$D_diff,
                        model = list(past_obs = 1, past_mean = 1), link = "log", distr = "poisson",
                        xreg = S_tsglm_reg_mat)
summary(mod_S_D_tsglm)
#store 
S_tsglm_train_dat$resid_D <- residuals(mod_S_D_tsglm)
S_tsglm_train_dat$fitted_D <- fitted(mod_S_D_tsglm)
plot(mod_S_D_tsglm)
scoring(mod_S_D_tsglm)
# logarithmic   quadratic   spherical    rankprob      dawseb      normsq     sqerror 
# 0.8758707  -0.5762365  -0.7453765   0.3452923   0.1002411   1.3527198   0.8066783 

#____________Plotting, Predicting, and Evaluation_________________
#plotting model fit vs observed              
ggplot(S_tsglm_train_dat , aes(x=Date, y = Values)) + 
  geom_line(aes(y=fitted_I, color="fitted_I")) +
  geom_line(aes(y=I_diff, color="actual_I")) +
geom_line(aes(y=fitted_D, color="fitted_D")) + 
  geom_line(aes(y=D_diff, color="actual_D"))


#predict D
newdata_D2 <- subset(S_tsglm_test_dat, select = c(boxcox_aqi, mean_Ozone_ppm, 
                                                   boxcox_mean_PM2_5_µg_m3_LC,
                                                   parks_percent_change_from_baseline))
S_tsglm_test_dat$D_pred <- predict(mod_S_D_tsglm, n.ahead = 7, newxreg = newdata_D2)$pred
df1 <-  subset(S_tsglm_train_dat)
df1a <-  subset(S_tsglm_test_dat)

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



w2 + theme(legend.title=element_blank()) +ggtitle("TSGLM Death_diff Model, Sacramento, CA")+ theme_bw() + scale_fill_manual(values = c("lightskyblue", "darkolivegreen2")) + scale_color_manual(values = c("mediumvioletred", "royalblue", "orangered3")) 



##########
## MSPE ##
##########

mspe_val(S_tsglm_test_dat$I_diff, S_tsglm_test_dat$I_pred)
#178.8306
mspe_val(S_tsglm_test_dat$D_diff, S_tsglm_test_dat$D_pred)
#0.07466773








#___________TSGLM y = Death_diff with another covariate, I_diff with lag (TBD in thr function)___________
#This model is similar to the previous TSGLM model, but now we focus only on death_diff as the response
#variable and also introduced a new covariate, I_diff with a number of days as a lag (from 1-15 days).
#I will loop through model generation, altering the column associated with lagged I_diff so that for each
#iteration, the lag value will increase by 1 day (starting from -15). The mspe will be saved to a table,
#and then it will be clear which lag time produces the most effective model.

#new data frame to alter
S_tot_dat<- Sacramento_t

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
    tmp_reg_mat <-  data.matrix(subset(tmp_train, select = c(I_diff_t_n,boxcox_aqi, mean_Ozone_ppm, 
                                                             boxcox_mean_PM2_5_µg_m3_LC,
                                                             parks_percent_change_from_baseline)))
    #MODEL  
    tmp_mod_tsglm <- tsglm(tmp_train$D_diff,
                           model = list(past_obs = 1, past_mean = 1), link = "log", distr = "poisson",
                           xreg = tmp_reg_mat)
    
    #predictive/test component
    newdat <- subset(tmp_test, select = c(I_diff_t_n, boxcox_aqi, mean_Ozone_ppm, 
                                          boxcox_mean_PM2_5_µg_m3_LC,
                                          parks_percent_change_from_baseline))
    tmp_test$D_pred <- predict(tmp_mod_tsglm, n.ahead = 7, newxreg = newdat)$pred
    tmp_mspe <- mspe_val(tmp_test$D_diff, tmp_test$D_pred)
    
    total_table <- rbind(total_table, c(lag_amt, tmp_mspe))   
    
  }
  
  
  return(total_table)
}

model_compare(S_tot_dat)


#It is evident that the lowest mspe (0.002893892) is that of model with lag of 14 days for infection. Therefore, I will
#use this as a covariate in the model.
#create lag
#create covariate column for lagged infection_diff
S_tot_dat$I_diff_t_n <- S_tot_dat$I_diff
S_tot_dat$I_diff_t_n <- shift(S_tot_dat$I_diff_t_n, n=14, fill=NA, type="lag")

S_tot_train <- S_tot_dat[c(1:98),]
S_tot_test <- S_tot_dat[c(99:105),]

#get rid of first rows with NAs due to the lag created
S_tot_train<- na.omit(S_tot_train)

#covariates matrix
tmp_reg_mat <-  data.matrix(subset(S_tot_train, select = c(I_diff_t_n,boxcox_aqi, mean_Ozone_ppm, 
                                                            boxcox_mean_PM2_5_µg_m3_LC,
                                                            parks_percent_change_from_baseline)))
#MODEL  
tmp_mod_tsglm <- tsglm(S_tot_train$D_diff,
                       model = list(past_obs = 1, past_mean = 1), link = "log", distr = "poisson",
                       xreg = tmp_reg_mat)
S_tot_train$fitted_D <- fitted(tmp_mod_tsglm)
#predictive/test component
newdat <- subset(S_tot_test, select = c(I_diff_t_n,boxcox_aqi, mean_Ozone_ppm, 
                                         boxcox_mean_PM2_5_µg_m3_LC,
                                         parks_percent_change_from_baseline))
S_tot_test$D_pred <- predict(tmp_mod_tsglm, n.ahead = 7, newxreg = newdat)$pred
tmp_mspe <- mspe_val(S_tot_test$D_diff, S_tot_test$D_pred)
tmp_mspe
#0.002893892

df1 <-  subset(S_tot_train)
df1a <-  subset(S_tot_test)
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



w3 + theme(legend.title=element_blank())+ggtitle("TSGLM Death_diff Model, full- Sacramento, CA") + theme_bw() + scale_fill_manual(values = c("lightskyblue", "darkolivegreen2")) + scale_color_manual(values = c("mediumvioletred", "royalblue", "orangered3")) 



