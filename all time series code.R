##########################################################################################################
#----------------------------------- DATA CLEANING AND PREPARATION --------------------------------------#
##########################################################################################################

# load packages
library(readr)
library(dplyr)
library(tseries)
library(forecast)
library(haven)
library(fma)
library(lmtest)
library(zoo)
library(seasonal)
library(lubridate)
library(imputeTS)
library(psych,ggplot2)
library(caschrono)
library(TSA)
library(car)
library(quantmod)
library(ggplot2)
library(astsa)
library(ggthemes)
library(RColorBrewer)
library(tidyverse)
library(readxl)

# set work directory 
setwd('C:/Users/chely/Documents/Time Series and Forecasting/Class Data')

# import data sets and convert Date column from character to date type
ozone <- read_csv('Ozone_Raleigh2.csv', col_types = cols(Date = col_date(format = '%m/%d/%Y')))
NO_data <- read_csv('NO_Raleigh.csv', col_types = cols(Date = col_date(format = '%m/%d/%Y')))
so2_data <- read_csv('SO2_Raleigh.csv', col_types = cols(Date = col_date(format = '%m/%d/%Y')))
co_data <- read_csv('CO_Raleigh.csv', col_types = cols(Date = col_date(format = '%m/%d/%Y')))
weather <- read_csv('Raleigh_weather.csv', col_types = cols(DATE = col_date(format = '%m/%d/%Y')))

# Change column names of pollution conc so it is shorter and easier to handle
names(ozone)[5] <- 'maxOzoneC'
names(NO_data)[5] <- 'maxNOconc'
names(so2_data)[5] <- 'maxSO2conc'
names(co_data)[5] <- 'maxCOconc'
names(weather)[3] <- 'Date'

#----------------------------------- Check and Impute missing values-------------------------------------#

# All data sets have start date of 2014-01-01 and end date of 2020-05-31
days_seq <- seq(as.Date("2014/01/01"), as.Date("2020/05/31"), by="days") # returns date sequence
days_seq <- as.data.frame(days_seq) # convert date sequence to dataframe
colnames(days_seq) <- "Date" # change col name to "Date"

# Ozone data with missing values denoted
ozone <- left_join(days_seq, ozone, by ='Date') # 60 missing values
NO_data <- left_join(days_seq, NO_data, by ='Date') # 80 missing values
so2_data <- left_join(days_seq, so2_data, by ='Date') # 111 missing values
co_data <- left_join(days_seq, co_data, by = 'Date') # 104 missing values
weather <- left_join(days_seq, weather, by='Date') # 0 missing values

# Impute numeric missing values with either linear, spline, or stineman separately BEFORE joining
# Warning messages show for categorical variables
ozone <- na_interpolation(ozone, option = "spline")
NO_data <- na_interpolation(NO_data, option = "spline")
so2_data <- na_interpolation(so2_data, option = "spline")
co_data <- na_interpolation(co_data, option = "spline")

# Combine all 4 data sets into 1 total
all_data  <- ozone %>% full_join(NO_data,  by = 'Date') %>%
  full_join(so2_data, by = 'Date') %>%
  full_join(co_data, by = 'Date') %>%
  full_join(weather, by = 'Date')

# Remove for leap days: extract the month and day and then use negation to exclude the leap dates
all_data <- all_data[!(format(all_data$Date,"%m") == "02" & format(all_data$Date, "%d") == "29"), , drop = FALSE]

#------------------------------------- Clean data and make ts object ------------------------------------#

# Columns of interest
cols <- c('Date', 'maxOzoneC', 'maxNOconc', 'maxSO2conc', 'maxCOconc', 'PRCP', 'SNWD', 'TMAX', 'WSF2', 'WSF5', 'SNOW', 'AWND', 'TAVG', 'TMIN')

# Data with only cols of interest
all_data_clean <- all_data %>%
  select(all_of(cols))

#------------------------- EDA: check for multicollinearity among variables ------------------------------#

# Before we build a model, we need to get rid of redundant variables and also check which lags are informative for predicting our target, ozone conc. Think model selection using proc reg in SAS.

predictors <- c('maxNOconc', 'maxSO2conc', 'maxCOconc', 'PRCP', 'SNWD', 'TMAX', 'WSF2', 'WSF5', 'SNOW', 'AWND', 'TAVG', 'TMIN')

# Check pair plots of pollutants
pairs.panels(all_data_clean[,predictors[1:3]], method = 'pearson', density=T, ellipses = F)

# Check pair plots of weather metrics
pairs.panels(all_data_clean[,predictors[4:12]], method = 'pearson', density=T, ellipses = F)

# Slight correlation between NO and CO pollutants
# We only need 1 indicator for temperature --> use TAVG (good summary of min/max)
# WSF2 and WSF5 (and AWND) are correlated --> use WSF5 (better distribution)
# SNOW and SNWD are weakly correlated --> use SNOW (random pick)

updated_pred <- c('maxNOconc', 'maxSO2conc', 'maxCOconc', 'PRCP', 'WSF5', 'SNOW', 'TAVG')

# Look at correlations of pollutants and weather metrics
pairs.panels(all_data_clean[,updated_pred], method = 'pearson', density=T, ellipses = F)


#-------------------------------------------- Partition Data ---------------------------------------------#

# client requests that you have a test data set of the last 14 days and validation data set of 28 days

### Data Splits ###
training_data <- all_data_clean[1:2299,]
valid_data <- all_data_clean[2300:2327,]
test_data <- all_data_clean[2328:2341,]

#---------------------------------- Get time series objects ---------------------------------------------#

# Get time series objects for ozone
ozone_train_ts <- ts(training_data$maxOzoneC, frequency = 365)
ozone_valid_ts <- ts(valid_data$maxOzoneC, frequency = 365)
ozone_test_ts <- ts(test_data$maxOzoneC, frequency = 365)

#########################################################################################################
#---------------------------------------- SEASONAL ARIMA MODELING --------------------------------------#
#########################################################################################################

# Fit Fourier Series
reg.trig <- Arima(ozone_train_ts, xreg = fourier(ozone_train_ts, K=4))
summary(reg.trig)
# MAE = 0.006850434807, MAPE = 21.14232063

# Check residuals on trig model
# checkresiduals(arima.trig); fancier function to display residuals
plot(reg.trig$residuals, main = 'Residuals from Trig Method')

# Ndiffs() tells you how many differences are required to make the series stationarity
# Zero means you are already stationary
ndiffs(reg.trig$residuals)
# Looks like we are stationary!

trig_auto_arima <- auto.arima(reg.trig$residuals)
trig_auto_arima
# Recommends 2 nonseasonal AR terms (zero mean)

# Models with AR and MA terms
arima.trig <- Arima(ozone_train_ts, order=c(2,0,0), xreg = fourier(ozone_train_ts, K=4))
summary(arima.trig)
# MAE = 0.00580, MAPE = 17.74998

checkresiduals(arima.trig)

#----------------------------------------- White Noise for Trig Model ---------------------------------- #
White.LB_resid <- rep(NA, 20)

for(i in 1:20){
  White.LB_resid[i] <- Box.test(arima.trig$residuals, lag = i, type = "Lj", fitdf = 2)$p.value
}

# p-values >= 0.2 are recorded as 0.2 (for plotting purposes)
White.LB_resid <- pmin(White.LB_resid, 0.2)

# Let's look at a plot of these p-values (lags 1,2,...,10)
# The horizontal lines let us see which lags have p-values <0.05 and <0.01
barplot(White.LB_resid, 
        main = "Ljung-Box Test for Trig Residuals", 
        ylab = "Probabilities", 
        xlab = "Lags", 
        ylim = c(0, 0.2),
        names.arg = seq(0,19))

abline(h = 0.01, lty = "dashed", col = "black")
abline(h = 0.05, lty = "dashed", col = "black")

#------------------------------ Forecast on validation data --------------------------------------------#

# Validation forecast
validation_forecast <- forecast(arima.trig, xreg = fourier(ozone_train_ts, K=4, h=28)) # Changed this to arima.trig to use (2,0,0)

# Double-checking the length of the forecast to make sure it's 28
length(validation_forecast$mean)

# Actual minus predicted values = Error
validation_pred_error = valid_data$maxOzoneC - validation_forecast$mean

# Validation MAE: 0.004878189
MAE_valid = mean(abs(validation_pred_error))

# Validation MAPE: 12.68688
MAPE_valid = mean(abs(validation_pred_error)/abs(valid_data$maxOzoneC))*100

#------------------------ Plot actual vs predicted for the validation data ---------------------------#

# Create df with predicted values, actual values (validation data)
results_df <- data.frame(valid_pred = validation_forecast$mean[1:28],
                         actual = valid_data$maxOzoneC,
                         date = seq(as.Date('2020-04-20'), as.Date('2020-05-17'), by="days"))


# Title and legend needs to be added along with correct labels!!!
ggplot(results_df, aes(x=date)) + 
  geom_line(aes(y = valid_pred), color = "darkred") + 
  geom_line(aes(y = actual), color="steelblue", linetype="twodash")


#--------------------------------- Forecast on Testing Data -----------------------------------------#

# Testing forecast
validtesting_forecast = forecast(arima.trig, xreg = fourier(ozone_train_ts, K=4, h=42))
testing_forecast = validtesting_forecast$mean[29:42]

# Double-checking the forecast length
length(testing_forecast)

# Actual minus predicted values = Error
testing_pred_error = test_data$maxOzoneC - testing_forecast

# Testing MAE: 0.01436417
MAE_testing = mean(abs(testing_pred_error))

# Testing MAPE: 50.44742
MAPE_testing = mean(abs(testing_pred_error)/abs(test_data$maxOzoneC))*100

######################################################################################################
#------------------------------------ ARIMAX MODELING: PART ONE -------------------------------------#
######################################################################################################

# Step 1: Get preliminary regression with appropriate x-vars (only those that passed EDA screening)
ozone_model1 <- lm(maxOzoneC ~ maxNOconc + maxCOconc + PRCP + WSF5 + TMAX, 
                   data = training_data)

# Step 2: See VIFs to ensure no multicollinearity --> all are less than 10
car::vif(ozone_model1)

# Step 3: Plot the residuals from the regression model --> oscillates at 0 mean, no trend
plot(ozone_model1$residuals)

# Step 4: Check that residuals are stationary --> 0 difference required, residuals are stationary
ndiffs(ozone_model1$residuals, test = 'adf')

# Step 5: Look for potential AR, MA terms --> 4 AR, 0 diff, 1 MA terms
auto.arima(ozone_model1$residuals)

# Step 6: Put predictors in a numeric matrix w/o col names (needed for next step)
predictors <- training_data %>%
  select(maxNOconc, maxCOconc, PRCP, WSF5, TMAX) %>%
  as.matrix() %>%
  unname()

# Step 7: Build ARIMA model with error terms and predictors matrix
ozone_model2 <- Arima(ozone_train_ts, order=c(4,0,1), seasonal = c(0,0,0), xreg=predictors)
summary(ozone_model2)

### A more comprehensive look at ARIMA model ###
astsa::sarima(ozone_train_ts, p=4, d=0, q=1, P = 0, D = 0, Q = 0, 
              details = TRUE, xreg=predictors)

# Step 8: plot/compute correlation function on the residuals: both are relatively withing C.I.
Acf(ozone_model2$residuals)$acf
Pacf(ozone_model2$residuals)$acf

# Step 9: Check WN-- Good!
White.LB <- rep(NA, 20)
for(i in 1:20){
  White.LB[i] <- Box.test(ozone_model2$residuals, lag = i, type = "Lj", fitdf = 5)$p.value
}
White.LB <- pmin(White.LB, 0.2)
barplot(White.LB, 
        main = "Ljung-Box Test for Series P-values", 
        ylab = "Probabilities", 
        xlab = "Lags", 
        ylim = c(0, 0.2),
        names.arg = seq(0,19))

abline(h = 0.01, lty = "dashed", col = "black")
abline(h = 0.05, lty = "dashed", col = "black")

# want all lags to be insignificant and fail to reject OVERALL --> no autocorrelation of residuals!

#------------------------------- Lags of x-vars can also predict target -------------------------------#

# Step 10: Get lags of each x-var (going back 7 lags is a good start for daily data)
no1 <- lag(training_data$maxNOconc,1)
no2 <- lag(training_data$maxNOconc,2)
no3 <- lag(training_data$maxNOconc,3)
no4 <- lag(training_data$maxNOconc,4)
no5 <- lag(training_data$maxNOconc,5)
no6 <- lag(training_data$maxNOconc,6)
no7 <- lag(training_data$maxNOconc,7)

co1 <- lag(training_data$maxCOconc,1)
co2 <- lag(training_data$maxCOconc,2)
co3 <- lag(training_data$maxCOconc,3)
co4 <- lag(training_data$maxCOconc,4)
co5 <- lag(training_data$maxCOconc,5)
co6 <- lag(training_data$maxCOconc,6)
co7 <- lag(training_data$maxCOconc,7)

prcp1 <- lag(training_data$PRCP,1)
prcp2 <- lag(training_data$PRCP,2)
prcp3 <- lag(training_data$PRCP,3)
prcp4 <- lag(training_data$PRCP,4)
prcp5 <- lag(training_data$PRCP,5)
prcp6 <- lag(training_data$PRCP,6)
prcp7 <- lag(training_data$PRCP,7)

wsf5_1 <- lag(training_data$WSF5,1)
wsf5_2 <- lag(training_data$WSF5,2)
wsf5_3 <- lag(training_data$WSF5,3)
wsf5_4 <- lag(training_data$WSF5,4)
wsf5_5 <- lag(training_data$WSF5,5)
wsf5_6 <- lag(training_data$WSF5,6)
wsf5_7 <- lag(training_data$WSF5,7)

tmax1 <- lag(training_data$TMAX,1)
tmax2 <- lag(training_data$TMAX,2)
tmax3 <- lag(training_data$TMAX,3)
tmax4 <- lag(training_data$TMAX,4)
tmax5 <- lag(training_data$TMAX,5)
tmax6 <- lag(training_data$TMAX,6)
tmax7 <- lag(training_data$TMAX,7)

#----------------------------- Use backwards selection to get final model -------------------------#

# Step 11: Add all x-vars and their lags along with target to a new dataset
all.x.reg <- training_data %>%
  select(maxOzoneC, maxNOconc, maxCOconc, PRCP, WSF5, TMAX) %>%
  cbind(no1, no2, no3, no4, no5, no6, no7, 
        co1, co2, co3, co4, co5, co6, co7, 
        prcp1, prcp2, prcp3, prcp4, prcp5, prcp6, prcp7, 
        wsf5_1, wsf5_2, wsf5_3, wsf5_4, wsf5_5, wsf5_6, wsf5_7, 
        tmax1, tmax2, tmax3, tmax4, tmax5, tmax6, tmax7)

# Remove missing values in order to do stepwise regression
xreg_nonulls <- na.omit(all.x.reg)
model.all <- lm(maxOzoneC ~ ., data = xreg_nonulls)
step(model.all, direction = "both")

# Get col indices of selected x-vars using data.frame(colnames(xreg_nonulls))
choose.x <- match(c('maxNOconc', 'maxCOconc', 'PRCP', 'WSF5', 'TMAX', 'no1', 'no6', 'no7', 'co1', 'co3', 
                    'co5', 'co6', 'co7', 'prcp1', 'prcp2', 'prcp4', 'prcp6', 'prcp7', 'wsf5_1', 'wsf5_2',
                    'wsf5_4', 'wsf5_5', 'wsf5_6', 'wsf5_7', 'tmax3', 'tmax5'), names(xreg_nonulls))

# Subset x-vars from df and put into numeric matrix
x.arima <- as.matrix(all.x.reg[,choose.x])

# Final regression model
xreg.final <- lm(maxOzoneC ~ maxNOconc + maxCOconc + PRCP + WSF5 + 
                   TMAX + no1 + no6 + no7 + co1 + co3 + co5 + co6 + co7 + prcp1 + 
                   prcp2 + prcp4 + prcp6 + prcp7 + wsf5_1 + wsf5_2 + wsf5_4 + 
                   wsf5_5 + wsf5_6 + wsf5_7 + tmax3 + tmax5, data = xreg_nonulls)

plot(xreg.final$residuals) # plot residuals
ndiffs(xreg.final$residuals, test = 'adf') # residuals are stationary
auto.arima(xreg.final$residuals) # suggests 2 AR, 1 MA terms

# Final ARIMAX model
final.arimax <- Arima(ozone_train_ts, order=c(2,0,1), xreg=x.arima)
# fixed arg only needed for skipping ar or ma terms
# use plots to know which terms to skip over as an option

# Check stationarity: stationary residuals!
ndiffs(final.arimax$residuals, test = 'adf')
Acf(final.arimax$residuals)$acf
Pacf(final.arimax$residuals)$acf

# Check WN-- Good!
White.LB <- rep(NA, 20)
for(i in 1:20){
  White.LB[i] <- Box.test(final.arimax$residuals, lag = i, type = "Lj", fitdf = 3)$p.value
}
White.LB <- pmin(White.LB, 0.2)
barplot(White.LB, 
        main = "Ljung-Box Test for Series P-values", 
        ylab = "Probabilities", 
        xlab = "Lags", 
        ylim = c(0, 0.2),
        names.arg = seq(0,19))
abline(h = 0.01, lty = "dashed", col = "black")
abline(h = 0.05, lty = "dashed", col = "black")

######################################################################################################
#------------------------------------ ARIMAX MODELING: PART TWO -------------------------------------#
######################################################################################################

# Get time series objects of all relevant x-vars

# NO ts objects
NO_train_ts <- ts(training_data$maxNOconc, frequency = 365)
NO_valid_ts <- ts(valid_data$maxNOconc, frequency = 365)
NO_test_ts <- ts(test_data$maxNOconc, frequency = 365)

# CO ts objects
co_train_ts <- ts(training_data$maxCOconc, frequency = 365)
co_valid_ts <- ts(valid_data$maxCOconc, frequency = 365)
co_test_ts <- ts(test_data$maxCOconc, frequency = 365)

# PRCP ts objects
prcp_train_ts <- ts(training_data$PRCP, frequency = 365)
prcp_valid_ts <- ts(valid_data$PRCP, frequency = 365)
prcp_test_ts <- ts(test_data$PRCP, frequency = 365)

# WSF5 ts objects
wsf5_train_ts <- ts(training_data$WSF5, frequency = 365)
wsf5_valid_ts <- ts(valid_data$WSF5, frequency = 365)
wsf5_test_ts <- ts(test_data$WSF5, frequency = 365)

# SNOW ts objects
snow_train_ts <- ts(training_data$SNOW, frequency = 365)
snow_valid_ts <- ts(valid_data$SNOW, frequency = 365)
snow_test_ts <- ts(test_data$SNOW, frequency = 365)

# TMAX ts objects
tmax_train_ts <- ts(training_data$TMAX, frequency = 365)
tmax_valid_ts <- ts(valid_data$TMAX, frequency = 365)
tmax_test_ts <- ts(test_data$TMAX, frequency = 365)

# We will use the actual validation and test data for weather provided by the client because in the creation of our models because meteorologists will have the best predictions.

# Thus, find models for NO and CO to be plugged into the final ARIMAX model.

#----------------------------------------------- NO modeling ----------------------------------------#

#Test Arima model with seasonality
reg.trig <- Arima(NO_train_ts, xreg = fourier(NO_train_ts, K=4))
summary(reg.trig)

#Auto.arima to get PDQ terms
auto.arima(NO_train_ts) #(3,1,1)

#Creating Sin/Cosine Values
index.ts=seq(1,length(NO_train_ts))
x1.sin=sin(2*pi*index.ts*1/365)
x1.cos=cos(2*pi*index.ts*1/365)
x2.sin=sin(2*pi*index.ts*2/365)
x2.cos=cos(2*pi*index.ts*2/365)
x3.sin=sin(2*pi*index.ts*3/365)
x3.cos=cos(2*pi*index.ts*3/365)
x4.sin=sin(2*pi*index.ts*4/365)
x4.cos=cos(2*pi*index.ts*4/365)
x.reg=cbind(x1.sin,x1.cos,x2.sin,x2.cos,x3.sin,x3.cos,x4.sin,x4.cos)

#Trig Model, based off AR/MA terms from Auto.Arima
trig.1<-Arima(NO_train_ts,order=c(3,1,1),xreg=x.reg)
summary(trig.1) #67.47 #white noise worked!

#Check for white noise!
White.LB <- rep(NA, 20)
for(i in 1:20){
  White.LB[i] <- Box.test(trig.1$residuals, lag = i, type = "Lj", fitdf = 5)$p.value
}
White.LB <- pmin(White.LB, 0.2)
barplot(White.LB, 
        main = "Ljung-Box Test for Series P-values", 
        ylab = "Probabilities", 
        xlab = "Lags", 
        ylim = c(0, 0.2),
        names.arg = seq(0,19))
abline(h = 0.01, lty = "dashed", col = "black")
abline(h = 0.05, lty = "dashed", col = "black")
#passed the test


NO_training_forecast <- forecast(trig.1, xreg = x.reg, h = 56)
no_pred_v <- NO_training_forecast$mean[1:28] # predictions for validation
no_pred_test <- NO_training_forecast$mean[29:42] # predictions for testing data
no_pred_june <- NO_training_forecast$mean[43:56] # predictions for June


train_pred_error = training_data$maxNOconc - training_forecast$mean
MAE_train = mean(abs(train_pred_error))
print(MAE_train)


#Validation MAPE
MAPE_train = mean(abs(train_pred_error)/abs(training_data$maxNOconc))*100
print(MAPE_train)
#------------------------------------------- CO Modeling ---------------------------------------------#

plot(co_train_ts) # deterministic pulse intervention around 4

Acf(co_train_ts, lag = 35)$acf
Pacf(co_train_ts, lag = 21)$pacf

decomp_co <- stl(co_train_ts, s.window = 20)
plot(decomp_co, main = 'Time Plot of Daily Average CO (ppm)')

# First model
# auto.arima(co_train_ts) # suggests 2 AR, 0 diff, 3 MA
# co_model1 <- Arima(co_train_ts, order=c(2,0,3))
# checkresiduals(co_model1$residuals)
# ndiffs(co_model1$residuals) # stationary residuals
# summary(co_model1) # MAE = 0.0995242, MAPE = 35.02221

training_data %>%
  filter(training_data$co_i==1)

# Get column for coding intervention var
training_data$co_i <- ifelse(training_data$maxCOconc>2,1,0)
# Get model with intervention col
co_model2 <- forecast::Arima(co_train_ts, xreg = training_data$co_i, method = "ML")
checkresiduals(co_model2$residuals)
ndiffs(co_model2$residuals) # stationary residuals
auto.arima(co_model2$residuals) # 4 AR, 0 diff, 1 MA
co_model3 <- Arima(co_train_ts, order=c(4,0,1), xreg = training_data$co_i, method = "ML")
summary(co_model3) # MAE = 0.09868918 , MAPE = 34.44208 

# Check WN for first model-- good!
# White.LB <- rep(NA, 20)
# for(i in 1:20){
#   White.LB[i] <- Box.test(co_model1$residuals, lag = i, type = "Lj", fitdf = 5)$p.value
# }
# White.LB <- pmin(White.LB, 0.2)
# barplot(White.LB, 
#         main = "Ljung-Box Test for Series P-values", 
#         ylab = "Probabilities", 
#         xlab = "Lags", 
#         ylim = c(0, 0.2),
#         names.arg = seq(0,19))
# abline(h = 0.01, lty = "dashed", col = "black")
# abline(h = 0.05, lty = "dashed", col = "black")


# Check WN for second model-- good!
White.LB <- rep(NA, 20)
for(i in 1:20){
  White.LB[i] <- Box.test(co_model3$residuals, lag = i, type = "Lj", fitdf = 5)$p.value
}
White.LB <- pmin(White.LB, 0.2)
barplot(White.LB, 
        main = "Ljung-Box Test for Series P-values", 
        ylab = "Probabilities", 
        xlab = "Lags", 
        ylim = c(0, 0.2),
        names.arg = seq(0,19))
abline(h = 0.01, lty = "dashed", col = "black")
abline(h = 0.05, lty = "dashed", col = "black")

# Forecast validation and test data
co_pred <- forecast(co_model3, xreg = training_data$co_i, h = 56)
co_pred_v <- co_pred$mean[1:28] # predictions for validation
co_pred_test <- co_pred$mean[29:42] # predictions for testing data
co_pred_june <- co_pred$mean[43:56] # predictions for June

######################################################################################################
#------------------------------------ ARIMAX MODELING: PART THREE -----------------------------------#
######################################################################################################

#---------------------------- FORECAST OZONE VALIDATION WITH PREDICTED X VARS -----------------------#
# Significant lags using validation data from stepwise model selection
no1 <- lag(no_pred_v,1)
no6 <- lag(no_pred_v,6)
no7 <- lag(no_pred_v,7)

co1 <- lag(co_pred_v,1)
co3 <- lag(co_pred_v,3)
co5 <- lag(co_pred_v,5)
co6 <- lag(co_pred_v,6)
co7 <- lag(co_pred_v,7)

prcp1 <- lag(valid_data$PRCP,1)
prcp2 <- lag(valid_data$PRCP,2)
prcp4 <- lag(valid_data$PRCP,4)
prcp6 <- lag(valid_data$PRCP,6)
prcp7 <- lag(valid_data$PRCP,7)

wsf5_1 <- lag(valid_data$WSF5,1)
wsf5_2 <- lag(valid_data$WSF5,2)
wsf5_4 <- lag(valid_data$WSF5,4)
wsf5_5 <- lag(valid_data$WSF5,5)
wsf5_6 <- lag(valid_data$WSF5,6)
wsf5_7 <- lag(valid_data$WSF5,7)

tmax3 <- lag(valid_data$TMAX,3)
tmax5 <- lag(valid_data$TMAX,5)

#----------------------------------------------------------------------------------------------------#

# X vars to include in x.arima
maxNOconc <- no_pred_v # predicted
maxCOconc <- co_pred_v # predicted
PRCP <- valid_data$PRCP
WSF5 <- valid_data$WSF5
TMAX <- valid_data$TMAX

# get external regressors with actual weather validation and predicted pollutants validation along with lags using validation data
x.arima_v <- cbind(maxNOconc, maxCOconc, PRCP, WSF5, TMAX, no1, no6, no7, co1, co3, co5, co6, co7, 
                   prcp1, prcp2, prcp4, prcp6, prcp7, wsf5_1, wsf5_2, wsf5_4, wsf5_5, 
                   wsf5_6, wsf5_7, tmax3, tmax5)

# Need to fill missing values with previous data from training in order to use forecast()

# Forecast validation and test data
arimax_pred <- forecast(final.arimax, xreg = x.arima_v, h = 28)
# have same regressors in same order and title

arimax_pred_v <- arimax_pred$mean # predictions for validation

# Calculate validation data error
valid_error <- valid_data$maxOzoneC - arimax_pred_v
# MAE = 0.003837417
MAE_v <- mean(abs(valid_error))
# MAPE = 9.54276
MAPE_v <- mean(abs(valid_error)/abs(valid_data$maxOzoneC))*100

#------------------------------- FORECAST OZONE TEST WITH PREDICTED X VARS --------------------------#
# Significant lags using test data from stepwise model selection
no1 <- lag(no_pred_test,1)
no6 <- lag(no_pred_test,6)
no7 <- lag(no_pred_test,7)

co1 <- lag(co_pred_test,1)
co3 <- lag(co_pred_test,3)
co5 <- lag(co_pred_test,5)
co6 <- lag(co_pred_test,6)
co7 <- lag(co_pred_test,7)

prcp1 <- lag(test_data$PRCP,1)
prcp2 <- lag(test_data$PRCP,2)
prcp4 <- lag(test_data$PRCP,4)
prcp6 <- lag(test_data$PRCP,6)
prcp7 <- lag(test_data$PRCP,7)

wsf5_1 <- lag(test_data$WSF5,1)
wsf5_2 <- lag(test_data$WSF5,2)
wsf5_4 <- lag(test_data$WSF5,4)
wsf5_5 <- lag(test_data$WSF5,5)
wsf5_6 <- lag(test_data$WSF5,6)
wsf5_7 <- lag(test_data$WSF5,7)

tmax3 <- lag(test_data$TMAX,3)
tmax5 <- lag(test_data$TMAX,5)


# X vars to include in x.arima
maxNOconc <- no_pred_test # predicted
maxCOconc <- co_pred_test # predicted

maxNOconc <- no_pred_june # predicted
maxCOconc <- co_pred_june # predicted
PRCP <- test_data$PRCP
WSF5 <- test_data$WSF5
TMAX <- test_data$TMAX

x.arima_test <- cbind(maxNOconc, maxCOconc, PRCP, WSF5, TMAX, no1, no6, no7, co1, co3, co5, co6, co7, 
                      prcp1, prcp2, prcp4, prcp6, prcp7, wsf5_1, wsf5_2, wsf5_4, wsf5_5, 
                      wsf5_6, wsf5_7, tmax3, tmax5)

# Need to fill missing values with data from validation in order to use forecast()

# Forecast validation and test data
arimax_pred <- forecast(final.arimax, xreg = x.arima_test, h = 14)
arimax_june <- arimax_pred$mean
arimax_pred_test <- arimax_pred$mean # predictions for test

# Calculate test data error
test_error <- test_data$maxOzoneC - arimax_pred_test
# MAE = 0.00883002
MAE_test <- mean(abs(test_error))
# MAPE = 32.07289
MAPE_test <- mean(abs(test_error)/abs(test_data$maxOzoneC))*100

######################################################################################################
#----------------------------------------- UCM MODELING IN SAS  -------------------------------------#
######################################################################################################

#   
#   *Create library;
# libname time "/home/u37604878/Time Series";
# 
# *Append training/validation;
# proc append base=time.training data=time.valid force;
# run;
# 
# /* VISUALIZE OZONE CONCENTRATION */
#   %macro ts(s=maxOzoneC);
# proc sgplot data=time.training;
# series x=date y=&s;
# run;
# %mend;
# quit;
# %ts;
# 
# 
# /* FITTING FOR OZONE CONCENTRATION */
#   *Start with level model;
# proc ucm data=time.training plots=all;
# level;
# irregular;
# model maxOzonec;
# forecast back=28 lead=28 plot=forecasts;
# run;
# 
# *Check potential trig modeling for seasonality;
# proc ucm data=time.training plots=all;
# level;
# irregular;
# model maxOzonec;
# season length=365 type=trig keeph=1 to 7;
# estimate plot=wn;
# forecast back=28 lead=28 plot=forecasts;
# run;
# 
# proc ucm data=time.training plots=all;
# level;
# irregular;
# model maxOzonec;
# estimate plot=wn;
# splineseason length=365 knots=2, 3;
# forecast back=28 lead=28 plot=forecasts;
# run;
# 
# *Evaluate trend component;
# proc ucm data=time.training plots=residuals(loess);
# level;
# irregular;
# model maxOzonec;
# estimate plot=wn;
# slope;
# forecast back=28 lead=28 plot=forecasts;
# run;
# 
# *See plots of potential predictors;
# %ts;
# %ts(s=maxNOconc);
# %ts(s=maxSO2conc);
# %ts(s=maxCOconc);
# %ts(s=TMAX);
# %ts(s=PRCP);
# %ts(s=WSF5);
# %ts(s=AWND);
# 
# *Regression model for predictors of ozone. Use regressors from ARIMA section;
# %let vars = maxNOconc maxSO2conc maxCOconc TMAX TAVG PRCP WSF2 WSF5 AWND;
# proc glmselect data=time.training plots=all;
# Backward: model maxOzonec=&vars/selection=backward select=SL slstay=0.01;
# run;
# 
# *Model Xs from model;
# proc ucm data=time.training plots=residuals(loess);
# level;
# irregular;
# model maxOzonec = maxNOconc maxCOconc TMAX PRCP WSF5;
# estimate plot=wn;
# forecast back=28 lead=28 plot=forecasts;
# run;
# 
# *Test various AR, MA, seasonal terms based on autocorrelation plots for white noise;
# proc ucm data=time.training plots=residuals(acf pacf loess);
# level;
# irregular p=12; * Tested multiples of AR, MA, and seasonal terms: p = x q = x sp = x sq = x;
# model maxOzonec=maxNOconc maxCOconc TMAX PRCP WSF5;
# estimate plot=wn;
# forecast back=28 lead=28 plot=forecasts;
# run;
# 
# 
# *Append training/test;
# proc append base=time.training data=time.test force;
# run;
# 
# *Forecast and export validation/test;
# ods excel file="/home/u37604878/Time Series/UCM_VT.xlsx";
# proc ucm data=time.training plots=none;
# level;
# irregular p=12;
# model maxOzonec=maxNOconc maxCOconc TMAX PRCP WSF5;
# estimate plot=wn;
# forecast back=42 lead=42 plot=forecasts;
# run;
# ods excel close;
# 
# ods excel file="/home/u37604878/Time Series/UCM_T.xlsx";
# proc ucm data=time.training plots=none;
# level;
# irregular p=12;
# model maxOzonec=maxNOconc maxCOconc TMAX PRCP WSF5;
# estimate plot=wn;
# forecast back=14 lead=14 plot=forecasts;
# run;
# ods excel close;

######################################################################################################
#---------------------------------------- ENSEMBLE MODELING -----------------------------------------#
######################################################################################################
# Ensemble Model
# cbind validation sets from all models to get ensemble validation data
# cbind test sets from all models to get ensemble test data
# use rowMeans(df) to get average daily validation data

# import UCM and Seasonal ARIMA validation and test data sets
UCM_v <- read_excel('C:/Users/chely/Documents/Time Series and Forecasting/FINAL PROJECT/UCM_Validation.xlsx')
UCM_test <- read_excel('C:/Users/chely/Documents/Time Series and Forecasting/FINAL PROJECT/UCM_Test.xlsx')
sarima_v <- read_csv('C:/Users/chely/Documents/Time Series and Forecasting/FINAL PROJECT/validation_forecast_SARIMA.csv')
sarima_test <- read_csv('C:/Users/chely/Documents/Time Series and Forecasting/FINAL PROJECT/testing_forecast_SARIMA.csv')

# Validation data sets
UCM_v <- UCM_v$`Forecast (Validation`
sarima_v <- sarima_v$x
arimax_v <- arimax_pred_v # run arimax code

# Test data sets
UCM_test <- UCM_test$`Forecast (test)`
UCM_test <- UCM_test[1:14]
sarima_test <- sarima_test$x
arimax_test <- arimax_pred_test # run arimax code

# Combine all validation data sets into a df
df_v <- cbind(UCM_v, sarima_v, arimax_v)

# Average the validation data sets to get ensemble validation data
ensemble_v <- rowMeans(df_v)

# Combine all test data sets into a df
df_test <- cbind(UCM_test, sarima_test, arimax_test)

# Average the test data sets to get ensemble test data
ensemble_test <- rowMeans(df_test)

#-------------------------------------- MAE and MAPE for validation ---------------------------------#

# Calculate validation data error
en_valid_error <- valid_data$maxOzoneC - ensemble_v
# MAE = 0.00369162
MAE_en_valid <- mean(abs(en_valid_error))
# MAPE = 9.387289
MAPE_en_valid <- mean(abs(en_valid_error)/abs(valid_data$maxOzoneC))*100

#----------------------------------------- MAE and MAPE for test ------------------------------------#
# Calculate test data error
en_test_error <- test_data$maxOzoneC - ensemble_test
# MAE = 0.01097605
MAE_en_test <- mean(abs(en_test_error))
# MAPE = 38.93559
MAPE_en_test <- mean(abs(en_test_error)/abs(test_data$maxOzoneC))*100

######################################################################################################
#----------------------------------------- PLOTTING RESULTS -----------------------------------------#
######################################################################################################

# test validation and test data sets from all models
test <- read_csv('C:/Users/chely/Documents/Time Series and Forecasting/FINAL PROJECT/all_test_data.csv')
valid <- read_csv('C:/Users/chely/Documents/Time Series and Forecasting/FINAL PROJECT/all_validation_data.csv')

valid$Date <- seq(as.Date('2020-4-20'), as.Date('2020-5-17'), by = 'days')
test$Date <- seq(as.Date('2020-5-18'), as.Date('2020-5-31'), by = 'days')

valid %>% select(2:7) %>% 
  rename(Actual = actual,
         `ARIMAX` = arimax_v,
         `Ensemble` = ensemble_v, 
         `SARIMA` = sarima_v,
         UCM = UCM_v)%>% 
  pivot_longer(-Date, 'Model', 'value') %>% 
  ggplot(aes(x = Date, y = value, color = Model, size = Model, linetype = Model)) +
  geom_line() +
  #ylim(c(.025,.055))+
  scale_size_manual(values = c(1,1,.5,.5,.5)) +
  scale_color_manual(values = c('#e94f37', '#385a7c', 'black', 'black', 'black'))+
  scale_linetype_manual(values = c(1,1,2,3,4))+
  labs(title = 'Predicted versus Actual Daily Ozone Concentrations: Validation Data',
       #subtitle = 'Validation Data',
       x = 'Data', 
       y = 'Daily Ozone Concentration')+
  theme_excel_new() 

#ggsave('~/MSA 21/AA 502/Time Series/Final/valid.png', width = 9, height = 4)


test %>% select(c(4,6,7)) %>%
  rename(ARIMAX = arimax_test, 
         Actual = actual) %>% 
  pivot_longer(-Date, 'Model', 'value') %>% 
  ggplot(aes(x = Date, y = value, color = Model)) + 
  #ylim(c(.02, .06))+
  geom_line(size = 1) +
  labs(title = 'Predicted versus Actual Daily Ozone Concentrations: Test Data',
       #subtitle = 'Validation Data',
       x = 'Data', 
       y = 'Daily Ozone Concentration')+
  scale_color_manual(values = c('#e94f37', '#385a7c'))+
  theme_excel_new()

#ggsave('~/MSA 21/AA 502/Time Series/Final/test.png', width = 9, height = 4)



tmp <-  valid %>% select(c(4,6,7)) %>% rbind(select(test, c(4,6,7)), use.names = F) %>% 
  rename(Actual = actual,
         `ARIMAX` = arimax_v)
tmp %>% 
  pivot_longer(-Date, 'Model', 'value') %>% 
  ggplot(aes(x = Date, y = value, color = Model)) +
  geom_line(size = 1) + 
  labs(title = 'Predicted versus Actual Daily Ozone Concentrations: Validation and Test Data',
       #subtitle = 'Validation Data',
       x = 'Data', 
       y = 'Daily Ozone Concentration')+
  scale_color_manual(values = c('#e94f37', '#385a7c')) +
  geom_vline(xintercept = as.numeric(tmp$Date[29]), linetype = 4) + 
  theme_excel_new()

#ggsave('~/MSA 21/AA 502/Time Series/Final/valid_test.png', width = 9, height = 4)
