library(readxl)
library(astsa)
library(dplyr)


### load data set
fatalities <- read_excel("C:/Users/eberl/OneDrive/Desktop/Semester I/Time Series/fatalities.xlsx")
fatalities$seatbelt <- ifelse(fatalities$Year>1967, 1 , 0)
fatalities$drinking <- ifelse(fatalities$Year>1984, 1 , 0)

### Looking at main TS
par(mfrow=c(1,1))
plot(fatalities$Year, fatalities$Deaths, xlab='Year',ylab='Fatalities', main='Motor Vehicle Fatalities', 
     type='l',col='royalblue')
plot( fatalities$Year, fatalities$Population, xlab='Year',ylab='Population',type="l", 
      main = 'US Population', col="royalblue" )
plot( fatalities$Year, fatalities$`Vehicle miles traveled (billions)`,xlab='Year',ylab='VMT (Billions)',
      main='Vehicle Miles Traveled (Billions)', type="l", col="royalblue" )
plot( fatalities$Year, fatalities$Motorcycle,xlab='Year',ylab='Fatalities', 
      main = 'Motorcycle Fatalities', type="l", col="royalblue" )
plot( fatalities$Year, fatalities$Bicycle,xlab='Year',ylab='Fatalities', 
      main = 'Bicycle Fatalities', type="l", col="royalblue" )

# Graphs with perspective
plot(fatalities$Year, fatalities$`Fatalities per 100 million VMT`, xlab='Year', 
     ylab='Fatalities per 100 mil VMT', main='Fatalities per 100 million VMT', type='l', col='royalblue')
plot(fatalities$Year, fatalities$`Fatalities per 100,000 population`, xlab='Year', 
     ylab='Fatalities per 100,000 Population',main='Fatalities per 100,000 population', type='l',col='royalblue')
plot(fatalities$Year, fatalities$`Change in per capita fatalities from previous year`, xlab='Year',
     ylab="Change in per capita fatalities", main='Change in per capita Fatalities from Previous Year',type='l', col='royalblue')


### Looking at trends

plot(fatalities$Year, fatalities$Population)
plot(fatalities$Year, fatalities$`Vehicle miles traveled (billions)`) #driving way more
plot(fatalities$Year, fatalities$`Vehicle miles traveled (billions)`/fatalities$Population)
plot(fatalities$Year, fatalities$Population) #very linear, increasing
plot(fatalities$Year, fatalities$`Fatalities per 100 million VMT`) #VMD increases so much that fataliteis drop
plot(fatalities$Year, fatalities$`Fatalities per 100,000 population`) #been decreasing for a bit

plot(fatalities$Deaths, fatalities$Population)
plot(fatalities$Deaths, fatalities$`Vehicle miles traveled (billions)`)
plot(fatalities$Deaths, fatalities$`Fatalities per 100 million VMT`) ## Negative trend
plot(fatalities$Deaths, fatalities$`Fatalities per 100,000 population`) ## Positive trend
plot(fatalities$Deaths, fatalities$Motorcycle) 
plot(fatalities$Deaths, fatalities$Bicycle) # A little bit of a positive trend

### Checking out stationarity
acf2(fatalities$Deaths) # possible AR(1) process?
acf2(diff(d)) # stationary
acf2(log(fatalities$Deaths)) # don't need to log

diff<- diff(fatalities$Deaths)
plot(diff, ylab='Differenced Once')
diff2<- diff(diff(fatalities$Deaths))
plot(diff2, ylab='Differenced Twice')
diff3<- diff(diff(diff(fatalities$Deaths)))
plot(diff3, ylab='Differenced Thrice')
# conclude - differencing once is best

### Finding a sarima model for motorized vehicular deaths
d<-fatalities$Deaths
plot(d)
# No seasonal, so D=0
# d=2 and d=3 make the plot look better
# d=1 has the best looking pacf. 
#[1,2] means P=0, Q=1 is recommended by AIC.

# pq
uplim=4
aicmat=matrix(double((uplim+1)^2),uplim+1,uplim+1)
for (i in 0:uplim){
  for (j in 0:uplim){
    aicmat[i+1,j+1]<-sarima(d,i,1,j,0,0,0,0,details=F)$AIC
    print(aicmat)}} #[1,2] means p=0, q=1 is recommended by AIC.

# Final sarima model
sarima(d,0,1,1,0,0,0,0) # AIC = 18.33253

# If you difference twice, AIC = 18.37316 [1,5] but a smaller model is AIC = 18.37554 for [1,3]
sarima(d,0,2,2,0,0,0,0)
# If you difference thrice, AIC = 18.47656 [1,4]
sarima(d,0,2,3,0,0,0,0)


### Trying another sarima model - reverse order finding and using seasonal
#Looking for periodic behavior.
deaths<-fatalities$Deaths
x=deaths-mean(deaths) #Centering.
plot(x)
n=length(x)
I=abs(fft(x))^2/n
P=(4/n)*I[1:(n/2)]
plot((0:(n/2-1))/n,P,type='o',xlab='Frequency',ylab='Scaled Periodogram')

which.max(P) # first fundamental freq 1/n
period=n/1
period 

# seasonality?
plot(d)
Dbirth<-diff(d,2)
Dbirth<-diff(d,3)
Dbirth<-diff(d,4)
Dbirth<-diff(d,5)
Dbirth<-diff(d,6)
Dbirth<-diff(d,7)
Dbirth<-diff(d,8)
Dbirth<-diff(d,9)
Dbirth<-diff(d,10)
Dbirth<-diff(d,11) 
Dbirth<-diff(d,12) #12 looks ok, but 7 looks better I think
acf2(Dbirth)

diffd<-diff(Dbirth) #d=1 looks best
plot(diffd)  
acf2(diffd) 

# PQ RUN AGAIN
uplim=4
aicmat2=matrix(double((uplim+1)^2),uplim+1,uplim+1)
for (i in 0:uplim){
  for (j in 0:uplim){
    aicmat2[i+1,j+1]=sarima(d,0,1,0,i,1,j,7,details=F, tol=.001)$AIC
    print(aicmat2)}} #[1,2] means P=0, Q=1 is recommended by AIC.

#pq
uplim=4
aicmat=matrix(double((uplim+1)^2),uplim+1,uplim+1)
for (i in 0:uplim){
  for (j in 0:uplim){
    aicmat[i+1,j+1]<-sarima(d,i,1,j,0,1,1,7,details=F)$AIC
    print(aicmat)}} #[1,2] means p=0, q=1 is recommended by AIC.

# Final sarima model
sarima(d,0,1,1,0,1,1,7) # AIC = 17.57777 
# Very low p-values. Doesn't make sense. Dismiss seasonality. 

### CCF's - consider lags in liner model
plot(diff(fatalities$`Vehicle miles traveled (billions)`))
ccf(diff(newfatalities$Deaths),diff(newfatalities$`Vehicle miles traveled (billions)`))
ccf(diff(newfatalities$Deaths),diff(diff(newfatalities$`Vehicle miles traveled (billions)`)))

fullset <- fatalities %>% filter(Year > 1979)
ccf(fullset$Deaths,fullset$Bicycle, main='CCF of Motor Vehicle and Bicycle Deaths')
ccf(fullset$Deaths,fullset$Year, main='CCF of Motor Vehicle and Year')
ccf(fullset$Deaths,fullset$Motorcycle, main='CCF of Motor Vehicle and Motorcycle Deaths')
ccf(fullset$Deaths,fullset$Population, main='CCF of Motor Vehicle and Population')
ccf(fatalities$Deaths,fatalities$seatbelt, main='CCF of Motor Vehicle and Seatbelt') 
#Use deaths to predict future seatbelt values - USELESS. Can vaguely use seatbelt values to predict deaths at close lags

# differenced (stationary) ccfs 
ccffatalities<-fatalities %>% filter(Year > 1899)
diffmot<-diff(fityears$Motorcycle)
ccf(diffd,diffmot, main='CCF of Fatalities and Motorcycle Deaths - Both Differenced')
diffbike<-diff(fityears$Bicycle)
ccf(diffd,diffbike, main='CCF of Fatalities and Bicyle Deaths - Both Differenced')
#population doesn't look like white noise after any amount of differencing

cor(newfatalities$Year, newfatalities$Population)
#ccf(fatalities$Deaths,fatalities$Population)
ccf(fatalities$Deaths,fatalities$seatbelt,main='CCF of Fatalities and Seatbelt Indicator')
ccf(fatalities$Deaths,fatalities$drinking,main='CCF of Fatalities and Drinking Indicator')

ccffatalities<-fatalities %>% filter(Year > 1920)
ccf(ccffatalities$Population,ccffatalities$`Vehicle miles traveled (billions)`, main='CCF of Population and VMT')
cor(ccffatalities$Population,ccffatalities$`Vehicle miles traveled (billions)`)
# Looking at year a predictor - especially with lags
# Include year. No lagged variables. 


### Linear Model with ARMA Errors
cor(fatalities$Population,fatalities$`Vehicle miles traveled (billions)`, use="complete.obs") # won't go in together
cor(fatalities$Motorcycle,fatalities$Bicycle, use="complete.obs") # they're okay together but drinking and motorcycle can't
cor(fatalities$drinking,fatalities$`Vehicle miles traveled (billions)`, use="complete.obs")

out1=lm(d~fatalities$Year+fatalities$`Vehicle miles traveled (billions)`+ fatalities$Bicycle  + fatalities$Motorcycle)
#What I learned: Year is a super, super important predictor
summary(out1)

#Now let's do regression with ARMA errors.
#First let's find the residuals.

LM_Residuals=out1$residuals
plot(LM_Residuals,type='l')
acf2(res) # Looks like white noise.
acf2(res^2)
acf2(abs(res)) # These do too, so we can say iid

predict(out1)
predict(out1, newdata = subset(fatalities, Year = 2019), interval = "prediction")
plot(out1$fitted.values)
fityears<- fatalities %>% filter(Year > 1979)
plot(fityears$Year, out1$fitted.values, xlab='Year',ylab='Fitted Value', main = 'Fitted Values')

par(mfrow=c(2,2))
plot(out1)
par(mfrow=c(1,1))

xreg = cbind(fityears$Year,fityears$`Vehicle miles traveled (billions)`,fityears$Bicycle,fityears$Motorcycle)
out2=arima(fityears$Deaths,xreg=xreg,order=c(0,0,0))
summary(out2$aic) # can use sarima.for to make predictions if we have future values
# without future population/bicycle/motorcycle data OR lagged variables, we can't make predictions
# aic is much higher than other models

### Looking at shortened data (After 1945)
newfatalities <- fatalities %>% filter(Year > 1945)
plot(fatalities$`Change in per capita fatalities from previous year`) # settles around 1945
plot(fatalities$Year, log(fatalities$Deaths), type='o') # confirms that after 1945 makes sense. 
# also note residuals in main arima model settle after 1945. 

#Looking for periodic behavior.
newdeaths<-newfatalities$Deaths
x=newdeaths-mean(newfatalities$Deaths) #Centering.
plot(x, type='o')
n=length(x)
I=abs(fft(x))^2/n
P=(4/n)*I[1:(n/2)]
plot((0:(n/2-1))/n,P,type='o',xlab='Frequency',ylab='Scaled Periodogram')
which.max(P)
period=n/4
period # Period of 18.25 dominates - seems coincidental. 
# No periodic or seasonal behavior. 

### Modeling deaths since 1945 



#Looking for periodic behavior.
newdeaths<-newfatalities$Deaths
x=newdeaths-mean(newfatalities$Deaths) #Centering.
plot(x, type='o')
n=length(x)
I=abs(fft(x))^2/n
P=(4/n)*I[1:(n/2)]
plot((0:(n/2-1))/n,P,type='o',xlab='Frequency',ylab='Scaled Periodogram')
which.max(P)
period=n/4
period # Period of 18.25 dominates - seems coincidental. NO seasonal components.
    
# Data was looked at since 1920, 1930, 1945
newd<-newfatalities$Deaths
plot(newfatalities$Year,newfatalities$Deaths, xlab='Year',ylab='Fatalities', main='Motor Vehicle Fatalities (Since 1945)', type='o')
plot(log(newd)) # doesn't do a whole lot
plot(diff(newd)) # d=1 is enough
acf2(diff(newd))

# pq
uplim=4
aicmat=matrix(double((uplim+1)^2),uplim+1,uplim+1)
for (i in 0:uplim){
  for (j in 0:uplim){
    aicmat[i+1,j+1]<-sarima(newd,i,1,j,0,0,0,0,details=F)$AIC
    print(aicmat)}} #[1,2] means p=0, q=1 is recommended by AIC.

# Final sarima model
sarima(newd,0,1,1,0,0,0,0) # AIC = 18.01259 fits better than on whole data set

## Forcasting
newd<-ts(newd, start=1946, end=2018)
sarima.for(d,5,1,1,0) #Forecasts for 5 years.
sarima.for(newd,5,1,1,0)





