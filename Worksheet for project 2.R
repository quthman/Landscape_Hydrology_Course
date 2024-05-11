####Fitting poisson
##Get rainny days only
dt <- read.csv("F:/Dropbox (Personal)/Landscape Hydrology/PST.csv")
#utils::View(dt)
str(dt)
library(dplyr)
library(lubridate)
dtn <- dt %>% filter(STATION == "USC00082834" )
dtn$date <- as.Date(as.character(dtn$date), format = "%Y-%m-%d")
day <- day(dtn$date)
month <- month(dtn$date)
year <- year(dtn$date)
dtn <- cbind(dtn,month,day,year)
yearly_rainfall <- aggregate(alpha_p ~ year, dtn, FUN = sum)
daily_365_rainfall <- aggregate(alpha_p ~ day+month, dtn, FUN = mean)

Rain_SC <- dtn %>% filter(alpha_p>0)

library(seas)
library(MASS)
pois.SCppt<-fitdistr(Rain_SC$alpha_p, densfun="poisson")

inter_arr_time <- interarrival(dtn,var = "alpha_p")
pois.SCti<-fitdistr(na.exclude(inter_arr_time$dry), densfun="poisson")

Sim_SCrainfall <- data.frame(Date <- seq(min(as.Date(0, origin="2016-01-01")), max(as.Date(365, origin="2016-01-01")), by=1),
                             Rainfall <- rep(0,366))
colnames(Sim_SCrainfall) <- c("Date", "Rainfall")
i=1
while(Sim_SCrainfall$Date[i] < as.Date(365, origin="2016-01-01")){
  n <-  round(rexp(n=1, rate=1/pois.SCti$estimate))
  if(n+1+i < 365){
    Sim_SCrainfall$Rainfall[n+1+i] <- rpois(1, lambda=pois.SCppt$estimate)
    i= n+1+i}
  else {
    break
  }
}
Sim_SCrainfall <- Sim_SCrainfall[-366,]
Sim_SCrainfall <- cbind(Sim_SCrainfall, daily_365_rainfall$alpha_p)
plot(Sim_SCrainfall$Rainfall, type="l")
lines(Sim_Scirainfall$alpha_p, col="red")
sum(Sim_SCrainfall$Rainfall, na.rm = T)
sum(Sim_SCrainfall$`daily_365_rainfall$alpha_p`, na.rm = T)


Sim_SCrainfall$Rainfall[which(Sim_SCrainfall$Rainfall==0)] <- NA
Sim_SCrainfall$`daily_365_rainfall$alpha_p`[which(Sim_SCrainfall$`daily_365_rainfall$alpha_p`<=12)] <- NA

plot(ecdf(Sim_SCrainfall$Rainfall))
lines(ecdf(Sim_SCrainfall$`daily_365_rainfall$alpha_p`))

ecdfplot(~ Rainfall + `daily_365_rainfall$alpha_p`, 
         data=Sim_SCrainfall, auto.key=list(space='right'),
         xlab = "P (cm)", ylab = "Quantile")

ks.test(Sim_SCrainfall$Rainfall, Sim_SCrainfall$`daily_365_rainfall$alpha_p`, alternative = "two.sided", exact = T)


Sim_SCrainfall <- Sim_SCrainfall %>% mutate(Month= month(Date))


RainPoisson = function(ndays, lambda, alpha){
  
  rain = rep(NA, ndays)
  for(i in 1:ndays){
    n1 = n=rexp(1, rate=
    if(n1 < lambda){
      n2 = runif(1, min = 0, max = 1)
      rain[i] = (1/alpha)*log(1/(1-n2))
    }
    else{rain[i] = 0}
  }
  
  return(rain)
}
Sim_SC.monthly_sum <- aggregate(Rainfall ~  Month, Sim_SCrainfall, FUN = sum)
Sim_SC.monthly_mean <- aggregate(Rainfall ~ Month, Sim_SCrainfall, FUN = mean)
Sim_SC.monthly_se <- aggregate(Rainfall ~ Month, Sim_SCrainfall, FUN =  std.error)


ks.test(SC.monthly_mean$PRCP, Sim_SC.monthly_sum$Rainfall, alternative = "two.sided", exact = T)

test<-fitdistr(SC.monthly_mean$PRCP, densfun="poisson")
plot(ecdf(Sim_SC.monthly_sum$Rainfall))
lines(ecdf(SC.monthly_mean$PRCP), col="red")
lines(ecdf(rpois(5000, lambda=test$estimate)), col="blue")
