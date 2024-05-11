
##Site A Withlacochee river using daily rainfall between 1961-1970
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

library(MASS)
library(seas)
inter_arr_time <- interarrival(dtn,var = "alpha_p")
b <- ecdf(inter_arr_time$dry)
plot(b)
abc <- mean(inter_arr_time$dry, na.rm = T)
lamda_p <- 1/abc
alpha_p_mean <- (mean(dt$alpha_p, na.rm = T))/10 #to convert to cm
lamda_intensity <- 1/alpha_p_mean

#Dr. Jawitz
tmax.yrs <- 1
tmax <- tmax.yrs*365

daily.rain <- numeric(tmax)
event.count <- numeric(1)

event.time <- 1
nonzero <- 1

while(event.time <= tmax)
{interval <- round(rexp(1, lamda_p), digits = 0)
event.time <- event.time + interval
amount <- rexp(1, alpha_p_mean)
daily.rain[event.time] <- amount
event.count[nonzero] <- amount
nonzero < nonzero + 1
  
}
daily.rain <- daily.rain[-(366:379)]
sum(daily.rain, na.rm = T)



#library(Ecohydmod)
#rainsimulation <- RainPoisson(365,lambda = lamda_p, alpha = alpha_p_mean)
#sum(rainsimulation)


#theoretical expected rainfall
alpha_t <- rexp(365, rate = lamda_p)#using rexp does not change the graph

library(ggplot2)
daily_365_rainfall <- aggregate(alpha_p ~ day+month, dtn, FUN = mean)
dt_daily_rainfall <- cbind(daily_365_rainfall, daily.rain, alpha_t)
dt_daily_rainfall$Date <- as.Date( paste( daily_365_rainfall$month , daily_365_rainfall$day , sep = "-" )  , format = "%m-%d" )

plot(ecdf())
lines(ecdf())

ggplot(dt_daily_rainfall, aes(x=Date, y=alpha_p))+
  geom_bar(stat = "identity", position = "identity")+
  theme_bw()+
  scale_x_date(date_labels = ("%m/%d"))

ggplot(dt_daily_rainfall, aes(x=Date, y=rainsimulation))+
  geom_bar(stat = "identity", position = "identity")+
  theme_bw()+
  scale_x_date(date_labels = ("%m/%d"))

#Use excel to plot the cdf. Adding legend to ggplot is freaking me out.
#ggplot(dt_daily_rainfall)+stat_ecdf(aes(alpha_p),color="red")+
#  stat_ecdf(aes(rainsimulation),color="blue")+
#  stat_ecdf(aes(alpha_t),color="green")+
#  theme_bw()+
#  labs(y = "Quantile", x="Precipitation")


#cdf plot  
library(lattice)
library(latticeExtra)


vals <- data.frame(dt_daily_rainfall)
colnames(vals)[colnames(vals)=="alpha_p"] <- "Observed_P"
colnames(vals)[colnames(vals)=="alpha_t"] <- "Theoretical_P"

ecdfplot(~ Observed_P + rainsimulation + Theoretical_P, 
         data=vals, auto.key=list(space='right'),
         xlab = "Precipitation (cm)", ylab = "Quantile")



#komolorov-smirnov test
ks.test(dt_daily_rainfall$alpha_p,dt_daily_rainfall$rainsimulation, alternative = "less")
ks.test(dt_daily_rainfall$alpha_p,dt_daily_rainfall$alpha_t, alternative = "less")
ks.test(dt_daily_rainfall$alpha_t,dt_daily_rainfall$rainsimulation, alternative = "less")

# Can you get lamda_s? Yes I can. It is the same value with lamda_p
colnames(dt_daily_rainfall)[colnames(dt_daily_rainfall)=="Date"] <- "date"
#inter_arr_time_s <- interarrival(dt_daily_rainfall,var = "rainsimulation")
#d <- ecdf(inter_arr_time_s$dry)
#plot(d)
#abs <- mean(inter_arr_time$dry, na.rm = T)
#lamda_s <- 1/abs



#Stochastic recharge. calculate for recharge to get run off
alpha_s_mean <- mean(rainsimulation)
Zr <- 40 # units is in cm values are from class slides
theta_fc <- 0.3
theta_pwp <- 0.15
ETmax <- 0.9*alpha_s_mean#Look at Sanford
ETmax
phi <- ETmax/(alpha_s_mean*lamda_p)
gamma <- Zr*(theta_fc-theta_pwp)/alpha_s_mean
k <- exp(-gamma)*gamma^(gamma/phi)*(phi/gamma); 
library(zipfR)
l <- gamma(gamma/phi)-Igamma(gamma/phi, gamma, lower = TRUE)
m <- k/l
lamda_d <- lamda_p*m
lamda_d
Recharge <- RainPoisson(365,lambda = lamda_d, alpha = alpha_s_mean)
sum(Recharge)

plot(dt_daily_rainfall$date, Recharge, type = "l")

#Convolution integral for runoff
k <- 0.42 #from project 1, plot a declining Q against a time to have an idea of k
f_tau <- dexp(seq(1,365), rate = k)
#Recharge[is.na(Recharge)] <- 0
Q <- convolve(f_tau, Recharge)

plot(dt_daily_rainfall$date, Q, type = "l")

#theoretical expected runoff
beta_scale <- alpha_s_mean*k #beta
alpha_shape <- lamda_d/k #alpha
Q_t <- rgamma(365,shape = alpha_shape,scale = beta_scale)
plot(dt_daily_rainfall$date, Q_t, type = "l")

#visualization of Q, Recharge and alpha
dt_daily_rainfall_v <- cbind(dt_daily_rainfall, Q, Recharge)
g1 <- ggplot(dt_daily_rainfall_v, aes(date, rainsimulation)) +
  geom_bar(stat = 'identity', fill = "blue") +
  theme_bw() +
  ylab("Precip.") +
  labs(title = "Site A: Dobes Hole Lake-Withlacoochee River, HUC-10 Watershed ID 0310020802 675 km²") +
  scale_y_reverse()+
  theme(axis.title.x    = element_blank(),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank())

g2 <- ggplot(dt_daily_rainfall_v,aes(date, Recharge))+
  geom_bar(stat = 'identity', fill = "light blue") +
  ylab("Recharge")+
  theme_bw()+
  scale_y_reverse()+
  theme(axis.title.x    = element_blank(),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank())

g3 <- ggplot(dt_daily_rainfall_v)+
  geom_line(aes(date, Q)) +
    ylab("Streamflow")+
  scale_x_date(date_labels = "%b")+
  theme_bw()
  

library(grid)
library(gridExtra)
g1 <- ggplot_gtable(ggplot_build(g1))
g2 <- ggplot_gtable(ggplot_build(g2))
g3 <- ggplot_gtable(ggplot_build(g3))

maxWidth = unit.pmax(g1$widths[2:3], g2$widths[2:3], g3$widths[2:3])

g1$widths[2:3] <- maxWidth
g2$widths[2:3] <- maxWidth
g3$widths[2:3] <- maxWidth
grid.arrange(g1, g2, g3, ncol = 1, heights = c(2, 2, 3))

#Use excel to plot the cdf. Adding legend to ggplot is freaking me out.
runoff <- data.frame(cbind(Q, Q_t))
#ggplot(runoff)+stat_ecdf(aes(Q),color="red")+
#  stat_ecdf(aes(Q_t),color="dark red")+
#  theme_bw()+
#  labs(y = "Quantile", x="Runoff")

colnames(runoff)[colnames(runoff)=="Q_t"] <- "Theoretical_Q"
ecdfplot(~ Q + Theoretical_Q, 
         data=runoff, auto.key=list(space='right'),
         xlab = "Runoff (cm)", ylab = "Quantile")



#Comparison of simulated runoff to theoretical runoff
ks.test(Q,Q_t, alternative = "greater")



