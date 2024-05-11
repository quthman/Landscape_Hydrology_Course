#### INPUT PARAMETERS ####
##
#### heterogeneity and reaction parameters ####
gamma<-0 #structured heterogeneity(equals 0.05)/ correlation of solute sources to travel time
#for homogenous solute equals 0 or 0.01
sigma.w <- 0 #random heterogeneity equals 0.05; homogeneous equals 0.01
damkohler.interarrival <- 1 #Da1
damkohler.transport <- 1 #Da2
damkohler.longterm <- 1 #Da3

# simulation time
tmax.yrs <- 10
tmax <- tmax.yrs*365

#### climate parameters ####
rain.per.year <- 133          # [cm/yr] mean rainfall amount 
interarrival.time.mean <- 6   # [d] mean rainfall interarrival time
lambda.time <- 1/interarrival.time.mean   # [1/days] mean of the rainfall arrival time exponential PDF
rainfall.mean <- rain.per.year/365/lambda.time # [cm/day] # it was divided by lamda.time because it does not rain all the days in the year.
lambda.intensity <- 1/rainfall.mean       # [1/cm] mean of the rainfall depth exponential PDF 
ET.max <- 0.9*rain.per.year/365                    # [cm/day] potential evapotranspiration 

#### soil and hydrology paramaters ####
z.r <- 40            # [cm], root zone thickness
z.vz <- 40           #[cm], thickness of the entire vadose zone including root zone # how can root zone thickness and vadose zone thickness be the same number.
theta.fc <- 0.3    #[-] water content at field capacity
theta.wp <- 0.15     #[-] water content at permanent wilting point
theta.res <- 0.15    # [-] residual water content
soil.storage <- z.r*(theta.fc-theta.wp)  # [cm] soil water storage root zone
soil.storage.vz <- z.vz*(theta.fc-theta.wp) #[cm] soil water storage entire vadose zone # no difference between vadose zone and root zone
R <- 1               # retardation factor in Harman et al. [2011]
mean.tr <- 5       # hydraulic response time [days]

theta.sat <- 0.4     # [-] saturated water content aquifer
aquifer.z <- 20      # [cm] aquifer thickness

#### HYDRAULIC RESPONSE - UNIT HYDROGRAPH ####
##
# dimensionless ratios
gamma.harman <- soil.storage/rainfall.mean   #ratio in Harman et al. [2011] for root zone
gamma.harman.vz <- soil.storage.vz/rainfall.mean #ratio in Harman et al. [2011] for entire vadose zone
bigF.harman <- (theta.fc-theta.res)/(theta.fc-theta.wp) # Harman et al. [2011]
phi <- ET.max/rainfall.mean/lambda.time         # Harman et al. [2011] dryness index

# [1/days] mean of the drainage arrival time exponential PDF root zone
lambda.d <- lambda.time/
  (gamma(gamma.harman/phi)-pgamma(gamma.harman,gamma.harman/phi, lower=FALSE)*gamma(gamma.harman/phi))*
  exp(-gamma.harman)*gamma.harman^(gamma.harman/phi)*phi/gamma.harman

# [1/days] mean of the drainage arrival time exponential PDF entire vadose zone
lambda.d.vz <- lambda.time/
  (gamma(gamma.harman.vz/phi)-pgamma(gamma.harman.vz,gamma.harman.vz/phi, lower=FALSE)*gamma(gamma.harman.vz/phi))*
  exp(-gamma.harman.vz)*gamma.harman.vz^(gamma.harman.vz/phi)*phi/gamma.harman.vz

# "effective gamma" from Harman et al. [2011] entire vadose zone
lambda.e.vz <- lambda.d.vz
gamma.harman.e.vz <- gamma.harman.vz*(lambda.d.vz/lambda.time)^.5

# travel time parameters vadose zone (=waiting time)
# these use "effective gamma" from Harman et al. [2011]
mean.waiting <- R*bigF.harman*gamma.harman.e.vz/lambda.e.vz  # equation (6), mean travel time vadose zone
sd.waiting <- (2*R*bigF.harman*gamma.harman.e.vz)^.5/lambda.e.vz # standard deviation travel time vadose zone

second.moment <- sd.waiting^2+mean.waiting^2
mu.ln.waiting <- 2*log(mean.waiting)-log(second.moment)/2
sigma.ln.waiting <- (log(second.moment)-2*log(mean.waiting))^.5

# linear reservoir recession constant, tr (k)
lambda.tr <- 1/mean.tr  # [1/day]

hydraulic.res.time <- rexp(3650, lambda.tr)# hydraulic response time distribution
CVh <- sd(hydraulic.res.time, na.rm = T)/mean(hydraulic.res.time, na.rm = T)

#### SOLUTE RESPONSE - TRAVEL TIME ####
##
# travel time as sum of vadose zone travel (waiting) time and groundwater travel time
# lambda.intensity does not change from vadose filtering
# total drainage does change = lambda.d/lambda.intensity 

drainage.mean <- lambda.d/lambda.intensity   # [cm/day] mean drainage from root zone
mean.g <- aquifer.z*theta.sat/drainage.mean  #  [days], equation (7)

# vadose and groundwater travel times: exponential
lambda.waiting <- 1/mean.waiting
lambda.g <- 1/mean.g

# m1 and m2 based on lambda time vadose zone and lambda time groundwater
m1.tp <- 1/lambda.g+1/lambda.waiting  #equation (9)
m2.tp <- 2/lambda.g^2+2/lambda.waiting^2+2/(lambda.g*lambda.waiting) #equation (10)

# hyperexponential fit of travel time with LN
mu.ln.time <- 2*log(m1.tp)-log(m2.tp)/2
sigma.ln.time <- (log(m2.tp)-2*log(m1.tp))^.5
mean.travel.time <- exp(mu.ln.time+sigma.ln.time^2/2) #mean total travel time
#sigma.travel.time <- sqrt(log((0.125^2+1),10)/log(exp(1), 10))# I assume 0.78 to be CV of log normal distribution



travel.time <- rlnorm(3650, meanlog = mu.ln.time, sdlog = sigma.ln.time)
CVt <- sd(travel.time, na.rm = T)/mean(travel.time, na.rm = T)

# times mean drainage rate, "mean log(streamtube length)": needed for redistribution of mobilized masses to following events
mu.ln.tube <- mu.ln.time+log(drainage.mean)
sigma.ln.tube <- sigma.ln.time 


#### SOURCE CONCENTRATION PARAMETERS ####
## 
# immobile concentration correlation to travel time - random heterogeneity
c.c <- 1                  # see equation (14), 1 for simplicity
mu.w<-sigma.w^2/2         # unit-mean log-normally distributed random variable, see equation (14)

# mean immobile concentration
mean.c.im.number <- 1 # [mg/L]
mean.c.im <- matrix(data=mean.c.im.number,nrow=tmax,ncol=1) # matrix to account for temporal variance

#Damkoehler 3: logistic long-term change of source concentration
max.c.im<-10*mean.c.im.number
longterm.change.time <- mean.travel.time/damkohler.longterm
logi.shape <- 2/(longterm.change.time)
mean.c.im<-matrix(data=(max.c.im/(1+exp(-logi.shape*(seq(1,tmax,1)-longterm.change.time)))),nrow=tmax,ncol=1)

# rate constants for reactions
# Damkoehler 1: drainage event inter-arrival, initial mobile concentration
rate.interarrival <- damkohler.interarrival*lambda.d #equation (19)

#  Damkoehler 2: mobile concentration (low value = little degradation)
rate.transport <- damkohler.transport/mean.travel.time #equation (21)

#### DIMENSION VARIABLES ####
## 
# create vectors
#set.seed(1)
#  rain events - effective (drainage)
daily.rain <- numeric(tmax)
interval <- numeric(1)
nonzero.rain <- numeric(1)
nonzero.times <- numeric(1)

# discharge and flux
discharge <- numeric(1)
mass.flux <- numeric(1)

# immobile and flux-avg concentrations
event.c.m <- numeric(1)

### EFFECTIVE RAINFALL EVENTS - TIMING AND AMOUNT
## 
N<-100*tmax*lambda.d #this version avoids a time-consuming loop [help from Gavan McGrath]
interval.l<-round(rexp(N,lambda.d),digits=0)
realtimes<-cumsum(interval.l) # random integers from exponential distribution
realtimes<-unique(realtimes) # get rid of two events on one day
realtimes<-realtimes[2:length(realtimes)] #get rid of first day as this can be zero
number.events<-length(which(realtimes<tmax)) #cut vector to correct length
interval<-interval.l[1:number.events]
nonzero.times<-realtimes[1:number.events]
nonzero.rain<-rexp(number.events,lambda.intensity) # lambda.intensity does not change from filtering (Botter 2007)

# vector of running sums of nonzero rain (effective)
rain.sum <- cumsum(nonzero.rain)

#### DISCHARGE ####
##
#definition of matries and vectors for discharge
event.discharge <- matrix(nrow = tmax, ncol = number.events)
event.discharge[is.na(event.discharge)] <- 0 
event.flux <- matrix(nrow = tmax, ncol = number.events)
event.flux[is.na(event.flux)] <- 0 

for (events.c in 1:number.events) {
  #discharge speeded up version
  startingtime=nonzero.times[events.c]
  if (startingtime==tmax) {
    #do nothing
  } else {
    event.discharge[(startingtime+1):(tmax),events.c] <- nonzero.rain[events.c]*dexp(1:((tmax)-startingtime), rate = lambda.tr)
  }
}

#### SOLUTE FLUX ####
##
# distribute solute among events
# initialize with zero values so can add to each other
event.flux[is.na(event.flux)] <- 0  

# k.count counter starting with each event
# k running count of events
for (events.c in 1:number.events) {
  time.ev <- nonzero.times[events.c]
  
  # event initial mobile concentration after immoble-mobile exchange
  event.c.m[events.c] <- mean.c.im[time.ev]*(1-exp(-rate.interarrival*interval[events.c])) #equation (18)
  # mean mobile concentration scaling parameter
  a.c <- event.c.m[events.c]*exp(-0.5*(gamma^2*sigma.ln.time^2+c.c^2*sigma.w^2)-gamma*mu.ln.time-c.c*mu.w) #equation (17)
  
  k.count <- 0
  
  # Track allocation of "events.c" to later events
  for (k in events.c:number.events) {
    k.count <- k.count + 1
    
    ## fraction of solute mass in this event
    event.timer <- nonzero.times[k]+1
    if (k.count==1) {
      event.fraction <- plnorm(rain.sum[k.count],mu.ln.tube,sigma.ln.tube) 
    }
    else {
      event.fraction <- plnorm(rain.sum[k.count],mu.ln.tube,sigma.ln.tube)-plnorm(rain.sum[(k.count-1)],mu.ln.tube,sigma.ln.tube)       
    }  
    
    # add 1 because zero Q on day of event
    if (nonzero.times[k]<tmax) {
      event.timer <- nonzero.times[k]+1
    } else event.timer <- nonzero.times[k]
    
    # From time of event k to tmax
    tau.real <- seq(event.timer,tmax)
    tau <- seq(1,(length(tau.real)))
    c.m.tau.vector <- a.c*tau^gamma*(rlnorm((length(tau)), mean = mu.w, sd = sigma.w))^c.c #equation (14)
    ### MOBILE CONCENTRATION
    # exponential decay along streamtube
    c.m.decay.tau.vector <- c.m.tau.vector*exp(-rate.transport*tau) #equation (20)
    # adding flux to the matrix
    event.flux[event.timer:tmax,events.c] <- event.flux[event.timer:tmax,events.c] + c.m.decay.tau.vector*event.discharge[event.timer:tmax,k]*event.fraction            
  }   
}

#### OUTPUT DISCHARGE, LOAD and CONCENTRATION [daily] ####
##  
# sum of rows 
discharge <- rowSums(event.discharge, na.rm = TRUE)
mass.flux <- rowSums(event.flux, na.rm = TRUE)
concentration=mass.flux/discharge

library(ggplot2)
library(lubridate)
Date <- seq(from = as.Date("2001-01-01"), 
            by = "day", length.out = 3650)
Days <- seq(1,3650, by = 1)

dts <- data.frame(cbind(discharge, mass.flux, concentration))

#mass flux - discharge time series
ggplot(dts, aes(x=Days))+
  geom_line(aes(y = discharge, colour = "Discharge"))+
  geom_line(aes(y = mass.flux/10, colour = "Mass flux"))+
  scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Mass flux [mg/day]"))+
  scale_colour_manual(values = c("blue", "red"))+
  labs(y = "Discharge [cm/day]",
       x = "Days",
       colour = "Parameter")+
  theme(legend.position = c(0.8, 0.9))+theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#solute concentration - discharge pattern
plot(dts$discharge, type="l")
lines(dts$concentration, col="red")

ggplot(dts, aes(x=Days))+
  geom_line(aes(y = discharge, colour = "Discharge"))+
  geom_line(aes(y = concentration/10, colour = "Solute concentration"))+
  scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Solute concentration [mg/L]"))+
  scale_colour_manual(values = c("blue", "red"))+
  labs(y = "Discharge [cm/day]",
       x = "Days",
       colour = "Parameter")+
  theme(legend.position = c(0.8, 0.9))+theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#Question 1b
plot(ecdf(travel.time))
lines(ecdf(hydraulic.res.time), col="red")
plot(ecdf(hydraulic.res.time))

library(lattice)
library(latticeExtra)
vals <- cbind(hydraulic.res.time, travel.time)
vals <- data.frame(vals)
ecdfplot(~hydraulic.res.time+travel.time, 
         data=vals, auto.key=list(space='right'),
         xlab = "Days", ylab = "Quantile")
ecdfplot(vals$hydraulic.res.time, xlab = "Days", ylab = "Quantile")
ecdfplot(vals$travel.time, xlab = "Days", ylab = "Quantile")
#Question 1a
#Mass flux discharge simulation time series 
ggplot(dts, aes(x=Days, y = mass.flux, colour="Mass discharge flux [mg/day]"))+
  geom_line()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  labs(x="Days", y="Mass discharge flux [mg/day]")

#Question 1c
#Solute concentration simulation time series
ggplot(dts, aes(x=Days, y = concentration, colour="Solute concentration [mg/L]"))+
  geom_line()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  labs(x="Days", y="Solute concentration [mg/L]")

#### PART 2 Project 3 #####
#Do fourier transformation of Discharge, Mass.flux and Concentration to get travel.time distribution
#Let me simulate recharge/drainage from the root zone using drainage.mean

drainage <- numeric(tmax)
event.count <- numeric(1)

event.time <- 1
nonzero <- 1

while(event.time <= tmax)
{interval <- round(rexp(1, lambda.d), digits = 0)
event.time <- event.time + interval
amount <- rexp(1, lambda.intensity)
drainage[event.time] <- amount
event.count[nonzero] <- amount
nonzero < nonzero + 1

}
drainage <- drainage[1:3650]

fft_mass.flux <- fft(mass.flux)
concentration <- ifelse(is.nan(concentration),0,concentration)
fft_concentration <- fft(concentration)
fft_drainage <- fft(drainage)
fft_travel.time <- fft_mass.flux/fft_drainage/fft_concentration
inverse.travel.time <- fft(fft_travel.time, inverse = T)/length(fft_travel.time)
plot(Re(inverse.travel.time), type = "l")
lines(travel.time, col="red")

plot(travel.time, type = "l")
lines(Re(inverse.travel.time), col="red")

ks.test(travel.time,Re(inverse.travel.time))
Summary(Re(fft_travel.time))

spectrum(fft_travel.time)
spectrum(fft_drainage)
spectrum(fft_concentration)
spectrum(travel.time)
spec.pgram(fft_travel.time, xlab = "Frequency [1/days]", ylab = "Spectrum power")
spec.pgram(travel.time, xlab = "Frequency [1/days]", ylab = "Spectrum power")# spec.pgram and spectrum does the same thing
spec.pgram(fft_drainage, xlab = "Frequency [1/days]", ylab = "Spectrum power")
spec.pgram(fft_concentration, xlab = "Frequency [1/days]", ylab = "Spectrum power")
spec.pgram(fft_mass.flux, xlab = "Frequency [1/days]", ylab = "Spectrum power")

CVt <- sd(travel.time, na.rm = T)/mean(travel.time, na.rm = T)
CVft <- sd(Re(fft_travel.time), na.rm = T)/mean(Re(fft_travel.time), na.rm = T)

### new idea

con <- convolve(drainage, travel.time)
