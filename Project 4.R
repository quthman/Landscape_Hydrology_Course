
#### Project 4 ####
#### Pervious layer ####
alpha_p_mean <- 1.33
lamda_p <- 0.1699387 #per day
lamda_intensity <- 0.7518797 #per cm
tmax.yrs <- 1
tmax <- tmax.yrs*365

daily.rain <- numeric(tmax)
event.count <- numeric(1)

event.time <- 1
nonzero <- 1
set.seed(10)
while(event.time <= tmax)
{interval <- round(rexp(1, lamda_p), digits = 0)
event.time <- event.time + interval
amount <- rexp(1, lamda_intensity)
daily.rain[event.time] <- amount
event.count[nonzero] <- amount
nonzero < nonzero + 1

}
sum(daily.rain[1:365], na.rm = T)
rainsimulation <- daily.rain[1:365]
#sum(rainsimulation, na.rm = T)

#alpha_t <- rexp(365, lamda_p)#using rexp does not change the graph
library(lubridate)
Date <- seq(from = as.Date("2020-01-01"), 
            by = "day", length.out = 365)
#dt_daily_rainfall <- data.frame(cbind(rainsimulation, alpha_t))
#dt_daily_rainfall$Date <- as.Date(dt_daily_rainfall$Date)

library(ggplot2)
#ggplot(dt_daily_rainfall, aes(x=Date, y=rainsimulation))+
 # geom_bar(stat = "identity", position = "identity")+
  #theme_bw()+
  #scale_x_date(date_labels = ("%m/%d"))

#rainsimulation[which(rainsimulation==0)] <- NA
#plot(ecdf(rainsimulation))
#lines(ecdf(alpha_t), col = "red")

#Stochastic recharge. calculate for recharge to get run off
alpha_s_mean <- mean(rainsimulation, na.rm = T)
Zr <- 40 # units is in cm values are from class slides
theta_fc <- 0.3
theta_pwp <- 0.15
ETmax <- 0.9*alpha_s_mean #Check Sanford
ETmax
phi <- ETmax/(alpha_p_mean*lamda_p)
gamma <- Zr*(theta_fc-theta_pwp)/alpha_p_mean
k <- exp(-gamma)*gamma^(gamma/phi)*(phi/gamma); 
library(zipfR)
l <- gamma(gamma/phi)-Igamma(gamma/phi, gamma, lower = TRUE)
m <- k/l; m
lamda_d <- lamda_p*m
lamda_d
lamda_d_intensity <- 1/alpha_p_mean
1/lamda_d;1/lamda_p

#Dr. Jawitz
tmax.yrs <- 1
tmax <- tmax.yrs*365

daily.recharge <- numeric(tmax)
event.count <- numeric(1)

event.time <- 1
nonzero <- 1
set.seed(10)
while(event.time <= tmax)
{interval <- round(rexp(1, lamda_d), digits = 0)
event.time <- event.time + interval
amount <- rexp(1, lamda_d_intensity)
daily.recharge[event.time] <- amount
event.count[nonzero] <- amount
nonzero < nonzero + 1

}

Recharge <- daily.recharge[1:365]
sum(Recharge, na.rm = T)
#Recharge <- Recharge*0.7 # 0.7 is 1-fraction of impervious surface
dp <- mean(Recharge)
plot(Date, Recharge, type = "l")

#Convolution integral for runoff
k. <- 25 #from project 1, plot a declining Q against a time to have an idea of k
k_ <- 1/k.
f_tau <- dexp(seq(1,365), rate = k_)
Recharge[is.na(Recharge)] <- 0
Qp <- convolve(f_tau, Recharge)
plot(Date, Qp)

#theoretical expected runoff
#beta_scale <- alpha_s_mean*k_ #beta
#alpha_shape <- lamda_d/k_ #alpha
#Q_t <- rgamma(366,shape = alpha_shape,scale = beta_scale)
#plot(Date, Q_t, type = "l")


#dt_daily_rainfall_v <- data.frame(cbind(rainsimulation, Q, Recharge))
#g1 <- ggplot(dt_daily_rainfall_v, aes(Date, rainsimulation)) +
#  geom_bar(stat = 'identity', fill = "blue") +
 # theme_bw() +
  #ylab("Precip.") +
  #labs(title = "Lagos, Nigeria 1,200 km²") +
  #scale_y_reverse()+
  #theme(axis.title.x    = element_blank(),
   #     axis.text.x     = element_blank(),
    #    axis.ticks.x    = element_blank())

#g2 <- ggplot(dt_daily_rainfall_v,aes(Date, Recharge))+
  #geom_bar(stat = 'identity', fill = "light blue") +
  #ylab("Recharge")+
  #theme_bw()+
  #scale_y_reverse()+
  #theme(axis.title.x    = element_blank(),
   #     axis.text.x     = element_blank(),
    #    axis.ticks.x    = element_blank())

#g3 <- ggplot(dt_daily_rainfall_v)+
  #geom_line(aes(Date, Q)) +
  #ylab("Streamflow")+
  #scale_x_date(date_labels = "%b")+
  #theme_bw()


library(grid)
library(gridExtra)
#g1 <- ggplot_gtable(ggplot_build(g1))
#g2 <- ggplot_gtable(ggplot_build(g2))
#g3 <- ggplot_gtable(ggplot_build(g3))

#maxWidth = unit.pmax(g1$widths[2:3], g2$widths[2:3], g3$widths[2:3])

#g1$widths[2:3] <- maxWidth
#g2$widths[2:3] <- maxWidth
#g3$widths[2:3] <- maxWidth
#grid.arrange(g1, g2, g3, ncol = 1, heights = c(2, 2, 3))

#plot(ecdf(Q))
#lines(ecdf(Q_t), col = "red")


#### Impervious layer ####
#Stochastic recharge. calculate for recharge to get run off
#alpha_s_mean <- mean(rainsimulation, na.rm = T)
#Zr1 <- 0.1 # units is in cm values are from class slides
#theta_fc1 <- 0.15
#theta_pwp1 <- 0.13
#ETmax1 <- 0.4*alpha_s_mean#Look at Sanford
#ETmax1
#phi1 <- ETmax1/(alpha_s_mean*lamda_p)
#gamma1 <- Zr1*(theta_fc1-theta_pwp1)/alpha_s_mean
#k1 <- exp(-gamma1)*gamma1^(gamma1/phi1)*(phi1/gamma1); 
#library(zipfR)
#l1 <- gamma(gamma1/phi1)-Igamma(gamma1/phi1, gamma1, lower = TRUE)
#m1 <- k1/l1
#lamda_d1 <- lamda_p*m1
#lamda_d1
#lamda_d_intensity <- 1/alpha_s_mean
  
#Dr. Jawitz
#tmax.yrs <- 1
#tmax <- tmax.yrs*365
  
#daily.recharge1 <- numeric(tmax)
#event.count <- numeric(1)
  
#event.time <- 1
#nonzero <- 1
#set.seed(10000)
#while(event.time <= tmax)
#{interval <- round(rexp(1, lamda_p), digits = 0)
#event.time <- event.time + interval
#amount <- rexp(1, lamda_intensity)
#daily.recharge1[event.time] <- amount
#event.count[nonzero] <- amount
#nonzero < nonzero + 1
  
#}
  
Recharge1 <- rainsimulation
sum(Recharge1, na.rm = T)
#Recharge1 <- Recharge1*0.3 # 0.3 is fraction of impervious surface
plot(Date, Recharge1, type = "l")
  
#Convolution integral for runoff
k_. <- 5 #from project 1, plot a declining Q against a time to have an idea of k
k_1 <- 1/k_.
f_tau1 <- dexp(seq(1,365), rate = k_1)
Recharge1[is.na(Recharge1)] <- 0
Qi <- convolve(f_tau1, Recharge1)


#### Multiply by fraction of pervious and impervious surface
ui <- seq(0.05, 1, by = 0.05)
up <- 1-ui
Qpf <- up*mean(Qp) 
Qif <- ui*mean(Qi)
Q <- Qpf + Qif
Q

Qfl <- c(0.012432734,0.100999516,0.090764837,0.009121503,0.132522726,0.011135437,
         0.147629435,0.001857744,7.71763E-05,0.042055569,0.119292504,0.008246839,
         0.055881151,0.021609363,0.055164514,0.180200037,0.039310298,0.000694587,
         0.009161929,0.08957963,0.092098884,0.151684498,0.196125182)

plot(Qfl)

#### Moments ####
m.term.1 <- lamda_d/lamda_d_intensity*up
m.term.2 <- lamda_p/lamda_intensity*ui
moment <- m.term.1 + m.term.2
moment; m.term.1; m.term.2 
####Variance####
v.term.1 <- lamda_d*k_*up^2/lamda_d_intensity^2
v.term.2 <- lamda_p*k_1*ui^2/lamda_d_intensity^2
vterm.3 <- 2*lamda_d*k_*k_1*ui*up/lamda_d_intensity^2
vterm.4 <- (dp*lamda_d_intensity +2)/(ui+up)
var. <- v.term.1 + v.term.2 + vterm.3*vterm.4
var.
v.term.1
v.term.2
####skewness####
term.1 <- 2*lamda_d*k_^2*up^3/lamda_d_intensity^3
term.2 <- 2*lamda_p*k_1^2*ui^3/lamda_d_intensity^3
term.3 <- 3*lamda_d*k_*k_1*up*ui/lamda_d_intensity^3
term.4 <- 2*k_*up*(3+dp*lamda_d_intensity)/(k_1+2*k_)
term.5 <- k_1*ui*(6+4*dp*lamda_d_intensity+dp^2*lamda_d_intensity^2)/(2*k_1+k_)
s_k <- term.1 + term.2 + term.3*(term.4 + term.5)
skewness <- s_k/var.
s.term.1 <- term.1/var.
s.term.2 <- term.2/var.
skewness; s.term.1; s.term.2
CVp <- sqrt(k_/lamda_d)
CVi <- sqrt(k_1/lamda_p)
CVp; CVi
plot(ecdf(Qif))
lines(ecdf(Qpf), col = "red")
new_data <- cbind(ui, Q, Qpf, Qif, moment, m.term.1, m.term.2, v.term.1, v.term.2, var., 
                  skewness, s.term.1, s.term.2 )
new_data <- data.frame(new_data)
plot(ui, log10(s.term.1), type = "l")

####Visualization of results#####
library(ggplot2)
ggplot(new_data, aes(x=ui))+
  geom_line(aes(y=Q, colour = "Total Flow"))+
  geom_line(aes(y=Qpf, colour = "Pervious Flow"))+
  geom_line(aes(y=Qif, colour = "Impervious Flow"))+
  scale_colour_manual(values = c("blue", "green", "purple"))+
  labs(y = "Discharge [cm/day]",
       x = "Fraction of impervious surface",
       colour = "")+
  theme(legend.position = c(0.8, 0.9))+theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggplot(new_data, aes(x=ui))+
  geom_line(aes(y=log10(var.), colour = "var {Q}/A^2"))+
  geom_line(aes(y=log10(Qpf), colour = "var {Qp}/A^2"))+
  geom_line(aes(y=log10(Qif), colour = "var {Qi}/A^2"))+
  scale_colour_manual(values = c("blue", "green", "purple"))+
  labs(y = "var{Q}/A^2 [cm/day]^2",
       x = "Fraction of impervious surface",
       colour = "")+
  theme(legend.position = c(0.8, 0.9))+theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggplot(new_data, aes(x=ui))+
  geom_line(aes(y=log10(skewness), colour = "S{Q}/A^3"))+
  geom_line(aes(y=log10(Qpf), colour = "S{Qp}/A^3"))+
  geom_line(aes(y=log10(Qif), colour = "S{Qi}/A^3"))+
  scale_colour_manual(values = c("blue", "green", "purple"))+
  labs(y = "S{Q}/A^3 [cm/day]^3",
       x = "Fraction of impervious surface",
       colour = "")+
  theme(legend.position = c(0.8, 0.9))+theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
library(lattice)
library(latticeExtra)

ecdfplot(~ Q + Qpf + Qif, 
         data=new_data, auto.key=list(space='right'),
         xlab = "Discharge", ylab = "Quantile")
plot(dgamma(Qif, shape = lamda_p/k_.), type = "l")
plot(dgamma(Qpf, shape = lamda_d/k_), type = "l")


Qo <- c(0.942763585,0.764608843,2.16175951,0.010584178,0.100999516,0.150136247,0.94185707,0.39295965,0.466198853,0.397665199,
        0.904506193,0.699441883,0.503634134,0.161439949,0.088818892,0.006642674,0.283255383,0.167812139,0.123553738,0.001253196,
        0.326128653,1.001096462,0.000176403,0.03005098,0.295091287,0.480682037,0.235404241,0.969797339,0.246523141,0.571104593,
        0.717249547,0.06988038,0.045081982,0.898788532,0.014523844,0.054881167,0.558559507,0.777263152,0.136370516,0.995684474,
        1.511841749,0.000694587,0,1.345991358,1.01373454,0.648104486,1.393528281,0.301263185,0,0,
        0.004586477,0.08957963,0.147680781,0.158929269,0.733836326,0.151684498,0.21970107,0.543088219,0.32593755,1.04096491,
        0.208091795)
Year <- seq(1958,2018, by = 1)
Year[1:2] <- c(1930, 1931)
plot(Year, Qo, type="l", ylab = "Observed Q")
library(hydroTSM)
fdc(Q, ylab="Q [cm/day]", main="Flow Duration Curve - Model")
fdc(Qo, ylab="Q [cm/day]", main="Flow Duration Curve - Observed Q")

cdf_Qif <- 1 - (1:length(Qif)/length(Qif))
Qif_s <- sort(Qif)
plot(cdf_Qif, log10(Qif_s))
cdf_Qp <-1- (1:length(Qp)/length(Qp))
Qp_s <- sort(Qp)
plot(cdf_Qp, log10(Qp_s))
cdf_Q <- 1 - (1:length(Q)/length(Q))
Q_s <- sort(Q)
plot(cdf_Q, log10(Q_s))
cdf_Qo <- 1 - (1:length(Qo)/length(Qo))
Qo_s <- sort(Qo)
plot(cdf_Qo, log10(Qo_s))
lines(cdf_Q, log10(Q_s), type = "l")


####Plot of Qp, Qi and P####

dat <- data.frame(cbind(rainsimulation, Qp, Qi, Recharge, Recharge1))
g1 <- ggplot(dat, aes(Date, rainsimulation)) +
  geom_bar(stat = 'identity', fill = "blue") +
  theme_bw() +
  ylab("Precip.") +
  labs(title = "Withlacoochee River Basin (Pervious Surface)") +
  scale_y_reverse()+
  theme(axis.title.x    = element_blank(),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

g2 <- ggplot(dat, aes(Date, Recharge))+
  geom_bar(stat = 'identity', fill = "light blue") +
  ylab("Recharge [cm/day]")+
  theme_bw()+
  scale_y_reverse()+
  theme(axis.title.x    = element_blank(),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

g3 <- ggplot(dat)+
  geom_line(aes(Date, Qp)) +
  ylab("Qp")+
  scale_x_date(date_labels = "%b")+labs(x="Month of the year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


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

g4 <- ggplot(dat, aes(Date, rainsimulation)) +
  geom_bar(stat = 'identity', fill = "blue") +
  theme_bw() +
  ylab("Precip.") +
  labs(title = "Withlacoochee River Basin (Impervious Surface)") +
  scale_y_reverse()+
  theme(axis.title.x    = element_blank(),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

g5 <- ggplot(dat, aes(Date, Recharge1))+
  geom_bar(stat = 'identity', fill = "light blue") +
  ylab("Recharge [cm/day]")+
  theme_bw()+
  scale_y_reverse()+
  theme(axis.title.x    = element_blank(),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

g6 <- ggplot(dat)+
  geom_line(aes(Date, Qi)) +
  ylab("Qi")+
  scale_x_date(date_labels = "%b")+labs(x="Month of the year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


library(grid)
library(gridExtra)
g4 <- ggplot_gtable(ggplot_build(g4))
g5 <- ggplot_gtable(ggplot_build(g5))
g6 <- ggplot_gtable(ggplot_build(g6))

maxWidth = unit.pmax(g4$widths[2:3], g5$widths[2:3], g6$widths[2:3])

g4$widths[2:3] <- maxWidth
g5$widths[2:3] <- maxWidth
g6$widths[2:3] <- maxWidth
grid.arrange(g4, g5, g6, ncol = 1, heights = c(2, 2, 3))

