dta <- read.csv("F:/Dropbox (Personal)/Landscape Hydrology/KS.csv")
str(dta)
ks.test(dta$A, dta$B)
