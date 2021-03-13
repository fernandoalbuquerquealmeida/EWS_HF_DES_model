# Survival analyses for generating survival models used in simulation
# Text to the right of the "#" sign are comments and are not processed
# Comments are an essential part of model code and should be used to describe what the code is doing
# This promotes transparency and helps avoid and identify errors

# Remove previous contents
rm(list=ls())

# Set the working directory by replacing the text "C:/DES R model" with your directory where you have the model files (copy address as text)
# NOTE: Use / rather than \\ symbol
# NOTE: Put the working directory between ""
setwd("C:/Users/Fernando/Documents/PhD/TEN-HMS (Fernando)/Discrete event simulation model")

#### Options ####
# Generate plots? (T/F)
generate.plots <- T

# Save models? (T/F)
save.models <- T

#### Instalation of required packages ####
# Function for installing and loading R Packages
InstallAndLoad <- function(required.packages) {
  
    remaining.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])];
    
    if(length(remaining.packages)) 
    {
        install.packages(remaining.packages)
    }
    for(package_name in required.packages)
    {
        library(package_name, character.only = TRUE, quietly = TRUE)
    }
}

# Specify the list of required packages to be installed and loaded
required.packages <- c("survival", "flexsurv", "lme4", "ggplot2", "survminer", "gridExtra")

# Call the Function
InstallAndLoad(required.packages)

# create custom colours (for parametric survival fitting)
BM.Dblue <- rgb(0,45,92,max=255)
BM.red <- rgb(225, 55, 60, max=255)
BM.yellow <- rgb(238, 224, 30, max=255)
BM.pink <- rgb(211,78,147,max=255)
BM.blue <- rgb(69, 185, 209, max=255)
BM.orange <- rgb(255,128,0,max=255)
BM.green <- rgb(0,153,0,max=255)
BM.Dyellow <- rgb(214, 200, 16, max=255)

#### Data importing, cleaning and formatting ####
# Import data from .txt files
pat.char <- read.csv("./Data/baseline_population.csv")
surv.dat <- read.csv("./Data/survival_data.txt", sep = ';', stringsAsFactors = FALSE, na.strings = "")


# Convert dates to a date vectors readable in R
for(i in 2:ncol(surv.dat)) {
  surv.dat[,i] <- as.Date(surv.dat[,i], format = "%d-%m-%Y")
}

# Find first date of contact for each patient
surv.dat$first.contact.date <- as.Date(apply(surv.dat[2:3], 1, min, na.rm = TRUE))

# Find lastest date of contact for each patient
surv.dat$last.contact.death <- ifelse(!is.na(surv.dat$death.date), as.Date(surv.dat$death.date),
                                      as.Date(apply(surv.dat[3:11], 1, max, na.rm = TRUE)))

surv.dat$last.contact.hosp <- ifelse(!is.na(surv.dat$hosp.date), as.Date(surv.dat$hosp.date),
                                     as.Date(apply(surv.dat[3:11], 1, max, na.rm = TRUE)))

# Create a new column with time-to-event (death or censoring)
surv.dat$time.to.death <- (as.numeric(surv.dat$last.contact.death) - as.numeric(surv.dat$first.contact.date))

# Create a new column with time-to-event (hospitalisation or censoring)
surv.dat$time.to.hosp <- (as.numeric(surv.dat$last.contact.hosp) - as.numeric(surv.dat$first.contact.date))

# Create new column with binomial data for event/no event (death)
surv.dat$death <- ifelse(!is.na(surv.dat$death.date), 1, 0)

# Create new column with binomial data for event/no event (hospitalisation)
surv.dat$hosp <- ifelse(!is.na(surv.dat$hosp.date), 1, 0)

# Add survival object for death
surv.dat$surv.obj.death <- Surv(time = surv.dat$time.to.death, event = surv.dat$death)

# Add survival object for hospitalisation
surv.dat$surv.obj.hosp <- Surv(time = surv.dat$time.to.hosp, event = surv.dat$hosp)

# Merge survival data and baseline characteristics data
names(surv.dat)[names(surv.dat) == "patient.study.number"] <- "Patient.Study.Number"
surv.dat <- merge(pat.char, surv.dat, by = "Patient.Study.Number")

#### Save patients in a table ####
write.csv(surv.dat,"./Data/patient_list.csv")

# Remove patients with negative time-to-event: assumed error in data management
surv.dat <- subset(surv.dat, surv.dat$time.to.death > 0)

# Subset patients in UC and EWS groups
surv.dat <- subset(surv.dat, surv.dat$intervention %in% c(0, 1))

#### Death ####
#### Kaplan Meier curves for death ####
# Fit Kaplan-Meier curve
death.fit <- survfit(formula = surv.obj.death ~ factor(intervention), data = surv.dat)

# Show object
summary(death.fit)

#log rank test can be applied to test the difference between two survival curves
survdiff(surv.dat$surv.obj.death ~ factor(surv.dat$intervention), rho=1)
survdiff(surv.dat$surv.obj.death ~ factor(surv.dat$intervention), rho=0)

#plot Kaplan-Meier curves
if (generate.plots==T) {
tiff('./Graphs/Kaplan-Meier estimate (death).tiff',
     units="cm", width=23, height=15, res=600, compression = 'lzw')

km.death.plot <- ggsurvplot(fit = death.fit, censor = TRUE, risk.table = TRUE, risk.table.y.text.col = FALSE,
                 title = 'Kaplan-Meier estimate (death)',
                 legend.title = '', legend.labs =  c('UC', 'EWS'),
                 palette = c(1, 2)) #use this to control colour of lines

print(km.death.plot)

dev.off()
}

#### Regression formula for death ####
reg.formula.death <- surv.obj.death ~ factor(intervention) + scale(ejection.fraction) + scale(age) + scale(sbp) +
  scale(bmi) + scale(creatinine) + factor(nyha.class) + factor(gender) + factor(smoker) + factor(diabetes) +
  factor(copd) + factor(recent.diagnosis) + factor(no.beta.blocker) + factor(no.ace) + scale(age.ef) + scale(sbp.ef)

#### Cox regression for death ####
# Intervention only
surv.death.cox.int.only <- coxph(surv.obj.death ~ factor(intervention), data = surv.dat)

summary(surv.death.cox.int.only)

res.death.cox.int.only <- round(coef(summary(surv.death.cox.int.only))[,c(2,3,5)],3)
res.death.cox.int.only <- as.data.frame(res.death.cox.int.only)
write.csv(res.death.cox.int.only, file = "./Results/Survival analyses/res.death.cox.int.only.csv")

# Test proportional hazards assumption
cox.zph(surv.death.cox.int.only)

res.ph.test.death.int.only <- cox.zph(surv.death.cox.int.only)
res.ph.test.death.int.only <- as.data.frame(res.ph.test.death.int.only$table)
write.csv(res.ph.test.death.int.only, file = "./Results/Survival analyses/res.ph.test.death.int.only.csv")


# Test influential observations
if (generate.plots==T) {
  tiff('./Graphs/Influential observations test - intervention only (death).tiff',
       units = "cm", width = 23, height = 15, res = 800, compression = 'lzw')
  
  io.test.death.plot.in.only <- ggcoxdiagnostics(surv.death.cox.int.only, xlab = "The index number of observations",
                                type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_bw())
  
  print(io.test.death.plot.in.only)
  
  dev.off()
}

# With covariates
surv.death.cox <- coxph(reg.formula.death, data = surv.dat)

summary(surv.death.cox)

res.death.cox <- round(coef(summary(surv.death.cox))[,c(2,3,5)],3)
res.death.cox <- as.data.frame(res.death.cox)
write.csv(res.death.cox, file = "./Results/Survival analyses/res.death.cox.csv")

# Test proportional hazards assumption
cox.zph(surv.death.cox)

res.ph.test.death <- cox.zph(surv.death.cox)
res.ph.test.death <- as.data.frame(res.ph.test.death$table)
write.csv(res.ph.test.death, file = "./Results/Survival analyses/res.ph.test.death.csv")


# Test influential observations
if (generate.plots==T) {
tiff('./Graphs/Influential observations test (death).tiff',
     units = "cm", width = 23, height = 15, res = 800, compression = 'lzw')

io.test.death.plot <- ggcoxdiagnostics(surv.death.cox, xlab = "The index number of observations",
                      type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_bw())

print(io.test.death.plot)

dev.off()
}

# Plot log cumulative hazards
if (generate.plots==T) {
  tiff('./Graphs/log cumulative hazards, by treatment (death).tiff',
       units = "cm", width = 23, height = 15, res = 800, compression = 'lzw')
  
plot(death.fit, fun = "cloglog", main = "log cumulative hazards, by treatment (death)", xlab = "ln(time)", ylab = "log cumulative hazards", col = 1:2)
legend(x = "topleft", col = 1:2, lwd = 2, bty = "n", legend = c("UC", "EWS"))
  
  dev.off()
}


#### Parametric regression for death ####
# Intervention only
# Using an Exponential distribution
surv.death.exp.int.only <- flexsurvreg(surv.obj.death ~ factor(intervention), data = surv.dat, dist = "exp")

# Using a Weibull distribution
surv.death.weibull.int.only <- flexsurvreg(surv.obj.death ~ factor(intervention), data = surv.dat, dist = "weibull")

# Using a Log-normal distribution
surv.death.lnorm.int.only <- flexsurvreg(surv.obj.death ~ factor(intervention), data = surv.dat, dist = "lnorm")

# Using a Log-logistic distribution
surv.death.llogis.int.only <- flexsurvreg(surv.obj.death ~ factor(intervention), data = surv.dat, dist = "llogis")

# Using a Gompertz distribution
surv.death.gompertz.int.only <- flexsurvreg(surv.obj.death ~ factor(intervention), data = surv.dat, dist = "gompertz")

# Using a Generalised Gamma distribution
surv.death.gengamma.int.only <- flexsurvreg(surv.obj.death ~ factor(intervention), data = surv.dat, dist = "gengamma")


# All covariates
# Using an Exponential distribution
surv.death.exp <- flexsurvreg(reg.formula.death, data = surv.dat, dist = "exp")

# Using a Weibull distribution
surv.death.weibull <- flexsurvreg(reg.formula.death, data = surv.dat, dist = "weibull")

# Using a Log-normal distribution
surv.death.lnorm <- flexsurvreg(reg.formula.death, data = surv.dat, dist = "lnorm")

# Using a Log-logistic distribution
surv.death.llogis <- flexsurvreg(reg.formula.death, data = surv.dat, dist = "llogis")

# Using a Gompertz distribution
surv.death.gompertz <- flexsurvreg(reg.formula.death, data = surv.dat, dist = "gompertz")

# Using a Generalised Gamma distribution
surv.death.gengamma <- flexsurvreg(reg.formula.death, data = surv.dat, dist = "gengamma")


# Export results of parametric survival model (all covariates)
res.death.exp <- as.data.frame(round(surv.death.exp$res.t,3))
write.csv(res.death.exp, file = "./Results/Survival analyses/res.death.exp.csv")

res.death.weibull <- as.data.frame(round(surv.death.weibull$res.t,3))
write.csv(res.death.weibull, file = "./Results/Survival analyses/res.death.weibull.csv")

res.death.lnorm <- as.data.frame(round(surv.death.lnorm$res.t,3))
write.csv(res.death.lnorm, file = "./Results/Survival analyses/res.death.lnorm.csv")

res.death.llogis <- as.data.frame(round(surv.death.llogis$res.t,3))
write.csv(res.death.llogis, file = "./Results/Survival analyses/res.death.llogis.csv")

res.death.gompertz <- as.data.frame(round(surv.death.gompertz$res.t,3))
write.csv(res.death.gompertz, file = "./Results/Survival analyses/res.death.gompertz.csv")

res.death.gengamma <- as.data.frame(round(surv.death.gengamma$res.t,3))
write.csv(res.death.gengamma, file = "./Results/Survival analyses/res.death.gengamma.csv")


#### Plots for visual inspection for death ####
# Time horizon for plots (in days)
time.horizon <- 3652.5 # 10 years

# Define average patient
# UC
average.patient.uc <- as.numeric(c(0,
                                   mean(surv.dat$ejection.fraction.scale, na.rm = TRUE),
                                   mean(surv.dat$age.scale, na.rm = TRUE),
                                   mean(surv.dat$sbp.scale, na.rm = TRUE),
                                   mean(surv.dat$bmi.scale, na.rm = TRUE),
                                   mean(surv.dat$creatinine.scale, na.rm = TRUE),
                                   mean(surv.dat$nyha.class.2, na.rm = TRUE),
                                   mean(surv.dat$nyha.class.3, na.rm = TRUE),
                                   mean(surv.dat$nyha.class.4, na.rm = TRUE),
                                   mean(surv.dat$gender, na.rm = TRUE),
                                   mean(surv.dat$smoker, na.rm = TRUE),
                                   mean(surv.dat$diabetes, na.rm = TRUE),
                                   mean(surv.dat$copd, na.rm = TRUE),
                                   mean(surv.dat$recent.diagnosis, na.rm = TRUE),
                                   mean(surv.dat$no.beta.blocker, na.rm = TRUE),
                                   mean(surv.dat$no.ace, na.rm = TRUE),
                                   mean(surv.dat$age.ef.scale, na.rm = TRUE),
                                   mean(surv.dat$sbp.ef.scale, na.rm = TRUE)))

# EWS
average.patient.ews <- as.numeric(c(1,
                                    mean(surv.dat$ejection.fraction.scale, na.rm = TRUE),
                                    mean(surv.dat$age.scale, na.rm = TRUE),
                                    mean(surv.dat$sbp.scale, na.rm = TRUE),
                                    mean(surv.dat$bmi.scale, na.rm = TRUE),
                                    mean(surv.dat$creatinine.scale, na.rm = TRUE),
                                    mean(surv.dat$nyha.class.2, na.rm = TRUE),
                                    mean(surv.dat$nyha.class.3, na.rm = TRUE),
                                    mean(surv.dat$nyha.class.4, na.rm = TRUE),
                                    mean(surv.dat$gender, na.rm = TRUE),
                                    mean(surv.dat$smoker, na.rm = TRUE),
                                    mean(surv.dat$diabetes, na.rm = TRUE),
                                    mean(surv.dat$copd, na.rm = TRUE),
                                    mean(surv.dat$recent.diagnosis, na.rm = TRUE),
                                    mean(surv.dat$no.beta.blocker, na.rm = TRUE),
                                    mean(surv.dat$no.ace, na.rm = TRUE),
                                    mean(surv.dat$age.ef.scale, na.rm = TRUE),
                                    mean(surv.dat$sbp.ef.scale, na.rm = TRUE)))

# Calculate parameters for the average patient in UC
# Exponential UC
reg.coef.death.exp.uc  <- as.numeric(coef(surv.death.exp))
rate.death.exp.uc <- exp(reg.coef.death.exp.uc[1]+sum(tail(reg.coef.death.exp.uc,n=-1)*average.patient.uc))

# Weibull UC
reg.coef.death.weibull.uc  <- as.numeric(coef(surv.death.weibull))
shape.death.weibull.uc <- exp(reg.coef.death.weibull.uc[1])
scale.death.weibull.uc <- exp(sum(tail(reg.coef.death.weibull.uc,n=-1)*c(1,average.patient.uc)))

# Log-normal UC
reg.coef.death.lnorm.uc  <- as.numeric(coef(surv.death.lnorm))
meanlog.death.lnorm.uc <- reg.coef.death.lnorm.uc[1]+sum(tail(reg.coef.death.lnorm.uc,n=-2)*average.patient.uc)
sdlog.death.lnorm.uc <- exp(reg.coef.death.lnorm.uc[2])

# Log-logistic UC
reg.coef.death.llogis.uc  <- as.numeric(coef(surv.death.llogis))
shape.death.llogis.uc <- exp(reg.coef.death.llogis.uc[1])
scale.death.llogis.uc <- exp(reg.coef.death.llogis.uc[2]+sum(tail(reg.coef.death.llogis.uc,n=-2)*average.patient.uc))

# Gompertz UC
reg.coef.death.gompertz.uc  <- as.numeric(coef(surv.death.gompertz))
shape.death.gompertz.uc <- reg.coef.death.gompertz.uc[1]
rate.death.gompertz.uc <- exp(reg.coef.death.gompertz.uc[2]+sum(tail(reg.coef.death.gompertz.uc,n=-2)*average.patient.uc))  

# Generalised gamma UC
reg.coef.death.gengamma.uc  <- as.numeric(coef(surv.death.gengamma))
mu.death.gengamma.uc <- reg.coef.death.gengamma.uc[1]+sum(tail(reg.coef.death.gengamma.uc,n=-3)*average.patient.uc)
sigma.death.gengamma.uc <- exp(reg.coef.death.gengamma.uc[2])
Q.death.gengamma.uc <- reg.coef.death.gengamma.uc[3]


# Calculate parameters for the average patient in EWS
# Exponential EWS
reg.coef.death.exp.ews  <- as.numeric(coef(surv.death.exp))
rate.death.exp.ews <- exp(reg.coef.death.exp.ews[1]+sum(tail(reg.coef.death.exp.ews,n=-1)*average.patient.ews))

# Weibull EWS
reg.coef.death.weibull.ews  <- as.numeric(coef(surv.death.weibull))
shape.death.weibull.ews <- exp(reg.coef.death.weibull.ews[1])
scale.death.weibull.ews <- exp(sum(tail(reg.coef.death.weibull.ews,n=-1)*c(1,average.patient.ews)))

# Log-normal EWS
reg.coef.death.lnorm.ews  <- as.numeric(coef(surv.death.lnorm))
meanlog.death.lnorm.ews <- reg.coef.death.lnorm.ews[1]+sum(tail(reg.coef.death.lnorm.ews,n=-2)*average.patient.ews)
sdlog.death.lnorm.ews <- exp(reg.coef.death.lnorm.ews[2])

# Log-logistic EWS
reg.coef.death.llogis.ews  <- as.numeric(coef(surv.death.llogis))
shape.death.llogis.ews <- exp(reg.coef.death.llogis.ews[1])
scale.death.llogis.ews <- exp(reg.coef.death.llogis.ews[2]+sum(tail(reg.coef.death.llogis.ews,n=-2)*average.patient.ews))

# Gompertz EWS
reg.coef.death.gompertz.ews  <- as.numeric(coef(surv.death.gompertz))
shape.death.gompertz.ews <- reg.coef.death.gompertz.ews[1]
rate.death.gompertz.ews <- exp(reg.coef.death.gompertz.ews[2]+sum(tail(reg.coef.death.gompertz.ews,n=-2)*average.patient.ews))  

# Generalised gamma EWS
reg.coef.death.gengamma.ews  <- as.numeric(coef(surv.death.gengamma))
mu.death.gengamma.ews <- reg.coef.death.gengamma.ews[1]+sum(tail(reg.coef.death.gengamma.ews,n=-3)*average.patient.ews)
sigma.death.gengamma.ews <- exp(reg.coef.death.gengamma.ews[2])
Q.death.gengamma.ews <- reg.coef.death.gengamma.ews[3]


# Plot KM curves and parametric survival model fitting curves (death) - intervention as the only variable
if (generate.plots==T) {
tiff('./Graphs/Parametric survival model fitting (death).tiff',
     units="cm", width=23, height=15, res=1200, compression = 'lzw')

plot(death.fit, main = "Parametric survival model fitting (death)", xlab = "time (days)", ylab = "survival probability", xlim = c(0, time.horizon), conf.int	= F, lwd = 2, lty = c(1,2))
lines(surv.death.exp.int.only, ci =F, t = c(0:time.horizon), lwd = 1, lty = 3, col = BM.orange)
lines(surv.death.weibull.int.only, ci =F, t = c(0:time.horizon), lwd = 1, lty = 4, col = BM.Dyellow)
lines(surv.death.lnorm.int.only, ci =F,t = c(0:time.horizon), lwd = 1, lty = 5, col = BM.pink)
lines(surv.death.llogis.int.only, ci =F, t = c(0:time.horizon), lwd = 1, lty = 6, col = BM.blue)
lines(surv.death.gompertz.int.only, ci =F, t = c(0:time.horizon), lwd = 1, lty = 7, col = BM.Dblue)
lines(surv.death.gengamma.int.only, ci =F, t = c(0:time.horizon), lwd = 1, lty = 8, col = BM.yellow)
legend(x = "bottomleft", col = c(1, 1, BM.orange, BM.Dyellow, BM.pink, BM.blue, BM.Dblue, BM.yellow),
       lwd = 2, lty = 1:8, bty = "n", cex = 0.8,
       legend = c("UC", "EWS", "exponential", "weibull", "log-normal", "log-logistic", "gompertz", "gengamma"))

dev.off()
}


# Plot KM curves and parametric survival model fitting curves (death) - all covariates
if (generate.plots==T) {
  tiff('./Graphs/Parametric survival model fitting, all covariates (death).tiff',
       units="cm", width=23, height=15, res=1200, compression = 'lzw')
  
  plot(death.fit, main = "Parametric survival model fitting, all covariates (death)", xlab = "time (days)", ylab = "survival probability", xlim = c(0, time.horizon), conf.int	= F, lwd = 2, lty = c(1,2))
  curve((1-pexp(x, rate = rate.death.exp.uc)), xlim = c(0,time.horizon), lwd = 1, lty = 3, col = BM.orange, add = TRUE)
  curve((1-pweibull(x, shape = shape.death.weibull.uc, scale = scale.death.weibull.uc)), xlim = c(0,time.horizon), lwd = 1, lty = 4, col = BM.Dyellow, add = TRUE)
  curve((1-plnorm(x, meanlog = meanlog.death.lnorm.uc, sdlog = sdlog.death.lnorm.uc)), xlim = c(0,time.horizon), lwd = 1, lty = 5, col = BM.pink, add = TRUE)
  curve((1-pllogis(x, shape = shape.death.llogis.uc, scale = scale.death.llogis.uc)), xlim = c(0,time.horizon), lwd = 1, lty = 6, col = BM.blue, add = TRUE)
  curve((1-pgompertz(x, shape = shape.death.gompertz.uc, rate = rate.death.gompertz.uc)), xlim = c(0,time.horizon), lwd = 1, lty = 7, col = BM.Dblue, add = TRUE)
  curve((1-pgengamma(x, mu = mu.death.gengamma.uc, sigma = sigma.death.gengamma.uc, Q = Q.death.gengamma.uc)), xlim = c(0,time.horizon), lwd = 1, lty = 8, col = BM.yellow, add = TRUE)
  curve((1-pexp(x, rate = rate.death.exp.ews)), xlim = c(0,time.horizon), lwd = 1, lty = 3, col = BM.orange, add = TRUE)
  curve((1-pweibull(x, shape = shape.death.weibull.ews, scale = scale.death.weibull.ews)), xlim = c(0,time.horizon), lwd = 1, lty = 4, col = BM.Dyellow, add = TRUE)
  curve((1-plnorm(x, meanlog = meanlog.death.lnorm.ews, sdlog = sdlog.death.lnorm.ews)), xlim = c(0,time.horizon), lwd = 1, lty = 5, col = BM.pink, add = TRUE)
  curve((1-pllogis(x, shape = shape.death.llogis.ews, scale = scale.death.llogis.ews)), xlim = c(0,time.horizon), lwd = 1, lty = 6, col = BM.blue, add = TRUE)
  curve((1-pgompertz(x, shape = shape.death.gompertz.ews, rate = rate.death.gompertz.ews)), xlim = c(0,time.horizon), lwd = 1, lty = 7, col = BM.Dblue, add = TRUE)
  curve((1-pgengamma(x, mu = mu.death.gengamma.ews, sigma = sigma.death.gengamma.ews, Q = Q.death.gengamma.ews)), xlim = c(0,time.horizon), lwd = 1, lty = 8, col = BM.yellow, add = TRUE)
  legend(x = "bottomleft", col = c(1, 1, BM.orange, BM.Dyellow, BM.pink, BM.blue, BM.Dblue, BM.yellow),
         lwd = 2, lty = 1:8, bty = "n", cex = 0.8,
         legend = c("UC", "EWS", "exponential", "weibull", "log-normal", "log-logistic", "gompertz", "gengamma"))
  
  dev.off()
}


# Compare goodness-of-fit
# Using AIC
surv.death.exp.int.only$AIC
surv.death.weibull.int.only$AIC
surv.death.lnorm.int.only$AIC
surv.death.llogis.int.only$AIC
surv.death.gompertz.int.only$AIC
surv.death.gengamma.int.only$AIC

# Using BIC
BIC(surv.death.exp.int.only)
BIC(surv.death.weibull.int.only)
BIC(surv.death.lnorm.int.only)
BIC(surv.death.llogis.int.only)
BIC(surv.death.gompertz.int.only)
BIC(surv.death.gengamma.int.only)

# Create and export goodness-of-fit table
parametric.distributions <- c("Exponential", "Weibull", "log-normal", "log-logistic", "Gompertz", "Gengamma")
aic.surv.death <- c(surv.death.exp.int.only$AIC, surv.death.weibull.int.only$AIC, surv.death.lnorm.int.only$AIC, surv.death.llogis.int.only$AIC, surv.death.gompertz.int.only$AIC, surv.death.gengamma.int.only$AIC)
bic.surv.death <- c(BIC(surv.death.exp.int.only), BIC(surv.death.weibull.int.only), BIC(surv.death.lnorm.int.only), BIC(surv.death.llogis.int.only), BIC(surv.death.gompertz.int.only), BIC(surv.death.gengamma.int.only))

goodness.fit.death <- data.frame(parametric.distributions, aic.surv.death, bic.surv.death)
names(goodness.fit.death) <- c("Distribution", "AIC", "BIC")
write.csv(goodness.fit.death, file = "./Results/Survival analyses/goodness.fit.death.csv")


# Extract 3-year survival for each of the different survival models
# Intervention only
summary(surv.death.exp.int.only, type = "survival", t = 1095.75)
summary(surv.death.weibull.int.only, type = "survival", t = 1095.75)
summary(surv.death.lnorm.int.only, type = "survival", t = 1095.75)
summary(surv.death.llogis.int.only, type = "survival", t = 1095.75)
summary(surv.death.gompertz.int.only, type = "survival", t = 1095.75)
summary(surv.death.gengamma.int.only, type = "survival", t = 1095.75)

# All covariates
# For UC
pexp(1095.75, rate = rate.death.exp.uc)
pweibull(1095.75, shape = shape.death.weibull.uc, scale = scale.death.weibull.uc)
plnorm(1095.75, meanlog = meanlog.death.lnorm.uc, sdlog = sdlog.death.lnorm.uc)
pllogis(1095.75, shape = shape.death.llogis.uc, scale = scale.death.llogis.uc)
pgompertz(1095.75, shape = shape.death.gompertz.uc, rate = rate.death.gompertz.uc)
pgengamma(1095.75, mu = mu.death.gengamma.uc, sigma = sigma.death.gengamma.uc, Q = Q.death.gengamma.uc)

# For EWS
pexp(1095.75, rate = rate.death.exp.ews)
pweibull(1095.75, shape = shape.death.weibull.ews, scale = scale.death.weibull.ews)
plnorm(1095.75, meanlog = meanlog.death.lnorm.ews, sdlog = sdlog.death.lnorm.ews)
pllogis(1095.75, shape = shape.death.llogis.ews, scale = scale.death.llogis.ews)
pgompertz(1095.75, shape = shape.death.gompertz.ews, rate = rate.death.gompertz.ews)
pgengamma(1095.75, mu = mu.death.gengamma.ews, sigma = sigma.death.gengamma.ews, Q = Q.death.gengamma.ews)


#### Save model for death ####
if (save.models==T) {
save(surv.death.exp, surv.death.weibull, surv.death.lnorm, surv.death.llogis, surv.death.gompertz, surv.death.gengamma,
     file = "./Regression models/surv_death_model.rda")
}

#### Hospitalisation ####
#### Kaplan Meier curves for hospitalisation ####
# Fit Kaplan-Meier curve
hosp.fit <- survfit(formula = surv.obj.hosp ~ factor(intervention), data = surv.dat)

# Show object
summary(hosp.fit)

#log rank test can be applied to test the difference between two survival curves
survdiff(surv.dat$surv.obj.hosp ~ factor(surv.dat$intervention), rho=1)
survdiff(surv.dat$surv.obj.hosp ~ factor(surv.dat$intervention), rho=0)

#plot Kaplan-Meier curves
if (generate.plots==T) {
tiff('./Graphs/Kaplan-Meier estimate (hospitalisation).tiff',
     units="cm", width=23, height=15, res=1200, compression = 'lzw')

km.hosp.plot <- ggsurvplot(fit = hosp.fit, censor = TRUE, risk.table = TRUE, risk.table.y.text.col = FALSE,
                title = 'Kaplan-Meier estimate (hospitalisation)',
                legend.title = '', legend.labs =  c('UC', 'EWS'),
                palette = c(1, 2)) #use this to control colour of lines

print(km.hosp.plot)

dev.off()
}

#### Regression formula for hospitalisation ####
reg.formula.hosp <- surv.obj.hosp ~ factor(intervention) + factor(intervention) + scale(ejection.fraction) + scale(age) + scale(sbp) +
  scale(bmi) + scale(creatinine) + factor(nyha.class) + factor(gender) + factor(smoker) + factor(diabetes) +
  factor(copd) + factor(recent.diagnosis) + factor(no.beta.blocker) + factor(no.ace) + scale(age.ef) + scale(sbp.ef)

#### Cox regression for hospitalisation ####
# Intervention only
surv.hosp.cox.int.only <- coxph(surv.obj.hosp ~ factor(intervention), data = surv.dat)

summary(surv.hosp.cox.int.only)

res.hosp.cox.int.only <- round(coef(summary(surv.hosp.cox.int.only))[,c(2,3,5)],3)
res.hosp.cox.int.only <- as.data.frame(res.hosp.cox.int.only)
write.csv(res.hosp.cox.int.only, file = "./Results/Survival analyses/res.hosp.cox.int.only.csv")

# Test proportional hazards assumption
cox.zph(surv.hosp.cox.int.only)

res.ph.test.hosp.int.only <- cox.zph(surv.hosp.cox.int.only)
res.ph.test.hosp.int.only <- as.data.frame(res.ph.test.hosp.int.only$table)
write.csv(res.ph.test.hosp.int.only, file = "./Results/Survival analyses/res.ph.test.hosp.int.only.csv")


# Test influential observations
if (generate.plots==T) {
  tiff('./Graphs/Influential observations test - intervention only (hospitalisation).tiff',
       units = "cm", width = 23, height = 15, res = 800, compression = 'lzw')
  
  io.test.hosp.plot.in.only <- ggcoxdiagnostics(surv.hosp.cox.int.only, xlab = "The index number of observations",
                               type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_bw())
  
  print(io.test.hosp.plot.in.only)
  
  dev.off()
}

# With covariates
surv.hosp.cox <- coxph(reg.formula.hosp, data = surv.dat)

summary(surv.hosp.cox)

res.hosp.cox <- round(coef(summary(surv.hosp.cox))[,c(2,3,5)],3)
res.hosp.cox <- as.data.frame(res.hosp.cox)
write.csv(res.hosp.cox, file = "./Results/Survival analyses/res.hosp.cox.csv")

# Test proportional hazards assumption
cox.zph(surv.hosp.cox)

res.ph.test.hosp <- cox.zph(surv.hosp.cox)
res.ph.test.hosp <- as.data.frame(res.ph.test.hosp$table)
write.csv(res.ph.test.hosp, file = "./Results/Survival analyses/res.ph.test.hosp.csv")

# Test influential observations
if (generate.plots==T) {
tiff('./Graphs/Influential observations test (hospitalisation).tiff',
     units = "cm", width = 23, height = 15, res = 800, compression = 'lzw')

io.test.hosp.plot <- ggcoxdiagnostics(surv.hosp.cox, xlab = "The index number of observations",
                      type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_bw())

print(io.test.hosp.plot)

dev.off()
}

# Plot log cumulative hazards
if (generate.plots==T) {
  tiff('./Graphs/log cumulative hazards, by treatment (hospitalisation).tiff',
       units = "cm", width = 23, height = 15, res = 800, compression = 'lzw')
  
  plot(hosp.fit, fun = "cloglog", main = "log cumulative hazards, by treatment (hospitalisation)", xlab = "ln(time)", ylab = "log cumulative hazards", col = 1:2)
  legend(x = "topleft", col = 1:2, lwd = 2, bty = "n", legend = c("UC", "EWS"))
  
  dev.off()
}

#### Parametric regression for hospitalisation ####
# Intervention only
# Using an Exponential distribution
surv.hosp.exp.int.only <- flexsurvreg(surv.obj.hosp ~ factor(intervention), data = surv.dat, dist = "exp")

# Using a Weibull distribution
surv.hosp.weibull.int.only <- flexsurvreg(surv.obj.hosp ~ factor(intervention), data = surv.dat, dist = "weibull")

# Using a Log-normal distribution
surv.hosp.lnorm.int.only <- flexsurvreg(surv.obj.hosp ~ factor(intervention), data = surv.dat, dist = "lnorm")

# Using a Log-logistic distribution
surv.hosp.llogis.int.only <- flexsurvreg(surv.obj.hosp ~ factor(intervention), data = surv.dat, dist = "llogis")

# Using a Gompertz distribution
surv.hosp.gompertz.int.only <- flexsurvreg(surv.obj.hosp ~ factor(intervention), data = surv.dat, dist = "gompertz")

# Using a Generalised Gamma distribution
surv.hosp.gengamma.int.only <- flexsurvreg(surv.obj.hosp ~ factor(intervention), data = surv.dat, dist = "gengamma")


# Using an Exponential distribution
surv.hosp.exp <- flexsurvreg(reg.formula.hosp, data = surv.dat, dist = "exp")

# Using a Weibull distribution
surv.hosp.weibull <- flexsurvreg(reg.formula.hosp, data = surv.dat, dist = "weibull")

# Using a Log-normal distribution
surv.hosp.lnorm <- flexsurvreg(reg.formula.hosp, data = surv.dat, dist = "lnorm")

# Using a Log-logistic distribution
surv.hosp.llogis <- flexsurvreg(reg.formula.hosp, data = surv.dat, dist = "llogis")

# Using a Gompertz distribution
surv.hosp.gompertz <- flexsurvreg(reg.formula.hosp, data = surv.dat, dist = "gompertz")

# Using a generalised gamma distribution
surv.hosp.gengamma <- flexsurvreg(reg.formula.hosp, data = surv.dat, dist = "gengamma")

# Export results of parametric survival models
res.hosp.exp <- as.data.frame(round(surv.hosp.exp$res.t,3))
write.csv(res.hosp.exp, file = "./Results/Survival analyses/res.hosp.exp.csv")

res.hosp.weibull <- as.data.frame(round(surv.hosp.weibull$res.t,3))
write.csv(res.hosp.weibull, file = "./Results/Survival analyses/res.hosp.weibull.csv")

res.hosp.lnorm <- as.data.frame(round(surv.hosp.lnorm$res.t,3))
write.csv(res.hosp.lnorm, file = "./Results/Survival analyses/res.hosp.lnorm.csv")

res.hosp.llogis <- as.data.frame(round(surv.hosp.llogis$res.t,3))
write.csv(res.hosp.llogis, file = "./Results/Survival analyses/res.hosp.llogis.csv")

res.hosp.gompertz <- as.data.frame(round(surv.hosp.gompertz$res.t,3))
write.csv(res.hosp.gompertz, file = "./Results/Survival analyses/res.hosp.gompertz.csv")

res.hosp.gengamma <- as.data.frame(round(surv.hosp.gengamma$res.t,3))
write.csv(res.hosp.gengamma, file = "./Results/Survival analyses/res.hosp.gengamma.csv")

#### Plots for visual inspection for hospitalisation ####
# Time horizon for plots (in days)
time.horizon <- 3652.5 # 10 years

# Define average patient
# UC
average.patient.uc <- as.numeric(c(0,
                                   mean(surv.dat$ejection.fraction.scale, na.rm = TRUE),
                                   mean(surv.dat$age.scale, na.rm = TRUE),
                                   mean(surv.dat$sbp.scale, na.rm = TRUE),
                                   mean(surv.dat$bmi.scale, na.rm = TRUE),
                                   mean(surv.dat$creatinine.scale, na.rm = TRUE),
                                   mean(surv.dat$nyha.class.2, na.rm = TRUE),
                                   mean(surv.dat$nyha.class.3, na.rm = TRUE),
                                   mean(surv.dat$nyha.class.4, na.rm = TRUE),
                                   mean(surv.dat$gender, na.rm = TRUE),
                                   mean(surv.dat$smoker, na.rm = TRUE),
                                   mean(surv.dat$diabetes, na.rm = TRUE),
                                   mean(surv.dat$copd, na.rm = TRUE),
                                   mean(surv.dat$recent.diagnosis, na.rm = TRUE),
                                   mean(surv.dat$no.beta.blocker, na.rm = TRUE),
                                   mean(surv.dat$no.ace, na.rm = TRUE),
                                   mean(surv.dat$age.ef.scale, na.rm = TRUE),
                                   mean(surv.dat$sbp.ef.scale, na.rm = TRUE)))

# EWS
average.patient.ews <- as.numeric(c(1,
                                    mean(surv.dat$ejection.fraction.scale, na.rm = TRUE),
                                    mean(surv.dat$age.scale, na.rm = TRUE),
                                    mean(surv.dat$sbp.scale, na.rm = TRUE),
                                    mean(surv.dat$bmi.scale, na.rm = TRUE),
                                    mean(surv.dat$creatinine.scale, na.rm = TRUE),
                                    mean(surv.dat$nyha.class.2, na.rm = TRUE),
                                    mean(surv.dat$nyha.class.3, na.rm = TRUE),
                                    mean(surv.dat$nyha.class.4, na.rm = TRUE),
                                    mean(surv.dat$gender, na.rm = TRUE),
                                    mean(surv.dat$smoker, na.rm = TRUE),
                                    mean(surv.dat$diabetes, na.rm = TRUE),
                                    mean(surv.dat$copd, na.rm = TRUE),
                                    mean(surv.dat$recent.diagnosis, na.rm = TRUE),
                                    mean(surv.dat$no.beta.blocker, na.rm = TRUE),
                                    mean(surv.dat$no.ace, na.rm = TRUE),
                                    mean(surv.dat$age.ef.scale, na.rm = TRUE),
                                    mean(surv.dat$sbp.ef.scale, na.rm = TRUE)))

# Calculate parameters for the average patient in UC
# Exponential UC
reg.coef.hosp.exp.uc  <- as.numeric(coef(surv.hosp.exp))
rate.hosp.exp.uc <- exp(reg.coef.hosp.exp.uc[1]+sum(tail(reg.coef.hosp.exp.uc,n=-1)*average.patient.uc))

# Weibull UC
reg.coef.hosp.weibull.uc  <- as.numeric(coef(surv.hosp.weibull))
shape.hosp.weibull.uc <- exp(reg.coef.hosp.weibull.uc[1])
scale.hosp.weibull.uc <- exp(sum(tail(reg.coef.hosp.weibull.uc,n=-1)*c(1,average.patient.uc)))

# Log-normal UC
reg.coef.hosp.lnorm.uc  <- as.numeric(coef(surv.hosp.lnorm))
meanlog.hosp.lnorm.uc <- reg.coef.hosp.lnorm.uc[1]+sum(tail(reg.coef.hosp.lnorm.uc,n=-2)*average.patient.uc)
sdlog.hosp.lnorm.uc <- exp(reg.coef.hosp.lnorm.uc[2])

# Log-logistic UC
reg.coef.hosp.llogis.uc  <- as.numeric(coef(surv.hosp.llogis))
shape.hosp.llogis.uc <- exp(reg.coef.hosp.llogis.uc[1])
scale.hosp.llogis.uc <- exp(reg.coef.hosp.llogis.uc[2]+sum(tail(reg.coef.hosp.llogis.uc,n=-2)*average.patient.uc))

# Gompertz UC
reg.coef.hosp.gompertz.uc  <- as.numeric(coef(surv.hosp.gompertz))
shape.hosp.gompertz.uc <- reg.coef.hosp.gompertz.uc[1]
rate.hosp.gompertz.uc <- exp(reg.coef.hosp.gompertz.uc[2]+sum(tail(reg.coef.hosp.gompertz.uc,n=-2)*average.patient.uc))  

# Generalised gamma UC
reg.coef.hosp.gengamma.uc  <- as.numeric(coef(surv.hosp.gengamma))
mu.hosp.gengamma.uc <- reg.coef.hosp.gengamma.uc[1]+sum(tail(reg.coef.hosp.gengamma.uc,n=-3)*average.patient.uc)
sigma.hosp.gengamma.uc <- exp(reg.coef.hosp.gengamma.uc[2])
Q.hosp.gengamma.uc <- reg.coef.hosp.gengamma.uc[3]


# Calculate parameters for the average patient in EWS
# Exponential EWS
reg.coef.hosp.exp.ews  <- as.numeric(coef(surv.hosp.exp))
rate.hosp.exp.ews <- exp(reg.coef.hosp.exp.ews[1]+sum(tail(reg.coef.hosp.exp.ews,n=-1)*average.patient.ews))

# Weibull EWS
reg.coef.hosp.weibull.ews  <- as.numeric(coef(surv.hosp.weibull))
shape.hosp.weibull.ews <- exp(reg.coef.hosp.weibull.ews[1])
scale.hosp.weibull.ews <- exp(sum(tail(reg.coef.hosp.weibull.ews,n=-1)*c(1,average.patient.ews)))

# Log-normal EWS
reg.coef.hosp.lnorm.ews  <- as.numeric(coef(surv.hosp.lnorm))
meanlog.hosp.lnorm.ews <- reg.coef.hosp.lnorm.ews[1]+sum(tail(reg.coef.hosp.lnorm.ews,n=-2)*average.patient.ews)
sdlog.hosp.lnorm.ews <- exp(reg.coef.hosp.lnorm.ews[2])

# Log-logistic EWS
reg.coef.hosp.llogis.ews  <- as.numeric(coef(surv.hosp.llogis))
shape.hosp.llogis.ews <- exp(reg.coef.hosp.llogis.ews[1])
scale.hosp.llogis.ews <- exp(reg.coef.hosp.llogis.ews[2]+sum(tail(reg.coef.hosp.llogis.ews,n=-2)*average.patient.ews))

# Gompertz EWS
reg.coef.hosp.gompertz.ews  <- as.numeric(coef(surv.hosp.gompertz))
shape.hosp.gompertz.ews <- reg.coef.hosp.gompertz.ews[1]
rate.hosp.gompertz.ews <- exp(reg.coef.hosp.gompertz.ews[2]+sum(tail(reg.coef.hosp.gompertz.ews,n=-2)*average.patient.ews))  

# Generalised gamma EWS
reg.coef.hosp.gengamma.ews  <- as.numeric(coef(surv.hosp.gengamma))
mu.hosp.gengamma.ews <- reg.coef.hosp.gengamma.ews[1]+sum(tail(reg.coef.hosp.gengamma.ews,n=-3)*average.patient.ews)
sigma.hosp.gengamma.ews <- exp(reg.coef.hosp.gengamma.ews[2])
Q.hosp.gengamma.ews <- reg.coef.hosp.gengamma.ews[3]



# Plot KM curves and parametric survival model fitting curves (hospitalisation) - intervention as the only covariate
if (generate.plots==T) {
tiff('./Graphs/Parametric survival model fitting (hospitalisation).tiff',
     units="cm", width=23, height=15, res=1200, compression = 'lzw')

  plot(hosp.fit, main = "Parametric survival model fitting (hospitalisation)", xlab = "time (days)", ylab = "survival probability", xlim = c(0, time.horizon), conf.int	= F, lwd = 2, lty = c(1,2))
  lines(surv.hosp.exp.int.only, ci =F, t = c(0:time.horizon), lwd = 1, lty = 3, col = BM.orange)
  lines(surv.hosp.weibull.int.only, ci =F, t = c(0:time.horizon), lwd = 1, lty = 4, col = BM.Dyellow)
  lines(surv.hosp.lnorm.int.only, ci =F,t = c(0:time.horizon), lwd = 1, lty = 5, col = BM.pink)
  lines(surv.hosp.llogis.int.only, ci =F, t = c(0:time.horizon), lwd = 1, lty = 6, col = BM.blue)
  lines(surv.hosp.gompertz.int.only, ci =F, t = c(0:time.horizon), lwd = 1, lty = 7, col = BM.Dblue)
  lines(surv.hosp.gengamma.int.only, ci =F, t = c(0:time.horizon), lwd = 1, lty = 8, col = BM.yellow)
  legend(x = "bottomleft", col = c(1, 1, BM.orange, BM.Dyellow, BM.pink, BM.blue, BM.Dblue, BM.yellow),
         lwd = 2, lty = 1:8, bty = "n", cex = 0.8,
         legend = c("UC", "EWS", "exponential", "weibull", "log-normal", "log-logistic", "gompertz", "gengamma"))
  

dev.off()
}


# Plot KM curves and parametric survival model fitting curves (hospitalisation) - all covariates
# UC
if (generate.plots==T) {
  tiff('./Graphs/Parametric survival model fitting, all covariates, UC (hospitalisation).tiff',
       units="cm", width=23, height=15, res=1200, compression = 'lzw')
  
  plot(survfit(formula = surv.obj.hosp ~ 1, data = subset(surv.dat, surv.dat$intervention==0)), main = "Parametric survival model fitting, UC (hospitalisation)", xlab = "time (days)", ylab = "survival probability", xlim = c(0, 3652.5), conf.int	= F)
  curve((1-pexp(x, rate = rate.hosp.exp.uc)), xlim = c(0,time.horizon), lwd = 1, lty = 3, col = BM.orange, add = TRUE)
  curve((1-pweibull(x, shape = shape.hosp.weibull.uc, scale = scale.hosp.weibull.uc)), xlim = c(0,time.horizon), lwd = 1, lty = 4, col = BM.Dyellow, add = TRUE)
  curve((1-plnorm(x, meanlog = meanlog.hosp.lnorm.uc, sdlog = sdlog.hosp.lnorm.uc)), xlim = c(0,time.horizon), lwd = 1, lty = 5, col = BM.pink, add = TRUE)
  curve((1-pllogis(x, shape = shape.hosp.llogis.uc, scale = scale.hosp.llogis.uc)), xlim = c(0,time.horizon), lwd = 1, lty = 6, col = BM.blue, add = TRUE)
  curve((1-pgompertz(x, shape = shape.hosp.gompertz.uc, rate = rate.hosp.gompertz.uc)), xlim = c(0,time.horizon), lwd = 1, lty = 7, col = BM.Dblue, add = TRUE)
  curve((1-pgengamma(x, mu = mu.hosp.gengamma.uc, sigma = sigma.hosp.gengamma.uc, Q = Q.hosp.gengamma.uc)), xlim = c(0,time.horizon), lwd = 1, lty = 8, col = BM.yellow, add = TRUE)
  legend(x = "bottomleft", col = c(1, BM.orange, BM.Dyellow, BM.pink, BM.blue, BM.Dblue, BM.yellow),
         lwd = 2, lty = 1:8, bty = "n", cex = 0.8,
         legend = c("UC", "exponential", "weibull", "log-normal", "log-logistic", "gompertz", "gengamma"))
  
  dev.off()
}

# EWS
if (generate.plots==T) {
  tiff('./Graphs/Parametric survival model fitting, all covariates, EWS (hospitalisation).tiff',
       units="cm", width=23, height=15, res=1200, compression = 'lzw')
  
  plot(survfit(formula = surv.obj.hosp ~ 1, data = subset(surv.dat, surv.dat$intervention==1)), main = "Parametric survival model fitting, EWS (hospitalisation)", xlab = "time (days)", ylab = "survival probability", xlim = c(0, 3652.5), conf.int	= F)
  curve((1-pexp(x, rate = rate.hosp.exp.ews)), xlim = c(0,time.horizon), lwd = 1, lty = 3, col = BM.orange, add = TRUE)
  curve((1-pweibull(x, shape = shape.hosp.weibull.ews, scale = scale.hosp.weibull.ews)), xlim = c(0,time.horizon), lwd = 1, lty = 4, col = BM.Dyellow, add = TRUE)
  curve((1-plnorm(x, meanlog = meanlog.hosp.lnorm.ews, sdlog = sdlog.hosp.lnorm.ews)), xlim = c(0,time.horizon), lwd = 1, lty = 5, col = BM.pink, add = TRUE)
  curve((1-pllogis(x, shape = shape.hosp.llogis.ews, scale = scale.hosp.llogis.ews)), xlim = c(0,time.horizon), lwd = 1, lty = 6, col = BM.blue, add = TRUE)
  curve((1-pgompertz(x, shape = shape.hosp.gompertz.ews, rate = rate.hosp.gompertz.ews)), xlim = c(0,time.horizon), lwd = 1, lty = 7, col = BM.Dblue, add = TRUE)
  curve((1-pgengamma(x, mu = mu.hosp.gengamma.ews, sigma = sigma.hosp.gengamma.ews, Q = Q.hosp.gengamma.ews)), xlim = c(0,time.horizon), lwd = 1, lty = 8, col = BM.yellow, add = TRUE)
  legend(x = "bottomleft", col = c(1, BM.orange, BM.Dyellow, BM.pink, BM.blue, BM.Dblue, BM.yellow),
         lwd = 2, lty = 1:8, bty = "n", cex = 0.8,
         legend = c("EWS", "exponential", "weibull", "log-normal", "log-logistic", "gompertz", "gengamma"))
  
  dev.off()
}


# Compare goodness-of-fit
# Using AIC
surv.hosp.exp$AIC
surv.hosp.weibull$AIC
surv.hosp.lnorm$AIC
surv.hosp.llogis$AIC
surv.hosp.gompertz$AIC
surv.hosp.gengamma$AIC

# Using BIC
BIC(surv.hosp.exp)
BIC(surv.hosp.weibull)
BIC(surv.hosp.lnorm)
BIC(surv.hosp.llogis)
BIC(surv.hosp.gompertz)
BIC(surv.hosp.gengamma)

# Create and export goodness-of-fit table
parametric.distributions <- c("Exponential", "Weibull", "log-normal", "log-logistic", "Gompertz", "Gengamma")
aic.surv.hosp <- c(surv.hosp.exp$AIC, surv.hosp.weibull$AIC, surv.hosp.lnorm$AIC, surv.hosp.llogis$AIC, surv.hosp.gompertz$AIC, surv.hosp.gengamma$AIC)
bic.surv.hosp <- c(BIC(surv.hosp.exp), BIC(surv.hosp.weibull), BIC(surv.hosp.lnorm), BIC(surv.hosp.llogis), BIC(surv.hosp.gompertz), BIC(surv.hosp.gengamma))

goodness.fit.hosp <- data.frame(parametric.distributions, aic.surv.hosp, bic.surv.hosp)
names(goodness.fit.hosp) <- c("Distribution", "AIC", "BIC")
write.csv(goodness.fit.hosp, file = "./Results/Survival analyses/goodness.fit.hosp.csv")


#### Save model for hospitalisation ####
if (save.models==T) {
save(surv.hosp.exp, surv.hosp.weibull, surv.hosp.lnorm, surv.hosp.llogis, surv.hosp.gompertz, surv.hosp.gengamma,
     file = "./Regression models/surv_hosp_model.rda")
}
