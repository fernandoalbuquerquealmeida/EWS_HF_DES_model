# Regression analysis for generating model used for calculating probability of dying in hospital
# Text to the right of the "#" sign are comments and are not processed
# Comments are an essential part of model code and should be used to describe what the code is doing
# This promotes transparency and helps avoid and identify errors

#### Instalation of required packages ####
# Remove previous contents
rm(list=ls())

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
required.packages <- c("lme4")

# Call the Function
InstallAndLoad(required.packages)

#### Data import, cleaning and formatting ####
# Set the working directory by replacing the text "C:/DES R model" with your directory where you have the model files (copy address as text)
# NOTE: Use / rather than \\ symbol
# NOTE: Put the working directory between ""
setwd("C:/Users/Fernando/Documents/PhD/TEN-HMS (Fernando)/Discrete event simulation model")

# Import data from .txt files
pat.char <- read.csv("./Data/baseline_population.csv")
hosp.dat <- read.csv("./Data/hospitalisation_data.txt", sep = ';', stringsAsFactors = FALSE, na.strings = "")

# Convert dates to a date vectors readable in R
for(i in c(2,3,4,6)) {
    hosp.dat[,i] <- as.Date(hosp.dat[,i], format = "%d-%m-%Y")
}

# Merge hospitalisation data with baseline characteristics data
names(hosp.dat)[names(hosp.dat) == "patient.study.number"] <- "Patient.Study.Number"
hosp.dat <- merge(pat.char, hosp.dat, by = "Patient.Study.Number")

# Replace values for date of birth that are incorrect and/or missing with the average of the population
hosp.dat$date.of.birth[is.na(hosp.dat$date.of.birth)] <- mean(hosp.dat$date.of.birth, na.rm = TRUE)
hosp.dat$date.of.birth[hosp.dat$age < 10] <- mean(hosp.dat$date.of.birth, na.rm = TRUE)

hosp.dat$age <- (hosp.dat$admit.date - hosp.dat$date.of.birth) / 365

hosp.dat$previous.hosp <- ifelse(hosp.dat$hospitalisation.number == 1, 0, 1)

max.hosp <- aggregate(hospitalisation.number ~ Patient.Study.Number, hosp.dat, max)

colnames(max.hosp) <- c("Patient.Study.Number", "last.hosp.number")

hosp.dat <- merge(hosp.dat, max.hosp, by = "Patient.Study.Number", all.x = TRUE)

hosp.dat$death <- ifelse(is.na(hosp.dat$death.in.hospital), 0,ifelse(hosp.dat$death.in.hospital == "yes"
                         & hosp.dat$hospitalisation.number == hosp.dat$last.hosp.number, 1, 0))

hosp.dat$age <- as.numeric(hosp.dat$age)

# Fit regression model
p.death.hosp <- glm(factor(death) ~ age + factor(gender) + factor(myocardial.infarction) +
                    factor(chronic.atrial.fibrillation) + factor(diabetes) + factor(copd) + previous.hosp,
                    data = hosp.dat, family = binomial('logit'))

# Show model results
summary(p.death.hosp)

res.p.death.hosp <- as.data.frame(round(coef(summary(p.death.hosp)),3))
write.csv(res.p.death.hosp, file = "./Results/Survival analyses/res.p.death.hosp.csv")


# Save model
save(p.death.hosp, file = "./Regression models/death_in_hospital_model.rda")

# Return the estimated probabilities for the data used in the regression equation
predict(p.death.hosp, type="response")

# Predict probability of dying for a hypothetical patient
predict(p.death.hosp, newdata = data.frame(age = 60, gender = 1, myocardial.infarction = 0,
        chronic.atrial.fibrillation = 0, diabetes = 1, copd = 0, previous.hosp = 0), type="response")
