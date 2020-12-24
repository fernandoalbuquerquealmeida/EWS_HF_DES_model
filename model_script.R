#### Start of script ####
# Main script for running model and outputting results
# Text to the right of the "#" sign are comments and are not processed
# Comments are an essential part of model code and should be used to describe what the code is doing
# This promotes transparency and helps avoid and identify errors

# Remove previous contents
rm(list=ls())

# Set the working directory by filling in the text between quotation marks in line 13 of the code with the directory where you have the model files (copy address as text)
# Use / rather than \\ symbol
# Delete # before line 13
# setwd("")

# Load functions from model_functions.R
# NOTE: Files must be in the same directory (working directory)
source("model_functions.R")

# Load regression models for survival and death in hospital (Regression models folder should be in working directory)
load("./Regression models/surv_death_model.rda")
load("./Regression models/surv_hosp_model.rda")
load("./Regression models/death_in_hospital_model.rda")

# Required pacakages for running model
required.packages <- c("survival", "flexsurv", "lme4", "ggplot2", "survminer", "gridExtra")


#### Set options ####
dbg.mode     <- T   # Debug mode: prints the simulated events list and accrued costs, life years and qalys as the simulation is being run (T/F)
ind.pat.data <- T   # Determines whether individual level patient details are stored and output (T/F)
export.excel <- T   # Set to T for exporting results to Excel
new.patients <- T   # Set to T for creating a new set of simulated patients and F for keeping current set of simulated patients

# Deterministic analysis
deterministic <- T    # Perform deterministic analysis? (T/F)
npats <- 1000         # Number of simulated patients in the deterministic analysis

# Probabilistic sensitivity analysis
probabilistic <- F    # Perform probabilistic analysis?
npats.psa <- 10     # Number of patients in the probabilistic sensitivity analysis simulations
nloops.psa <- 5      # Number of loops in the probabilistic sensitivity analysis simulations


#### Define parameters ####
# Interventions: usual care (uc), early warning system (ews), or early warning system + diagnostic algorithm (ewsda)
# Global variables
# Discount rates
drly <- 0.00000001   # Discount rate for life years
drq  <- 0.015   # Discount rate for QALYs
drc  <- 0.040   # Discount rate for costs

# Parametric survival distributions: "exponential", "weibull", "lognormal", "loglogistic", "gompertz", "gengamma"
dist.death.uc  <- "weibull"     # Death UC
dist.death.ews <- "weibull"     # Death EWS
dist.hosp.uc   <- "lognormal"     # Hospitalisation UC
dist.hosp.ews  <- "lognormal"     # Hospitalisation EWS

# Diagnostic algorithm characteristics
sensitivity         <- 0.96
false.positive.rate <- 0.54                 # 1 - specificity
avoid.prob          <- 0.50                 # Probability of avoiding hospitalisation given alarm
se.avoid.prob       <- avoid.prob * 0.1     # Standard error for probability of avoiding hospitalisation given alarm

# Utility multipliers
umult.outpat    <- 1                        # Utility multiplier for outpatient visit
se.umult.outpat <- umult.outpat * 0.1       # Standard error for utility multiplier for outpatient visit
umult.hosp 	    <- 0.82                     # Utility multiplier for hospitalisation
se.umult.hosp   <- umult.hosp * 0.1         # Standard error for utility multiplier for hospitalisation

# Time to outpatient visit (years)
time.outpat.visit.uc     <- 0.234           # For UC
se.time.outpat.visit.uc  <- 0.017           # Standard error for UC
time.outpat.visit.ews    <- 0.141           # For EWS
se.time.outpat.visit.ews <- 0.005           # Standard error for EWS

# Costs
# Event costs (per event)
cost.outpat.uc	   <- 46.33                     # Cost of outpatient visit in usual care arm
se.cost.outpat.uc  <- cost.outpat.uc * 0.1      # Standard error for cost of outpatient visit in usual care arm
cost.outpat.ews    <- 44.63                     # Cost of outpatient visit in EWS arm
se.cost.outpat.ews <- cost.outpat.ews * 0.1     # Standard error for cost of outpatient visit in EWS arm
cost.hosp          <- 4937.36                   # Cost of hospitalisation
se.cost.hosp       <- cost.hosp * 0.1           # Standard error for cost of hospitalisation
cost.death         <- 1                         # Cost of death
se.cost.death      <- cost.death * 0.1          # Standard error for cost of death

# Upkeep costs (per patient per year)
cost.uc     <- 705.71                           # UC costs
se.cost.uc  <- cost.uc * 0.1                    # Standard error for UC costs
cost.ews    <- 2621.70                          # EWS costs
se.cost.ews <- cost.ews * 0.1                   # Standard error for EWS costs

# Alarm costs (per alarm)
cost.fpalarm    <- 18                           # Cost of managing  a false positive alarm
se.cost.fpalarm <- cost.fpalarm * 0.1           # Standard error for cost of managing  a false positive alarm


#### Run the simulation ####
set.seed(9)                           # In order to fix simulation results, use set.seed(). Useful for debugging as removes random variation; define any integer
RunSim(deterministic, probabilistic)  # Start simulation (run this line together with the code from the start, after setting desired options and variable inputs)











#### Examine results ####
#### Summary statistics (per patient) ####
# Total (discounted) costs for ewsda arm
tot.costs.ewsda/npats

# Total (discounted) qalys for ewsda arm
tot.qalys.ewsda/npats

# Total (discounted) life years for ewsda arm
tot.lifeyears.ewsda/npats

# Total (discounted) costs for ews arm
tot.costs.ews/npats

# Total (discounted) qalys for ews arm
tot.qalys.ews/npats

# Total (discounted) life years for ews arm
tot.lifeyears.ews/npats

# Total (discounted) costs for uc arm
tot.costs.uc/npats

# Total (discounted) qalys for uc arm
tot.qalys.uc/npats

# Total (discounted) life years for uc arm
tot.lifeyears.uc/npats



# Total difference in costs between ewsda and ews
tot.dcosts.ews.ewsda/npats

# Total difference in qalys between ewsda and ews
tot.dqalys.ews.ewsda/npats

# Total difference in life years between ewsda and ews
tot.dlifeyears.ews.ewsda/npats

# Total difference in costs between ewsda and uc
tot.dcosts.uc.ewsda/npats

# Total difference in qalys between ewsda and uc
tot.dqalys.uc.ewsda/npats

# Total difference in life years between ewsda and uc
tot.dlifeyears.uc.ewsda/npats

# Total difference in costs between ews and uc
tot.dcosts.uc.ews/npats

# Total difference in qalys between ews and uc
tot.dqalys.uc.ews/npats

# Total difference in life years between ews and uc
tot.dlifeyears.uc.ews/npats


#### Number of events (per patient) ####
# Number of deaths not related to hospitalisation in ewsda arm
sum(sapply(PatData, function(x) x$ewsda$deathother))/npats

# Number of deaths from hospitalisation in ewsda arm
sum(sapply(PatData, function(x) x$ewsda$deathhosp))/npats

# Total number of hospitalisations in ewsda arm
sum(sapply(PatData, function(x) x$ewsda$nhosp))/npats

# Total number of avoided hospitalisations in ewsda arm
sum(sapply(PatData, function(x) x$ewsda$navoidhosp))/npats

# Total number of outpatient visits in ewsda arm
sum(sapply(PatData, function(x) x$ewsda$noutpat))/npats


# Number of deaths not related to hospitalisation in ews arm
sum(sapply(PatData, function(x) x$ews$deathother))/npats

# Number of deaths from hospitalisation in ews arm
sum(sapply(PatData, function(x) x$ews$deathhosp))/npats

# Total number of hospitalisations in ews arm
sum(sapply(PatData, function(x) x$ews$nhosp))/npats

# Total number of outpatient visits in ews arm
sum(sapply(PatData, function(x) x$ews$noutpat))/npats


# Number of deaths not related to hospitalisation in uc arm
sum(sapply(PatData, function(x) x$uc$deathother))/npats

# Number of deaths from hospitalisation in uc arm
sum(sapply(PatData, function(x) x$uc$deathhosp))/npats

# Total number of hospitalisations in uc arm
sum(sapply(PatData, function(x) x$uc$nhosp))/npats

# Total number of outpatient visits in uc arm
sum(sapply(PatData, function(x) x$uc$noutpat))/npats

#### Individual patient information ####
# For example, the 9th patient (change number inside [[]] for different patient)
PatData[[9]]$ewsda   # ewsda copy of the patient

PatData[[9]]$ews     # ews copy of the patient

PatData[[9]]$uc      # uc intervention copy of the patient

# To access the initial characteristics of the patient
PatData[[9]]$ewsda$patchar
PatData[[9]]$ews$patchar
PatData[[9]]$uc$patchar

# To access the patient history
PatData[[9]]$ewsda$pathist
PatData[[9]]$ews$pathist
PatData[[9]]$uc$pathist

#### ICERs ####
# ICER (cost/LY) for ewsda vs ews
tot.dcosts.ews.ewsda / tot.dlifeyears.ews.ewsda

# ICER (cost/QALY) for ewsda vs ews
tot.dcosts.ews.ewsda / tot.dqalys.ews.ewsda

# ICER (cost/LY) for ewsda vs UC
tot.dcosts.uc.ewsda / tot.dlifeyears.uc.ewsda

# ICER (cost/QALY) for ewsda vs UC
tot.dcosts.uc.ewsda / tot.dqalys.uc.ewsda

# ICER (cost/LY) for ews vs UC
tot.dcosts.uc.ews / tot.dlifeyears.uc.ews

# ICER (cost/QALY) for ews vs UC
tot.dcosts.uc.ews / tot.dqalys.uc.ews
