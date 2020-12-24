# Functions called by model_script.R that represent the "model engine"
# Text to the right of the "#" sign are comments that are not processed by R
# Comments are an essential part of model code and should be used to describe what the code is doing
# This promotes transparency and helps avoid and identify errors

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

# IndividualPatientData: store individual patient data

IndividualPatientData <- function(ind.pat.data,npats){

if (ind.pat.data==T){
    PatData <<- vector("list", length = npats) # empty list with number of patients
}
    
}

# CreatePatients: create individual patient list
CreatePatients <- function(n, new.patients){
    
    if (new.patients == T){
        
    pat.list <- read.csv("./Data/patient_list.csv", sep = ",")
    
    output   <- pat.list[sample(nrow(pat.list), size = n, replace = TRUE),]
        
    write.csv(output, file = "./Data/simulated_patients.csv", row.names = FALSE)
    
    return(output)
    
    }
}

# SurvivalModels: determine survival models used in the simulation

SurvivalModels <- function(dist.death.uc, dist.death.ews, dist.hosp.uc, dist.hosp.ews){
    
    if (dist.death.uc=="exponential") {
        surv.death.uc <<- surv.death.exp
        
    } else if (dist.death.uc=="weibull") {
        surv.death.uc <<- surv.death.weibull
        
    } else if (dist.death.uc=="lognormal") {
        surv.death.uc <<- surv.death.lnorm
        
    } else if (dist.death.uc=="loglogistic") {
        surv.death.uc <<- surv.death.llogis
        
    } else if (dist.death.uc=="gompertz") {
        surv.death.uc <<- surv.death.gompertz
        
    } else if (dist.death.uc=="gengamma") {
        surv.death.uc <<- surv.death.gengamma}
    
    if (dist.death.ews=="exponential") {
        surv.death.ews <<- surv.death.exp
        
    } else if (dist.death.ews=="weibull") {
        surv.death.ews <<- surv.death.weibull
        
    } else if (dist.death.ews=="lognormal") {
        surv.death.ews <<- surv.death.lnorm
        
    } else if (dist.death.ews=="loglogistic") {
        surv.death.ews <<- surv.death.llogis
        
    } else if (dist.death.ews=="gompertz") {
        surv.death.ews <<- surv.death.gompertz
        
    } else if (dist.death.ews=="gengamma") {
        surv.death.ews <<- surv.death.gengamma}
    
    if (dist.hosp.uc=="exponential") {
        surv.hosp.uc <<- surv.hosp.exp
        
    } else if (dist.hosp.uc=="weibull") {
        surv.hosp.uc <<- surv.hosp.weibull
        
    } else if (dist.hosp.uc=="lognormal") {
        surv.hosp.uc <<- surv.hosp.lnorm
        
    } else if (dist.hosp.uc=="loglogistic") {
        surv.hosp.uc <<- surv.hosp.llogis
        
    } else if (dist.hosp.uc=="gompertz") {
        surv.hosp.uc <<- surv.hosp.gompertz
        
    } else if (dist.hosp.uc=="gengamma") {
        surv.hosp.uc <<- surv.hosp.gengamma}
        
    if (dist.hosp.ews=="exponential") {
        surv.hosp.ews <<- surv.hosp.exp
        
    } else if (dist.hosp.ews=="weibull") {
        surv.hosp.ews <<- surv.hosp.weibull
        
    } else if (dist.hosp.ews=="lognormal") {
        surv.hosp.ews <<- surv.hosp.lnorm
        
    } else if (dist.hosp.ews=="loglogistic") {
        surv.hosp.ews <<- surv.hosp.llogis
        
    } else if (dist.hosp.ews=="gompertz") {
        surv.hosp.ews <<- surv.hosp.gompertz
        
    } else if (dist.hosp.ews=="gengamma") {
        surv.hosp.ews <<- surv.hosp.gengamma}
    
}

# Interventions: usual care (uc), early warning system (ews), or early warning system + diagnostic algorithm (ewsda)

# GetNextEvent: initialise the event list and identify which of the events will be processed by the ReactEvt function
GetNextEvent <- function(intervention){

    pat.char.tte <- as.numeric(c(curpat$intervention,
                                 (curpat$ejection.fraction - mean(pat.list$ejection.fraction))/sd(pat.list$ejection.fraction),
                                 (curpat$age - mean(pat.list$age))/sd(pat.list$age),
                                 (curpat$sbp - mean(pat.list$sbp))/sd(pat.list$sbp),
                                 (curpat$bmi - mean(pat.list$bmi))/sd(pat.list$bmi),
                                 (curpat$creatinine - mean(pat.list$creatinine))/sd(pat.list$creatinine),
                                 curpat$nyha.class.2,
                                 curpat$nyha.class.3,
                                 curpat$nyha.class.4,
                                 curpat$gender,
                                 curpat$smoker,
                                 curpat$diabetes,
                                 curpat$copd,
                                 curpat$recent.diagnosis,
                                 curpat$no.beta.blocker,
                                 curpat$no.ace,
                                 (curpat$age.ef - mean(pat.list$age.ef))/sd(pat.list$age.ef),
                                 (curpat$sbp.ef - mean(pat.list$sbp.ef))/sd(pat.list$sbp.ef)))
    
    if (intervention=="uc"){

    ## Time-to-outpatient_visit
    
    outpat      <- parameters$time.outpat.visit.uc
    
    
    ## Time-to-hospitalisation
    
    if (dist.hosp.uc=="exponential") {
        
        rate.hosp.exp <- exp(reg.coef.hosp.uc[1]+sum(tail(reg.coef.hosp.uc,n=-1)*pat.char.tte))
        
        hosp       <- qexp(simulation.draws[simulation.index,patient.index], rate = rate.hosp.exp)/365.25
        
    } else if (dist.hosp.uc=="weibull") {
        
        shape.hosp.weibull <- exp(reg.coef.hosp.uc[1])
        
        scale.hosp.weibull <- exp(sum(tail(reg.coef.hosp.uc,n=-1)*c(1,pat.char.tte)))
        
        hosp       <- qweibull(simulation.draws[simulation.index,patient.index], shape = shape.hosp.weibull, scale = scale.hosp.weibull)/365.25
        
    } else if (dist.hosp.uc=="lognormal") {
        
        meanlog.hosp.lnorm <- reg.coef.hosp.uc[1]+sum(tail(reg.coef.hosp.uc,n=-2)*pat.char.tte)
        
        sdlog.hosp.lnorm <- exp(reg.coef.hosp.uc[2])
        
        hosp       <- qlnorm(simulation.draws[simulation.index,patient.index], meanlog = meanlog.hosp.lnorm, sdlog = sdlog.hosp.lnorm)/365.25
        
    } else if (dist.hosp.uc=="loglogistic") {
        
        shape.hosp.llogis <- exp(reg.coef.hosp.uc[1])
        
        scale.hosp.llogis <- exp(reg.coef.hosp.uc[2]+sum(tail(reg.coef.hosp.uc,n=-2)*pat.char.tte))
        
        hosp       <- qllogis(simulation.draws[simulation.index,patient.index], shape = shape.hosp.llogis, scale = scale.hosp.llogis)/365.25
        
    } else if (dist.hosp.uc=="gompertz") {
        
        shape.hosp.gompertz <- reg.coef.hosp.uc[1]
        
        rate.hosp.gompertz <- exp(reg.coef.hosp.uc[2]+sum(tail(reg.coef.hosp.uc,n=-2)*pat.char.tte))  
        
        hosp       <- qgompertz(simulation.draws[simulation.index,patient.index], shape = shape.hosp.gompertz, rate = rate.hosp.gompertz)/365.25
        
    } else if (dist.hosp.uc=="gengamma") {
        
        mu.hosp.gengamma <- reg.coef.hosp.uc[1]+sum(tail(reg.coef.hosp.uc,n=-3)*pat.char.tte)
        
        sigma.hosp.gengamma <- exp(reg.coef.hosp.uc[2])
        
        Q.hosp.gengamma <- reg.coef.hosp.uc[3]
        
        hosp       <- qgengamma(simulation.draws[simulation.index,patient.index], mu = mu.hosp.gengamma, sigma = sigma.hosp.gengamma, Q = Q.hosp.gengamma)/365.25
        
    }
    
    ## Time-to-death
        
    if (dist.death.uc=="exponential") {
        
        rate.death.exp <- exp(reg.coef.death.uc[1]+sum(tail(reg.coef.death.uc,n=-1)*pat.char.tte))
        
        death       <- qexp(simulation.draws[simulation.index+1,patient.index], rate = rate.death.exp)/365.25

    } else if (dist.death.uc=="weibull") {
        
        shape.death.weibull <- exp(reg.coef.death.uc[1])
        
        scale.death.weibull <- exp(sum(tail(reg.coef.death.uc,n=-1)*c(1,pat.char.tte)))
        
        death       <- qweibull(simulation.draws[simulation.index+1,patient.index], shape = shape.death.weibull, scale = scale.death.weibull)/365.25
        
    } else if (dist.death.uc=="lognormal") {
        
        meanlog.death.lnorm <- reg.coef.death.uc[1]+sum(tail(reg.coef.death.uc,n=-2)*pat.char.tte)
        
        sdlog.death.lnorm <- exp(reg.coef.death.uc[2])
        
        death       <- qlnorm(simulation.draws[simulation.index+1,patient.index], meanlog = meanlog.death.lnorm, sdlog = sdlog.death.lnorm)/365.25
        
    } else if (dist.death.uc=="loglogistic") {
        
        shape.death.llogis <- exp(reg.coef.death.uc[1])
        
        scale.death.llogis <- exp(reg.coef.death.uc[2]+sum(tail(reg.coef.death.uc,n=-2)*pat.char.tte))
        
        death       <- qllogis(simulation.draws[simulation.index+1,patient.index], shape = shape.death.llogis, scale = scale.death.llogis)/365.25

    } else if (dist.death.uc=="gompertz") {
        
        shape.death.gompertz <- reg.coef.death.uc[1]
        
        rate.death.gompertz <- exp(reg.coef.death.uc[2]+sum(tail(reg.coef.death.uc,n=-2)*pat.char.tte))  
        
        death       <- qgompertz(simulation.draws[simulation.index+1,patient.index], shape = shape.death.gompertz, rate = rate.death.gompertz)/365.25
        
    } else if (dist.death.uc=="gengamma") {
        
        mu.death.gengamma <- reg.coef.death.uc[1]+sum(tail(reg.coef.death.uc,n=-3)*pat.char.tte)
        
        sigma.death.gengamma <- exp(reg.coef.death.uc[2])
        
        Q.death.gengamma <- reg.coef.death.uc[3]
        
        death       <- qgengamma(simulation.draws[simulation.index+1,patient.index], mu = mu.death.gengamma, sigma = sigma.death.gengamma, Q = Q.death.gengamma)/365.25
        
        }
    

    evtlist     <- data.frame(
        
        # This code assigns a vector of names stored as strings and allows the user to refer to the corresponding event times by name
        
        evtname = c("outpat", "hosp", "death"),
        
        evttime = c(outpat, hosp, death)

    )  
    
    # Sort event list by event time
    evtlist <- evtlist[order(evtlist$evttime),]
    
    # Identify next event
    nextevt <- evtlist$evtname[1]
    nextevttime <- evtlist$evttime[1] 
    
    # Create a list containing name and time of next event
    output <- list(evt=nextevt, evttime=nextevttime)
    
    return(output)
    
} else if (intervention=="ews"){
   
    ## Time-to-outpatient_visit
    
    outpat      <- parameters$time.outpat.visit.ews
    
    
    ## Time-to-hospitalisation
    
    if (dist.hosp.ews=="exponential") {
        
        rate.hosp.exp <- exp(reg.coef.hosp.ews[1]+sum(tail(reg.coef.hosp.ews,n=-1)*pat.char.tte))
        
        hosp       <- qexp(simulation.draws[simulation.index,patient.index], rate = rate.hosp.exp)/365.25
        
    } else if (dist.hosp.ews=="weibull") {
        
        shape.hosp.weibull <- exp(reg.coef.hosp.ews[1])
        
        scale.hosp.weibull <- exp(sum(tail(reg.coef.hosp.ews,n=-1)*c(1,pat.char.tte)))
        
        hosp       <- qweibull(simulation.draws[simulation.index,patient.index], shape = shape.hosp.weibull, scale = scale.hosp.weibull)/365.25
        
    } else if (dist.hosp.ews=="lognormal") {
        
        meanlog.hosp.lnorm <- reg.coef.hosp.ews[1]+sum(tail(reg.coef.hosp.ews,n=-2)*pat.char.tte)
        
        sdlog.hosp.lnorm <- exp(reg.coef.hosp.ews[2])
        
        hosp       <- qlnorm(simulation.draws[simulation.index,patient.index], meanlog = meanlog.hosp.lnorm, sdlog = sdlog.hosp.lnorm)/365.25
        
    } else if (dist.hosp.ews=="loglogistic") {
        
        shape.hosp.llogis <- exp(reg.coef.hosp.ews[1])
        
        scale.hosp.llogis <- exp(reg.coef.hosp.ews[2]+sum(tail(reg.coef.hosp.ews,n=-2)*pat.char.tte))
        
        hosp       <- qllogis(simulation.draws[simulation.index,patient.index], shape = shape.hosp.llogis, scale = scale.hosp.llogis)/365.25
        
    } else if (dist.hosp.ews=="gompertz") {
        
        shape.hosp.gompertz <- reg.coef.hosp.ews[1]
        
        rate.hosp.gompertz <- exp(reg.coef.hosp.ews[2]+sum(tail(reg.coef.hosp.ews,n=-2)*pat.char.tte))  
        
        hosp       <- qgompertz(simulation.draws[simulation.index,patient.index], shape = shape.hosp.gompertz, rate = rate.hosp.gompertz)/365.25
        
    } else if (dist.hosp.ews=="gengamma") {
        
        mu.hosp.gengamma <- reg.coef.hosp.ews[1]+sum(tail(reg.coef.hosp.ews,n=-3)*pat.char.tte)
        
        sigma.hosp.gengamma <- exp(reg.coef.hosp.ews[2])
        
        Q.hosp.gengamma <- reg.coef.hosp.ews[3]
        
        hosp       <- qgengamma(simulation.draws[simulation.index,patient.index], mu = mu.hosp.gengamma, sigma = sigma.hosp.gengamma, Q = Q.hosp.gengamma)/365.25
        
    }
    
    ## Time-to-death
    
    if (dist.death.ews=="exponential") {
        
        rate.death.exp <- exp(reg.coef.death.ews[1]+sum(tail(reg.coef.death.ews,n=-1)*pat.char.tte))
        
        death       <- qexp(simulation.draws[simulation.index+1,patient.index], rate = rate.death.exp)/365.25
        
    } else if (dist.death.ews=="weibull") {
        
        shape.death.weibull <- exp(reg.coef.death.ews[1])
        
        scale.death.weibull <- exp(sum(tail(reg.coef.death.ews,n=-1)*c(1,pat.char.tte)))
        
        death       <- qweibull(simulation.draws[simulation.index+1,patient.index], shape = shape.death.weibull, scale = scale.death.weibull)/365.25
        
    } else if (dist.death.ews=="lognormal") {
        
        meanlog.death.lnorm <- reg.coef.death.ews[1]+sum(tail(reg.coef.death.ews,n=-2)*pat.char.tte)
        
        sdlog.death.lnorm <- exp(reg.coef.death.ews[2])
        
        death       <- qlnorm(simulation.draws[simulation.index+1,patient.index], meanlog = meanlog.death.lnorm, sdlog = sdlog.death.lnorm)/365.25
        
    } else if (dist.death.ews=="loglogistic") {
        
        shape.death.llogis <- exp(reg.coef.death.ews[1])
        
        scale.death.llogis <- exp(reg.coef.death.ews[2]+sum(tail(reg.coef.death.ews,n=-2)*pat.char.tte))
        
        death       <- qllogis(simulation.draws[simulation.index+1,patient.index], shape = shape.death.llogis, scale = scale.death.llogis)/365.25
        
    } else if (dist.death.ews=="gompertz") {
        
        shape.death.gompertz <- reg.coef.death.ews[1]
        
        rate.death.gompertz <- exp(reg.coef.death.ews[2]+sum(tail(reg.coef.death.ews,n=-2)*pat.char.tte))  
        
        death       <- qgompertz(simulation.draws[simulation.index+1,patient.index], shape = shape.death.gompertz, rate = rate.death.gompertz)/365.25
        
    } else if (dist.death.ews=="gengamma") {
        
        mu.death.gengamma <- reg.coef.death.ews[1]+sum(tail(reg.coef.death.ews,n=-3)*pat.char.tte)
        
        sigma.death.gengamma <- exp(reg.coef.death.ews[2])
        
        Q.death.gengamma <- reg.coef.death.ews[3]
        
        death       <- qgengamma(simulation.draws[simulation.index+1,patient.index], mu = mu.death.gengamma, sigma = sigma.death.gengamma, Q = Q.death.gengamma)/365.25
        
    }
    
    evtlist     <- data.frame(
        
        # This code assigns a vector of names stored as strings and allows the user to refer to the corresponding event times by name
        
        evtname = c("outpat", "hosp", "death"),
        
        evttime = c(outpat, hosp, death)
        
    )  
    
    # Sort event list by event time
    evtlist <- evtlist[order(evtlist$evttime, evtlist$evtname),]
    
    # Identify next event
    nextevt <- evtlist$evtname[1]
    nextevttime <- evtlist$evttime[1] 
    
    # Create a list containing name and time of next event
    output <- list(evt=nextevt, evttime=nextevttime)
    
    return(output)

} else if (intervention=="ewsda"){
    
    ## Time-to-outpatient_visit
    
    outpat      <- parameters$time.outpat.visit.ews
    
    
    ## Time-to-hospitalisation
    
    if (dist.hosp.ews=="exponential") {
        
        rate.hosp.exp <- exp(reg.coef.hosp.ews[1]+sum(tail(reg.coef.hosp.ews,n=-1)*pat.char.tte))
        
        hosp       <- qexp(simulation.draws[simulation.index,patient.index], rate = rate.hosp.exp)/365.25
        
    } else if (dist.hosp.ews=="weibull") {
        
        shape.hosp.weibull <- exp(reg.coef.hosp.ews[1])
        
        scale.hosp.weibull <- exp(sum(tail(reg.coef.hosp.ews,n=-1)*c(1,pat.char.tte)))
        
        hosp       <- qweibull(simulation.draws[simulation.index,patient.index], shape = shape.hosp.weibull, scale = scale.hosp.weibull)/365.25
        
    } else if (dist.hosp.ews=="lognormal") {
        
        meanlog.hosp.lnorm <- reg.coef.hosp.ews[1]+sum(tail(reg.coef.hosp.ews,n=-2)*pat.char.tte)
        
        sdlog.hosp.lnorm <- exp(reg.coef.hosp.ews[2])
        
        hosp       <- qlnorm(simulation.draws[simulation.index,patient.index], meanlog = meanlog.hosp.lnorm, sdlog = sdlog.hosp.lnorm)/365.25
        
    } else if (dist.hosp.ews=="loglogistic") {
        
        shape.hosp.llogis <- exp(reg.coef.hosp.ews[1])
        
        scale.hosp.llogis <- exp(reg.coef.hosp.ews[2]+sum(tail(reg.coef.hosp.ews,n=-2)*pat.char.tte))
        
        hosp       <- qllogis(simulation.draws[simulation.index,patient.index], shape = shape.hosp.llogis, scale = scale.hosp.llogis)/365.25
        
    } else if (dist.hosp.ews=="gompertz") {
        
        shape.hosp.gompertz <- reg.coef.hosp.ews[1]
        
        rate.hosp.gompertz <- exp(reg.coef.hosp.ews[2]+sum(tail(reg.coef.hosp.ews,n=-2)*pat.char.tte))  
        
        hosp       <- qgompertz(simulation.draws[simulation.index,patient.index], shape = shape.hosp.gompertz, rate = rate.hosp.gompertz)/365.25
        
    } else if (dist.hosp.ews=="gengamma") {
        
        mu.hosp.gengamma <- reg.coef.hosp.ews[1]+sum(tail(reg.coef.hosp.ews,n=-3)*pat.char.tte)
        
        sigma.hosp.gengamma <- exp(reg.coef.hosp.ews[2])
        
        Q.hosp.gengamma <- reg.coef.hosp.ews[3]
        
        hosp       <- qgengamma(simulation.draws[simulation.index,patient.index], mu = mu.hosp.gengamma, sigma = sigma.hosp.gengamma, Q = Q.hosp.gengamma)/365.25
        
    }
    
    
    ## Avoid hospitalisation?
    
    avoidhosp   <- ifelse(runif(1) < parameters$sensitivity * parameters$avoid.prob, hosp, Inf)
    
    
    ## Time-to-death
    
    if (dist.death.ews=="exponential") {
        
        rate.death.exp <- exp(reg.coef.death.ews[1]+sum(tail(reg.coef.death.ews,n=-1)*pat.char.tte))
        
        death       <- qexp(simulation.draws[simulation.index+1,patient.index], rate = rate.death.exp)/365.25
        
    } else if (dist.death.ews=="weibull") {
        
        shape.death.weibull <- exp(reg.coef.death.ews[1])
        
        scale.death.weibull <- exp(sum(tail(reg.coef.death.ews,n=-1)*c(1,pat.char.tte)))
        
        death       <- qweibull(simulation.draws[simulation.index+1,patient.index], shape = shape.death.weibull, scale = scale.death.weibull)/365.25
        
    } else if (dist.death.ews=="lognormal") {
        
        meanlog.death.lnorm <- reg.coef.death.ews[1]+sum(tail(reg.coef.death.ews,n=-2)*pat.char.tte)
        
        sdlog.death.lnorm <- exp(reg.coef.death.ews[2])
        
        death       <- qlnorm(simulation.draws[simulation.index+1,patient.index], meanlog = meanlog.death.lnorm, sdlog = sdlog.death.lnorm)/365.25
        
    } else if (dist.death.ews=="loglogistic") {
        
        shape.death.llogis <- exp(reg.coef.death.ews[1])
        
        scale.death.llogis <- exp(reg.coef.death.ews[2]+sum(tail(reg.coef.death.ews,n=-2)*pat.char.tte))
        
        death       <- qllogis(simulation.draws[simulation.index+1,patient.index], shape = shape.death.llogis, scale = scale.death.llogis)/365.25
        
    } else if (dist.death.ews=="gompertz") {
        
        shape.death.gompertz <- reg.coef.death.ews[1]
        
        rate.death.gompertz <- exp(reg.coef.death.ews[2]+sum(tail(reg.coef.death.ews,n=-2)*pat.char.tte))  
        
        death       <- qgompertz(simulation.draws[simulation.index+1,patient.index], shape = shape.death.gompertz, rate = rate.death.gompertz)/365.25
        
    } else if (dist.death.ews=="gengamma") {
        
        mu.death.gengamma <- reg.coef.death.ews[1]+sum(tail(reg.coef.death.ews,n=-3)*pat.char.tte)
        
        sigma.death.gengamma <- exp(reg.coef.death.ews[2])
        
        Q.death.gengamma <- reg.coef.death.ews[3]
        
        death       <- qgengamma(simulation.draws[simulation.index+1,patient.index], mu = mu.death.gengamma, sigma = sigma.death.gengamma, Q = Q.death.gengamma)/365.25
        
    }
    

    evtlist     <- data.frame(
        
        # This code assigns a vector of names stored as strings and allows the user to refer to the corresponding event times by name
    
        evtname = c("outpat", "hosp", "avoidhosp", "death"),
        
        evttime = c(outpat, hosp, avoidhosp, death)
        
    )  
        
    # Sort event list by event time
    evtlist <- evtlist[order(evtlist$evttime, evtlist$evtname),]
    
    # Identify next event
    nextevt <- evtlist$evtname[1]
    nextevttime <- evtlist$evttime[1] 
        
    # Create a list containing name and time of next event
    output <- list(evt=nextevt, evttime=nextevttime)
        
    return(output)
    
}
    }

# AddOngoing: calculate additional outcomes accrued from previous to current event
# Inputs:
# lcldrly: discount rate for life years
# lcldrq: discount rate for qalys
# lcldrc: discount rate for costs
# lclprevtime: time of previous event
# lclcurtime: time of current event
# lclvally: fixed value of life years
# lclvalq: fixed value of qalys
# lclvalc: fixed value of costs

# For compound continuous discounting we need to use instantaneous rate
InstantDiscount <- function(rate){     
    log(1+rate)                               
}

AddOngoing <- function(lcldrly=parameters$drly, lcldrq=parameters$drq, lcldrc=parameters$drc, lclprevtime, lclcurtime, lclvally, lclvalq, lclvalc){ 
    
    Instantdrly <- InstantDiscount(lcldrly)
    Instantdrq <- InstantDiscount(lcldrq)
    Instantdrc <- InstantDiscount(lcldrc)
    
    # calculate additional life years
    addlifeyears <- ((lclvally)/(0 - Instantdrly)) * (exp(lclcurtime * ( 0 - Instantdrly)) - exp(lclprevtime * (0 - Instantdrly)))
    
    # calculate additional qalys
    addqalys <- ((lclvalq)/(0 - Instantdrq)) * (exp(lclcurtime * ( 0 - Instantdrq)) - exp(lclprevtime * (0 - Instantdrq)))
    
    # calculate additional costs
    addcosts <- ((lclvalc)/(0 - Instantdrc)) * (exp(lclcurtime * ( 0 - Instantdrc)) - exp(lclprevtime * (0 - Instantdrc)))
    
    # combine additional costs, qalys, and life years in a list
    output <- list(addlifeyears=addlifeyears, addqalys=addqalys, addcosts=addcosts)
    
    return(output)
}

# AddInstant: add instantaneous outcomes
# Inputs:
# lcldrly: discount rate for life years
# lcldrq: discount rate for qalys
# lcldrc: discount rate for costs
# lclcurtime: time of current event
# lclvally: fixed value of life years
# lclvalq: fixed value of qalys 
# lclvalc: fixed value of costs

AddInstant <- function(lcldrly=parameters$drly, lcldrq=parameters$drq, lcldrc=parameters$drc, lclcurtime, lclvally, lclvalq, lclvalc){
    
    addinstlifeyears <- lclvally * ((1+lcldrly)^(-lclcurtime))    # Note use of DISCRETE TIME discounting for instantaneous costs and benefits
    addinstqalys <- lclvalq * ((1+lcldrq)^(-lclcurtime))
    addinstcosts <- lclvalc * ((1+lcldrc)^(-lclcurtime))
    
    # combine additional costs, qalys, and life years in a list
    output <- list(addinstlifeyears=addinstlifeyears, addinstqalys=addinstqalys, addinstcosts=addinstcosts)
    
    return(output)
}

# ReactEvt: react to the next event (as identified in the GetNextEvt function)
# thisevt: a two element list containing the output from GetNextEvt
# $evt: event 
# $evttime: event time

ReactEvt <- function(thisevt, intervention){
    
    evt <- thisevt$evt                       # Identify event type
    prevtime <- curtime                      # Identify time of previous event
    curtime  <<- thisevt$evttime + prevtime  # Identify time of current event
    
    pat.char.dth.hosp <- as.numeric(c(curpat$age,
                                      curpat$gender,
                                      curpat$myocardial.infarction,
                                      curpat$chronic.atrial.fibrillation,
                                      curpat$diabetes,
                                      curpat$copd,
                                      curpat$previous.hosp))
    
    log.odds.dth.hosp <- sum(reg.coef.dth.hosp*c(1, pat.char.dth.hosp)) # log-odds
    p.dth.hosp <- exp(log.odds.dth.hosp)/(1+exp(log.odds.dth.hosp))     # Probability of dying in hospital
    
    
    if (intervention=="uc"){
        
        # Debugging line: if debugging is enabled in model_script, R will print out event that happened so far
        if (dbg.mode==T){
            print(data.frame(evtname=evt,evttime=curtime), row.names = FALSE)
        }
        
        # Usual care logic
        
        if (evt=="outpat"){
            if (ind.pat.data==T){                                                  # Code chunks like this are processed only if individual patient info is being saved (ind=TRUE) in model_script
                this.PatData$uc$noutpat <<- this.PatData$uc$noutpat + 1   # This records the number of outpatient visits
            }
            
            additionals     <- AddOngoing(lclprevtime=prevtime, lclcurtime=curtime, lclvally=1, lclvalq=utilmlt, lclvalc=parameters$cost.uc)
            instadditionals <- AddInstant(lclcurtime=curtime, lclvally=0, lclvalq=0, lclvalc=parameters$cost.outpat.uc)
            
            thslifeyears <<- thslifeyears + additionals$addlifeyears + instadditionals$addinstlifeyears
            thsqalys <<- thsqalys + additionals$addqalys + instadditionals$addinstqalys
            thscosts <<- thscosts + additionals$addcosts + instadditionals$addinstcosts 
            
            utilmlt <<- ifelse(curpat$nyha.class < 4, utilmlt * parameters$umult.outpat, utilmlt)
            
        } else if (evt=="hosp"){  
            if (ind.pat.data==T){
                this.PatData$uc$nhosp <<- this.PatData$uc$nhosp + 1     # Records that patient was hospitalised
            }
            
            # Hospitalisation logic
            
            additionals     <- AddOngoing(lclprevtime=prevtime, lclcurtime=curtime, lclvally=1, lclvalq=utilmlt, lclvalc=parameters$cost.uc)
            instadditionals <- AddInstant(lclcurtime=curtime, lclvally=0, lclvalq=0, lclvalc=parameters$cost.hosp)
            
            thslifeyears <<- thslifeyears + additionals$addlifeyears + instadditionals$addinstlifeyears
            thsqalys <<- thsqalys + additionals$addqalys + instadditionals$addinstqalys
            thscosts <<- thscosts + additionals$addcosts + instadditionals$addinstcosts 
            
            utilmlt <<- ifelse(curpat$nyha.class < 4, utilmlt * parameters$umult.hosp, utilmlt)
            
            patdies <- simulation.draws[simulation.index+2,patient.index] < p.dth.hosp   # If the random number between 0 and 1 is lower than the predicted probability of dying in hospital (p.dth.hosp), "patdies" is set to TRUE

            if (patdies==T) { 
                if (ind.pat.data==T){
                    this.PatData$uc$deathhosp <<- this.PatData$uc$deathhosp + 1     # This indicates that the patient died in hospital
                }
                
                simulation.clock <<- Inf
            }
            
        } else if (evt=="death"){
            if(ind.pat.data==T){
                this.PatData$uc$deathother <<- 1     # This indicates that the patient died outside of the hospital
            }
            
            additionals     <- AddOngoing(lclprevtime=prevtime, lclcurtime=curtime, lclvally=1, lclvalq=utilmlt, lclvalc=parameters$cost.uc)
            instadditionals <- AddInstant(lclcurtime=curtime, lclvally=0, lclvalq=0, lclvalc=parameters$cost.death)
            
            thslifeyears <<- thslifeyears + additionals$addlifeyears + instadditionals$addinstlifeyears
            thsqalys <<- thsqalys + additionals$addqalys + instadditionals$addinstqalys
            thscosts <<- thscosts + additionals$addcosts + instadditionals$addinstcosts 
            
            simulation.clock <<- Inf # Set current time to Infinity so that patient level loop stops
            
        } 
        
    } else if (intervention=="ews"){ 
        
        # ews logic - follows same format as uc logic
        
        if (dbg.mode==T){
            print(data.frame(evtname=evt,evttime=curtime), row.names = FALSE)
        }
        
        if (evt=="outpat"){     
            if (ind.pat.data==T){
                this.PatData$ews$noutpat <<- this.PatData$ews$noutpat + 1 
            }      
            
            additionals     <- AddOngoing(lclprevtime=prevtime, lclcurtime=curtime, lclvally=1, lclvalq=utilmlt, lclvalc=parameters$cost.ews)
            instadditionals <- AddInstant(lclcurtime=curtime, lclvally=0, lclvalq=0, lclvalc=parameters$cost.outpat.ews)
            
            thslifeyears <<- thslifeyears + additionals$addlifeyears + instadditionals$addinstlifeyears
            thsqalys <<- thsqalys + additionals$addqalys + instadditionals$addinstqalys
            thscosts <<- thscosts + additionals$addcosts + instadditionals$addinstcosts 
            
            utilmlt <<- ifelse(curpat$nyha.class < 4, utilmlt * parameters$umult.outpat, utilmlt)
            
        } else if (evt=="hosp"){
            if (ind.pat.data==T){
                this.PatData$ews$nhosp <<- this.PatData$ews$nhosp + 1
            }
            
            additionals     <- AddOngoing(lclprevtime=prevtime, lclcurtime=curtime, lclvally=1, lclvalq=utilmlt, lclvalc=parameters$cost.ews)
            instadditionals <- AddInstant(lclcurtime=curtime, lclvally=0, lclvalq=0, lclvalc=parameters$cost.hosp)
            
            thslifeyears <<- thslifeyears + additionals$addlifeyears + instadditionals$addinstlifeyears
            thsqalys <<- thsqalys + additionals$addqalys + instadditionals$addinstqalys
            thscosts <<- thscosts + additionals$addcosts + instadditionals$addinstcosts 
            
            utilmlt <<- ifelse(curpat$nyha.class < 4, utilmlt * parameters$umult.hosp, utilmlt)
            
            patdies <- simulation.draws[simulation.index+2,patient.index] < p.dth.hosp   # If the random number between 0 and 1 is lower than the predicted probability of dying in hospital (p.dth.hosp), "patdies" is set to TRUE
            
            if (patdies==T) { 
                if (ind.pat.data==T){
                    this.PatData$ews$deathhosp <<- this.PatData$ews$deathhosp + 1
                    
                }
                
                simulation.clock <<- Inf
            }
            
        } else if (evt=="death"){     
            if (ind.pat.data==T){
                this.PatData$ews$deathother <<- this.PatData$ews$deathother + 1
            }
            
            additionals     <- AddOngoing(lclprevtime=prevtime, lclcurtime=curtime, lclvally=1, lclvalq=utilmlt, lclvalc=parameters$cost.ews)
            instadditionals <- AddInstant(lclcurtime=curtime, lclvally=0, lclvalq=0, lclvalc=parameters$cost.death)
            
            thslifeyears <<- thslifeyears + additionals$addlifeyears + instadditionals$addinstlifeyears
            thsqalys <<- thsqalys + additionals$addqalys + instadditionals$addinstqalys
            thscosts <<- thscosts + additionals$addcosts + instadditionals$addinstcosts 
            
            simulation.clock <<- Inf # Set current time to Infinity so patient level loop stops
            
        }
        
    } else if (intervention=="ewsda"){ 
        
        # ewsda logic - follows same format as ews logic, adding the possibility of avoiding hospitalisation and also accounting for false positive costs
        
        if (dbg.mode==T){
            print(data.frame(evtname=evt,evttime=curtime), row.names = FALSE)
        }
        
        if (evt=="outpat"){     
            if (ind.pat.data==T){
                this.PatData$ewsda$noutpat <<- this.PatData$ewsda$noutpat + 1 
            }      
            
            additionals     <- AddOngoing(lclprevtime=prevtime, lclcurtime=curtime, lclvally=1, lclvalq=utilmlt, lclvalc=parameters$cost.ews + 365 * parameters$cost.fpalarm * parameters$false.positive.rate)
            instadditionals <- AddInstant(lclcurtime=curtime, lclvally=0, lclvalq=0, lclvalc=parameters$cost.outpat.ews)
            
            thslifeyears <<- thslifeyears + additionals$addlifeyears + instadditionals$addinstlifeyears
            thsqalys <<- thsqalys + additionals$addqalys + instadditionals$addinstqalys
            thscosts <<- thscosts + additionals$addcosts + instadditionals$addinstcosts 
            
            utilmlt <<- ifelse(curpat$nyha.class < 4, utilmlt * parameters$umult.outpat, utilmlt)
            
        } else if (evt=="avoidhosp"){     
            if (ind.pat.data==T){
                this.PatData$ewsda$navoidhosp <<- this.PatData$ewsda$navoidhosp + 1 
            }      
            
            additionals     <- AddOngoing(lclprevtime=prevtime, lclcurtime=curtime, lclvally=1, lclvalq=utilmlt, lclvalc=parameters$cost.ews + 365 * parameters$cost.fpalarm * parameters$false.positive.rate)
            instadditionals <- AddInstant(lclcurtime=curtime, lclvally=0, lclvalq=0, lclvalc=parameters$cost.outpat.ews)
            
            thslifeyears <<- thslifeyears + additionals$addlifeyears + instadditionals$addinstlifeyears
            thsqalys <<- thsqalys + additionals$addqalys + instadditionals$addinstqalys
            thscosts <<- thscosts + additionals$addcosts + instadditionals$addinstcosts 
            
            utilmlt <<- ifelse(curpat$nyha.class < 4, utilmlt * parameters$umult.outpat, utilmlt)
            
            
        } else if (evt=="hosp"){
            if (ind.pat.data==T){
                this.PatData$ewsda$nhosp <<- this.PatData$ewsda$nhosp + 1
            }
            
            additionals     <- AddOngoing(lclprevtime=prevtime, lclcurtime=curtime, lclvally=1, lclvalq=utilmlt, lclvalc=parameters$cost.ews + 365 * parameters$cost.fpalarm * parameters$false.positive.rate)
            instadditionals <- AddInstant(lclcurtime=curtime, lclvally=0, lclvalq=0, lclvalc=parameters$cost.hosp)
            
            thslifeyears <<- thslifeyears + additionals$addlifeyears + instadditionals$addinstlifeyears
            thsqalys <<- thsqalys + additionals$addqalys + instadditionals$addinstqalys
            thscosts <<- thscosts + additionals$addcosts + instadditionals$addinstcosts 
            
            utilmlt <<- ifelse(curpat$nyha.class < 4, utilmlt * parameters$umult.hosp, utilmlt)
            
            patdies <- simulation.draws[simulation.index+2,patient.index] < p.dth.hosp   # If the random number between 0 and 1 is lower than the predicted probability of dying in hospital (p.dth.hosp), "patdies" is set to TRUE
            
            if (patdies==T) { 
                if (ind.pat.data==T){
                    this.PatData$ewsda$deathhosp <<- this.PatData$ewsda$deathhosp + 1
                    
                }
                
                simulation.clock <<- Inf
            }
            
        } else if (evt=="death"){     
            if (ind.pat.data==T){
                this.PatData$ewsda$deathother <<- this.PatData$ewsda$deathother + 1
            }
            
            additionals     <- AddOngoing(lclprevtime=prevtime, lclcurtime=curtime, lclvally=1, lclvalq=utilmlt, lclvalc=parameters$cost.ews + 365 * parameters$cost.fpalarm * parameters$false.positive.rate)
            instadditionals <- AddInstant(lclcurtime=curtime, lclvally=0, lclvalq=0, lclvalc=parameters$cost.death)
            
            thslifeyears <<- thslifeyears + additionals$addlifeyears + instadditionals$addinstlifeyears
            thsqalys <<- thsqalys + additionals$addqalys + instadditionals$addinstqalys
            thscosts <<- thscosts + additionals$addcosts + instadditionals$addinstcosts 
            
            simulation.clock <<- Inf # Set current time to Infinity so patient level loop stops
            
        }
        
        else {
            # if this code is reached, something has gone wrong
            stop("Event type not recognised")
        } 
    }
    # set to return nothing (NULL object) so earlier operations are not returned by accident
    return(NULL)
    
}

# UpdatePatientCharacteristics: update characteristics of the patient based on event type and elapsed time between previous and current event
# thisevt: a two element list containing the output from GetNextEvt
# $evt: event type 
# $updtime: elapsed time between the previous and current event (i.e. update time)

UpdatePatientCharacteristics <- function(thisevt, intervention){      # This function updates patient characteristics based on event, as identified in the GetNextEvt function
    
    evt      <- thisevt$evt                 # Identify event type
    updtime  <- thisevt$evttime             # Identify time elapsed between previous and current events
    
    if (intervention=="uc"){
        
        curpat$age <<- curpat$age + updtime
        
        curpat$previous.hosp <<- ifelse(evt == "hosp", 1, curpat$previous.hosp)
        
        if(evt == "hosp")curpat$hospitalisation.count <<-  curpat$hospitalisation.count + 1
        
        
        if (ind.pat.data==T){
            this.PatData$uc$pathist <<- rbind(this.PatData$uc$pathist, data.frame(age = round(curpat$age, 2), evtname=evt,evttime=round(curtime, 2)))
        } 
    }
    
    if (intervention=="ews"){
        
        curpat$age <<- curpat$age + updtime
        
        curpat$previous.hosp <<- ifelse(evt == "hosp", 1, curpat$previous.hosp)
        
        if(evt == "hosp")curpat$hospitalisation.count <<-  curpat$hospitalisation.count + 1
        
        
        if (ind.pat.data==T){
            this.PatData$ews$pathist <<- rbind(this.PatData$ews$pathist, data.frame(age = round(curpat$age, 2), evtname=evt,evttime=round(curtime, 2)))
        } 
    }
    
    if (intervention=="ewsda"){
        
        curpat$age <<- curpat$age + updtime
        
        curpat$previous.hosp <<- ifelse(evt == "hosp", 1, curpat$previous.hosp)
        
        if(evt == "hosp")curpat$hospitalisation.count <<-  curpat$hospitalisation.count + 1
        
        
        if (ind.pat.data==T){
            this.PatData$ewsda$pathist <<- rbind(this.PatData$ewsda$pathist, data.frame(age = round(curpat$age, 2), evtname=evt,evttime=round(curtime, 2)))
        } 
    }
    
    simulation.index <<- simulation.index + 3  # Used for updating random draws within same patient and across patients 
}

# RunSim: Enclose simulation within a function
RunSim <- function(deterministic, probabilistic){
    
    InstallAndLoad(required.packages)
    
    simulation.draws <<- data.frame(replicate(max(npats, npats.psa),runif(100000)))
    
    CreatePatients(max(npats, npats.psa), new.patients)
    
    pat.list <<- read.csv("./Data/patient_list.csv", sep = ",")
    
    sim.pats <<- read.csv("./Data/simulated_patients.csv", sep = ",")

    
    #### Deterministic analysis ####
    
    if (deterministic == T) {
    
    IndividualPatientData(ind.pat.data,npats)
    
    SurvivalModels(dist.death.uc, dist.death.ews, dist.hosp.uc, dist.hosp.ews)
    
    parameters <<- list(drly = drly,
                        drq = drq,
                        drc = drc,
                        dist.death.uc = dist.death.uc,
                        dist.death.ews = dist.death.ews,
                        dist.hosp.uc = dist.hosp.uc,
                        dist.hosp.ews = dist.hosp.ews,
                        sensitivity = sensitivity,
                        false.positive.rate = false.positive.rate,
                        avoid.prob = avoid.prob,
                        umult.outpat = umult.outpat,
                        umult.hosp = umult.hosp,
                        time.outpat.visit.uc = time.outpat.visit.uc,
                        time.outpat.visit.ews = time.outpat.visit.ews,
                        cost.outpat.uc = cost.outpat.uc,
                        cost.outpat.ews = cost.outpat.ews,
                        cost.hosp = cost.hosp,
                        cost.death = cost.death,
                        cost.uc = cost.uc,
                        cost.ews = cost.ews,
                        cost.fpalarm = cost.fpalarm
    )
    
    reg.coef.death.uc  <<- as.numeric(coef(surv.death.uc))
    reg.coef.death.ews <<- as.numeric(coef(surv.death.ews))
    reg.coef.hosp.uc   <<- as.numeric(coef(surv.hosp.uc))
    reg.coef.hosp.ews  <<- as.numeric(coef(surv.hosp.ews))
    
    reg.coef.dth.hosp  <<- as.numeric(coef(p.death.hosp))
    
    
    # initialise variable of total costs and QALYs (note use of superassignment operator (<<-))
    tot.lifeyears.uc <<- 0 # total life years accrued by all patients - UC
    tot.qalys.uc <<- 0 # total QALYs accrued by all patients - UC
    tot.costs.uc <<- 0 # total costs accrued by all patients - UC
    
    tot.lifeyears.ews <<- 0 # total life years accrued by all patients - EWS
    tot.qalys.ews <<- 0 # total QALYs accrued by all patients - EWS
    tot.costs.ews <<- 0 # total costs accrued by all patients - EWS
    
    tot.lifeyears.ewsda <<- 0 # total life years accrued by all patients - EWS+DA
    tot.qalys.ewsda <<- 0 # total QALYs accrued by all patients - EWS+DA
    tot.costs.ewsda <<- 0 # total costs accrued by all patients - EWS+DA
    
    tot.dlifeyears.uc.ews <<- 0 # difference in life years between UC and EWS
    tot.dqalys.uc.ews <<- 0 # difference in QALYs between UC and EWS
    tot.dcosts.uc.ews <<- 0 # difference in costs between UC and EWS
    
    tot.dlifeyears.uc.ewsda <<- 0 # difference in life years between UC and EWS+DA
    tot.dqalys.uc.ewsda <<- 0 # difference in QALYs between UC and EWS+DA
    tot.dcosts.uc.ewsda <<- 0 # difference in costs between UC and EWS+DA
    
    tot.dlifeyears.ews.ewsda <<- 0 # difference in life years between EWS and EWS+DA
    tot.dqalys.ews.ewsda <<- 0 # difference in QALYs between EWS and EWS+DA
    tot.dcosts.ews.ewsda <<- 0 # difference in costs between EWS and EWS+DA
    
    # Inner loop, repeat for each patient
    for (i in 1:npats){
        
        patient.index <<- i
        
        if (ind.pat.data==T){
            this.PatData <<- list(
                uc=list(
                    noutpat=0,
                    nhosp=0,
                    deathhosp=0,
                    deathother=0
                ),
                ews=list(
                    noutpat=0,
                    nhosp=0,
                    deathhosp=0,
                    deathother=0
                ),
                ewsda=list(
                    noutpat=0,
                    navoidhosp=0,
                    nhosp=0,
                    deathhosp=0,
                    deathother=0
                )
            )
        }
        
        # Debugging line
        if (dbg.mode==T){
            cat(paste("\n######\n[", i, "]\n######\n"))    # The cat function prints text to console  
            cat(paste("UC Patient\n\n"))
        }
        
        # For the uc patient:
        # Draw random patient from created patient list
        # Set current time, set as global variable to 0
        curpat <<- sim.pats[i,]
        
        # Define uc group 
        curpat$intervention <<- 0
        
        if (ind.pat.data==T){
            this.PatData$uc$patchar <<- data.frame(pat.id = curpat$Patient.Study.Number,
                                                   ejection.fraction = curpat$ejection.fraction,
                                                   age = curpat$age,
                                                   sbp = curpat$sbp,
                                                   bmi = curpat$bmi,
                                                   creatinine = curpat$creatinine,
                                                   nyha.class = curpat$nyha.class,
                                                   gender = curpat$gender,
                                                   smoker = curpat$smoker,
                                                   diabetes = curpat$diabetes,
                                                   copd = curpat$copd,
                                                   recent.diagnosis = curpat$recent.diagnosis,
                                                   no.beta.blocker = curpat$no.beta.blocker,
                                                   no.ace = curpat$no.ace,
                                                   age.ef = curpat$age.ef,
                                                   sbp.ef = curpat$sbp.ef,
                                                   myocardial.infarction = curpat$myocardial.infarction,
                                                   atrial.fibrillation = curpat$chronic.atrial.fibrillation,
                                                   hypertension = curpat$hypertension,
                                                   utility = curpat$utility)
        }
        
        if (dbg.mode==T){
            print(data.frame(pat.id = curpat$Patient.Study.Number,
                             ejection.fraction = curpat$ejection.fraction,
                             age = curpat$age,
                             sbp = curpat$sbp,
                             bmi = curpat$bmi,
                             creatinine = curpat$creatinine,
                             nyha.class = curpat$nyha.class,
                             gender = curpat$gender,
                             smoker = curpat$smoker,
                             diabetes = curpat$diabetes,
                             copd = curpat$copd,
                             recent.diagnosis = curpat$recent.diagnosis,
                             no.beta.blocker = curpat$no.beta.blocker,
                             no.ace = curpat$no.ace,
                             age.ef = curpat$age.ef,
                             sbp.ef = curpat$sbp.ef,
                             myocardial.infarction = curpat$myocardial.infarction,
                             atrial.fibrillation = curpat$chronic.atrial.fibrillation,
                             hypertension = curpat$hypertension,
                             utility = curpat$utility), row.names = FALSE)
        }
        
        simulation.clock <<- 0
        curtime <<- 0
        utilmlt <<- curpat$utility
        
        simulation.index <<- 1
        
        # Life years, QALYs, and costs for this patient
        thslifeyears <<- 0
        thsqalys <<- 0
        thscosts <<- 0
        
        while(simulation.clock < Inf){
            
            # Identify next event
            Evt <- GetNextEvent(intervention="uc")
            
            # Process event
            ReactEvt(Evt, intervention="uc")
            
            # Update patient characteristics
            UpdatePatientCharacteristics(Evt, intervention="uc")
            
            if(dbg.mode==T){        
                print(paste("Life years:", round(thslifeyears, 2), "; Qalys:", round(thsqalys, 2), "; Costs:", round(thscosts,0)))
            }
            if (ind.pat.data==T){
                this.PatData$uc$thslifeyears <<- thslifeyears
                this.PatData$uc$thsqalys <<- thsqalys
                this.PatData$uc$thscosts <<- thscosts        
            }
        }
        
        # Record outcomes
        tot.lifeyears.uc <<- tot.lifeyears.uc + thslifeyears
        tot.qalys.uc <<- tot.qalys.uc + thsqalys
        tot.costs.uc <<- tot.costs.uc + thscosts
        
        tot.dlifeyears.uc.ews <<- tot.dlifeyears.uc.ews - thslifeyears 
        tot.dqalys.uc.ews <<- tot.dqalys.uc.ews - thsqalys 
        tot.dcosts.uc.ews <<- tot.dcosts.uc.ews - thscosts
        
        tot.dlifeyears.uc.ewsda <<- tot.dlifeyears.uc.ewsda - thslifeyears 
        tot.dqalys.uc.ewsda <<- tot.dqalys.uc.ewsda - thsqalys 
        tot.dcosts.uc.ewsda <<- tot.dcosts.uc.ewsda - thscosts
        
    # For the ews patient:
    if (dbg.mode==T){
        cat(paste("\nEWS Patient\n\n"))}
    
    # Reset patient characteristics, current time, and initial utility
    curpat <<- sim.pats[i,]
    
    # Define ews group 
    curpat$intervention <<- 1
        
        if (ind.pat.data==T){
            this.PatData$ews$patchar <<- data.frame(pat.id = curpat$Patient.Study.Number,
                                                    ejection.fraction = curpat$ejection.fraction,
                                                    age = curpat$age,
                                                    sbp = curpat$sbp,
                                                    bmi = curpat$bmi,
                                                    creatinine = curpat$creatinine,
                                                    nyha.class = curpat$nyha.class,
                                                    gender = curpat$gender,
                                                    smoker = curpat$smoker,
                                                    diabetes = curpat$diabetes,
                                                    copd = curpat$copd,
                                                    recent.diagnosis = curpat$recent.diagnosis,
                                                    no.beta.blocker = curpat$no.beta.blocker,
                                                    no.ace = curpat$no.ace,
                                                    age.ef = curpat$age.ef,
                                                    sbp.ef = curpat$sbp.ef,
                                                    myocardial.infarction = curpat$myocardial.infarction,
                                                    atrial.fibrillation = curpat$chronic.atrial.fibrillation,
                                                    hypertension = curpat$hypertension,
                                                    utility = curpat$utility)
        }
        
        if (dbg.mode==T){
            print(data.frame(pat.id = curpat$Patient.Study.Number,
                             ejection.fraction = curpat$ejection.fraction,
                             age = curpat$age,
                             sbp = curpat$sbp,
                             bmi = curpat$bmi,
                             creatinine = curpat$creatinine,
                             nyha.class = curpat$nyha.class,
                             gender = curpat$gender,
                             smoker = curpat$smoker,
                             diabetes = curpat$diabetes,
                             copd = curpat$copd,
                             recent.diagnosis = curpat$recent.diagnosis,
                             no.beta.blocker = curpat$no.beta.blocker,
                             no.ace = curpat$no.ace,
                             age.ef = curpat$age.ef,
                             sbp.ef = curpat$sbp.ef,
                             myocardial.infarction = curpat$myocardial.infarction,
                             atrial.fibrillation = curpat$chronic.atrial.fibrillation,
                             hypertension = curpat$hypertension,
                             utility = curpat$utility), row.names = FALSE)
        }
        
        simulation.clock <<- 0
        curtime <<- 0
        utilmlt <<- curpat$utility
        
        simulation.index <<- 1
        
        # Life years, QALYs, and costs for this patient
        thslifeyears <<- 0
        thsqalys <<- 0
        thscosts <<- 0
        
        while(simulation.clock < Inf){
            
            # Identify next event
            Evt <- GetNextEvent(intervention="ews")
            
            # Process event
            ReactEvt(Evt, intervention="ews")
            
            # Update patient characteristics
            UpdatePatientCharacteristics(Evt, intervention="ews")
            
            if(dbg.mode==T){        
                print(paste("Life years:", round(thslifeyears, 2), "; Qalys:", round(thsqalys, 2), "; Costs:", round(thscosts,0)))            }
            if (ind.pat.data==T){
                this.PatData$ews$thslifeyears <<- thslifeyears
                this.PatData$ews$thsqalys <<- thsqalys
                this.PatData$ews$thscosts <<- thscosts        
            }
            
        } 
        
        if (ind.pat.data==T){
            PatData[[i]] <<- this.PatData
        }
        
        # Record outcomes
        tot.lifeyears.ews <<- tot.lifeyears.ews + thslifeyears
        tot.qalys.ews <<- tot.qalys.ews + thsqalys
        tot.costs.ews <<- tot.costs.ews + thscosts
        
        tot.dlifeyears.uc.ews <<- tot.dlifeyears.uc.ews + thslifeyears
        tot.dqalys.uc.ews <<- tot.dqalys.uc.ews + thsqalys
        tot.dcosts.uc.ews <<- tot.dcosts.uc.ews + thscosts
        
        tot.dlifeyears.ews.ewsda <<- tot.dlifeyears.ews.ewsda - thslifeyears 
        tot.dqalys.ews.ewsda <<- tot.dqalys.ews.ewsda - thsqalys 
        tot.dcosts.ews.ewsda <<- tot.dcosts.ews.ewsda - thscosts
    
    # For the ewsda patient:
    if (dbg.mode==T){
        cat(paste("\nEWS+DA Patient\n\n"))}
    
    # Reset patient characteristics, current time, and initial utility
    curpat <<- sim.pats[i,]
    
    # Define ewsda group (same as ews group)
    curpat$intervention <<- 1
        
        if (ind.pat.data==T){
            this.PatData$ewsda$patchar <<- data.frame(pat.id = curpat$Patient.Study.Number,
                                                      ejection.fraction = curpat$ejection.fraction,
                                                      age = curpat$age,
                                                      sbp = curpat$sbp,
                                                      bmi = curpat$bmi,
                                                      creatinine = curpat$creatinine,
                                                      nyha.class = curpat$nyha.class,
                                                      gender = curpat$gender,
                                                      smoker = curpat$smoker,
                                                      diabetes = curpat$diabetes,
                                                      copd = curpat$copd,
                                                      recent.diagnosis = curpat$recent.diagnosis,
                                                      no.beta.blocker = curpat$no.beta.blocker,
                                                      no.ace = curpat$no.ace,
                                                      age.ef = curpat$age.ef,
                                                      sbp.ef = curpat$sbp.ef,
                                                      myocardial.infarction = curpat$myocardial.infarction,
                                                      atrial.fibrillation = curpat$chronic.atrial.fibrillation,
                                                      hypertension = curpat$hypertension,
                                                      utility = curpat$utility)
        }
        
        if (dbg.mode==T){
            print(data.frame(pat.id = curpat$Patient.Study.Number,
                             ejection.fraction = curpat$ejection.fraction,
                             age = curpat$age,
                             sbp = curpat$sbp,
                             bmi = curpat$bmi,
                             creatinine = curpat$creatinine,
                             nyha.class = curpat$nyha.class,
                             gender = curpat$gender,
                             smoker = curpat$smoker,
                             diabetes = curpat$diabetes,
                             copd = curpat$copd,
                             recent.diagnosis = curpat$recent.diagnosis,
                             no.beta.blocker = curpat$no.beta.blocker,
                             no.ace = curpat$no.ace,
                             age.ef = curpat$age.ef,
                             sbp.ef = curpat$sbp.ef,
                             myocardial.infarction = curpat$myocardial.infarction,
                             atrial.fibrillation = curpat$chronic.atrial.fibrillation,
                             hypertension = curpat$hypertension,
                             utility = curpat$utility), row.names = FALSE)
        }
        
        simulation.clock <<- 0
        curtime <<- 0
        utilmlt <<- curpat$utility
        
        simulation.index <<- 1
        
        # QALYs and costs for this patient
        thslifeyears <<- 0
        thsqalys <<- 0
        thscosts <<- 0
        
        while(simulation.clock < Inf){
            
            # Identify next event
            Evt <- GetNextEvent(intervention="ewsda")
            
            # Process event
            ReactEvt(Evt, intervention="ewsda")
            
            # Update patient characteristics
            UpdatePatientCharacteristics(Evt, intervention="ewsda")
            
            if(dbg.mode==T){        
                print(paste("Life years:", round(thslifeyears, 2), "; Qalys:", round(thsqalys, 2), "; Costs:", round(thscosts,0)))            }
            if (ind.pat.data==T){
                this.PatData$ewsda$thslifeyears <<- thslifeyears
                this.PatData$ewsda$thsqalys <<- thsqalys
                this.PatData$ewsda$thscosts <<- thscosts        
            }
            
        } 
        
        if (ind.pat.data==T){
            PatData[[i]] <<- this.PatData
        }
        
        
        # Record outcomes
        tot.lifeyears.ewsda <<- tot.lifeyears.ewsda + thslifeyears
        tot.qalys.ewsda <<- tot.qalys.ewsda + thsqalys
        tot.costs.ewsda <<- tot.costs.ewsda + thscosts
        
        tot.dlifeyears.uc.ewsda <<- tot.dlifeyears.uc.ewsda + thslifeyears
        tot.dqalys.uc.ewsda <<- tot.dqalys.uc.ewsda + thsqalys
        tot.dcosts.uc.ewsda <<- tot.dcosts.uc.ewsda + thscosts
        
        tot.dlifeyears.ews.ewsda <<- tot.dlifeyears.ews.ewsda + thslifeyears 
        tot.dqalys.ews.ewsda <<- tot.dqalys.ews.ewsda + thsqalys 
        tot.dcosts.ews.ewsda <<- tot.dcosts.ews.ewsda + thscosts
        
    }
    
    if (export.excel==T) {
    
    individual.patient.data <- data.frame(
        
        pat.id=rep(NA, npats),
        ejection.fraction=rep(NA, npats),
        age=rep(NA, npats),
        sbp=rep(NA, npats),
        bmi=rep(NA, npats),
        creatinine=rep(NA, npats),
        nyha.class=rep(NA, npats),
        gender=rep(NA, npats),
        smoker=rep(NA, npats),
        diabetes=rep(NA, npats),
        copd=rep(NA, npats),
        recent.diagnosis=rep(NA, npats),
        no.beta.blocker=rep(NA, npats),
        no.ace=rep(NA, npats),
        age.ef=rep(NA, npats),
        sbp.ef=rep(NA, npats),
        myocardial.infarction=rep(NA, npats),
        atrial.fibrillation=rep(NA, npats),
        hypertension=rep(NA, npats),
        utility=rep(NA, npats),
        age.death.uc=rep(NA, npats),
        age.death.ews=rep(NA, npats),
        age.death.ewsda=rep(NA, npats),
        
        uc.outpatcount=rep(NA, npats),
        uc.hospcount=rep(NA, npats),
        uc.deathhosp=rep(NA,npats),
        uc.deathother=rep(NA, npats),
        uc.costs=rep(NA, npats),
        uc.lifeyears=rep(NA, npats),
        uc.qalys=rep(NA, npats),
        
        ews.outpatcount=rep(NA, npats),
        ews.hospcount=rep(NA, npats),
        ews.deathhosp=rep(NA,npats),
        ews.deathother=rep(NA, npats),
        ews.costs=rep(NA, npats),
        ews.lifeyears=rep(NA, npats),
        ews.qalys=rep(NA, npats),
        
        ewsda.outpatcount=rep(NA, npats),
        ewsda.hospcount=rep(NA, npats),
        ewsda.avoidhospcount=rep(NA, npats),
        ewsda.deathhosp=rep(NA,npats),
        ewsda.deathother=rep(NA, npats),
        ewsda.costs=rep(NA, npats),
        ewsda.lifeyears=rep(NA, npats),
        ewsda.qalys=rep(NA, npats)
    )
    
    # Fill in with data generated from the simulation
    for (i in 1:npats){
        
        individual.patient.data[i, "pat.id"]                  <- PatData[[i]]$uc$patchar$pat.id
        individual.patient.data[i, "ejection.fraction"]       <- PatData[[i]]$uc$patchar$ejection.fraction
        individual.patient.data[i, "age"]                     <- PatData[[i]]$uc$patchar$age
        individual.patient.data[i, "sbp"]                     <- PatData[[i]]$uc$patchar$sbp
        individual.patient.data[i, "bmi"]                     <- PatData[[i]]$uc$patchar$bmi
        individual.patient.data[i, "creatinine"]              <- PatData[[i]]$uc$patchar$creatinine
        individual.patient.data[i, "nyha.class"]              <- PatData[[i]]$uc$patchar$nyha.class
        individual.patient.data[i, "gender"]                  <- PatData[[i]]$uc$patchar$gender
        individual.patient.data[i, "smoker"]                  <- PatData[[i]]$uc$patchar$smoker
        individual.patient.data[i, "diabetes"]                <- PatData[[i]]$uc$patchar$diabetes
        individual.patient.data[i, "copd"]                    <- PatData[[i]]$uc$patchar$copd
        individual.patient.data[i, "recent.diagnosis"]        <- PatData[[i]]$uc$patchar$recent.diagnosis
        individual.patient.data[i, "no.beta.blocker"]         <- PatData[[i]]$uc$patchar$no.beta.blocker
        individual.patient.data[i, "no.ace"]                  <- PatData[[i]]$uc$patchar$no.ace
        individual.patient.data[i, "age.ef"]                  <- PatData[[i]]$uc$patchar$age.ef
        individual.patient.data[i, "sbp.ef"]                  <- PatData[[i]]$uc$patchar$sbp.ef
        individual.patient.data[i, "myocardial.infarction"]   <- PatData[[i]]$uc$patchar$myocardial.infarction
        individual.patient.data[i, "atrial.fibrillation"]     <- PatData[[i]]$uc$patchar$atrial.fibrillation
        individual.patient.data[i, "hypertension"]            <- PatData[[i]]$uc$patchar$hypertension
        individual.patient.data[i, "utility"]                 <- PatData[[i]]$uc$patchar$utility
        individual.patient.data[i, "age.death.uc"]            <- PatData[[i]]$uc$patchar$age+PatData[[i]]$uc$thslifeyears
        individual.patient.data[i, "age.death.ews"]           <- PatData[[i]]$uc$patchar$age+PatData[[i]]$ews$thslifeyears
        individual.patient.data[i, "age.death.ewsda"]         <- PatData[[i]]$uc$patchar$age+PatData[[i]]$ewsda$thslifeyears
        
        individual.patient.data[i, "uc.outpatcount"] <- PatData[[i]]$uc$noutpat
        individual.patient.data[i, "uc.hospcount"]   <- PatData[[i]]$uc$nhosp
        individual.patient.data[i, "uc.deathhosp"]   <- PatData[[i]]$uc$deathhosp
        individual.patient.data[i, "uc.deathother"]  <- PatData[[i]]$uc$deathother
        individual.patient.data[i, "uc.costs"]       <- PatData[[i]]$uc$thscosts
        individual.patient.data[i, "uc.lifeyears"]   <- PatData[[i]]$uc$thslifeyears
        individual.patient.data[i, "uc.qalys"]       <- PatData[[i]]$uc$thsqalys
        
        individual.patient.data[i, "ews.outpatcount"] <- PatData[[i]]$ews$noutpat
        individual.patient.data[i, "ews.hospcount"]   <- PatData[[i]]$ews$nhosp
        individual.patient.data[i, "ews.deathhosp"]   <- PatData[[i]]$ews$deathhosp
        individual.patient.data[i, "ews.deathother"]  <- PatData[[i]]$ews$deathother
        individual.patient.data[i, "ews.costs"]       <- PatData[[i]]$ews$thscosts
        individual.patient.data[i, "ews.lifeyears"]   <- PatData[[i]]$ews$thslifeyears
        individual.patient.data[i, "ews.qalys"]       <- PatData[[i]]$ews$thsqalys
        
        individual.patient.data[i, "ewsda.outpatcount"]     <- PatData[[i]]$ewsda$noutpat
        individual.patient.data[i, "ewsda.hospcount"]       <- PatData[[i]]$ewsda$nhosp
        individual.patient.data[i, "ewsda.avoidhospcount"]  <- PatData[[i]]$ewsda$navoidhosp
        individual.patient.data[i, "ewsda.deathhosp"]       <- PatData[[i]]$ewsda$deathhosp
        individual.patient.data[i, "ewsda.deathother"]      <- PatData[[i]]$ewsda$deathother
        individual.patient.data[i, "ewsda.costs"]           <- PatData[[i]]$ewsda$thscosts
        individual.patient.data[i, "ewsda.lifeyears"]       <- PatData[[i]]$ewsda$thslifeyears
        individual.patient.data[i, "ewsda.qalys"]           <- PatData[[i]]$ewsda$thsqalys
    }
   
    # Calculate total and average for simulated populations
    individual.patient.data <- rbind(individual.patient.data, c(NA, colSums(individual.patient.data[,-1])/npats))
    
    # Create patient number column and average
    individual.patient.data <- data.frame(patient.number=c(1:npats, "Average"), individual.patient.data)
    
    # Create.csv file
    write.csv(individual.patient.data, file = "./Results/Individual_patient_data.csv")
    
    # Create deterministic table of results
    deterministic.results <- data.frame(
        
        Comparison = c("ewsda vs ews", "ewsda vs uc", "ews vs uc"),
        difference.costs = c(tot.dcosts.ews.ewsda, tot.dcosts.uc.ewsda, tot.dcosts.uc.ews),
        difference.lifeyears = c(tot.dlifeyears.ews.ewsda, tot.dlifeyears.uc.ewsda, tot.dlifeyears.uc.ews),
        difference.qalys = c(tot.dqalys.ews.ewsda, tot.dqalys.uc.ewsda, tot.dqalys.uc.ews),
        ICER.euro.lifeyear = c(tot.dcosts.ews.ewsda / tot.dlifeyears.ews.ewsda, tot.dcosts.uc.ewsda / tot.dlifeyears.uc.ewsda, tot.dcosts.uc.ews / tot.dlifeyears.uc.ews),
        ICER.euro.QALY = c(tot.dcosts.ews.ewsda / tot.dqalys.ews.ewsda, tot.dcosts.uc.ewsda / tot.dqalys.uc.ewsda, tot.dcosts.uc.ews / tot.dqalys.uc.ews)
        
    )
    
    write.csv(deterministic.results, file = "./Results/Deterministic_results.csv")
    
    }
    }
    
    
    #### Probabilistic sensitivity analysis ####
    
    if (probabilistic == T){
        
        ind.pat.data <<- F
        dbg.mode     <<- F
        
        SurvivalModels(dist.death.uc, dist.death.ews, dist.hosp.uc, dist.hosp.ews)
        
        PSA.sim <<- vector("list", length = nloops.psa)
        
        parameters.psa <- list(drly = rep(drly, nloops.psa),
                               drq = rep(drq, nloops.psa),
                               drc = rep(drc, nloops.psa),
                               dist.death.uc = rep(dist.death.uc,nloops.psa),
                               dist.death.ews = rep(dist.death.ews,nloops.psa),
                               dist.hosp.uc = rep(dist.hosp.uc,nloops.psa),
                               dist.hosp.ews = rep(dist.hosp.ews,nloops.psa),
                               sensitivity = rep(sensitivity, nloops.psa),
                               false.positive.rate = rep(false.positive.rate, nloops.psa),
                               avoid.prob = pmax(rep(0,nloops.psa), rnorm(nloops.psa, mean = avoid.prob, sd = se.avoid.prob)),
                               umult.outpat = pmax(rep(0,nloops.psa), rnorm(nloops.psa, mean = umult.outpat, sd = se.umult.outpat)),
                               umult.hosp = pmax(rep(0,nloops.psa), rnorm(nloops.psa, mean = umult.hosp, sd = se.umult.hosp)),
                               time.outpat.visit.uc = pmax(rep(0,nloops.psa), rnorm(nloops.psa, mean = time.outpat.visit.uc, sd = se.time.outpat.visit.uc)),
                               time.outpat.visit.ews = pmax(rep(0,nloops.psa), rnorm(nloops.psa, mean = time.outpat.visit.ews, sd = se.time.outpat.visit.ews)),
                               cost.outpat.uc = rgamma(nloops.psa, shape = (cost.outpat.uc/se.cost.outpat.uc)^2 , scale = se.cost.outpat.uc^2/cost.outpat.uc),
                               cost.outpat.ews = rgamma(nloops.psa, shape = (cost.outpat.ews/se.cost.outpat.ews)^2 , scale = se.cost.outpat.ews^2/cost.outpat.ews),
                               cost.hosp = rgamma(nloops.psa, shape = (cost.hosp/se.cost.hosp)^2 , scale = se.cost.hosp^2/cost.hosp),
                               cost.death = rgamma(nloops.psa, shape = (cost.death/se.cost.death)^2 , scale = se.cost.death^2/cost.death),
                               cost.uc = rgamma(nloops.psa, shape = (cost.uc/se.cost.uc)^2 , scale = se.cost.uc^2/cost.uc),
                               cost.ews = rgamma(nloops.psa, shape = (cost.ews/se.cost.ews)^2 , scale = se.cost.ews^2/cost.ews),
                               cost.fpalarm = rgamma(nloops.psa, shape = (cost.fpalarm/se.cost.fpalarm)^2 , scale = se.cost.fpalarm^2/cost.fpalarm)
        )
        
        # Outer loop, repeat for each PSA simulation
        
        for (j in 1:nloops.psa) {
            
            parameters <<- lapply(parameters.psa, function(x) x[[j]])
            
            reg.coef.death.uc  <<- as.numeric(mvtnorm::rmvnorm(1,coef(surv.death.uc),vcov(surv.death.uc)))
            reg.coef.death.ews <<- as.numeric(mvtnorm::rmvnorm(1,coef(surv.death.ews),vcov(surv.death.ews)))
            reg.coef.hosp.uc   <<- as.numeric(mvtnorm::rmvnorm(1,coef(surv.hosp.uc),vcov(surv.hosp.uc)))
            reg.coef.hosp.ews  <<- as.numeric(mvtnorm::rmvnorm(1,coef(surv.hosp.ews),vcov(surv.hosp.ews)))
            
            reg.coef.dth.hosp  <<- as.numeric(mvtnorm::rmvnorm(1,coef(p.death.hosp),vcov(p.death.hosp)))
            
            set.seed(9)
            
            this.PSA.sim <<- list(dcosts.ews.ewsda=0,
                                   dcosts.uc.ewsda=0,
                                   dcosts.uc.ews=0,
                                   dlifeyears.ews.ewsda=0,
                                   dlifeyears.uc.ewsda=0,
                                   dlifeyears.uc.ews=0,
                                   dqalys.ews.ewsda=0,
                                   dqalys.uc.ewsda=0,
                                   dqalys.uc.ews=0,
                                   ICER.euro.ly.ews.ewsda=0,
                                   ICER.euro.ly.uc.ewsda=0,
                                   ICER.euro.ly.uc.ews=0,
                                   ICER.euro.qaly.ews.ewsda=0,
                                   ICER.euro.qaly.uc.ewsda=0,
                                   ICER.euro.qaly.uc.ews=0)
        
        # Trace PSA
        cat(paste("\n## PSA loop [", j, "] ##"))

        # initialise variable of total costs and QALYs (note use of superassignment operator (<<-))
        tot.lifeyears.uc <<- 0 # total life years accrued by all patients - UC
        tot.qalys.uc <<- 0 # total QALYs accrued by all patients - UC
        tot.costs.uc <<- 0 # total costs accrued by all patients - UC
        
        tot.lifeyears.ews <<- 0 # total life years accrued by all patients - EWS
        tot.qalys.ews <<- 0 # total QALYs accrued by all patients - EWS
        tot.costs.ews <<- 0 # total costs accrued by all patients - EWS
        
        tot.lifeyears.ewsda <<- 0 # total life years accrued by all patients - EWS+DA
        tot.qalys.ewsda <<- 0 # total QALYs accrued by all patients - EWS+DA
        tot.costs.ewsda <<- 0 # total costs accrued by all patients - EWS+DA
        
        tot.dlifeyears.uc.ews <<- 0 # difference in life years between UC and EWS
        tot.dqalys.uc.ews <<- 0 # difference in QALYs between UC and EWS
        tot.dcosts.uc.ews <<- 0 # difference in costs between UC and EWS
        
        tot.dlifeyears.uc.ewsda <<- 0 # difference in life years between UC and EWS+DA
        tot.dqalys.uc.ewsda <<- 0 # difference in QALYs between UC and EWS+DA
        tot.dcosts.uc.ewsda <<- 0 # difference in costs between UC and EWS+DA
        
        tot.dlifeyears.ews.ewsda <<- 0 # difference in life years between EWS and EWS+DA
        tot.dqalys.ews.ewsda <<- 0 # difference in QALYs between EWS and EWS+DA
        tot.dcosts.ews.ewsda <<- 0 # difference in costs between EWS and EWS+DA
        
        # Inner loop, repeat for each patient
        for (i in 1:npats.psa){
            
            patient.index <<- i
            
           # Trace PSA
            cat(paste("\nPatient ", i, ""))

            # For the uc patient:
            # Draw random patient from created patient list
            # Set current time, set as global variable to 0
            curpat <<- sim.pats[i,]
            
            # Define uc group 
            curpat$intervention <<- 0
            
            simulation.clock <<- 0
            curtime <<- 0
            utilmlt <<- curpat$utility
            
            simulation.index <<- 1
            
            # Life years, QALYs, and costs for this patient
            thslifeyears <<- 0
            thsqalys <<- 0
            thscosts <<- 0
            
            while(simulation.clock < Inf){
                
                # Identify next event
                Evt <- GetNextEvent(intervention="uc")
                
                # Process event
                ReactEvt(Evt, intervention="uc")
                
                # Update patient characteristics
                UpdatePatientCharacteristics(Evt, intervention="uc")

            }
            
            # Record outcomes
            tot.lifeyears.uc <<- tot.lifeyears.uc + thslifeyears
            tot.qalys.uc <<- tot.qalys.uc + thsqalys
            tot.costs.uc <<- tot.costs.uc + thscosts
            
            tot.dlifeyears.uc.ews <<- tot.dlifeyears.uc.ews - thslifeyears 
            tot.dqalys.uc.ews <<- tot.dqalys.uc.ews - thsqalys 
            tot.dcosts.uc.ews <<- tot.dcosts.uc.ews - thscosts
            
            tot.dlifeyears.uc.ewsda <<- tot.dlifeyears.uc.ewsda - thslifeyears 
            tot.dqalys.uc.ewsda <<- tot.dqalys.uc.ewsda - thsqalys 
            tot.dcosts.uc.ewsda <<- tot.dcosts.uc.ewsda - thscosts
            
            # For the ews patient:
            # Reset patient characteristics, current time, and initial utility
            curpat <<- sim.pats[i,]
            
            # Define ews group 
            curpat$intervention <<- 1
            
            simulation.clock <<- 0
            curtime <<- 0
            utilmlt <<- curpat$utility
            
            simulation.index <<- 1
            
            # Life years, QALYs, and costs for this patient
            thslifeyears <<- 0
            thsqalys <<- 0
            thscosts <<- 0
            
            while(simulation.clock < Inf){
                
                # Identify next event
                Evt <- GetNextEvent(intervention="ews")
                
                # Process event
                ReactEvt(Evt, intervention="ews")
                
                # Update patient characteristics
                UpdatePatientCharacteristics(Evt, intervention="ews")
                
            } 
            
            # Record outcomes
            tot.lifeyears.ews <<- tot.lifeyears.ews + thslifeyears
            tot.qalys.ews <<- tot.qalys.ews + thsqalys
            tot.costs.ews <<- tot.costs.ews + thscosts
            
            tot.dlifeyears.uc.ews <<- tot.dlifeyears.uc.ews + thslifeyears
            tot.dqalys.uc.ews <<- tot.dqalys.uc.ews + thsqalys
            tot.dcosts.uc.ews <<- tot.dcosts.uc.ews + thscosts
            
            tot.dlifeyears.ews.ewsda <<- tot.dlifeyears.ews.ewsda - thslifeyears 
            tot.dqalys.ews.ewsda <<- tot.dqalys.ews.ewsda - thsqalys 
            tot.dcosts.ews.ewsda <<- tot.dcosts.ews.ewsda - thscosts
            
            # For the ewsda patient:
            # Reset patient characteristics, current time, and initial utility
            curpat <<- sim.pats[i,]
            
            # Define ewsda group (same as ews group)
            curpat$intervention <<- 1
            
            simulation.clock <<- 0
            curtime <<- 0
            utilmlt <<- curpat$utility
            
            simulation.index <<- 1
            
            # QALYs and costs for this patient
            thslifeyears <<- 0
            thsqalys <<- 0
            thscosts <<- 0
            
            while(simulation.clock < Inf){
                
                # Identify next event
                Evt <- GetNextEvent(intervention="ewsda")
                
                # Process event
                ReactEvt(Evt, intervention="ewsda")
                
                # Update patient characteristics
                UpdatePatientCharacteristics(Evt, intervention="ewsda")

            } 

            # Record outcomes
            tot.lifeyears.ewsda <<- tot.lifeyears.ewsda + thslifeyears
            tot.qalys.ewsda <<- tot.qalys.ewsda + thsqalys
            tot.costs.ewsda <<- tot.costs.ewsda + thscosts
            
            tot.dlifeyears.uc.ewsda <<- tot.dlifeyears.uc.ewsda + thslifeyears
            tot.dqalys.uc.ewsda <<- tot.dqalys.uc.ewsda + thsqalys
            tot.dcosts.uc.ewsda <<- tot.dcosts.uc.ewsda + thscosts
            
            tot.dlifeyears.ews.ewsda <<- tot.dlifeyears.ews.ewsda + thslifeyears 
            tot.dqalys.ews.ewsda <<- tot.dqalys.ews.ewsda + thsqalys 
            tot.dcosts.ews.ewsda <<- tot.dcosts.ews.ewsda + thscosts
         
               
        }
        
        
        this.PSA.sim$tot.lifeyears.uc <<- tot.lifeyears.uc
        this.PSA.sim$tot.qalys.uc <<- tot.qalys.uc
        this.PSA.sim$tot.costs.uc <<- tot.costs.uc
        this.PSA.sim$tot.lifeyears.ews <<- tot.lifeyears.ews
        this.PSA.sim$tot.qalys.ews <<- tot.qalys.ews
        this.PSA.sim$tot.costs.ews <<- tot.costs.ews
        this.PSA.sim$tot.lifeyears.ewsda <<- tot.lifeyears.ewsda
        this.PSA.sim$tot.qalys.ewsda <<- tot.qalys.ewsda
        this.PSA.sim$tot.costs.ewsda <<- tot.costs.ewsda
        this.PSA.sim$ICER.euro.ly.ews.ewsda <<- tot.dcosts.ews.ewsda / tot.dlifeyears.ews.ewsda
        this.PSA.sim$ICER.euro.ly.uc.ewsda <<- tot.dcosts.uc.ewsda / tot.dlifeyears.uc.ewsda
        this.PSA.sim$ICER.euro.ly.uc.ews <<- tot.dcosts.uc.ews / tot.dlifeyears.uc.ews
        this.PSA.sim$ICER.euro.qaly.ews.ewsda <<- tot.dcosts.ews.ewsda / tot.dqalys.ews.ewsda
        this.PSA.sim$ICER.euro.qaly.uc.ewsda <<- tot.dcosts.uc.ewsda / tot.dqalys.uc.ewsda
        this.PSA.sim$ICER.euro.qaly.uc.ews <<- tot.dcosts.uc.ews / tot.dqalys.uc.ews
        
        PSA.sim[[j]] <<- this.PSA.sim
        
        }
        
        
        if (export.excel==T) {
        
        psa.results <- data.frame(
            
            tot.lifeyears.uc=rep(NA, nloops.psa),
            tot.qalys.uc=rep(NA, nloops.psa),
            tot.costs.uc=rep(NA, nloops.psa),
            tot.lifeyears.ews=rep(NA, nloops.psa),
            tot.qalys.ews=rep(NA, nloops.psa),
            tot.costs.ews=rep(NA, nloops.psa),
            tot.lifeyears.ewsda=rep(NA, nloops.psa),
            tot.qalys.ewsda=rep(NA, nloops.psa),
            tot.costs.ewsda=rep(NA, nloops.psa),
            ICER.euro.ly.ews.ewsda=rep(NA, nloops.psa),
            ICER.euro.ly.uc.ewsda=rep(NA, nloops.psa),
            ICER.euro.ly.uc.ews=rep(NA, nloops.psa),
            ICER.euro.qaly.ews.ewsda=rep(NA, nloops.psa),
            ICER.euro.qaly.uc.ewsda=rep(NA, nloops.psa),
            ICER.euro.qaly.uc.ews=rep(NA, nloops.psa)
            )
        
        # Fill in with data generated from the simulation
        for (j in 1:nloops.psa){
            
            psa.results[j, "tot.lifeyears.uc"]       <- PSA.sim[[j]]$tot.lifeyears.uc
            psa.results[j, "tot.qalys.uc"]       <- PSA.sim[[j]]$tot.qalys.uc
            psa.results[j, "tot.costs.uc"]       <- PSA.sim[[j]]$tot.costs.uc
            psa.results[j, "tot.lifeyears.ews"]       <- PSA.sim[[j]]$tot.lifeyears.ews
            psa.results[j, "tot.qalys.ews"]       <- PSA.sim[[j]]$tot.qalys.ews
            psa.results[j, "tot.costs.ews"]       <- PSA.sim[[j]]$tot.costs.ews
            psa.results[j, "tot.lifeyears.ewsda"]       <- PSA.sim[[j]]$tot.lifeyears.ewsda
            psa.results[j, "tot.qalys.ewsda"]       <- PSA.sim[[j]]$tot.qalys.ewsda
            psa.results[j, "tot.costs.ewsda"]       <- PSA.sim[[j]]$tot.costs.ewsda
            psa.results[j, "ICER.euro.ly.ews.ewsda"]       <- PSA.sim[[j]]$ICER.euro.ly.ews.ewsda
            psa.results[j, "ICER.euro.ly.uc.ewsda"]       <- PSA.sim[[j]]$ICER.euro.ly.uc.ewsda
            psa.results[j, "ICER.euro.ly.uc.ews"]       <- PSA.sim[[j]]$ICER.euro.ly.uc.ews
            psa.results[j, "ICER.euro.qaly.ews.ewsda"]       <- PSA.sim[[j]]$ICER.euro.qaly.ews.ewsda
            psa.results[j, "ICER.euro.qaly.uc.ewsda"]       <- PSA.sim[[j]]$ICER.euro.qaly.uc.ewsda
            psa.results[j, "ICER.euro.qaly.uc.ews"]       <- PSA.sim[[j]]$ICER.euro.qaly.uc.ews
            
        }
        
        # Calculate average of all PSA loops
        psa.results <- rbind(psa.results, c(colSums(psa.results)/nloops.psa))
        
        # Create loop number column and average
        psa.results <- data.frame(loop.number=c(1:nloops.psa, "Average"), psa.results)
        
        write.csv(psa.results, file = "./Results/PSA_results.csv")
        
        }
    }
    
}