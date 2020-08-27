## R script to reproduce the open population spatial capture recapture (OPSCR) model from 
## Milleret et al. Instrumented individuals are not representative of the population. Survival comparison of wolverines with and without GPS collars 

## LOAD LIBRARIES 
library(nimble)

## SET WORKING DIRECTORY 
setwd("C:/Personal_Cloud/OneDrive/Work/PublicCodesGit/Public/WolverineGPSSurvival")
## SOURCE NIMBLE CUSTOM FUNCTIONS 
source("pointProcess.R")
source("dbin_LESS_Cached_MultipleCovResponseGPS.R")

## LOAD DATA 
load("22.18.J_FaBeforeGPSSexCovIslandsALL_1.RData")


#RUN NIMBLE MODEL (demonstration with single chain and low number of iterations)
model <- nimbleModel( code = modelCode
                      , constants = nimConstants
                      , data = nimData
                      , inits = nimInits
                      , check = FALSE       
                      , calculate = FALSE)  
cmodel <- compileNimble(model)
cmodel$calculate()
MCMCconf <- configureMCMC(model = model, monitors = c(nimParams),
                          control = list(reflective = TRUE, adaptScaleOnly = TRUE),
                          useConjugacy = FALSE) 
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
Runtime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC
                                                  , nburnin = 0
                                                  , niter = 100
                                                  , nchains = 1
                                                  , samplesAsCodaMCMC = TRUE))
