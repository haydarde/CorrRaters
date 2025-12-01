
library(irrCAC)
library(irr)
library(rTableICC)
library(runjags)
library(GreyZones)
library(tictoc)
library(coda)

source("runModel.R")

table <- matrix(c(345,	35,	1,	0,
                  129,	79,	13,	1,
                  25,	39,	40,	8,
                  2,	7,	13,	30), nrow = 4, ncol = 4, byrow = TRUE)

classicalKappa <- kappa2.table(table,weights = linear.weights(1:ncol(table)))
classicalKappaEst <- classicalKappa$coeff.val
classicalKappaCIlimits <- as.numeric(unlist(strsplit(gsub("[()]", "", classicalKappa$coeff.ci),",")))

tic()
res <- runModel(table, n.chains = 2, adapt = 500, burnin = 150000, sample = ceiling(1000/2), thin = 31 )
toc()

codaSamples <- res$codaSamples

# If additional runs is needed:
# runJagsOutExtended <- extend.jags(runjags.object = res$runJagsOut, burnin = 300000, sample = 1000)
# 
# runJagsOutExtended <- extend.jags(runjags.object = runJagsOutExtended, burnin = 1000000, sample = 5000)
# 
# codaSamples <- as.mcmc.list(runJagsOutExtended)
# 
# diagMCMC( codaSamples , parName="kappa" )
# diagMCMC( codaSamples , parName="rho_rater1_ratter2" )
# diagMCMC( codaSamples , parName="rho_R12_T" )
# diagMCMC( codaSamples , parName="rho_R12_U" )
# graphics.off()

# save(codaSamples, file ="codaSamples_TTE_CMR.RData")

# To further apply burn-in and thinning if additional iterations are run
# nChains <- 2
# furtherThin <- 1
# thiningSequence <- seq(55000,nrow(codaSamples[[1]]), furtherThin) 
# newCodaSamples <- mcmc.list() 
# for ( i in 1:nChains){   
#   newCodaSamples[[i]] <- as.mcmc(codaSamples[[i]][thiningSequence,]) 
# }
# 
# diagMCMC( newCodaSamples , parName="kappa" )
# diagMCMC( newCodaSamples , parName="rho_rater1_ratter2" )
# diagMCMC( newCodaSamples , parName="rho_R12_T" )
# diagMCMC( newCodaSamples , parName="rho_R12_U" )

# summaryNewCodaSamples <- summary(newCodaSamples)



res$summaryInfo <- summary(runJagsOutExtended)

# Get the required point and credible interval estimates
kappaEst <- res$summaryInfo[1,2]
kappaEstMean <- res$summaryInfo[1,2]
kappaCIlimits <- c(res$summaryInfo[1,1], res$summaryInfo[1,3])
rhoR1R2Est <- res$summaryInfo[2,2]
rhoR1R2CIlimits <- c(res$summaryInfo[2,1], res$summaryInfo[2,3])

rhoR1R2EstT <- res$summaryInfo[3,2]
rhoR1R2CIlimitsT <- c(res$summaryInfo[3,1], res$summaryInfo[3,3])

rhoR1R2EstU <- res$summaryInfo[4,2]
rhoR1R2CIlimitsU <- c(res$summaryInfo[4,1], res$summaryInfo[4,3])

mean(c(rhoR1R2Est,rhoR1R2EstT,rhoR1R2EstU))

rhoR1R2Est <- res$summaryInfo[2,4]
rhoR1R2CIlimits <- c(res$summaryInfo[2,1], res$summaryInfo[2,3])

rhoR1R2EstT <- res$summaryInfo[3,4]
rhoR1R2CIlimitsT <- c(res$summaryInfo[3,1], res$summaryInfo[3,3])

rhoR1R2EstU <- res$summaryInfo[4,4]
rhoR1R2CIlimitsU <- c(res$summaryInfo[4,1], res$summaryInfo[4,3])

rHatKappa <- res$summaryInfo[1,11]
rHatRho <- res$summaryInfo[2,11]
rHatRhoT <- res$summaryInfo[3,11]
rHatRhoU <- res$summaryInfo[4,11]


