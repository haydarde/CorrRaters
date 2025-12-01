
library(irrCAC)
library(irr)
library(rTableICC)
library(runjags)
library(GreyZones)
library(tictoc)
library(coda)

source("runModel.R")

table <- matrix(c(264,21,4,0,0,
                  184,16,3,1,1,
                  5,1,0,0,0,
                  5,0,1,0,0,
                  0,0,0,0,0), nrow =5, ncol =5, byrow = TRUE)

classicalKappa <- kappa2.table(table,weights = linear.weights(1:ncol(table)))
classicalKappaEst <- classicalKappa$coeff.val
classicalKappaCIlimits <- as.numeric(unlist(strsplit(gsub("[()]", "", classicalKappa$coeff.ci),",")))

AC2 <- gwet.ac1.table(table, weights = quadratic.weights(1:5))
AC2Est <- AC2$coeff.val
AC2CIlimits <- as.numeric(unlist(strsplit(gsub("[()]", "", AC2$coeff.ci),",")))   

detectGreyZones(table)  

BP2 <- bp2.table(table, weights = quadratic.weights(1:5))
BP2Est <- BP2$coeff.val
BP2CIlimits <- as.numeric(unlist(strsplit(gsub("[()]", "", BP2$coeff.ci),",")))  

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
# diagMCMC( codaSamples , parName="GwetAC2" )
# diagMCMC( codaSamples , parName="BP" )
# diagMCMC( codaSamples , parName="kappa" )
# diagMCMC( codaSamples , parName="rho_rater1_ratter2" )
# diagMCMC( codaSamples , parName="rho_R12_T" )
# diagMCMC( codaSamples , parName="rho_R12_U" )
# # graphics.off()
# 
# 
# nChains <- 2
# furtherThin <- 1
# thiningSequence <- seq(500,nrow(codaSamples[[1]]), furtherThin) 
# newCodaSamples <- mcmc.list() 
# for ( i in 1:nChains){   
#   newCodaSamples[[i]] <- as.mcmc(codaSamples[[i]][thiningSequence,]) 
# }
# 
# diagMCMC( newCodaSamples , parName="kappa" )
# diagMCMC( newCodaSamples , parName="GwetAC2" )
# diagMCMC( newCodaSamples , parName="BP" )
# diagMCMC( newCodaSamples , parName="rho_rater1_ratter2" )
# diagMCMC( newCodaSamples , parName="rho_R12_T" )
# diagMCMC( newCodaSamples , parName="rho_R12_U" )
# 
# summaryNewCodaSamples <- summary(newCodaSamples)

res$summaryInfo <- summary(codaSamples)

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

mean(c(rhoR1R2Est,rhoR1R2EstT,rhoR1R2EstU))


rHatKappa <- res$summaryInfo[1,11]
rHatRho <- res$summaryInfo[2,11]
rHatRhoT <- res$summaryInfo[3,11]
rHatRhoU <- res$summaryInfo[4,11]


