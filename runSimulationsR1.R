setwd("~/Documents/makaleler/Agreement/toGitHub_R1")

library(irrCAC)
library(irr)
library(rTableICC)
library(runjags)
library(GreyZones)
library(tictoc)
library(coda)

source("runModelR1_V8.R")
source("DBDA2E-utilities.R")


simRuns <- read.csv("simulationDriver.csv")

simRep <- 200 
simStart <- 1
simEnd <- 300

nSimRun <- simEnd - simStart + 2
  
allResults <- NULL
saveRes <- data.frame()


for (i in simStart:simEnd){ 
  cellProbs <- read.csv(paste0("cellProbabilities/",simRuns$File[i]))
  
  p <- as.matrix(cellProbs[,2:(nrow(cellProbs)+1)])
  if (ncol(cellProbs) > nrow(cellProbs) ){# implies a GZ scenario
    kappaTrue <- cellProbs[1,(ncol(cellProbs)-1)] 
  } else {
    AC2True <- true_AC2(probs = p, weights = linear.weights(1:ncol(table)))
    kappaTrue <- true_kappa(probs = p, weights = linear.weights(1:ncol(table)))
  }
  rhoR1R2True <- simRuns$rhoR1R1True[i]
  j <- 0
  while (j < simRep){
    if (any(p == 0)){
      p[which(p == 0)] <- 10^-6
    }
    
    table <- rTable.RxC(p = p, N = simRuns$N[i])$rTable
    classicalKappa <- kappa2.table(table, weights = linear.weights(1:ncol(table)))
    classicalKappaEst <- classicalKappa$coeff.val
    classicalKappaCIlimits <- as.numeric(unlist(strsplit(gsub("[()]", "", classicalKappa$coeff.ci),",")))
    
    classicalAC2 <-  gwet.ac1.table(table, weights = linear.weights(1:ncol(table)))
    classicalAC2Est <- classicalAC2$coeff.val
    classicalAC2CIlimits <- as.numeric(unlist(strsplit(gsub("[()]", "", classicalAC2$coeff.ci),",")))

    res1 <- tryCatch({
      tic()
        res <- runModel(table, n.chains = 3, burnin = 10000, sample = 1000, thin = 3)
      toc()
       if ((res$summaryInfo[1, 11] < 1.1) & (res$summaryInfo[2, 11] < 1.1)){ # Accept the run based on Rhat
        
        j <- j+1
        
        allResults <-rbind(allResults,cbind(rep(i,46),rep(j,46),res$summaryInfo))
       
        kappaEst <- res$summaryInfo[1,2] 
        kappaEstMean <- res$summaryInfo[1,2]  
        kappaCIlimits <- c(res$summaryInfo[1,1], res$summaryInfo[1,3])  

        AC2Est <- res$summaryInfo[6,2]
        AC2EstMean <- res$summaryInfo[6,2]  
        AC2CIlimits <- c(res$summaryInfo[6,1], res$summaryInfo[6,3]) 

        rhoR1R2Est <- res$summaryInfo[2,2]
        rhoR1R2CIlimits <- c(res$summaryInfo[2,1], res$summaryInfo[2,3])
        
        rhoR1R2EstT <- res$summaryInfo[3,2]
        rhoR1R2CIlimitsT <- c(res$summaryInfo[3,1], res$summaryInfo[3,3])
        
        rhoR1R2EstU <- res$summaryInfo[4,2]
        rhoR1R2CIlimitsU <- c(res$summaryInfo[4,1], res$summaryInfo[4,3])
        
        rhoR1R2EstR <- res$summaryInfo[5,2]
        rhoR1R2CIlimitsR <- c(res$summaryInfo[5,1], res$summaryInfo[5,3])
        
        rHatKappa <- res$summaryInfo[1,11]
        rHatRho <- res$summaryInfo[2,11]
        rHatRhoT <- res$summaryInfo[3,11]
        rHatRhoU <- res$summaryInfo[4,11]
        
        saveRes <- rbind(saveRes, data.frame(Scenario = i, Replication = j, 
                                      BayesianKappaMedian = kappaEst, BIsLL = kappaCIlimits[1] , BIsUL = kappaCIlimits[2], 
                                      rho = rhoR1R2Est, rhoBIsLL = rhoR1R2CIlimits[1], rhoBIsUL = rhoR1R2CIlimits[2],
                                      rhoT = rhoR1R2EstT, rhoBIsLLT = rhoR1R2CIlimitsT[1], rhoBIsULT = rhoR1R2CIlimitsT[2],
                                      rhoU = rhoR1R2EstU, rhoBIsLLU = rhoR1R2CIlimitsU[1], rhoBIsULU = rhoR1R2CIlimitsU[2],
                                      rhoR = rhoR1R2EstR, rhoBIsLLR = rhoR1R2CIlimitsR[1], rhoBIsULR = rhoR1R2CIlimitsR[2],
                                      ClassicalKappa = classicalKappaEst, CIsLL = classicalKappaCIlimits[1], CIsUL = classicalKappaCIlimits[2],
                                      TrueKappa = kappaTrue, rhoTrue = rhoR1R2True, TrueAC2 = AC2True,
                                      rHatKappa = rHatKappa, rHatRho = rHatRho, rHatRhoT = rHatRhoT, rHatRhoU = rHatRhoU,
                                      BayesianAC2Median = AC2Est, AC2BIsLL = AC2CIlimits[1] , AC2BIsUL = AC2CIlimits[2], BayesianAC2Mean = res$summaryInfo[6,4], 
                                      ClassicalAC2 = classicalAC2Est, AC2CIsLL = classicalAC2CIlimits[1], AC2CIsUL = classicalAC2CIlimits[2])
          )
       
      } else {
        cat("j:", j, " repeating ...")
      }
      
    }, error = function(err){
      print(paste0("ERROR:",err))
      j <- j-1 
    })
  }
}

write.csv(saveRes, paste0("AllEstimates_CohenKappa",simStart,"_",simEnd,".csv"))

write.csv(allResults, paste0("allResults_CohenKappa",simStart,"_",simEnd,".csv"))
 


