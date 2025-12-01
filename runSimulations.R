

library(irrCAC)
library(irr)
library(rTableICC)
library(runjags)
library(GreyZones)
library(tictoc)

source("runModel.R")

simRuns <- read.csv("simulationDriver.csv")

simRep <- 100 # 1000

simStart <- 1
simEnd <- 300 

nSimRun <- simEnd - simStart + 2
  
allResults <- NULL

MAEkappa <- matrix(NA, nrow = simRep, ncol = nSimRun)
MSEkappa <- matrix(NA, nrow = simRep, ncol = nSimRun)
MAPEkappa <- matrix(NA, nrow = simRep, ncol = nSimRun)
coverageKappa <- matrix(NA, nrow = simRep, ncol = nSimRun)

MAErhoR1R2a <- matrix(NA, nrow = simRep, ncol = nSimRun)
MSErhoR1R2 <- matrix(NA, nrow = simRep, ncol = nSimRun)
MAPErhoR1R2 <- matrix(NA, nrow = simRep, ncol = nSimRun)
coveragerhoR1R2 <- matrix(NA, nrow = simRep, ncol = nSimRun)

MAEclassicalKappa <- matrix(NA, nrow = simRep, ncol = nSimRun)
MSEclassicalKappa <- matrix(NA, nrow = simRep, ncol = nSimRun)
MAPEclassicalKappa <- matrix(NA, nrow = simRep, ncol = nSimRun)
coverageClassicalKappa <- matrix(NA, nrow = simRep, ncol = nSimRun)
saveRes <- data.frame()

for (i in simStart:simEnd){
  cellProbs <- read.csv(paste0("cellProbabilities/",simRuns$File[i]))
  
  p <- as.matrix(cellProbs[,2:(nrow(cellProbs)+1)])
  if (ncol(cellProbs) > nrow(cellProbs) ){# implies a GZ scenario
    kappaTrue <- cellProbs[1,(ncol(cellProbs)-1)] 
  } else {
    kappaTrue <- true_kappa(probs = p)
  }
  rhoR1R2True <- simRuns$rhoR1R1True[i]
  j <- 0
  while (j < simRep){
    if (any(p == 0)){
      p[which(p == 0)] <- 10^-6
    }
    j <- j+1
    table <- rTable.RxC(p = p, N = simRuns$N[i])$rTable
    classicalKappa <- kappa2.table(table,weights = linear.weights(1:ncol(table)))
    classicalKappaEst <- classicalKappa$coeff.val
    classicalKappaCIlimits <- as.numeric(unlist(strsplit(gsub("[()]", "", classicalKappa$coeff.ci),",")))

    res1 <- tryCatch({
      tic()
      res <- runModel(table)
      toc()
    
      allResults <-rbind(allResults,cbind(rep(i,46),rep(j,46),res$summaryInfo))
      
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
      
      rHatKappa <- res$summaryInfo[1,11]
      rHatRho <- res$summaryInfo[2,11]
      rHatRhoT <- res$summaryInfo[3,11]
      rHatRhoU <- res$summaryInfo[4,11]
      
      
      # Assess the goodness of fit
      MAEkappa[j, (i-simStart+1)] <- abs(kappaTrue - kappaEst)
      MSEkappa[j, (i-simStart+1)] <- (kappaTrue - kappaEst)^2
      MAPEkappa[j, (i-simStart+1)] <- abs(kappaTrue - kappaEst)/kappaTrue
      coverageKappa[j, (i-simStart+1)] <- ifelse(kappaEst >= kappaCIlimits[1] && kappaEst <= kappaCIlimits[2], 1, 0 )
      
      MAErhoR1R2a[j, (i-simStart+1)] <- abs(rhoR1R2True) - abs(rhoR1R2Est)
      MSErhoR1R2[j, (i-simStart+1)] <- (abs(rhoR1R2True) - abs(rhoR1R2Est))^2
      MAPErhoR1R2[j, (i-simStart+1)] <- abs(abs(rhoR1R2True) - abs(rhoR1R2Est))/abs(rhoR1R2True)
      coveragerhoR1R2[j, (i-simStart+1)] <- ifelse(rhoR1R2Est >= rhoR1R2CIlimits[1] && rhoR1R2Est <= rhoR1R2CIlimits[2], 1, 0 )
      
      MAEclassicalKappa[j, (i-simStart+1)] <- abs(kappaTrue - classicalKappaEst)
      MSEclassicalKappa[j, (i-simStart+1)] <- (kappaTrue - classicalKappaEst)^2
      MAPEclassicalKappa[j, (i-simStart+1)] <- abs(kappaTrue - classicalKappaEst)/kappaTrue
      coverageClassicalKappa[j, (i-simStart+1)] <- ifelse(classicalKappaEst >= classicalKappaCIlimits[1] && classicalKappaEst <= classicalKappaCIlimits[2], 1, 0 )
      
      saveRes <- rbind(saveRes, data.frame(Scenario = i, Replication = j, 
                                    BayesianKappaMedian = kappaEst, BIsLL = kappaCIlimits[1] , BIsUL = kappaCIlimits[2], 
                                    rho = rhoR1R2Est, rhoBIsLL = rhoR1R2CIlimits[1], rhoBIsUL = rhoR1R2CIlimits[2],
                                    rhoT = rhoR1R2EstT, rhoBIsLLT = rhoR1R2CIlimitsT[1], rhoBIsULT = rhoR1R2CIlimitsT[2],
                                    rhoU = rhoR1R2EstU, rhoBIsLLU = rhoR1R2CIlimitsU[1], rhoBIsULU = rhoR1R2CIlimitsU[2],
                                    BayesianKappaMean = res$summaryInfo[1,4], 
                                    ClassicalKappa = classicalKappaEst, CIsLL = classicalKappaCIlimits[1], CIsUL = classicalKappaCIlimits[2],
                                    TrueKappa = kappaTrue, rhoTrue = rhoR1R2True, TrueKappaGZ = ifelse((ncol(cellProbs) > nrow(cellProbs) ),cellProbs[1,(ncol(cellProbs)-2)], kappaTrue),
                                    rHatKappa = rHatKappa, rHatRho = rHatRho, rHatRhoT = rHatRhoT, rHatRhoU = rHatRhoU)
        )
       
      
      
    }, error = function(err){
      print(paste0("ERROR:",err))
      j <- j-1 # Do the same rep in case of an error
    })
  }
}

write.csv(saveRes, paste0("AllEstimates_",simStart,"_",simEnd,".csv"))

# write.csv(MAEkappa, paste0("MAEkappa_",simStart,"_",simEnd,".csv"))
# write.csv(MSEkappa, paste0("MSEkappa_",simStart,"_",simEnd,".csv"))
# write.csv(MAPEkappa, paste0("MAPEkappa_",simStart,"_",simEnd,".csv"))
# write.csv(coverageKappa, paste0("coverageKappa_",simStart,"_",simEnd,".csv"))
# 
# write.csv(MAErhoR1R2a, paste0("MAErhoR1R2a_",simStart,"_",simEnd,".csv"))
# write.csv(MSErhoR1R2, paste0("MSErhoR1R2_",simStart,"_",simEnd,".csv"))
# write.csv(MAPErhoR1R2, paste0("MAPErhoR1R2_",simStart,"_",simEnd,".csv"))
# write.csv(coveragerhoR1R2, paste0("coveragerhoR1R2_",simStart,"_",simEnd,".csv"))
# 
# write.csv(MAEclassicalKappa, paste0("MAEclassicalKappa_",simStart,"_",simEnd,".csv"))
# write.csv(MSEclassicalKappa, paste0("MSEclassicalKappa_",simStart,"_",simEnd,".csv"))
# write.csv(MAPEclassicalKappa, paste0("MAPEclassicalKappa_",simStart,"_",simEnd,".csv"))
# write.csv(coverageClassicalKappa, paste0("coverageClassicalKappa_",simStart,"_",simEnd,".csv"))

write.csv(allResults, paste0("allResults_",simStart,"_",simEnd,".csv"))


