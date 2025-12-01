true_kappa <- function(probs, weights = diag(ncol(probs))){
  q <- ncol(probs)
  pa <- sum(weights * probs) # percent agreement
  pe <- sum(weights*((probs %*% rep(1,q)) %*% t(t((t(rep(1,q)) %*% probs))))) # expected agreement
  kappa <- (pa - pe)/(1 - pe) # weighted kappa
  return(kappa)
}

creatDesignMatrix <- function(R = 3){
  N <- R^2
  k <- (R-1)*2+4 #4: u, tau, nu and tau-1
  X <- matrix(0, nrow = N, ncol = R)
  X[,1] <- 1
  for ( j in 1:(R-1)){
    X[((j-1)*R+1):((j-1)*R+R), (j+1)] <- rep(1, R)
    X[(N-R+1):N, (j+1)] <- rep(-1, R)
  }
  X <- cbind(X, apply(rbind(diag(1, (R-1)), rep(-1, (R-1))), 2, rep, R))
  X <- cbind(X,rep(1:R, each = R))
  X <- cbind(X,rep(1:R, R))
  X <- cbind(X,X[,(k-2)]-1)
  X[1:R,k] <- 1 
  return(list(des=X, k = k))
}

runModel <- function(table, maxN = 10, #xi1 = 0.2, xi2 = 0.3, eta1 = 0.5, eta2 = 0.6, alphaB = 4, alphaU = 0.5,                      
                     n.chains = 2, adapt = 500, burnin = 15000000, sample = ceiling(1000/2), thin = 11 ){

  
  R <- nrow(table)
  
  m <- R
  
  design <- creatDesignMatrix(R)
  Xall <- design$des
  
  desMatrixColNo <- design$k-3 # 9 #Column number where the design ends in the Xall
  X <- Xall[,1:desMatrixColNo] # Only main effects 
  
  Ncell <- nrow(X)#/2 # 2 raters
  
  greyZones <- detectGreyZones(table)
  greyZonesDesign <- as.integer(abs(greyZones$delta) > greyZones$tau_Delta) #Vector format suits better to the JAGS model
  
  indices <- matrix(NA, nrow = R, ncol = R)
  backIndices <- matrix(NA, nrow = R*R, ncol = 2)
  count <- 0
  for (i in 1:R){
    for (j in 1:R){
      count <- count + 1
      indices[i, j] <- count
      backIndices[count, 1] <- i
      backIndices[count, 2] <- j
    }
  }
  
  table <- array(t(table))
  
  data <- list( y = table, N = sum(table), X = as.matrix(X), k = ncol(X), R = R, 
                levelX = Xall[,(desMatrixColNo+1)], levelY = Xall[,(desMatrixColNo+2)], levelXeksiBir = Xall[,(desMatrixColNo+3)],
                indices = indices, backIndices = backIndices, Ncell = Ncell, designGZ = greyZonesDesign, m = m,  maxN = maxN)
                # xi1 = xi1, xi2 = xi2, eta1 = eta1, eta2 = eta2, alphaB = alphaB, alphaU = alphaU, )

  
  cat("
model {
  # y[1:9] ~ dmulti(theta[1:9], N) # theta is the probability vector! Multinomial sampling when total sample size is predefined.
  for ( i in 1:Ncell){
    y[i] ~ dpois(theta[i]*N) # Poisson sampling when total sample size is not predefined.
    theta[i] <- z[i]/sum(z[]) # This operation is to ensure sum to one constraint.
  }
  for ( i in 1:Ncell){
     w1[i] <-  1 - xi2 + alphaB*(xi1 - xi2) + alphaB
     w2[i] <-  xi2 - alphaB*(xi1 - xi2)
     w3[i] <-  1 - xi2 + alphaB*(xi1 - xi2)
     w4[i] <-  xi2 - alphaB*(xi1 - xi2) + alphaB
     # backIndices show the values that theta1 and theta2 takes for joint pdf calculations.
     for (j in 1:maxN){
      sum1[i,j] <- ifelse(j <= backIndices[i,1] , exp(logfact(backIndices[i,1]) - logfact(j) - logfact(abs(backIndices[i,1]-j)))*
                                           exp(logfact(m-backIndices[i,1]) - logfact(m-backIndices[i,1]) - logfact(abs(m-backIndices[i,1]-backIndices[i,2]+j)))*
                                           pow(w1[i],j)*pow(w2[i],backIndices[i,1]-j)*pow(w3[i],backIndices[i,2]-j)*pow(w4[i],m-backIndices[i,1]-backIndices[i,2]+j), 0)
     }
     sum2[i] <- sum(sum1[i,])
  
     B12[i] <- exp(loggam(m+1) - loggam(backIndices[i,1]+1) - loggam(m-backIndices[i,1]+1) ) * pow(1-xi1,backIndices[i,1]) * pow(xi1,m-backIndices[i,1]) * pow(1+alphaB,-m) * sum2[i]
  
     p_theta1_theta2[i] <- (1-eta1) * (1-eta2) * (m + ifelse(backIndices[i,1] == backIndices[i,2], m*alphaU, -alphaU ))/(m*pow(m+1,2)) 
                            + ((1-eta1)*eta2*pbin(backIndices[i,2],eta2,m))/(m+1) 
                            + ((1-eta2)*eta1*pbin(backIndices[i,1],eta1,m))/(m+1) + eta1*eta2*B12[i]
     
     for (j in 1:R){
      for (k in 1:R){
        sum3sum[i,j,k] <- exp( X[indices[k,j],] %*% beta + betaStar * score[2, indices[1,j]] * score[1, k] ) 
      }
      sum3[i,j] <- sum(sum3sum[i,j,])
     }
     
     z[i] <- ifelse(levelX[i] > 1, (exp( X[i,] %*% beta + betaStar * score[2, levelY[i]] * score[1, levelX[i]] + deltaGZ[i] * designGZ[i] ) / sum3[i,levelY[i]] )*
                                          (ilogit(score[2, levelX[i]] - alpha*backIndices[i,2])  - ilogit(score[2, levelXeksiBir[i]] - alpha*backIndices[i,2]) + 1E-6  ) * p_theta1_theta2[i], 
                                    (exp( X[i,] %*% beta + betaStar * score[2, levelY[i]] * score[1, levelX[i]] + deltaGZ[i] * designGZ[i]) / sum3[i,levelY[i]] )* 
                                          ilogit(score[2, levelX[i]] - alpha*backIndices[i,2]) * p_theta1_theta2[i] 
                                          
                    )
                    
                    
  }

  alpha ~ dnorm(0, 1E-6)T(0,)
  
  for (j in 1:(R*R)){
    deltaGZ[j] ~ dnorm(0, 1/4)
  }
  
  for (i in 1:2){
    for (j in 1:R){
      scoreStar[i,j] ~ dunif(0,R)
    }
    score[i, 1:R] <- sort(scoreStar[i, 1:R])
  }

  betaStar  ~ dnorm(0, 1E-6)
  for (i in 1:k){
    beta[i] ~ dnorm(0, 1E-6)
  }
  
  #scorediff <- abs(score[1,] - score[2,]) + 1E-6
  for ( i in 1:R){
    for ( j in 1:R){
      #weight[i,j] <- ifelse( i == j, 1, 1 - (abs(score[1,i] - score[2,j])/(max(score)-min(score)+1E-10))) # Linear weigths
      weight[i,j] <- ifelse( i == j, 1, 1 - (pow(score[1,i] - score[2,j], 2)/pow(max(score)-min(score)+1E-10, 2))) # Quadratic weights
    }
  }

  for ( i in 1:R){
    thetaRowSum[i] <- sum(theta[indices[i, ]]) # theta1., theta2., theta3., ...
    thetaColSum[i] <- sum(theta[indices[ ,i]]) # theta.1, theta.2, theta.3, ...
  }

  for ( i in 1:R){
    for ( j in 1:R){
      agreeSum[i,j] <- theta[indices[i,j]]*weight[i,j]+1E-10
      cagreeSum[i,j] <- weight[i,j]*thetaRowSum[i]*thetaColSum[j]+1E-10
    }
  }
  
  xi1 ~ dunif(0,1)# dbeta(1, 9)
  xi2 ~ dunif(0,1)#dbeta(1, 9)
  eta1 ~ dunif(0,1)#dbeta(1, 9)
  eta2 ~ dunif(0,1)#dbeta(1, 9)
  alphaU ~ dunif(-m, m)

  alphaB_LL <- ifelse(xi1 < xi2,  max(-(xi2/(1-xi1+xi2)), (xi2-1)/(1+xi1-xi2)), 
                      ifelse(xi1 > xi2, max(-(xi2/(1-xi1+xi2)), (xi2-1)/(1+xi1-xi2)),
                      max(-xi1, -(1-xi1))))
  alphaB_UL <- ifelse(abs(xi1-xi2) < 0.001, 500,
                      ifelse(xi1 < xi2,  (1-xi2)/(xi2-xi1+1E-5), 
                             ifelse(xi1 > xi2, xi2/(xi1-xi2+1E-5),
                      50)))
  alphaB ~ dunif(alphaB_LL, alphaB_UL)
  
  rho_R <- ((1-eta1) * (1-eta2) * (alphaU*(m+2)/12) + eta1*eta2*m*alphaB*xi1*(1-xi1)/(1+alphaB))/
    sqrt( ((1-eta1)*m*(((2*m+1)/6) - ((1-eta1)*m/4)) + eta1*m*(1-xi1)*xi1*(1-m*(1-eta1)))*
            ((1-eta2)*m*(((2*m+1)/6) - ((1-eta2)*m/4)) + eta2*m*(1-xi2)*xi2*(1-m*(1-eta2))))
  
  rho_R12_T <- (alphaB/(1+alphaB))*sqrt(xi1*(1-xi1)/xi2*(1-xi2))
  rho_R12_U <- alphaU/m
 
  rho_rater1_ratter2 <- (rho_R12_T+ rho_R12_U + rho_R12_R)/3
  
  agree <- sum(agreeSum)
  cagree <- sum(cagreeSum)
  kappa <- (agree-cagree)/(1-cagree)
  
  
  cagreeSumAC2 <- (thetaRowSum + thetaColSum + 1E-10)/2
  sumWeights <- sum(weight)
  
  cagreeAC2 <- sumWeights*sum(cagreeSumAC2*(1-cagreeSumAC2))/(R*(R-1))
  GwetAC2 <- (agree-cagreeAC2)/(1-cagreeAC2)
   
  cagreeBP <- (sumWeights + 1E-10)/pow(R,2)
  BP <- (agree-cagreeBP)/(1-cagreeBP)
  
  d <- agree-cagree
  
}", file="modelText.txt")
  
  params <- c("kappa","rho_rater1_ratter2", "rho_R12_T", "rho_R12_U", "rho_R12_R", "GwetAC2", "BP",  "betaStar", "score")
  # , "agree","cagree", "d",
  #             "theta", "beta", "betaStar", "score", 
  #             "p_theta1_theta2","deltaGZ", "GwetAC2", "BP")
  
  runJagsOut <- run.jags( method = "parallel",
                          model = "modelText.txt",
                          monitor = params,
                          data = data,
                          n.chains = n.chains,
                          adapt = adapt,
                          burnin = burnin, 
                          # inits = inits,
                          sample = sample, 
                          thin = thin,
                          summarise = FALSE,
                          plots = FALSE)
  
  summaryInfo <- summary(runJagsOut)#$summaries
  codaSamples = as.mcmc.list( runJagsOut )
  
  # diagMCMC( codaSamples , parName="rho_rater1_ratter2" )
  # diagMCMC( codaSamples , parName="rho_R12_T" )
  # diagMCMC( codaSamples , parName="rho_R12_U" )
  # graphics.off()
  return(list(summaryInfo = summaryInfo, codaSamples = codaSamples, runJagsOut = runJagsOut))
}