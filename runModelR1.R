true_kappa <- function(probs, weights = diag(ncol(probs))){
  q <- ncol(probs)
  pa <- sum(weights * probs) # percent agreement
  pe <- sum(weights*((probs %*% rep(1,q)) %*% t(t((t(rep(1,q)) %*% probs))))) # expected agreement
  kappa <- (pa - pe)/(1 - pe) # weighted kappa
  return(kappa)
}

true_AC2 <- function(probs, weights = diag(ncol(probs))){
  q <- ncol(probs)
  pa <- sum(weights * probs) # observed agreement
  # Marginal probabilities for each rater
  p_row <- rowSums(probs)  # rater 1
  p_col <- colSums(probs)  # rater 2
  # Expected agreement under Gwet's framework:
  p_avg <- 0.5 * (p_row + p_col)
  
  # Expected agreement (weighted)
  Pe_mat <- outer(p_avg, p_avg, "*")  
  pe  <- sum(weights * Pe_mat)

  AC2 <- (pa - pe)/(1 - pe) # weighted kappa
  return(AC2)
}

# P <- matrix(c(0.30, 0.05, 0.05,
#               0.04, 0.25, 0.06,
#               0.03, 0.07, 0.15),
#             nrow = 3, byrow = TRUE)
# 
# # Identity weights -> unweighted AC2
# true_AC2(P)#, weights = linear.weights(c(1,2,3)))
# 
# agreementToRatings <- function(A) {
#   if (is.null(rownames(A))) rownames(A) <- seq_len(nrow(A))
#   if (is.null(colnames(A))) colnames(A) <- seq_len(ncol(A))
#   
#   r1_cats <- rep(rownames(A), times = ncol(A))
#   r2_cats <- rep(colnames(A), each  = nrow(A))
#   counts  <- as.vector(A)
#   
#   r1 <- rep(r1_cats, times = counts)
#   r2 <- rep(r2_cats, times = counts)
#   
#   data.frame(rater1 = r1, rater2 = r2)
# }
# 
# library(irrCAC)
# ratings <- agreementToRatings(P*100)
# ac2_result <- gwet.ac1.raw(ratings, weights = "unweighted")# weights = linear.weights(c(1,2,3)))
# print(ac2_result)



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

runModel <- function(table, maxN = 10, inits = NULL, #xi1 = 0.2, xi2 = 0.3, eta1 = 0.5, eta2 = 0.6, alphaB = 4, alphaU = 0.5,                      
                     n.chains = 3, adapt = 500, burnin = 100000, sample = ceiling(1000/2), thin = 3 ){

  
  R <- nrow(table)
  
  m <- R
  
  design <- creatDesignMatrix(R)
  Xall <- design$des
  
  desMatrixColNo <- design$k-3 # 9 #Column number where the design ends in the Xall
  X <- Xall[,1:desMatrixColNo] # Only main effects 
  
  Ncell <- nrow(X)#/2 # 2 raters
  
  greyZones <- detectGreyZones(table)
  greyZonesDesign <- as.integer(abs(greyZones$delta) > greyZones$tau_Delta) #Vector format suits better to the JAGS model
  greyZoneExists <- as.integer(abs(greyZones$Delta) > greyZones$tau_Delta)
  
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
                indices = indices, backIndices = backIndices, Ncell = Ncell, designGZ = greyZonesDesign, greyZoneExists = greyZoneExists, m = m,  maxN = maxN)
                # xi1 = xi1, xi2 = xi2, eta1 = eta1, eta2 = eta2, alphaB = alphaB, alphaU = alphaU, )

  
  cat("
model {
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

  alpha ~ dnorm(0, 1/tau_alpha)T(0,)
  sigma_alpha ~ dnorm(0, 0.004)T(0,)
  tau_alpha <- 1 / pow(sigma_alpha, 2)

  for (j in 1:(R*R)){
    deltaGZ[j] ~ dnorm(0, 1/tau_delta)
  }
  sigma_delta ~ dnorm(0, 0.004)T(0,)
  tau_delta <- 1 / pow(sigma_delta, 2)
  
  for (i in 1:2){
    score_raw[i,1] <- 0
    for (j in 2:R){
      delta_score[i,j] ~ dnorm(0,0.004)T(0,) 
      score_raw[i,j] <- score_raw[i,j-1] + delta_score[i,j]
    }
    
    for (j in 1:R){
      score[i,j] <- score_raw[i,j] / score_raw[i,R]
    }
  }
  

  betaStar  ~ dnorm(0, 1/tau_betaStar)
  sigma_betaStar ~ dnorm(0, 0.004)T(0,)
  tau_betaStar <- 1 / pow(sigma_betaStar, 2)
  for (i in 1:k){
    beta[i] ~ dnorm(0, 1/tau_beta[i])
    sigma_beta[i] ~ dnorm(0, 0.004)T(0,)
    tau_beta[i] <- 1 / pow(sigma_beta[i], 2)
  }
  
  for ( i in 1:R){
    for ( j in 1:R){
      weight[i,j] <- ifelse( i == j, 1, 1 - (abs(score[1,i] - score[2,j])/(max(score)-min(score)+1E-10))) # Linear weigths
      # weightQuadratic[i,j] <- ifelse( i == j, 1, 1 - (pow(score[1,i] - score[2,j], 2)/pow(max(score)-min(score)+1E-10, 2))) # Quadratic weights
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
  
  # for ( i in 1:R){
  #   for ( j in 1:R){
  #     agreeSumQuadratic[i,j] <- theta[indices[i,j]]*weightQuadratic[i,j]+1E-10
  #     cagreeSumQuadratic[i,j] <- weightQuadratic[i,j]*thetaRowSum[i]*thetaColSum[j]+1E-10
  #   }
  # }
  
  
  logit_xi1  ~  dnorm(0, .5)      
  logit_xi2  ~ dnorm(0, .5)
  logit_eta1 ~ dnorm(0, .5) 
  logit_eta2 ~ dnorm(0, .5)
  logit_alphaU ~dnorm(0, .5)
  
  xi1    <- ilogit(logit_xi1)
  xi2    <- ilogit(logit_xi2)
  eta1   <- ilogit(logit_eta1)
  eta2   <- ilogit(logit_eta2)
  alphaU <- -1 + 2* ilogit(logit_alphaU)

  alphaB ~ dnorm(0, 0.004)
  
  rho_R12_R <- ((1-eta1) * (1-eta2) * (alphaU*(m+2)/12) + eta1*eta2*m*alphaB*xi1*(1-xi1)/(1+alphaB))/
    sqrt( ((1-eta1)*m*(((2*m+1)/6) - ((1-eta1)*m/4)) + eta1*m*(1-xi1)*xi1*(1-m*(1-eta1)))*
            ((1-eta2)*m*(((2*m+1)/6) - ((1-eta2)*m/4)) + eta2*m*(1-xi2)*xi2*(1-m*(1-eta2))))
  
  rho_R12_T <- (alphaB/(1+alphaB))*sqrt(xi1*(1-xi1)/xi2*(1-xi2))
  rho_R12_U <- alphaU/m

  rho_rater1_ratter2 <- (rho_R12_T+ rho_R12_U + rho_R12_R)/3

  agree <- sum(agreeSum)
  cagree <- sum(cagreeSum)
  kappa <- (agree-cagree)/(1-cagree)

  # agreeQuadratic <- sum(agreeSumQuadratic)
  # cagreeQuadratic <- sum(cagreeSumQuadratic)
  # kappaQuadratic <- (agreeQuadratic-cagreeQuadratic)/(1-cagreeQuadratic)  
  
  cagreeSumAC2 <- (thetaRowSum + thetaColSum + 1E-10)/2
  sumWeights <- sum(weight)
  
  # sumWeightsQuadratic <- sum(weightQuadratic)
  
  cagreeAC2 <- sumWeights*sum(cagreeSumAC2*(1-cagreeSumAC2))/(R*(R-1))
  GwetAC2 <- (agree-cagreeAC2)/(1-cagreeAC2)
  
  # cagreeAC2Quadratic <-  sumWeightsQuadratic*sum(cagreeSumAC2*(1-cagreeSumAC2))/(R*(R-1))
  # GwetAC2Quadratic <- (agree-cagreeAC2Quadratic)/(1-cagreeAC2Quadratic)
  
  cagreeBP <- (sumWeights + 1E-10)/pow(R,2)
  BP <- (agree-cagreeBP)/(1-cagreeBP)
  
  # cagreeBPQuadratic <- (sumWeightsQuadratic + 1E-10)/pow(R,2)
  # BPQuadratic <- (agree-cagreeBPQuadratic)/(1-cagreeBPQuadratic)
  
  d <- agree-cagree
  
}", file="modelText.txt")
  
  params <- c("kappa","rho_rater1_ratter2", "rho_R12_T", "rho_R12_U", "rho_R12_R", "GwetAC2", "BP",  
              "betaStar", "score", "xi1","xi2","eta1","eta2","alphaB","alphaU","alpha")#, 
              # "kappaQuadratic","GwetAC2Quadratic","BPQuadratic",)
  
  runJagsOut <- run.jags( method = "parallel",
                          model = "modelText.txt",
                          monitor = params,
                          data = data,
                          n.chains = n.chains,
                          adapt = adapt,
                          burnin = burnin, 
                          inits = ifelse(!is.null(inits), inits, NA),
                          sample = sample, 
                          thin = thin,
                          summarise = FALSE,
                          plots = FALSE)
  
  summaryInfo <- summary(runJagsOut)#$summaries
  codaSamples = as.mcmc.list( runJagsOut )

  return(list(summaryInfo = summaryInfo, codaSamples = codaSamples, runJagsOut = runJagsOut))
}