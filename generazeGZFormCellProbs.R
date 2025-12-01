library(GreyZones)

true_kappa <- function(probs, weights = diag(ncol(probs))){
  q <- ncol(probs)
  pa <- sum(weights * probs) # percent agreement
  pe <- sum(weights*((probs %*% rep(1,q)) %*% t(t((t(rep(1,q)) %*% probs))))) # expected agreement
  kappa <- (pa - pe)/(1 - pe) # weighted kappa
  return(kappa)
}

generateGZ_from_cell_probs <- function(cell_probs, r, gamma = 0.01, epsilon = 0.000001, maxPass = 10000){
  # cell_probs, an RxR matrix, brings in the table with no grey zones.
  # N is the total cell count.
  # r shows the first cell where a grey zone will be created.
  # gamma is the max step size in the search.
  # epsilon is the stopping criteria.
  # maxPass is the maximum number of search iterations.
  targetKappa <- true_kappa(probs = cell_probs)
  R <- nrow(cell_probs)
  if (r >= R){
    stop("A grey zone should be started from a cell less than the number of rows of the table.")
  }
  cell_probsRaw <- cell_probs
  cont <- TRUE
  direction <- "negative"
  passCount <- 0
  while (cont){
    passCount <- passCount + 1 
    # direction <- round(runif(1,0,1))
    step <- runif(1,0, gamma)
    # if ( direction == 0){
      cell_probs[r, r] = cell_probs[r, r] - step*cell_probs[r, r]
      cell_probs[(r), (r+1)] = cell_probs[(r), (r+1)] + step*cell_probs[r, r]
      # cell_probs[R, R] = cell_probs[R, R] + 0*step/4
      if (cell_probs[(r), (r+1)] < 0){
        cell_probs[(r), (r+1)] <- cell_probsRaw[(r), (r+1)]
      }
    # } else {
    #   cell_probs[r, r] = cell_probs[r, r] -  step*cell_probs[r, r]
    #   cell_probs[(r+1), r] = cell_probs[(r+1), (r)] + step*cell_probs[r, r]
    #   # cell_probs[R, R] = cell_probs[R, R] + 0*step/4
    #   if (cell_probs[(r+1), r] < 0){
    #     cell_probs[(r+1), r] <- cell_probsRaw[(r+1), r]
    #   }
    # }
    kappaNew <- true_kappa(probs = cell_probs)
    diff <- abs(kappaNew - targetKappa)
    # if (any(cell_probs<0)){
    #   print(cell_probs)
    # }
    # print(cont)
    GZres <- detectGreyZones(round(cell_probs*500))
    # print(passCount)
    # print(kappaNew)
    # print(diff)
    # print(GZres$Delta)
    # print(GZres$tau_Delta)
    # print(cell_probs)
    if ((diff <= epsilon) & (GZres$Delta > GZres$tau_Delta)){
      cont <- FALSE
      result <- 1
      # print(GZres$result)
    }
    if (passCount >= maxPass){
      cont <- FALSE
      result <- 0
    }
  }
  return(list(GZtableCellProbs = cell_probs, diff = diff, result = result, targetKappa = targetKappa, 
              GZkappa = kappaNew, numIter = passCount))
}

