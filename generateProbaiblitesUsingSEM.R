
# https://stats.oarc.ucla.edu/r/seminars/rsem/#s1a
#' ONLINE{newtest,
#'   author = {Bruin, J.},
#'   title = {newtest: command to compute new test {ONLINE}},
#'   month = FEB,
#'   year = {2011},
#'   url = {https://stats.oarc.ucla.edu/stata/ado/analysis/}
#' }

library(lavaan)
library(DescTools)
library(psych)
library(mnorm)


parValues <- read.csv("ParameterVauesForSims.csv")

for (i in 1:nrow(parValues)){
  filename <- paste0("cellProbabilities",parValues$Dimension[i],"x",parValues$Dimension[i],"_Group", parValues$group[i],".csv")
  beta_l1.x1 <- parValues$beta_l1.x1[i] # Regression coefficient between latent1 and Rater1
  beta_l2.x2 <- parValues$beta_l2.x2[i] # Regression coefficient between latent2 and Rater2
  beta_l1.l2 <- parValues$beta_l1.l2[i] # Regression coefficient between latent1 and latent2
  beta_x1.x2 <- parValues$beta_x1.x2[i] # Regression coefficient between Rater1 and Rater2
  if (parValues$Dimension[i] == 3){
    useBreaks <- qnormFast(c(10^-6, 1/3, 2/3, (1-10^-6)), mean = 0, sd = 1) # 3x3
    useLabels <- c("d", "o", "y") # 3x3
  } else if (parValues$Dimension[i] == 4){
    useBreaks <- qnormFast(c(10^-6, 1/4, 2/4, 3/4, (1-10^-6)), mean = 0, sd = 1) # 4x4
    useLabels <- c("d", "o", "y", "cy")
  } else if (parValues$Dimension[i] == 5){
    useBreaks <- qnormFast(c(10^-6, 1/5, 2/5, 3/5, 4/5, (1-10^-6)), mean = 0, sd = 1) # 5x5
    useLabels <- c("cd", "d", "o", "y", "cy")
  }
  popModel <- paste0("
          latent1 =~", beta_l1.x1, "*x1
          latent2 =~ ", beta_l2.x2, "*x2 
          latent1 ~~ ", beta_l1.l2,"*latent2
          x1 ~ ", beta_x1.x2, "*x2
          # Set all variances to 1
          # x1 ~~ 1*x1
          # x2 ~~ 1*x2
          x1~~x2
          # latent1 ~~ 1*latent1
          # latent2 ~~ 1*latent2
        ")
  # Population covariance
  fitPopModel <- fitted(sem(popModel, std.lv = TRUE))
 
  set.seed(1234)
  simData <- simulateData(popModel, sample.nobs=1000000, std.lv = TRUE, skewness = 0, kurtosis = 0, model.type = "sem")
 
  # We will get the breaks according to standard normal distribution.
  simData$x1_ord <- cut(simData$x1, breaks = useBreaks,
                        labels = useLabels, include.lowest = TRUE)
  simData$x2_ord <- cut(simData$x2, breaks = useBreaks,
                        labels = useLabels, include.lowest = TRUE)
  
  tablo <- table(simData$x1_ord, simData$x2_ord)
  probTablo <- tablo/sum(tablo)
  write.csv(probTablo, filename)
}

