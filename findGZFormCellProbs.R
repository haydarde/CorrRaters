
source("generazeGZFormCellProbs.R")

k <- 3 # Table size
g <- 1 # Group as described in Excel, LatenCor-AgreementSettings.xlsx and ParameterVauesForSims.xlsx/csv

for (k in 3:5){
  for (g in c(1,2,3,5,6,7,8,11,12,13)){
    
cellProbs <- read.csv(paste0("cellProbabilities/", "cellProbabilities", k, "x", k, "_Group",g,".csv"))

cellProbs <- as.matrix(cellProbs[,-1])

GZtableCellProbs <- generateGZ_from_cell_probs(cell_probs = cellProbs, r = 1, gamma = 0.001, epsilon = 0.035, maxPass = 1000)

print(paste0("cellProbabilities", k, "x", k, "_Group",g))
print(detectGreyZones(round(GZtableCellProbs$GZtableCellProbs*500))$result)

write.csv(GZtableCellProbs,paste0("cellProbabilities/", "GZ_cellProbabilities", k, "x", k, "_Group",g,".csv"))

  }
}



cellProbs <- read.csv(paste0("cellProbabilities/", "cellProbabilities", k, "x", k, "_Group",g,".csv"))

cellProbs <- as.matrix(cellProbs[,-1])

GZtableCellProbs <- generateGZ_from_cell_probs(cell_probs = cellProbs, r = 1, gamma = 0.05, epsilon = 0.1, maxPass = 1000)

GZtableCellProbs$GZtableCellProbs
