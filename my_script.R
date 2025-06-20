library("rstudioapi") # Load rstudioapi package
setwd(dirname(getActiveDocumentContext()$path)) # Set working directory to source file location
getwd() # Check updated working directory


source("RD_and_DT_Algorithm_copy.R")
#source("RD_and_DT_Algorithm_copy.R")


#Clutter_gen <- Clutter_gen(0.5, 3, 20)

#RD------------------------- 
results <- data.frame(
  Run = integer(),
  Lambda = numeric(),
  Length = numeric(),
  Cost = numeric(),
  NumDisambigs = integer()
)

lambda = 0

for (i in 1:50) {
  set.seed(i)
  obs_gen_para <- c(gamma = 0.3, d = 5, noPoints = 25)
  result <- ACS_Alg_C(obs_gen_para, k = 0, lambda)
  
  results[i, ] <- list(
    Run = i,
    Lambda = lambda,
    Length = result$Length_total,
    Cost = result$Cost_total,
    NumDisambigs = length(result$Disambiguate_state)
  )
}

saveRDS(results, file = "data_lambda_0.rds")
