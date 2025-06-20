# main_script.R

# Make sure working directory is the same as the script location (implicitly handled in GitHub Actions)
print(getwd())

# Load the helper script
source("RD_and_DT_Algorithm_copy.R")  # Ensure this file is in the same directory

results <- data.frame(
  Run = integer(),
  Lambda = numeric(),
  Length = numeric(),
  Cost = numeric(),
  NumDisambigs = integer()
)

lambda <- 0.5

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

saveRDS(results, file = "data_lambda_0_5.rds")
