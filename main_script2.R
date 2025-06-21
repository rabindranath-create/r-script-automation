# main_script2.R

# Make sure working directory is the same as the script location (implicitly handled in GitHub Actions)
print(getwd())

# Load the helper script
source("RD_and_DT_Algorithm_copy.R")  # Ensure this file is in the same directory

results_0 <- data.frame(
  Run = integer(),
  Lambda = numeric(),
  Length = numeric(),
  Cost = numeric(),
  NumDisambigs = integer()
)

lambda <- 0

for (i in 1:100) {
  set.seed(i)
  obs_gen_para <- c(gamma = 0.3, d = 5, noPoints = 25)
  result <- ACS_Alg_C(obs_gen_para, k = 2, lambda)
  
  results_0[i, ] <- list(
    Run = i,
    Lambda = lambda,
    Length = result$Length_total,
    Cost = result$Cost_total,
    NumDisambigs = length(result$Disambiguate_state)
  )
}

saveRDS(results_0, file = "data_25_2_0.rds")




results_05 <- data.frame(
  Run = integer(),
  Lambda = numeric(),
  Length = numeric(),
  Cost = numeric(),
  NumDisambigs = integer()
)

lambda <- 0.5

for (i in 1:100) {
  set.seed(100+i)
  obs_gen_para <- c(gamma = 0.3, d = 5, noPoints = 25)
  result <- ACS_Alg_C(obs_gen_para, k = 2, lambda)
  
  results_05[i, ] <- list(
    Run = i,
    Lambda = lambda,
    Length = result$Length_total,
    Cost = result$Cost_total,
    NumDisambigs = length(result$Disambiguate_state)
  )
}

saveRDS(results_05, file = "data_25_2_05.rds")




results_1 <- data.frame(
  Run = integer(),
  Lambda = numeric(),
  Length = numeric(),
  Cost = numeric(),
  NumDisambigs = integer()
)

lambda <- 1

for (i in 1:100) {
  set.seed(200+i)
  obs_gen_para <- c(gamma = 0.3, d = 5, noPoints = 25)
  result <- ACS_Alg_C(obs_gen_para, k = 2, lambda)
  
  results_1[i, ] <- list(
    Run = i,
    Lambda = lambda,
    Length = result$Length_total,
    Cost = result$Cost_total,
    NumDisambigs = length(result$Disambiguate_state)
  )
}

saveRDS(results_1, file = "data_25_2_1.rds")



results_15 <- data.frame(
  Run = integer(),
  Lambda = numeric(),
  Length = numeric(),
  Cost = numeric(),
  NumDisambigs = integer()
)

lambda <- 1.5

for (i in 1:100) {
  set.seed(300+i)
  obs_gen_para <- c(gamma = 0.3, d = 5, noPoints = 25)
  result <- ACS_Alg_C(obs_gen_para, k = 2, lambda)
  
  results_15[i, ] <- list(
    Run = i,
    Lambda = lambda,
    Length = result$Length_total,
    Cost = result$Cost_total,
    NumDisambigs = length(result$Disambiguate_state)
  )
}

saveRDS(results_15, file = "data_25_2_15.rds")



results_2 <- data.frame(
  Run = integer(),
  Lambda = numeric(),
  Length = numeric(),
  Cost = numeric(),
  NumDisambigs = integer()
)

lambda <- 2

for (i in 1:100) {
  set.seed(400+i)
  obs_gen_para <- c(gamma = 0.3, d = 5, noPoints = 25)
  result <- ACS_Alg_C(obs_gen_para, k = 2, lambda)
  
  results_2[i, ] <- list(
    Run = i,
    Lambda = lambda,
    Length = result$Length_total,
    Cost = result$Cost_total,
    NumDisambigs = length(result$Disambiguate_state)
  )
}

saveRDS(results_2, file = "data_25_2_2.rds")



results_25 <- data.frame(
  Run = integer(),
  Lambda = numeric(),
  Length = numeric(),
  Cost = numeric(),
  NumDisambigs = integer()
)

lambda <- 2.5

for (i in 1:100) {
  set.seed(500+i)
  obs_gen_para <- c(gamma = 0.3, d = 5, noPoints = 25)
  result <- ACS_Alg_C(obs_gen_para, k = 2, lambda)
  
  results_25[i, ] <- list(
    Run = i,
    Lambda = lambda,
    Length = result$Length_total,
    Cost = result$Cost_total,
    NumDisambigs = length(result$Disambiguate_state)
  )
}

saveRDS(results_25, file = "data_25_2_25.rds")





results_3 <- data.frame(
  Run = integer(),
  Lambda = numeric(),
  Length = numeric(),
  Cost = numeric(),
  NumDisambigs = integer()
)

lambda <- 3

for (i in 1:100) {
  set.seed(600+i)
  obs_gen_para <- c(gamma = 0.3, d = 5, noPoints = 25)
  result <- ACS_Alg_C(obs_gen_para, k = 2, lambda)
  
  results_3[i, ] <- list(
    Run = i,
    Lambda = lambda,
    Length = result$Length_total,
    Cost = result$Cost_total,
    NumDisambigs = length(result$Disambiguate_state)
  )
}

saveRDS(results_3, file = "data_25_2_3.rds")



results_35 <- data.frame(
  Run = integer(),
  Lambda = numeric(),
  Length = numeric(),
  Cost = numeric(),
  NumDisambigs = integer()
)

lambda <- 3.5

for (i in 1:100) {
  set.seed(700+i)
  obs_gen_para <- c(gamma = 0.3, d = 5, noPoints = 25)
  result <- ACS_Alg_C(obs_gen_para, k = 2, lambda)
  
  results_35[i, ] <- list(
    Run = i,
    Lambda = lambda,
    Length = result$Length_total,
    Cost = result$Cost_total,
    NumDisambigs = length(result$Disambiguate_state)
  )
}

saveRDS(results_35, file = "data_25_2_35.rds")




results_4 <- data.frame(
  Run = integer(),
  Lambda = numeric(),
  Length = numeric(),
  Cost = numeric(),
  NumDisambigs = integer()
)

lambda <- 4

for (i in 1:100) {
  set.seed(800+i)
  obs_gen_para <- c(gamma = 0.3, d = 5, noPoints = 25)
  result <- ACS_Alg_C(obs_gen_para, k = 2, lambda)
  
  results_4[i, ] <- list(
    Run = i,
    Lambda = lambda,
    Length = result$Length_total,
    Cost = result$Cost_total,
    NumDisambigs = length(result$Disambiguate_state)
  )
}

saveRDS(results_4, file = "data_25_2_4.rds")
