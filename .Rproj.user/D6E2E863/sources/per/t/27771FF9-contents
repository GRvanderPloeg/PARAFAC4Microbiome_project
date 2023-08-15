# Simulate microbiome count data with metaSPARSim

# From documentation
# N_sample_cond_A <- 5
# lib_size_cond_A <- 1000
# 
# cond_A_sim_param <- list()
# cond_A_sim_param$intensity <- c(100, 5000, 750, 70, 1000)
# cond_A_sim_param$variability <- c(0.1, 0.5, 0.3, 0.05, 0.3)
# cond_A_sim_param$lib_size <- rep (lib_size_cond_A, N_sample_cond_A)
# 
# fold_change_multiplier <- c(5, 0.04, 0.67, 1, 2)
# 
# N_sample_cond_B <- 5
# lib_size_cond_B <- 1000
# 
# cond_B_sim_param <- list()
# cond_B_sim_param$intensity <- cond_A_sim_param$intensity * fold_change_multiplier # apply fold-change
# cond_B_sim_param$variability <- cond_A_sim_param$variability
# cond_B_sim_param$lib_size <- rep (lib_size_cond_B, N_sample_cond_B)
# 
# metaSPARSim_param <- list(cond_A = cond_A_sim_param, cond_B = cond_B_sim_param)
# metaSPARSim_result <- metaSPARSim(metaSPARSim_param)

pattern1 = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
pattern2 = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)

noise = 0.0
averageLibSize = 10000
stdLibSize = 0
params = list()
simResults = list()

for(i in 1:9){
  for(j in 1:20){
    param = list()
    
    if(j <= 10){
      param$intensity = c(rep(pattern1[i], 50), rep(pattern2[i], 50))
    }
    else{
      param$intensity = c(rep(pattern2[i], 50), rep(pattern1[i], 50))
    }
    
    param$variability = rep(noise, 100)
    param$lib_size = rnorm(1, mean=averageLibSize, sd=stdLibSize)
    params[[j]] = param
  }
  names(params) = paste0("subject_", 1:20)
  simResults[[i]] = metaSPARSim(params)
}

result = data.frame(rbind(t(simResults[[1]]$counts), t(simResults[[2]]$counts), t(simResults[[3]]$counts), t(simResults[[4]]$counts), t(simResults[[5]]$counts),
                          t(simResults[[6]]$counts), t(simResults[[7]]$counts), t(simResults[[8]]$counts), t(simResults[[9]]$counts)))
result$timepoint = rep(1:9, each=20)
write.table(result, "./simData.csv", sep=",", row.names=FALSE, col.names=TRUE)
