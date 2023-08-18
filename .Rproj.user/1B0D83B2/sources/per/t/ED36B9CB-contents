# Simulate microbiome count data with metaSPARSim
library(metaSPARSim)
library(rTensor)

rlang::global_entrace()
options(rlang_backtrace_on_warning_report = "full")
options(rlang_backtrace_on_error_report = "full")

simulateCountData = function(subjectLoadings, featureLoadings, timeLoadings, numSubjectsPerGroup, numFeaturesPerGroup,
                             relativeNoise, avgLibSize, stdLibSize){
  
  set.seed(1)
  # Number of components is assumed to be the same for all modes, otherwise
  # PARAFAC would not be an appropriate model to use anyway.
  numTimepoints = nrow(timeLoadings)
  numComponents = nrow(subjectLoadings)
  numSubjects = sum(numSubjectsPerGroup)
  numFeatures = sum(numFeaturesPerGroup)
  numSubjectGroups = length(numSubjectsPerGroup)
  numFeatureGroups = length(numFeaturesPerGroup)
  
  # Create trilinear structure
  abundanceTable = list()
  coreArray = subjectLoadings %*% t(khatri_rao(timeLoadings, featureLoadings)) # matricized array
  
  for(i in 1:numTimepoints){
    abundanceTable[[i]] = coreArray[,(2*i-1):(2*i)]
  }
  
  # Correction for negativity
  # I think this is wrong, but I don't know how to fix it yet.
  # For now, avoid using negative loadings.
  
  # correctionTerm = abs(min(unlist(abundanceTable))) + 0.01 # to avoid adding zeroes
  # for(i in 1:numTimepoints){
  #   abundanceTable[[i]] = abundanceTable[[i]] + correctionTerm
  # }
  
  # Prepare modelling parameters per timepoint and simulate.
  params = list()
  paramIterator = 1
  featureIntensities = rnbinom(numFeatures, size=0.033, mu=15.75)
  featureVariability = rnorm(numFeatures, mean=3.6, sd=1.5)

  for(i in 1:numTimepoints){
    for(j in 1:numSubjectGroups){
      param = list()
      param$intensity = c()
      param$variability = c()
      for(k in 1:numFeatureGroups){
        param$intensity = c(param$intensity, rep(abundanceTable[[i]][j,k], numFeaturesPerGroup[k]))
        #param$variability = c(param$variability, rep(abundanceTable[[i]][j,k], numFeaturesPerGroup[k]))
      }
      param$intensity = featureIntensities * param$intensity
      #param$variability = param$variability * relativeNoise
      param$variability = featureVariability
      param$lib_size = round(rnorm(numSubjectsPerGroup[j], mean=avgLibSize, sd=stdLibSize))
      params[[paramIterator]] = param
      paramIterator = paramIterator + 1
    }
  }
  names(params) = paste0("situation_", 1:(paramIterator-1))
  simResults = metaSPARSim(params)
  
  return(list(abundanceTable, simResults))
}

# Rows: subject, features; columns: components
subjectLoadings = rbind(c(0.1, 1.0),
                        c(0.1, 0.5)) 

featureLoadings = rbind(c(0.8, 0.1),
                        c(0.1, 0.5))

timeLoadings = t(rbind(c(0.1, 0.11, 0.2, 0.55, 0.3, 0.25, 0.2, 0.11, 0.1),
                     c(0.9, 0.8, 0.65, 0.55, 0.5, 0.44, 0.3, 0.25, 0.65)))

numSubjectsPerGroup = c(10, 10)
numFeaturesPerGroup = c(20, 20)
relativeNoise = 0.0
avgLibSize = 10000
stdLibSize = 0

outcome = simulateCountData(subjectLoadings, featureLoadings, timeLoadings, numSubjectsPerGroup, numFeaturesPerGroup,
                            relativeNoise, avgLibSize, stdLibSize)

# Construct additional metadata (especially useful in unbalanced cases)
subjectMetadata = c()
for(i in 1:numTimepoints){
  for(j in 1:length(numSubjectsPerGroup)){
    subjectMetadata = c(subjectMetadata, rep(LETTERS[j], numSubjectsPerGroup[j]))
  }
}
subjectMetadata = cbind(1:sum(numSubjectsPerGroup), subjectMetadata)

featureMetadata = c()
for(i in 1:length(numFeaturesPerGroup)){
  featureMetadata = c(featureMetadata, rep(LETTERS[i], numFeaturesPerGroup[i]))
}
featureMetadata = cbind(1:sum(numFeaturesPerGroup), featureMetadata)

# Save data
write.table(subjectMetadata, "./subjectMetadata.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(featureMetadata, "./featureMetadata.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(subjectLoadings, "./subjectLoadings.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(featureLoadings, "./featureLoadings.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(timeLoadings, "./timeLoadings.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(t(outcome[[2]]$counts), "./simData.csv", sep=",", row.names=FALSE, col.names=TRUE)
write.table(do.call(cbind, outcome[[1]]), "./coreArray.csv", sep=",", row.names=FALSE, col.names=TRUE)

# Simulate based on real data
# data(HMP)
# params = estimate_parameter_from_data(data, data_norm, indpop,perc_not_zeros=.2)
# names(params)<-names(indpop)
# sim_data = metaSPARSim(params)
# write.table(sim_data$counts, "./simHMPData.csv", sep=",", row.names=FALSE, col.names=TRUE)
