# Summary functions for quick normalization checks
library(Polychrome)
library(lawstat)
library(tidyverse)

dataOverviewScatter <- function(df, row_lab, col_lab, titleA, titleB, titleC, titleD){
  
  row_num = 1:nrow(df)
  col_num = 1:ncol(df)
  robustSum = function(x) sum(x, na.rm=TRUE)
  robustStd = function(x) sd(x, na.rm=TRUE)
  
  # Row sum
  plotA = apply(df, 1, robustSum) %>% as_tibble() %>% mutate(row_number=row_num) %>% ggplot(aes(x=row_number,y=value)) + geom_point() + theme(axis.text.x=element_blank()) + ggtitle(titleA) + xlab(row_lab) + ylab("Row sum")
  
  # Col sum
  plotB = apply(df, 2, robustSum) %>% as_tibble() %>% mutate(col_number=col_num) %>% ggplot(aes(x=col_number,y=value)) + geom_point() + theme(axis.text.x=element_blank()) + ggtitle(titleB) + xlab(col_lab) + ylab("Column sum")
  
  # Row std
  plotC = apply(df, 1, robustStd) %>% as_tibble() %>% mutate(row_number=row_num) %>% ggplot(aes(x=row_number,y=value)) + geom_point() + theme(axis.text.x=element_blank()) + ggtitle(titleC) + xlab(row_lab) + ylab("Row std")
  
  # Col std
  plotD = apply(df, 2, robustStd) %>% as_tibble() %>% mutate(col_number=col_num) %>% ggplot(aes(x=col_number,y=value)) + geom_point() + theme(axis.text.x=element_blank()) + ggtitle(titleD) + xlab(col_lab) + ylab("Column std")
  
  ggarrange(plotA, plotB, plotC, plotD, ncol=2, nrow=2, labels=LETTERS[1:4])
}

dataOverviewHist <- function(df, row_lab, col_lab, titleA, titleB, titleC, titleD){
  
  row_num = 1:nrow(df)
  col_num = 1:ncol(df)
  robustSum = function(x) sum(x, na.rm=TRUE)
  robustStd = function(x) sd(x, na.rm=TRUE)
  
  # Row sum
  plotA = apply(df, 1, robustSum) %>% as_tibble() %>% ggplot(aes(x=value)) + geom_histogram(bins=30) + ggtitle(titleA) + xlab("Row sum") + ylab("Counts")
  
  # Col sum
  plotB = apply(df, 2, robustSum) %>% as_tibble() %>% ggplot(aes(x=value)) + geom_histogram(bins=30) + ggtitle(titleB) + xlab("Column sum") + ylab("Counts")
  
  # Row std
  plotC = apply(df, 1, robustStd) %>% as_tibble() %>% ggplot(aes(x=value)) + geom_histogram(bins=30) + ggtitle(titleC) + xlab("Row std") + ylab("Counts")
  
  # Col std
  plotD = apply(df, 2, robustStd) %>% as_tibble() %>% ggplot(aes(x=value)) + geom_histogram(bins=30) + ggtitle(titleD) + xlab("Row std") + ylab("Counts")
  
  ggarrange(plotA, plotB, plotC, plotD, ncol=2, nrow=2, labels=LETTERS[1:4])
}

heteroScedastic = function(df){
  colMeans = colMeans(df, na.rm=TRUE)
  colStds  = apply(df, 2, function(x) sd(x, na.rm=TRUE))
  
  result = cbind(colMeans, colStds)
  result %>% as_tibble() %>% ggplot(aes(x=colMeans, y=colStds)) + geom_point()
}

testNormality <- function(df){
  result = vector(length=ncol(df))
  
  for(i in 1:ncol(df)){
    testResult = shapiro.test(df[,i])
    result[i] = testResult$p.value
  }
  return(result)
}

testSymmetry <- function(df){
  result = vector(length=ncol(df))
  
  for(i in 1:ncol(df)){
    testResult = symmetry.test(df[,i], boot=FALSE)
    result[i] = testResult$p.value
  }
  return(result)
}

importPARAFAC = function(path, featureNames, dayVector, featureColumnOfInterest, numComponents){
  #  Output:
  #    [[1]]: Loadings in ID mode
  #    [[2]]: Loadings in feature mode
  #    [[3]]: Loadings in time mode
  #    [[4]]: Data as modeled, wide format
  #    [[5]]: Input data, wide format
  #    [[6]]: Data as modeled, long format
  #    [[7]]: Input data, long format
  #    [[8]]: Data as modeled, component 1, wide format
  #    [[9]]: Data as modeled, component 1, long format
  #   [[**]]: and so on until all components have a list item.
  
  id_mode = read.csv(paste0(path,"_mode1.csv"), header=FALSE)
  feature_mode = read.csv(paste0(path,"_mode2.csv"), header=FALSE)
  time_mode = read.csv(paste0(path,"_mode3.csv"), header=FALSE)
  
  if(dim(id_mode)[2] == (numComponents+1)){
    colnames(id_mode) = c(paste0("Component_", 1:numComponents), "subject")
    id_mode = id_mode %>% left_join(rf)
  }
  colnames(id_mode) = c(paste0("Component_", 1:numComponents), "subject", "RFgroup")
  id_mode = as_tibble(id_mode)
  
  colnames(feature_mode) = c(paste0("Component_", 1:numComponents), featureNames)
  feature_mode = as_tibble(feature_mode)
  
  colnames(time_mode) = c(paste0("Component_", 1:numComponents), "day")
  time_mode = as_tibble(time_mode) %>% mutate(days=dayVector)
  
  rawmodel = scan(paste0(path,"_model.csv"), sep=",")
  model_wide = matrix(0, nrow=nrow(id_mode), ncol=nrow(feature_mode)*nrow(time_mode))
  
  for(i in 1:nrow(id_mode)){
    model_wide[i,] = rawmodel[((nrow(feature_mode)*nrow(time_mode)*(i-1))+1):(nrow(feature_mode)*nrow(time_mode)*i)]
  }
  column_names = paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t1")
  
  for(i in 2:length(dayVector)){
    column_names = c(column_names, paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t", i))
  }
  colnames(model_wide) = column_names
  
  model_wide = as_tibble(model_wide)
  
  rawinput = scan(paste0(path,"_input.csv"), sep=",")
  input_wide = matrix(0, nrow=nrow(id_mode), ncol=nrow(feature_mode)*nrow(time_mode))
  
  for(i in 1:nrow(id_mode)){
    input_wide[i,] = rawinput[((nrow(feature_mode)*nrow(time_mode)*(i-1))+1):(nrow(feature_mode)*nrow(time_mode)*i)]
  }
  column_names = paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t1")
  
  for(i in 2:length(dayVector)){
    column_names = c(column_names, paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t", i))
  }
  colnames(input_wide) = column_names
  
  input_wide = as_tibble(input_wide)
  model_long = make_longer(model_wide, id_mode$subject)
  input_long = make_longer(input_wide, id_mode$subject)
  
  result = list(id_mode, feature_mode, time_mode, model_wide, input_wide, model_long, input_long)
  listIterator = 8
  
  for(i in 1:numComponents){
    modeled_component_raw = scan(paste0(path,"_component_", i, ".csv"), sep=",")
    modeled_component_wide = matrix(0, nrow=nrow(id_mode), ncol=nrow(feature_mode)*nrow(time_mode))
    
    for(i in 1:nrow(id_mode)){
      modeled_component_wide[i,] = modeled_component_raw[((nrow(feature_mode)*nrow(time_mode)*(i-1))+1):(nrow(feature_mode)*nrow(time_mode)*i)]
    }
    column_names = paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t1")
    
    for(i in 2:length(dayVector)){
      column_names = c(column_names, paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t", i))
    }
    colnames(modeled_component_wide) = column_names
    
    modeled_component_wide = as_tibble(modeled_component_wide)
    modeled_component_long = make_longer(modeled_component_wide, id_mode$subject)
    
    result[[listIterator]] = modeled_component_wide
    listIterator = listIterator + 1
    result[[listIterator]] = modeled_component_long
    listIterator = listIterator + 1
  }
  
  return(result)
}

importNPLS = function(path, featureNames, dayVector, featureColumnOfInterest){
  id_mode = read.csv(paste0(path,"_individual_mode.csv"), header=FALSE)
  feature_mode = read.csv(paste0(path,"_feature_mode.csv"), header=FALSE)
  time_mode = read.csv(paste0(path,"_time_mode.csv"), header=FALSE)
  numComponents = ncol(time_mode)
  
  colnames(id_mode) = c(paste0("Component_", 1:numComponents), "subject", "RFgroup")
  id_mode = as_tibble(id_mode)
  
  colnames(feature_mode) = c(paste0("Component_", 1:numComponents), featureNames)
  feature_mode = as_tibble(feature_mode)
  
  colnames(time_mode) = c(paste0("Component_", 1:numComponents))
  time_mode = as_tibble(time_mode) %>% mutate(days=dayVector)
  
  model = read.csv(paste0(path,"_model.csv"), header=FALSE)
  dim(model)
  
  ypred = read.csv(paste0(path, "_ypred.csv"), header=FALSE)
  colnames(ypred) = c(paste0("Component_", 1:numComponents), "subject", "y")
  ypred = as_tibble(ypred)
  
  
  column_names = paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t1")
  for(i in 2:length(dayVector)){
    column_names = c(column_names, paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t", i))
  }
  column_names
  colnames(model) = column_names
  
  model = as_tibble(model)
  
  return(list(id_mode, feature_mode, time_mode, model, ypred))
}

make_longer = function(data, subjects){
  df = data %>% mutate(subject=subjects) %>% pivot_longer(-subject, names_to=c("asv", "visit"), names_sep="_t") %>% as.data.frame()
  df$visit = as.integer(df$visit)
  df$value = as.double(df$value)
  df = df %>% as_tibble() %>% pivot_wider(id_cols=c("subject","visit"), names_from=asv, values_from=value) %>% select(-subject,-visit)
  return(df)
}