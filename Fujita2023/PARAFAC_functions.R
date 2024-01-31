# PARAFAC functions
# G.R. van der Ploeg
# University of Amsterdam, 2022

exclude_names = c("sample", "subject", "niche", "group", "RFgroup", "asv", "Species")
#exclude_names = c("subject", "RFgroup")

create_and_plot_parafac_model = function(df, num_components, title, featureName, featureMetaData, featureSorting, featureLabel){
  pfac = create_parafac_model(df, num_components)
  plot_parafac_model(pfac, title, featureName, featureMetaData, featureSorting, featureLabel)
}

create_parafac_model = function(df, num_components) {
  numeric_data = df %>% select(-any_of(exclude_names)) %>% select(-visit)
  
  if(min(df$visit) > 1){
    df$visit = df$visit - 1
  }
  
  universal_subjects = unique(df$subject)
  for(i in 1:max(unique(df$visit))){
    subset = df %>% filter(visit==i) %>% select(subject) %>% pull
    universal_subjects = intersect(universal_subjects, subset)
  }
  mydim <- c(length(universal_subjects), ncol(numeric_data), max(unique(df$visit)))
  # Number of individuals, number of features, number of visits
  
  preprocess_result = preprocess_pfac(df, mydim, num_components, universal_subjects, exclude_names)
  X = preprocess_result[[1]]
  grouped_data = preprocess_result[[2]]
  
  cl <- makeCluster(detectCores())
  ce <- clusterEvalQ(cl, library(multiway))
  clusterSetRNGStream(cl, 1)
  system.time({pfac <- parafac(X, nfac=num_components, nstart=1000, parallel=TRUE, cl=cl)})
  pfac
  stopCluster(cl)
  return(c(pfac, grouped_data))
}

plot_parafac_model = function(pfac, title, featureName, featureMetaData, featureSorting, featureLabel){
  num_components = ncol(pfac$A)
  pfac_output = postprocess_pfac(pfac, featureName, featureMetaData, featureSorting)
  feature_vector = pfac_output[[1]]
  individual_vector = pfac_output[[2]]
  time_vector = pfac_output[[3]]
  
  print(feature_vector)
  print(individual_vector)
  print(time_vector)
  
  plots = vector("list", length=3*num_components)
  
  for(i in 1:num_components){
    
    if (i == num_components){ # Handle bottom plots where x labels are shown
      plots[[3*(i-1)+1]] = feature_vector %>% ggplot(aes_string(x="number",y=paste0("value.",i),fill=featureSorting[1])) + geom_bar(stat="identity") + xlab(featureLabel) + ylab(paste("Comp. ", i, sep="")) + theme(legend.position="none")
      plots[[3*(i-1)+2]] = individual_vector %>% ggplot(aes_string(x="number",y=paste0("value.",i),fill="RFgroup")) + geom_bar(stat="identity") + xlab("individual") + theme(legend.position="none") + ylab("") + scale_fill_manual(values=RFcolours)
      plots[[3*(i-1)+3]] = time_vector %>% ggplot(aes_string(x="number",y=paste0("value.",i))) + geom_line() + xlab("visit") + ylab("")
    }
    else{
      plots[[3*(i-1)+1]] = feature_vector %>% ggplot(aes_string(x="number",y=paste0("value.",i),fill=featureSorting[1])) + geom_bar(stat="identity") + xlab("") + ylab(paste("Comp.",i)) + theme(legend.position="none")
      plots[[3*(i-1)+2]] = individual_vector %>% ggplot(aes_string(x="number",y=paste0("value.",i),fill="RFgroup")) + geom_bar(stat="identity") + xlab("") + theme(legend.position="none") + ylab("") + scale_fill_manual(values=RFcolours)
      plots[[3*(i-1)+3]] = time_vector %>% ggplot(aes_string(x="number",y=paste0("value.",i))) + geom_line() + xlab("") + ylab("")
    }
  }
  
  plot = ggarrange(plotlist=plots, nrow=num_components, ncol=3)
  annotate_figure(plot, top=text_grob(title))
}

preprocess_pfac = function(df, mydim, num_components, universal_subjects, exclude_names){
  Amat <- matrix(0, nrow = mydim[1], ncol = num_components)
  Bmat <- matrix(0, nrow = mydim[2], ncol = num_components)
  Cmat <- matrix(0, nrow = mydim[3], ncol = num_components)
  Xmat <- tcrossprod(Amat, krprod(Cmat, Bmat))
  Xmat <- array(Xmat, dim = mydim)
  X <- Xmat
  
  grouped_data = df %>% filter(subject %in% universal_subjects) %>% group_by(visit) %>% nest() %>% arrange(visit)
  
  for (i in 1:max(unique(df$visit))) {
    data = grouped_data$data[i][[1]] %>% arrange(subject)
    cube_slice = grouped_data$data[i][[1]] %>% arrange(subject) %>% select(-any_of(exclude_names))
    cube_slice = sweep(cube_slice, 2, colMeans(cube_slice), FUN="-")
    #cube_slice = sweep(cube_slice, 1, apply(cube_slice, 1, sd), FUN="/")
    grouped_data$data[i][[1]] = cbind(cube_slice, data %>% select(exclude_names)) %>% as_tibble()
    cube_slice = cube_slice %>% as_tibble()

    X[,,i] = cube_slice %>% data.matrix(.) %>% as.numeric(.) %>% array(., dim=c(mydim[1],mydim[2]))
  }
  
  return(list(X, grouped_data))
}

postprocess_pfac = function(pfac, featureName, featureMetaData, featureSorting){
  feature_vector = cbind(pfac$B, colnames(pfac$data[[1]] %>% select(-any_of(exclude_names))))
  colnames(feature_vector) = c(paste0("value.", 1:ncol(pfac$B)), featureName)
  feature_vector = as.data.frame(feature_vector)
  for(i in 1:ncol(pfac$B)){
    colname = paste("value.", i, sep="")
    feature_vector[colname] = as.double(unlist(feature_vector[colname])) 
    }
  feature_vector = feature_vector %>% as_tibble() %>% left_join(featureMetaData, by=featureName) %>% arrange_at(featureSorting) %>% mutate(number=1:nrow(.))
  
  individual_vector = cbind(pfac$A, as.factor(pfac$data[[1]]$RFgroup), pfac$data[[1]]$subject)
  colnames(individual_vector) = c(paste0("value.", 1:ncol(pfac$A)), "RFgroup", "subject")
  individual_vector = as.data.frame(individual_vector)
  for(i in 1:ncol(pfac$A)){
    colname = paste("value.", i, sep="")
    individual_vector[colname] = as.double(unlist(individual_vector[colname])) 
  }
  individual_vector = individual_vector %>% as_tibble() %>% arrange(RFgroup, subject) %>% mutate(number=1:nrow(.))
  
  time_vector = pfac$C
  colnames(time_vector) = paste0("value.", 1:ncol(pfac$C))
  time_vector = as.data.frame(time_vector)
  for(i in 1:ncol(pfac$C)){
    colname = paste("value.", i, sep="")
    time_vector[colname] = as.double(unlist(time_vector[colname])) 
  }
  time_vector = time_vector %>% as_tibble() %>% mutate(number=1:nrow(.))
  
  return(list(feature_vector, individual_vector, time_vector))
}

SSE_check = function(df, max_factors, title) {
  resultSSE = vector(mode="numeric", length=max_factors)
  for (i in 1:max_factors) {
    model = create_parafac_model(df, i)
    resultSSE[i] = model$SSE
  }
  resultSSE %>% as_tibble() %>% mutate(num_components=1:max_factors) %>% ggplot(aes(x=num_components,y=value)) + geom_point() + ylab("Sum squared error") + xlab("Number of components") + ggtitle(title)
}

var_explained_check = function(df, max_factors, title) {
  resultVarEx = vector(mode="numeric", length=max_factors)
  for (i in 1:max_factors) {
    model = create_parafac_model(df, i)
    resultVarEx[i] = model$Rsq
  }
  resultVarEx %>% as_tibble() %>% mutate(num_components=1:max_factors) %>% ggplot(aes(x=num_components,y=value)) + geom_point() + ylim(0,1) + ylab("Variance Explained (%)") + xlab("Number of components") + ggtitle(title)
}

centeredOutPlot = function(df, visits, excludeNames, title){
  # This is the data that is centered out of the input matrix of the PARAFACs by centering
  # per ASV per timepoint
  numFeatures = ncol(df %>% select(-any_of(excludeNames)))
  plotData = matrix(0, nrow=numFeatures, ncol=length(visits)) %>% as_tibble()
  
  for(i in 1:length(visits)){
    visitNumber = visits[i]
    plotData[,i] =  df %>% filter(visit==visitNumber) %>% select(-any_of(excludeNames)) %>% colMeans()
  }
  colnames(plotData) = paste0("visit_", visits)
  plotData$ASV = colnames(df %>% select(-any_of(excludeNames)))
  plotData = plotData %>% pivot_longer(-ASV, names_to="visitNumber")
  plotData$visit = rep(1:7, numFeatures)
  plotData$day = rep(c(-14,0,2,5,9,14,21), numFeatures)
  
  plotData %>% ggplot(aes(x=day,y=value,group=ASV,color=ASV)) + geom_line() + theme(legend.position = "none") + geom_vline(xintercept=0,color="red") + geom_vline(xintercept=14,color="red") + xlab("Time [day]") + ylab("Mean value (log transformed)") + ggtitle(title)
 
}