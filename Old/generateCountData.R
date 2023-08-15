# Generate count data
library(tidyverse)
library(patchwork)
library(timeOmics)
library(tscount)
set.seed(1234)

# Kodikara et al. - clustering simulated data based on time profiles
# RAW DATA
c1.0 <- c(0, 0.5,1,1.1,1.2,1.8,2.5,5,9)
c2.0<-c(-1, 8, -1, 0, 0.5, -2.5, -3, 2, 2)
c3.0 <-  c(-2,4, 8, 6,4.5,4,3.9, 3, 1)
c4.0 <- c(2, -5, 1, -4, 0.5, -3.5, 0, -3,0.5)

lst <-list()
lst[["c1.0"]]<-c1.0
lst[["c2.0"]]<-c2.0
lst[["c3.0"]]<-c3.0
lst[["c4.0"]]<-c4.0

for(i in 1:49){
  a<-round(rnorm(1,2,1),1)  # get random value with mean 2 and sd 1 and round it to one digit after the comma
  b<-round(rnorm(1,2,1),1)  # i.e. this is generally 2.0 but can be 1.0 or 3.0
  z1<-(c1.0+a)*abs(b+1)     # add it as offset + multiplicative factor
  z2<-(c2.0+a)*abs(b+1)
  z3<-(c3.0+a)*abs(b+1)
  z4<-(c4.0+a)*abs(b+1)
  lst[[paste0("c1.", i)]]<-assign(paste0("c1.", i),z1)
  lst[[paste0("c2.", i)]]<-assign(paste0("c2.", i),z2)
  lst[[paste0("c3.", i)]]<-assign(paste0("c3.", i),z3)
  lst[[paste0("c4.", i)]]<-assign(paste0("c4.", i),z4)
}
raw.data<-lst[order(names(lst))]

generate_LMMS_data <- function(raw_data, N_Ind, noise){
  data.gather <- raw_data %>% as.data.frame() %>% rownames_to_column("time") %>%
    gather(feature, value, -time)
  for(ind in 1:N_Ind){
    vect <- vector(length = nrow(data.gather), mode = "numeric")
    for(x in 1:length(vect)){
      vect[x] <- rnorm(1, mean = data.gather$value[x], sd = noise)
    }
    names.tmp <- colnames(data.gather)
    data.gather <- data.frame(data.gather, vect)
    colnames(data.gather) <- c(names.tmp, LETTERS[ind])
  }
  sim_data <- data.gather %>% dplyr::select(-c(value)) %>%
    gather(ind, value, -c(time, feature)) %>%
    mutate(sample = paste0(ind, "_", time)) %>%
    dplyr::select(feature, value, sample) %>%
    spread(feature, value) %>%
    column_to_rownames("sample") %>%
    as.matrix()
  return(sim_data)
}


sim_lmm<-function(rawData,noise, nInd){
  s1<-generate_LMMS_data(rawData,nInd,noise)
  time <- rep(1:9, nInd)
  
  lmms.output <- lmms::lmmSpline(data = s1, time = time,
                                 sampleID = rownames(s1), deri = FALSE,
                                 basis = "p-spline", numCores = 4, timePredict = 1:9,
                                 keepModels = TRUE)
  modelled.data <- t(slot(lmms.output, 'predSpline'))
  return(modelled.data)
}

rawData<-raw.data
nInd<-5
noise<-c(0.5,1.5,3)
n_iter <- 10 # Number of iterations of the loop

for (i in noise){
  list_of_frames <- replicate(n_iter, data.frame())
  for(j in 1:n_iter) {
    print(j)
    df_Scenario1<-sim_lmm(rawData,i,nInd)
    list_of_frames[[j]] <-df_Scenario1
  }
  #saveRDS(list_of_frames,file =  paste0("Data/clusData_",i,".Rdata")) 
}

# Kodikara et al. - differential abundance analysis based on simulated data
simulate<-function(mod,para, nIndiv, nTime,disper, nTaxa, nSc1, nSc2, nSc3, meta){
  
  TAXA<-setNames(data.frame(matrix(ncol = nTaxa, 
                                   nrow = nIndiv*nTime)),paste0("Taxa_",1:nTaxa))
  
  for(i in 1:nTaxa){
    if(i<=nSc1){#Time
      k=1
    }else if(i<=nSc1+nSc2){#Group
      k=2
    }else if(i<=nSc1+nSc2+nSc3){#Time+Group+Time*Group
      k=3
    }else {#No
      k=4
    }
    count<-c()
    for(j in 1:nIndiv){
      t<-c(tsglm.sim(n=nTime, param = para[[k]], model=model, 
                     xreg=matrix(c(1:nTime,
                                   rep(meta$Group[meta$Indiv==j][1],nTime),
                                   1:nTime*rep(meta$Group[meta$Indiv==j][1],nTime)),ncol=3), 
                     link="identity",
                     distr="nbinom", distrcoefs=c(size=1/disper))$ts)
      count<-c(count,t)
    }
    while(!any(count==0)){#At least one zero (Needed to run ZIGMM)
      count<-c()
      for(j in 1:nIndiv){
        t<-c(tsglm.sim(n=nTime, param = para[[k]], model=model, 
                       xreg=matrix(c(1:nTime,rep(meta$Group[meta$Indiv==j][1],nTime),
                                     1:nTime*rep(meta$Group[meta$Indiv==j][1],nTime)),ncol=3),
                       link="identity",
                       distr="nbinom", distrcoefs=c(size=1/disper))$ts)
        count<-c(count,t)
      } 
    }
    TAXA[,i]<-count
  }
  
  
  return(TAXA)
}

set.seed(1234)
model <- list(past_obs=1) #Only 1 AR parameter
disp_1<-0.1
param_1 <- list(list(intercept=runif(1,0,5), past_obs=0.04,  xreg=c(1.5,0,0)),
                list(intercept=runif(1,0,5), past_obs=0.04,  xreg=c(0,13,0)),
                list(intercept=runif(1,0,5), past_obs=0.04,  xreg=c(1.5,13,5)),
                list(intercept=runif(1,0,5), past_obs=0.04,  xreg=c(0,0,0)))

nIndiv=20;nTime=10
metaDF<-data.frame(Time=rep(c(1:nTime),nIndiv),
                   Indiv=rep(1:nIndiv,each=nTime),
                   Group=rep(0:1,each=nIndiv*nTime/2))

n_iter <- 10 # Number of iterations of the loop
list_of_frames <- replicate(n_iter, data.frame())
# Initializes the progress bar
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = n_iter, # Maximum value of the progress bar
                     style = 3,    # Progress bar style
                     width = 50,   # Progress bar width
                     char = "=")   # Character used to create the bar

for(i in 1:n_iter) {
  
  #---------------------
  # Code to be executed
  #---------------------
  df_Scenario1<-simulate(model,param_1,nIndiv,nTime, disp_1, 300, 10,10,10,metaDF)
  list_of_frames[[i]] <-df_Scenario1
  #---------------------
  
  # Sets the progress bar to the current state
  setTxtProgressBar(pb, i)
}

close(pb) # Close the connection

c_Sc1<-lapply(list_of_frames, function(x) {
  mutate(x,Library_size=rowSums(x))
})

ra_Sc1<-lapply(c_Sc1, function(x) {
  x[,-301]/x$Library_size
})

# My own work
realData = read.csv("./Martino2021/data/Halfvarson-IBD-Qiita-1629/table-matched.tsv", sep="\t", header=FALSE, skip=1)
hist(unlist(c(realData[,2:135])))

lambda_t = beta_0 + beta_1 * y_t_minus_1 + n_1*t + n_2*X + n_3*t*X

n_1 = 0 # time effect parameter
n_2 = 0 # group effect parameter
n_3 = 0 # time x group effect parameter

beta_0 = 1 # intercept parameter
beta_1 = 0.6 # auto-regression parameter
disp = 0.4 # dispersion parameter