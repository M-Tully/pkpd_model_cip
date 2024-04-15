rm(list=ls())

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("XXX")

library(gridExtra)
library(MASS)
library(extraDistr)
library(dplyr)
library(ggplot2)
library(dplyr)

# 1. PD model

# Initialisation function, creates a discretised normal distribution of parasite age
PKPD_model_initialisation = function(nStages, 
                                     ipl, 
                                     iplMu, 
                                     iplSigma,
                                     pmf ) {
  pAge = sapply(X = seq(1,nStages), mean = iplMu, sd = iplSigma, FUN = dnorm);
  
  x0 = ipl * (pAge)/sum(pAge);
  
  return(x0);
}

# Parasite model
PD_model = function(tMax, # Maximum simulation time
                    nStages, # number of hours of lifecycle of parasite (usually 48)
                    kill_window, # killing window of a drug; a vector with values of elements 0: no effect or 1: effect; e.g. effective at all ages: kill_window = rep(1, nStages)
                    ipl, # initial parasite load; e.g. 1e9
                    iplMu, # mean of the age distribution of parasites (between 1 and 48)
                    iplSigma, # spread of the age distribution of parasites (between 1 and 48)
                    pmf, # Parasite multiplication factor
                    Emax, # Maximum killing effect (between 0 and 1)
                    EC50, # Drug concentration at which E = Emax/2
                    gamma, # Sigmoidicity of concentration-effect relationship
                    conc, # drug concentration at effect site over hourly time points
                    seed,
                    tseq #what age sequestered (26)
                    
) {
  
  # Circulating parasites
  cirCount = rep(NA, tMax) %>% as.vector()
  # Total number of parasites
  NCount = rep(NA, tMax) %>% as.vector()
  
  totCountAgeVecHold = rep(NA, nStages) %>% as.vector()
  survProbHold = rep(1.0, nStages);
  NandM = list()
  
  # Initialising the age distribution of parasites
  totCountAgeVec = PKPD_model_initialisation(nStages, 
                                             ipl, 
                                             iplMu, 
                                             iplSigma,
                                             pmf);
  
  for (i in 1:tMax) {
    
    time = i - 1.0;
    if (time == 0) {
      NCount[i] = sum(totCountAgeVec);
      cirCount[i] = sum(totCountAgeVec[1:tseq]);
    }
    else {
      
      # Drug effect
      E = (Emax * conc[i]^gamma) / (conc[i]^gamma + EC50^gamma);
      # Probability of survival of parasites of age a (a in [1, 48]) to the next hour
      survProb = 1 - E * kill_window;
      totCountAgeVecHold[2:nStages] = totCountAgeVec[1:(nStages - 1)] * survProb[1:(nStages - 1)];
      totCountAgeVecHold[1] = pmf * totCountAgeVec[nStages] * survProb[nStages];
      totCountAgeVec = min(1, floor(sum(totCountAgeVecHold))) * totCountAgeVecHold;
      
      # Total number of parasites
      NCount[i] = sum(totCountAgeVec);
      # Circulating parasites
      cirCount[i] = sum(totCountAgeVec[1:tseq]);
      
    }
  }
  NandM[[1]] = cirCount;
  NandM[[2]] = NCount;
  return(NandM);
}


# 2. Function to simulate profiles & prepare for stan model
pd_input_fun2 <- function(doses, tmax, nPt, t, ipl, iplMu, iplSig, nStages, nObs,
                          Qd, pmf,Emax , EC50 , gamma, tLag,
                          tDose, concFull, 
                          CL_est, Vc_est, Q_est, Vp_est, ka_est, seed, tseq, mu_noise){  

  # check we have parameters for each individual (8 in our project)
  stopifnot(length(doses)==nPt, length(ipl)==nPt, length(iplMu)==nPt, length(iplSig)==nPt)   
  
  # create id vector where each id is repeated for each observation each patient has
  start = c()
  end = c()
  
  start[1] = 1
  end[1] = nObs[1]
  if (nPt >= 2){
    for(i in 2:nPt){
      start[i] = 1 + sum(nObs[1:(i-1)])
      end[i] = sum(nObs[1:(i)])} }
  
  id_vec = c()
  for(i in 1:nPt){id_vec[start[i]:end[i]] = i}
  
  # create data frame for each patient/obs time combo, where we have current conc, dose, tdose and then can add # parasites
  stan_input_df = list(id = id_vec,
                       time = t,
                       conc = rep(0,length(t)),
                       dose = rep(0,length(t)),
                       tDose = rep(0,length(t)))
  for (i in 1:nPt){
    stan_input_df$dose[start[i]:end[i]] = doses[i]
    stan_input_df$tDose[start[i]:end[i]] = tDose[i]} # match dose data by id (i.e. functions as stan_input_df$dose[id_vec==i])
  
  # implement the PD model for each participant and save into relevant cells of stan_input_df data.frame
  for(i in 1:nPt){   
    temp <- PD_model(
      tMax = tmax[i], # Maximum simulation time
      nStages = nStages, # number of hours of lifecycle of parasite (usually 48)
      kill_window = Qd, # killing window of a drug; a vector with values of elements 0: no effect or 1: effect; e.g. effective at all ages: kill_window = rep(1, nStages)
      ipl = ipl[i], # initial parasite load; e.g. 1e9
      iplMu = iplMu[i], # mean of the age distribution of parasites (between 1 and 48)
      iplSigma = iplSig[i], # spread of the age distribution of parasites (between 1 and 48)
      pmf = pmf[i], # Parasite multiplication factor
      Emax = Emax, # Maximum killing effect (between 0 and 1)
      EC50 = EC50, # Drug concentration at which E = Emax/2
      gamma = gamma, # Sigmoidicity of concentration-effect relationship
      conc = 1000*concFull[i,],# drug concentration at effect site over hourly time points
      tseq=tseq #age stop circulating
    )
    
    stan_input_df$nCir[start[i]:end[i]] <- temp[[1]][t[id_vec==i]]
    stan_input_df$nCount[start[i]:end[i]] <- temp[[2]][t[id_vec==i]]
  }
  
  cir_matrix <- stan_input_df$nCir
  count_matrix <- stan_input_df$nCount
  
  #sample for noise level
  set.seed(seed)
  noise <-rlnorm(1,log(mu_noise),0.01)
  
  cir_matrix_e <- stan_input_df$nCir*exp(rnorm(length(t),0,noise))
  count_matrix_e <- stan_input_df$nCount*exp(rnorm(length(t),0,noise))
  
  #introduce censoring here
  cir_matrix[cir_matrix < 50000] <- 1
  count_matrix[count_matrix < 50000] <- 1
  cir_matrix_e[cir_matrix_e < 50000] <- 1
  count_matrix_e[count_matrix_e < 50000] <- 1
  
  
  stan_input_list2 =list(
    nStages = nStages,
    nPt =nPt,
    nObs = as.array(nObs),
    N = sum(nObs),
    t = t,
    tmax= as.array(tmax),
    ID = id_vec,
    id = rep(1:nPt,each=tmax), #if tmax is not the same for every participant, this needs to be done like id_vec
    kill_window = Qd, 
    start = as.array(start),
    end = as.array(end),
    NCount_measured = log(count_matrix),
    cirCount_measured = log(cir_matrix),
    NCount_measured_e = log(count_matrix_e),
    cirCount_measured_e = log(cir_matrix_e),
    dose = as.array(doses),
    tLag = as.array(tLag),
    tDose = as.array(tDose),
    CL = as.array(CL_est),
    Vc = as.array(Vc_est),
    Q = as.array(Q_est),
    Vp = as.array(Vp_est),
    ka = as.array(ka_est),
    Low = c(1500,1,1,5, 0.05, 0.5, 1), # prior bounds
    Upp=c(2100,24,14,50, 1, 30, 10),
    tseq=tseq,
    ipl = as.array(ipl),
    iplMu = as.array(iplMu),
    iplSig = as.array(iplSig),
    Emax = as.array(Emax),
    EC50 = as.array(EC50),
    Gamma = as.array(gamma),
    PMF = pmf
  )
  
  listt <- list(pkpd_list = stan_input_list2,
                pkpd_dataframe = stan_input_df)
  
  
  return(stan_input_list2)
  
}


# 3. Function to use seed to create profile, run and save stan output
pd_sims_complete <- function(filename, 
                             nPt=8 # from 1 to 8 
){
  
  #First import corresponding data:
  
  # This has the simulated concentration profiles
  PK_profile <- readRDS(filename)
  
  #Now extract the fitted parameters from the stan outputs
  load(paste(substring(filename,6,18),".stan", sep=""))
  print(paste("stan file", paste(substring(filename,6,18),".stan", sep=""), "loading", sep=" "))
  
  #-----------------------------------------------------------------
  #conc_est <- summary(PKMVN_, pars=c("conc_hat"))$summary[,1]
  (pk_CL <- summary(PKMVN_, pars=c("thetaInd[1,1]","thetaInd[2,1]",
                                   "thetaInd[3,1]","thetaInd[4,1]",
                                   "thetaInd[5,1]", "thetaInd[6,1]",
                                   "thetaInd[7,1]","thetaInd[8,1]"))$summary[,1])
  (pk_Vc <-summary(PKMVN_, pars=c("thetaInd[1,2]","thetaInd[2,2]",
                                  "thetaInd[3,2]","thetaInd[4,2]",
                                  "thetaInd[5,2]", "thetaInd[6,2]",
                                  "thetaInd[7,2]","thetaInd[8,2]"))$summary[,1])
  (pk_Q <- summary(PKMVN_, pars=c("thetaInd[1,3]","thetaInd[2,3]",
                                  "thetaInd[3,3]","thetaInd[4,3]",
                                  "thetaInd[5,3]", "thetaInd[6,3]",
                                  "thetaInd[7,3]","thetaInd[8,3]"))$summary[,1])
  (pk_Vp <-summary(PKMVN_, pars=c("thetaInd[1,4]","thetaInd[2,4]",
                                  "thetaInd[3,4]","thetaInd[4,4]",
                                  "thetaInd[5,4]", "thetaInd[6,4]",
                                  "thetaInd[7,4]","thetaInd[8,4]"))$summary[,1])
  (pk_ka <- summary(PKMVN_, pars=c("thetaInd[1,5]","thetaInd[2,5]",
                                   "thetaInd[3,5]","thetaInd[4,5]",
                                   "thetaInd[5,5]", "thetaInd[6,5]",
                                   "thetaInd[7,5]","thetaInd[8,5]"))$summary[,1])
  
  
  #Second set up distribution for eta, patient variations from transformed pop avg
  
  #parameters order: N_0, mu_0, sig_0, PMF, EMax, Ec50, Gamma (Gamma fixed @5 in McCarthy)
  
  omega = diag((c(0.2, 0.2, 0.2, 0.0242, 0.2, 0.2, 0.2))) #4 values req
  
  #Use input to make a 7 digit seed (not the same as pk sims seed)
  set.seed(substring(filename,14,18))
  seed = (floor(runif(1,1000000,9999999)))
  
  #create correlation matrix from L:
  L = matrix(c((runif(49,0,1))), nrow=7)
  
  #make L a lower triangular matrix by removing all values above the diagonal
  L[1,2:7] <- L[2,3:7] <- L[3,4:7] <- L[4,5:7] <- L[5,6:7] <- L[6,7] <- 0
  
  #Eta will be MVN mean 0 and covariance Sigma:
  Sigma = omega%*%((L)%*%t(L))%*%omega
  
  eta = mvrnorm(8,rep(0,7),Sigma)
  
  # * 5000 changes from conc to total count approximately
  pop_avg = c(1800, 2, 3, 13,0.23,15.1,5) #N avg from McCarthy et al, mu_0, sig_0 & PMF ranges from Dini et al

  #time '0' is now at 72h after innoculation
  a = c(2100,24,14,50,1,30,10)
  b = c(1500,1,1,5,0.05,0.5,1)
  
  #if we want each row to be a patient, repeat the vector 8 times
  a_matrix = matrix(rep(a,8), nrow=8, byrow=T)
  b_matrix = matrix(rep(b,8), nrow=8, byrow=T)
  
  phi_i <- matrix(rep((log((pop_avg - a)/(b-pop_avg))),8),nrow=8,byrow=T)  + eta
  
  pd_paraMVN_ = (b_matrix*exp(phi_i)+a_matrix)/(exp(phi_i)+1)
  
  print(pd_paraMVN_)
  
  
  #-------------
  
  #Remove first 3 days as no data here for McCarthy et al, set t0 to 72h
  
  data_pdMVN_ <- pd_input_fun2(doses=rep(10,nPt),
                               tmax=rep(288,nPt),
                               nPt=nPt,
                               # t=floor(rep(seq(1,288,length.out=22),nPt)),
                               t=rep(c(72,96,108,120,132, 144, 156, 168, 168+c(4,8,12,16,24,30,36,48,60,72,96,120)),nPt),
                               ipl=pd_paraMVN_[1:nPt,1], #blood volume (eg, 5L)
                               iplMu=pd_paraMVN_[1:nPt,2],
                               iplSig=pd_paraMVN_[1:nPt,3],
                               nStages=40,
                               nObs=rep(20,nPt),
                               Qd = c(rep(1,40)), #killwindow assumed max (??)
                               pmf = pd_paraMVN_[1:nPt,4], # Parasite multiplication factor
                               Emax = pd_paraMVN_[1:nPt,5], # Maximum killing effect (between 0 and 1)
                               EC50 = pd_paraMVN_[1:nPt,6], # Drug concentration at which E = Emax/2 #prev 1.5
                               gamma = pd_paraMVN_[1:nPt,7],
                               tDose= rep(168,nPt),
                               tLag = rep(0,nPt),
                               concFull = matrix(PK_profile$concFull[1:nPt,],nrow=nPt), #output from pk sims!
                               CL_est = pk_CL[1:nPt],
                               Vc_est = pk_Vc[1:nPt],
                               Q_est = pk_Q[1:nPt],
                               Vp_est = pk_Vp[1:nPt],
                               ka_est = pk_ka[1:nPt],
                               seed=seed,
                               tseq=24,
                               mu_noise = 0.6)
  
  
  #Create an image of the profile:
  y <- c()
  for (i in 1: (20*nPt)){y[i] = max(log10(exp(data_pdMVN_$cirCount_measured_e[i])/5000),0)}
  plot_data <- data.frame("t"=data_pdMVN_$t, "y"=y, "ID"=as.factor(data_pdMVN_$ID))
  (pd_plot <- ggplot(data=plot_data, aes(x=t, y=y, group=ID, col=ID)) + theme_light() + 
      geom_line(alpha=0.5)+
      geom_point(size=0.2) + 
      labs(title=paste("seed =",seed, sep=" ")) +
      xlab("Time (days since innoculation)") +
      ylab("Parasites/mL") + 
      scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c("0","10","100","1,000", "10,000", "100,000" ), limits=c(0,5.5)) +
      # scale_x_continuous(breaks=seq(0,144,24), labels=seq(3,9,1), limits=c(0,145)))
      scale_x_continuous(breaks=seq(72,24*9,24), labels=seq(3,9,1), limits=c(72,216)))
  
  #profile plots for 1st 400 ish do not include fixed X scale & diff dim
  ggsave(paste("pd_plot_LLOQ",seed,".png",sep=""), pd_plot, width=7.5,height = 3)
  print(pd_plot)
  
  # save simulated profile 
  saveRDS(data_pdMVN_,file=paste("simu_MVN_pd_LLOQ_",seed,".rdata", sep=""))
  
  #Sample with Stan
  PDMVN_ <- stan("pd_hierarchy.stan", data=data_pdMVN_, chains=3, iter=800, warmup=300)
  # control = list(adapt_delta=0.99, stepsize_jitter=0.8, max_treedepth=10))
  # 
  save(PDMVN_,file=paste("MVN_pd_LLOQ_",seed,".stan", sep=""))
  
}

# Create list of names for files with data
filenames = list.files(pattern="simu_MVN_pk_*")

# do not include files less than 1.5 KB, as most are 6KB,
# and small size is a good sign importing was incomplete or interrupted
filenames <- filenames[(is.na(file.info(filenames)$size)==FALSE & file.info(filenames)$size >= 1500)]
# repeat, checking STAN output is there & of proper size
filenames <- filenames[(is.na(file.info(paste(substring(filenames,6,19),"stan", sep=""))$size)==FALSE & file.info(paste(substring(filenames,6,19),"stan", sep=""))$size >= 30000000)]
# repeat, checking R input data is saved
filenames <- filenames[(is.na(file.info(paste("simu_",substring(filenames,6,19),"rdata", sep=""))$size)==FALSE & file.info(paste("simu_",substring(filenames,6,19),"rdata", sep=""))$size >= 8000)]


seed_fun <- function(filename){ set.seed(substring(filename,14,18))
seed = (floor(runif(1,1000000,9999999)))
return(seed)}
  pd_seeds <- sapply(filenames,FUN=seed_fun)

sapply(filenames, FUN=pd_sims_complete)

