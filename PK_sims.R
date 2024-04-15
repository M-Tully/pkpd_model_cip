# install any not currently installed
#install.packages(c("gridExtra", "MASS", "extraDistr", "dplyr", "ggplot2", "rstan"))
library(gridExtra)
library(MASS)
library(extraDistr)
library(dplyr)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("XXXX")

#-------------------------------------------------------------------------------
# 1. Two compartment first-order absorption PK model (without ke0) source: Dini et al 2018
#-------------------------------------------------------------------------------

two_cpt_oral_firstabs = function(dose,  t,  tLag,  tDose,  CL,   Vc,  Q,  Vp,   ka) {
  
  beta = (((CL/Vc)+(Q/Vc)+(Q/Vp)) - sqrt((((CL/Vc)+(Q/Vc)+(Q/Vp))^2)-4*(CL/Vc)*(Q/Vp)))/2
  alpha = ((Q/Vp) * (CL/Vc)) / beta
  
  A = (ka/Vc) * (Q/Vp - alpha)/((ka - alpha) * (beta - alpha))
  B = (ka/Vc) * (Q/Vp - beta)/((ka - beta) * (alpha - beta))
  Ae = (A) 
  Be = (B) 
  Ce = - (Ae*(ka - alpha) + Be*(ka - beta)) / (ka)
  
  time = t - tLag - tDose
  if (time < 0){
    conc = 0
  } else {
    conc = dose * (Ae * exp(-alpha*time) + Be * exp(-beta*time) + Ce * exp(-time) - (Ae+Be+Ce) * exp(-ka*time))
  }
  
  return(conc)
}



#-------------------------------------------------------------------------------
# 2. Function to simulate a set of PK profiles & prepare/format for stan model
#-------------------------------------------------------------------------------

pk_input_fun <- function(doses, # length (nPt) vector of drug dose
                         tmax, # length (nPt) vector of total simulation length (inc growth)
                         nPt, # how many patients in the dataset
                         t, # a vector of element-wise sum:(nObs*nPt) length, each time of observation, for each pt
                         nObs, # length (nPt) vector of # observations
                         tLag,  # length (nPt) vector of drug delay/lagtime (set to 0 for Cipargamin)
                         tDose,  # length (nPt) vector of dose time (i.e. how long after innoculation)
                         CL, # length (nPt) vector of Clearance rate
                         Vc, # length (nPt) vector of Central compartment volume
                         Q, # length (nPt) vector of inter-compartmental clearance vol
                         Vp,  # length (nPt) vector of peripheral compartment volume
                         ka, # length (nPt) vector of absorption constant
                         seed=seed){  
  
  
  # we set up a vector to index where each participant's observations begin and end
  # (i.e flexible if 'nObs' differs by patient)
  start = c()
  end = c()
  
  # initiate patient 1 starts at 1 and ends at nObs[1]
  start[1] = 1
  end[1] = nObs[1]
  if (nPt >= 2){
    for(i in 2:nPt){
      start[i] = 1 + sum(nObs[1:(i-1)])
      end[i] = sum(nObs[1:(i)])} }
  
  #now using this 'start' and 'end' vector, create a vector of ID's, so each conc value w eproduce has a patient 
  id_vec = c()
  for(i in 1:nPt){id_vec[start[i]:end[i]] = i}
  
  # empty vector to store our output concentrations (note: at the observation times only)
  conc_vec <- c()
  # empty vector with 1 row per pt and 1 col per timepoint (for max possible tmax)
  # this will store the full hour-by-hour concentrations, that are mostly unobserved, for eahc patient
  conc_full <- matrix(0, ncol=max(tmax), nrow=nPt)
  
  # loop over each patient
  for(i in 1:nPt){   
    
    # run PK function for each time from 1:tmax for that patient, this makes a 1xtmax vector of conc
    Conc <- mapply(t=c(1:tmax[i]),function(t){two_cpt_oral_firstabs(dose = doses[i], 
                                                                    t = t, 
                                                                    tLag = tLag[i],
                                                                    tDose = tDose[i], #include growth phase for estimation of Ipl etc
                                                                    CL = CL[i], #clearance rate                  
                                                                    Vc = Vc[i], #central compartment volume       
                                                                    Q = Q[i], #inter-compartmental clearance     
                                                                    Vp = Vp[i], #peripheral compartment volume    
                                                                    ka = ka[i] #absorption duration   
    )})
    
    # save the concentrations from ONLY the observed timepoints (t). 
    # use 'start' and 'end' vectors to identify the correct indexes (e.g. same use as which(id_vec==i))
    conc_vec[start[i]:end[i]] = Conc[t[start[i]:end[i]]]
    # save the full concentration for use in generating the PD profile
    conc_full[i,1:tmax[i]] = Conc}
  
  # set seed as per this simualted dataset
  set.seed(seed)
  
  # create empty vector for log-concentration
  conc_log <- c()
  # create empty vector for concentration * multiplicative error
  conc_vec_er <- c()
  
  #for this simulated dataset, pick sigma^2 value for error dist
  er_lvl <- rlnorm(1,log(0.1),0.01) #inputs are on log scale - want expected value 0.1
  
  #simulate error at this level, multiply by 1000 to go from mg/L to ng/L
  conc_vec_er = conc_vec*1000*exp(rnorm(length(conc_vec),0,er_lvl))
  
  # take log10 of each concentration value, but set <1 to 0 to remove -ve concentrations
  for (j in 1:length(conc_vec)){
    conc_log[j] =  max(0,log10(conc_vec_er[j]))}
  
  # save all the inputs and the outputs in a list to return
  output_pk <- list(nPt = nPt,
                    N = sum(nObs),
                    t = as.array(t),
                    id = as.array(id_vec),
                    nObs = as.array(nObs),
                    tmax = as.array(tmax),
                    conc = conc_vec_er,    
                    log10_conc = conc_log,
                    start = as.array(start),
                    end = as.array(end),
                    tLag = as.array(tLag),
                    tDose = as.array(tDose),
                    concFull = conc_full, #every time, w/ no error - for sims only!
                    dose = as.array(doses),
                    Low = as.array(c(2.75,32.2,6.45,53.5,0.4595)), #limits to use in STAN model
                    Upp = as.array(c(11,128.8,25.8,214,1.838)))
  
  output <- c(output_pk)
  
  return(output)}




#-------------------------------------------------------------------------------
# 3. Function to use seed to create profile, run and save Stan output
#-------------------------------------------------------------------------------

pk_sims_complete <- function(seed){
  
  #First set up distribution for eta, patient variations from transformed pop avg
  omega = diag((c(0.325,0.227,0.1,1,0.1)))
  
  #Use input to set seed
  set.seed(seed)
  
  #create correlation matrix from L:
  L = matrix(c((runif(25,0,1))), nrow=5)
  
  #make L a lower triangular matrix by removing all values above the diagonal
  L[1,2:5] <- L[2,3:5] <- L[3,4:5] <- L[4,5] <- 0
  
  #Eta will be MVN w/ mean 0 and covariance Sigma:
  Sigma = omega%*%((L)%*%t(L))%*%omega
  
  eta = mvrnorm(8,rep(0,5),Sigma)
  
  # Input desired population average values
  pop_avg = c(5.5,64.4,12.9,107,0.919)
  
  a = 2*pop_avg #upper bounds double estimates
  b = 0.5*pop_avg #lower bounds half estimates
  
  #If we want each row to be 1 patient, repeat the vector 8 times
  a_matrix = matrix(rep(a,8), nrow=8, byrow=T)
  b_matrix = matrix(rep(b,8), nrow=8, byrow=T)
  
  # These calculations are essentially element-wise (cant divide vector by vector)
  phi_i <- matrix(rep((log((pop_avg - a)/(b-pop_avg))),8),nrow=8,byrow=T)  + eta
  
  # This is the inverse of the transformation in report section 2.3
  pk_paraMVN_ = (b_matrix*exp(phi_i)+a_matrix)/(exp(phi_i)+1)
  
  # Finally input all these values into the PK function
  data_pkMVN_ <- pk_input_fun(doses = rep(10,8),
                              tmax = rep(288,8),
                              nPt = 8,
                              t = rep(167+c(1,2,3,4,6,8,12,16,24,36,48,72,96,120),8), # here all patients have identical sampling schedules
                              nObs = c(rep(14,8)),
                              tLag = rep(0,8),
                              tDose = rep(168,8),
                              CL = pk_paraMVN_[,1], #clearance rate                 
                              Vc = pk_paraMVN_[,2], #central compartment volume     
                              Q = pk_paraMVN_[,3], #inter-compartmental clearance   
                              Vp = pk_paraMVN_[,4], #peripheral compartment volume    
                              ka = pk_paraMVN_[,5],
                              seed = seed)
  
  #Save simulated profile
  saveRDS(data_pkMVN_,file=paste("simu_MVN_pk_",seed,".rdata", sep=""))
  
  #Create a plot of the 8 profiles:
  pk_plot_data <- data.frame("t" = data_pkMVN_[["t"]], "conc" = data_pkMVN_[["conc"]], "id"=as.factor(data_pkMVN_[["id"]]))
  pk_plot <- ggplot(pk_plot_data, aes(x=t, y=log10(conc), group=id, col=id)) + 
    labs(title=paste("seed =",seed, sep=" ")) +
    #ylim(0,2) + 
    theme_light() + geom_line(cex=0.5) + geom_point(cex=0.5) + xlab("Time (hours)") +
    ylab("(ng/mL)") + scale_y_continuous(breaks=c(0,1,2,3), labels=c("0","10","100","1000")) +
    scale_x_continuous(breaks=seq(168,288,6), labels=seq(0,120,6))
  
  #Save plot
  ggsave(paste("pk_plot",seed,".png",sep=""), pk_plot, width=7.5,height = 3)
  
  #Sample with Stan - first is ones used for paper, second is for toy example
  #PKMVN_ <- stan("pk_hierarchy_cov_5para.stan", data=data_pkMVN_, chains=3, iter=2000, warmup=500,
  #               control = list(adapt_delta=0.99, stepsize_jitter=0.8, max_treedepth=15))
  
  PKMVN_ <- stan("pk_hierarchy.stan", data=data_pkMVN_, chains=3, iter=600, warmup=200,
                 control = list(adapt_delta=0.99, stepsize_jitter=0.8, max_treedepth=15))
  
  #Save Stan output
  save(PKMVN_,file=paste("MVN_pk_",seed,".stan", sep=""))
  
}

set.seed(123)
# HOW many seeds depends how many we want to run in this session, e.g. start with 5
seeds <- c(floor(runif(5,100000,999999)))
sapply(seeds, FUN=pk_sims_complete)
