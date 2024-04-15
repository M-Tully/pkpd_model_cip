# Need to take MCMC stan output, and extract summary stats
# which can be used to test for model accuracy.

# first import the stanfit files

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library(stringr)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Set wd
setwd("XXX")

## Create list of names of stan files
filenames = list.files(pattern="MVN_pd_LLOQ_*")

# do not load files less than 30 MB, as most are 31MB, (adjust this threshold for n chains & length!)
# and small size is a good sign importing was incomplete or interrupted
# this also clears the rdata files with similar names
filenames <- filenames[(is.na(file.info(filenames)$size)==FALSE &
                          file.info(filenames)$size >= 15000000) &
                         file.info(filenames)$mtime > "2024-01-01"]

n_datasets <- length(filenames)

# Create a vector of the seeds
seeds <- (substring(filenames,13,19))

# Create a function to extract key values
read.extract.pd <- function(name){
    load(name)

    stats <- (summary(PDMVN_, pars=c("thetaPop", "thetaInd")))[[1]][,c(1,4,6,8,9,10)]
    #extract mean, 2.5%, med, 97.5%, n_eff, rhat
    
    print(paste("loading simulated model fit",substring(name,8),sep=" "))
    return(stats[,c(3,5,6)]) # we have currently restricted this to only import posterior median, n_eff and rhat
}

test <- lapply(filenames, read.extract.pd)
saveRDS(test,file=paste("pd_extracted_n", n_datasets,".rdata", sep=""))

# if already saved can skip to this sep:
# test <- readRDS("pd_extracted.rdata")

# Create empty vectors to store results
ipl_med <- c()
iplMu_med <- c()
iplSd_med <- c()
PMF_med <- c()
Emax_med <- c()
Ec50_med <- c()
gamma_med <- c()
n_eff_min  <- c()
n_eff_mean  <- c()
r_hat_min  <- c()
r_hat_mean  <- c()


# Pull results out of function output 'test'
for (i in 1:(length(test))){ #number of sims
  ipl_med[i] <- test[[i]][1,1]
  iplMu_med[i] <- test[[i]][2,1]
  iplSd_med[i] <- test[[i]][3,1]
  PMF_med[i] <- test[[i]][4,1]
  Emax_med[i] <- test[[i]][5,1]
  Ec50_med[i] <- test[[i]][6,1]
  gamma_med[i] <- test[[i]][7,1]
  n_eff_min[i] <- min(test[[i]][,2], na.rm=T)     # note this is min within selected parameters only
  n_eff_mean[i] <- mean(test[[i]][,2], na.r=T)
  r_hat_min[i] <- min(test[[i]][,3],na.rm=T)  
  r_hat_mean[i] <- mean(test[[i]][,3],na.rm=T)}

# Move these vectors into one combined data frame
pd_extracted_pop <- data.frame(#row.names = filenames,
                               Ipl = ipl_med,
                               IplMu = iplMu_med,
                               IplSd= iplSd_med,
                               PMF = PMF_med,
                               Emax = Emax_med,
                               Ec50 = Ec50_med,
                               Gamma = gamma_med,
                               n_eff_min = n_eff_min,   # note this is min for selected parameters only
                               n_eff_mean = n_eff_mean,
                               r_hat_min = r_hat_min,
                               r_hat_mean = r_hat_mean)



# Create a list to store the plots
pd_pop_plots <- list()

# List actual population mean values
pd_pop_mean = c(1800, 2, 3, 13, 0.23, 15.1, 5)

# calculate mean of posterior medians
(sim_param <- as.numeric(colMeans(pd_extracted_pop)[1:7]))

# values for table:
sim_param
100*(1-sim_param/pd_pop_mean)
100*(sim_param/pd_pop_mean)
sim_param-pd_pop_mean

q2.5 <- c()
q97.5 <- c()
for (i in 1:7){
  q97.5[i] <- quantile(pd_extracted_pop[,i], probs=c(0.975))
  q2.5[i] <- quantile(pd_extracted_pop[,i], probs=c(0.025))}
q2.5
q97.5

#----------------------------------------------------------------------------
# Plots of distribution fo posterior medians
#----------------------------------------------------------------------------


# Base plot loads theme & data.frame
base = ggplot(pd_extracted_pop) + theme_minimal()

for (i in 1:7){
  pd_pop_plots[[i]] <- base +
    geom_histogram(aes_string(x=pd_extracted_pop[,i]), alpha=0.85,
                   fill="#83a2f7", binwidth=function(x){10 * IQR(x) / (length(x)^(1/1.5))}) +
    geom_vline(xintercept = pd_pop_mean[i], size=0.8, col="darkorange") +
    geom_vline(xintercept = mean(pd_extracted_pop[,i]), col = "#243c7d", linetype="dashed",size=0.7) +
    # geom_text(x=quantile(pk_extracted_pop[,i],0.8), y=3.3, label="Actual Mean") +
    # geom_text(x=quantile(pk_extracted_pop[,i],0.8), y=4.3, label="Sample Mean", col="darkblue") +
    xlab(paste(colnames(pd_extracted_pop)[i])) + ylab("Density") + 
    scale_y_continuous(breaks=c(5.88*(seq(0,10,2))), 
                       labels=c("0",  "0.02",  "0.04", "0.06", "0.08", "0.1")) +
    theme(text = element_text(size=20))
}


# Combined plot of histograms for pop values
meds <- do.call(grid.arrange,c(pd_pop_plots,nrow=2))

ggsave(paste("Median_PD_pop_plot_n",n_datasets,".png", sep=""), meds)

#----------------------------------------------------------------------------
# Posterior plots for an individual dataset
#----------------------------------------------------------------------------

read.pd.load <- function(name){
  if(is.na(file.info(name)$size)==FALSE & file.info(name)$size >= 15800000){
    load(name)
    print(paste("sim", name, "did load", sep=" "))
    return(PDMVN_)}
  else{
    print(paste("sim", name, "did not load", sep=" "))
  }
}

#pick randomly which n datasets we want to plot (which_files)
n = 2
which_files <- floor(sample(1:n_datasets, n))
which_ids <- seeds[which_files]

plot.pd.pos <- function(id){
  
  # import the full STAN output, which has the distribution of posterior samples
  stan <- read.pd.load(paste("MVN_pd_LLOQ_",id,".stan", sep=""))
  extr <- as.data.frame(extract(stan, pars=c("thetaPop"))[[1]])
  colnames(extr) <- c("Ipl ('000)", "IplMu", "IplSd", "PMF", "Emax", "Ec50","gamma")
  extr[,"Ipl ('000)"] <- extr[,"Ipl ('000)"] /1000
  
  # also import the simulated PD data that was used to fit the data, so we can compare
  rdata <- readRDS(paste("simu_MVN_pd_LLOQ_",id,".rdata", sep=""))
  
  # create nicely formatted ID labels
  ids <- c()
  for(i in 1:8){ids[i] <- paste("Patient", i,sep=" ")}
  
  # combine data for plotting, scaling from count on log scale into a concentration (per mL) on log10 scale
  rplot.data <- data.frame(t = rdata$t, 
                           Cir = log10(exp(rdata$cirCount_measured)/5000),
                           Cir_e = log10(exp(rdata$cirCount_measured_e)/5000),
                           id = as.factor(rep(ids, each = 20)) )
  
  # any values that were <1 became negative so change to 0
  rplot.data$Cir[rplot.data$Cir<0] <- 0
  rplot.data$Cir_e[rplot.data$Cir_e<0] <- 0
  
  # extract quantiles form full STAN posterior samples
  extr_plot <- extract(stan, pars=c("predict_cir"))[[1]]
  stanplot.data <- data.frame(t = rep(1:288, 8),
                              id = as.factor(rep(ids, each=288)),
                              med = NA,
                              p0.025 = NA, 
                              p0.975 = NA)
  
  for (i in 1:(288*8)){
    a <- quantile(extr_plot[,i], probs = c(0.025, 0.5, 0.975))
    # again, same conversion:
    # exp to go from log -> natural scale
    # /5000 to go from total to per mL
    # log10 to match scale in McCarthy
    for(j in 1:3){a[j] <- log10(exp(a[j])/5000)
    if(a[j]<0){a[j] <- 0}}
    stanplot.data[i,"med"] = a[2]
    stanplot.data[i,"p0.025"] = a[1]
    stanplot.data[i,"p0.975"] = a[3]
  }
  
  stanplot.data$Cir <- NA
  stanplot.data$Cir_e <- NA
  for (i in 1:160){
    stanplot.data$Cir[stanplot.data$t==rplot.data$t[i]] <- rplot.data$Cir[i]
    stanplot.data$Cir_e[stanplot.data$t==rplot.data$t[i]] <- rplot.data$Cir_e[i]
  }
  
  #plot without measurement error
  plot_pd <- ggplot(stanplot.data) + 
    geom_ribbon(aes(x=t, min=p0.025, max=p0.975), fill="lightblue") + 
    geom_line(aes(x=t, y=med)) +
    geom_point(aes(x=t, y=Cir)) + 
    geom_segment(aes(x=7*24, xend=7*24, y=0, yend=7.8), alpha=0.8, linetype="dashed",col="darkblue") +
    scale_y_continuous(labels=c("1", "100", "10,000", "100,000"),
                       breaks=seq(0,6,2), limits=c(0,7.8)) +
    ylab("Parasites/mL") +
    geom_hline(yintercept=log10(15), #linetype="dashed",
               alpha=0.5, lwd=0.4) +
    xlab("Days since innoculation (h)") +
    scale_x_continuous(labels=c("0","3", "6", "9", "12"),
                       breaks=seq(0,288,72), limits=c(0,288)) +
    theme_light() + theme(text = element_text(size = 15)) +
    facet_wrap(~ id, ncol=4) +
    annotate("rect", xmin=0, xmax=288, ymin=0, ymax=log10(15), alpha=0.15, fill="black") 
  
  ggsave(paste("pd_example_pos_",id,".png", sep=""), plot_pd,  width=10, height=6)
  
  #plto with measurement error
  plot_pd_e <- ggplot(stanplot.data) + 
    geom_ribbon(aes(x=t, min=p0.025, max=p0.975), fill="lightblue") + 
    geom_line(aes(x=t, y=med)) +
    geom_point(aes(x=t, y=Cir_e)) + 
    geom_segment(aes(x=7*24, xend=7*24, y=0, yend=7.8), alpha=0.8, linetype="dashed",col="darkblue") +
    scale_y_continuous(labels=c("1", "100", "10,000", "100,000"),
                       breaks=seq(0,6,2), limits=c(0,7.8)) +
    ylab("Parasites/mL") +
    geom_hline(yintercept=log10(15), #linetype="dashed",
               alpha=0.5, lwd=0.4) +
    xlab("Days since innoculation (h)") +
    scale_x_continuous(labels=c("0","3", "6", "9", "12"),
                       breaks=seq(0,288,72), limits=c(0,288)) +
    theme_light() + theme(text = element_text(size = 15)) +
    facet_wrap(~ id, ncol=4) +
    annotate("rect", xmin=0, xmax=288, ymin=0, ymax=log10(15), alpha=0.15, fill="black") 
  
  ggsave(paste("pd_example_error_pos_",id,".png", sep=""), plot_pd,  width=10, height=6)
  
  print(paste("Made plots for id:",id))
}

lapply(which_ids, plot.pd.pos)

#----------------------------------------------------------------------------
# Imaginary profile from the central input parameters vs the mean posterior median estimates
#----------------------------------------------------------------------------

# to create sample PD profiles, we need PK data to input for concentration levels
# we will base this off the input 'true' parameters used in the paper
# this required the two_cpt_oral_firstabs2, PD_model and  PKPD_model_initialisation functions found in PK_sims and PD_sims
# also at the bottom of this page, so can run from there

pk_pop_mean = c(5.5,64.6,12.9,107,0.919)
# create a PK profile to input
pk_profile <- c( mapply(t=c(1:100),function(t){
    two_cpt_oral_firstabs(10,t,0,0,pk_pop_mean[1],pk_pop_mean[2], 
                           pk_pop_mean[3], pk_pop_mean[4], pk_pop_mean[5])}))

# now we make a PD profile, firstly based on the input 'true' parameters, 
# and secondly using the mean of the posterior medians (estimated)

profiles <- c(
            #first the profile created by the 'true' input values
             PD_model(268, # Maximum simulation time
                        40, # number of hours of lifecycle of parasite (usually 48)
                        rep(1,40), # killing window of a drug; a vector with values of elements 0: no effect or 1: effect; e.g. effective at all ages: kill_window = rep(1, nStages)
                        sim_param[1], # initial parasite load; e.g. 1e9
                        sim_param[2], # mean of the age distribution of parasites (between 1 and 48)
                        sim_param[3], # spread of the age distribution of parasites (between 1 and 48)
                        sim_param[4], # Parasite multiplication factor
                        sim_param[5], # Maximum killing effect (between 0 and 1)
                        sim_param[6], # Drug concentration at which E = Emax/2
                        sim_param[7], # Sigmoidicity of concentration-effect relationship
                        1000*c(rep(0, 168), pk_profile), # drug concentration at effect site over hourly time points
                        123,
                        24)[[1]], 
              
             # repeat using the estimated parameters
             PD_model(268, # Maximum simulation time
                      40, # number of hours of lifecycle of parasite (usually 48)
                      rep(1,40), # killing window of a drug; a vector with values of elements 0: no effect or 1: effect; e.g. effective at all ages: kill_window = rep(1, nStages)
                      pd_pop_mean[1], # initial parasite load; e.g. 1e9
                      pd_pop_mean[2], # mean of the age distribution of parasites (between 1 and 48)
                      pd_pop_mean[3], # spread of the age distribution of parasites (between 1 and 48)
                      pd_pop_mean[4], # Parasite multiplication factor
                      pd_pop_mean[5], # Maximum killing effect (between 0 and 1)
                      pd_pop_mean[6], # Drug concentration at which E = Emax/2
                      pd_pop_mean[7], # Sigmoidicity of concentration-effect relationship
                      1000*c(rep(0, 168), pk_profile), # drug concentration at effect site over hourly time points
                      123,
                      24)[[1]])

# change to concentration (/5000) and log10 scale
profiles_log10 <- log10(profiles/5000)

#set anything below 0 to 0
profiles_log10[profiles_log10 <= 0] <- 0

# combine this all into 1 data frame
compare_data <- data.frame(Parameters=as.factor(rep(c("fitted", "'real'"), each=268)),
                           conc = profiles_log10, 
                           t=rep(1:268,2)) 

# create plto comparing these 2 profiles
ggplot(compare_data, aes(x=t, y=conc, col=Parameters, linetype=Parameters)) + theme_light() + 
  #geom_line(alpha=0.2,size=5)  +
  geom_line(size=1.3) + scale_color_manual(values=c("darkorange", "darkblue")) +
  ylab("Parasites/mL") + xlab("Time (hours)") +
  scale_x_continuous(breaks=c(24*(0:12))) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c("1", "10", "100", "1000", "10 000", "100 000"))


#----------------------------------------------------------------------------
# Repeat for 2.5% and 97.5% values (note these are the quantiles of the dist of MEDIANS, 
#                                   not the centre of the dist of 2.5% 97.5% posterior quantiles)
#----------------------------------------------------------------------------

profiles_q <- c(PD_model(268, # Maximum simulation time
                        40, # number of hours of lifecycle of parasite (usually 48)
                        rep(1,40), # killing window of a drug; a vector with values of elements 0: no effect or 1: effect; e.g. effective at all ages: kill_window = rep(1, nStages)
                        q97.5[1], # initial parasite load; e.g. 1e9
                        q97.5[2], # mean of the age distribution of parasites (between 1 and 48)
                        q97.5[3], # spread of the age distribution of parasites (between 1 and 48)
                        q97.5[4], # Parasite multiplication factor
                        q97.5[5], # Maximum killing effect (between 0 and 1)
                        q97.5[6], # Drug concentration at which E = Emax/2
                        q97.5[7], # Sigmoidicity of concentration-effect relationship
                        1000*c(rep(0, 168), pk_profile), # drug concentration at effect site over hourly time points
                        123,
                        24)[[1]], 
                PD_model(268, # Maximum simulation time
                         40, # number of hours of lifecycle of parasite (usually 48)
                         rep(1,40), # killing window of a drug; a vector with values of elements 0: no effect or 1: effect; e.g. effective at all ages: kill_window = rep(1, nStages)
                         q2.5[1], # initial parasite load; e.g. 1e9
                         q2.5[2], # mean of the age distribution of parasites (between 1 and 48)
                         q2.5[3], # spread of the age distribution of parasites (between 1 and 48)
                         q2.5[4], # Parasite multiplication factor
                         q2.5[5], # Maximum killing effect (between 0 and 1)
                         q2.5[6], # Drug concentration at which E = Emax/2
                         q2.5[7], # Sigmoidicity of concentration-effect relationship
                         1000*c(rep(0, 168), pk_profile), # drug concentration at effect site over hourly time points
                         123,
                         24)[[1]])

profiles_q_log10 <- log10(profiles_q/5000)
profiles_q_log10 [profiles_q_log10  <= 0] <- 0

compare_data_q <- data.frame(Parameters=as.factor(rep(c("97.5", "2.5"), each=268)), 
                            conc = profiles_q_log10, 
                            t=rep(1:268,2)) 

ggplot(compare_data_q, aes(x=t, y=conc, col=Parameters, linetype=Parameters)) + theme_light() + 
  geom_line(size=1.3) + 
  scale_color_manual(values=c("pink", "darkred")) +
  ylab("Parasites/mL") + xlab("Time (hours)") +
  scale_x_continuous(breaks=c(24*(0:12))) +
  geom_line(data=compare_data[compare_data$Parameters=="'actual'",], aes(x=t, y=conc), size=0.3, col="black") +
  scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c("1", "10", "100", "1000", "10 000", "100 000"))


#----------------------------------------------------------------------------
#create a profile for each set of parameters (fig 3 in the paper)
#----------------------------------------------------------------------------


Profiles_PD <- c(mapply(i=1:n_datasets, function(i){
  PD_model(268, # Maximum simulation time
           40, # number of hours of lifecycle of parasite (usually 48)
           rep(1,40), # killing window of a drug; a vector with values of elements 0: no effect or 1: effect; e.g. effective at all ages: kill_window = rep(1, nStages)
           pd_extracted_pop[i,1], # initial parasite load; e.g. 1e9
           pd_extracted_pop[i,2], # mean of the age distribution of parasites (between 1 and 48)
           pd_extracted_pop[i,3], # spread of the age distribution of parasites (between 1 and 48)
           pd_extracted_pop[i,4], # Parasite multiplication factor
           pd_extracted_pop[i,5], # Maximum killing effect (between 0 and 1)
           pd_extracted_pop[i,6], # Drug concentration at which E = Emax/2
           pd_extracted_pop[i,7], # Sigmoidicity of concentration-effect relationship
           1000*c(rep(0, 168), pk_profile), # drug concentration at effect site over hourly time points
           123,
           24 #what age sequestered (26)
  )[[1]]
}))

combo_PD <- data.frame(Conc = Profiles_PD,
                       ID = as.factor(rep(1:n_datasets, each=268)),
                       time = rep(1:268, n_datasets))
combo_PD$Conc = log10(combo_PD$Conc/5000)
combo_PD$Conc[combo_PD$Conc<0] <- 0


(fig <- ggplot() + theme_light() +
  geom_line(data=combo_PD, 
            aes(x=time, y=Conc, group=ID), lwd=0.5, col="grey40", alpha=0.07) +
  ylab("Parasites/mL") + xlab("Time (days)") +
  geom_line(data=compare_data[compare_data$Parameters=="'real'",], aes(x=t, y=conc), col="black", lwd=1.2) + 
  scale_x_continuous(breaks=48*c(0:6), labels=c("0", "2", "4", "6", "8", "10", "12")) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c("1", "10", "100", "1000", "10 000", "100 000")) +
  theme(text=element_text(size=15)) + 
  geom_vline(xintercept = 7*24, linetype="dashed"))
ggsave("00_pd_fig.png", fig, width=8, height=5)


#create the ribbon version of the above plot
ymin <- c()
ymax <- c()
ymid <- c()
for (i in 1:268){
  ymin[i] <- quantile(combo_PD$Conc[combo_PD$time==i], probs=c(0.025))
  ymax[i] <-quantile(combo_PD$Conc[combo_PD$time==i], probs=c(0.975))
  ymid[i] <-quantile(combo_PD$Conc[combo_PD$time==i], probs=c(0.5))
  }
ribbon <- data.frame(t = 1:268,
                     min = ymin,
                     max = ymax,
                     mid = ymid)


ggplot() + theme_light() +
  geom_ribbon(data = ribbon, aes(ymin =min, ymax =max, x=t), fill="chocolate1", col="chocolate1") +
   geom_line(data = ribbon, aes(y =mid, x=t), fill="chocolate1", col="chocolate3") +
  geom_line(data=compare_data[compare_data$Parameters=="'real'",], aes(x=t, y=conc), col="black", lwd=1.2) + 
  scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c("1", "10", "100", "1000", "10 000", "100 000")) +
  ylab("Parasites/mL") + xlab("Time (hours)") +
  scale_x_continuous(breaks=c(24*(1:12)), limits=c(12, 266)) +
  theme(text = element_text(size=25))





# extra functions required to make profile plots

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


