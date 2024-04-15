library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


library(stringr)
library(gridExtra)
library(ggplot2)
library(patchwork) # for plot inset

# Set wd
setwd("XXX")

# Create list of name files
filenames = list.files(pattern="MVN_pk_*")

# Create a vector of the seeds
seeds <- substring(filenames,8,13)

filenames <- filenames[(is.na(file.info(filenames)$size)==FALSE &
                          file.info(filenames)$size >= 300000) &
                         file.info(filenames)$mtime > "2024-01-01"]

# Create a function to extract key values, can also set size standard here instead. 
read.extract.pk <- function(name){
  load(name)
  stats <- (summary(PKMVN_, pars=c("thetaPop", "thetaInd")))[[1]][,c(1,4,6,8,9,10)]
  #extract mean, 2.5%, med, 97.5%, n_eff, rhat
  
  print(paste("loading simulated model fit",substring(name,8),sep=" "))
  return(stats)
  
}


test <- lapply(filenames, read.extract.pk)

# if this step has already been run, can skip all above and go to the commented readRDS step
saveRDS(test,file="pk_extracted.rdata")

#test <- readRDS("pk_extracted.rdata")

n_datasets <- length(test)

# Create empty vectors to store results
cl_med <- c()
vc_med <- c()
q_med <- c()
vp_med <- c()
ka_med <- c()
cl_2.5 <- c()
vc_2.5 <- c()
q_2.5 <- c()
vp_2.5 <- c()
ka_2.5 <- c()
cl_97.5 <- c()
vc_97.5 <- c()
q_97.5 <- c()
vp_97.5 <- c()
ka_97.5 <- c()
n_eff_min  <- c()
n_eff_mean  <- c()
r_hat_min  <- c()
r_hat_mean  <- c()

# Pull results out of function output 'test' list format
for (i in 1:(n_datasets)){ 
  cl_med[i] <- test[[i]][1,3]
  vc_med[i] <- test[[i]][2,3]
  q_med[i] <- test[[i]][3,3]
  vp_med[i] <- test[[i]][4,3]
  ka_med[i] <- test[[i]][5,3]
  n_eff_min[i] <- min(test[[i]][,2], na.rm=T)     # note this is min within selected parameters only
  n_eff_mean[i] <- mean(test[[i]][,2], na.r=T)
  r_hat_min[i] <- min(test[[i]][,3],na.rm=T)  
  r_hat_mean[i] <- mean(test[[i]][,3],na.rm=T)
  cl_2.5[i] <- test[[i]][1,2]
  vc_2.5[i] <- test[[i]][2,2]
  q_2.5[i] <- test[[i]][3,2]
  vp_2.5[i] <- test[[i]][4,2]
  ka_2.5[i] <- test[[i]][5,2]
  cl_97.5[i] <- test[[i]][1,4]
  vc_97.5[i] <- test[[i]][2,4]
  q_97.5[i] <- test[[i]][3,4]
  vp_97.5[i] <- test[[i]][4,4]
  ka_97.5[i] <- test[[i]][5,4]
  
  
  
  }

# Move these vectors into one combined data frame
pk_extracted_pop <- data.frame(row.names = filenames,
                               Cl = cl_med,
                               Vc = vc_med,
                               Q = q_med,
                               Vp = vp_med,
                               ka = ka_med,
                               n_eff_min = n_eff_min,   # note this is min for selected parameters only
                               n_eff_mean = n_eff_mean,
                               r_hat_min = r_hat_min,
                               r_hat_mean = r_hat_mean,
                               cl_2.5 = cl_2.5,
                               vc_2.5 =  vc_2.5,
                               q_2.5 = q_2.5,
                               vp_2.5 = vp_2.5,
                               ka_2.5 = ka_2.5,
                               cl_97.5 = cl_97.5,
                               vc_97.5 = vc_97.5,
                               q_97.5 =  q_97.5,
                               vp_97.5 = vp_97.5,
                               ka_97.5 = ka_97.5)


#-------------------------------------------------------------------------------
# Results calculations:
#-------------------------------------------------------------------------------

# List actual population mean values used to generate data
pk_pop_mean = c(5.5,64.6,12.9,107,0.919)

# extract means of the posterior medians
(sim_param <- as.numeric(colMeans(pk_extracted_pop)[1:5]))



#calculate bias: #incclude mean as well, 
(abs_bias = sim_param-pk_pop_mean[1:5]) 
(rel_bias = 100*(sim_param/pk_pop_mean[1:5]-1))


q2.5 <- c()
q97.5 <- c()
for (i in 1:5){
  q97.5[i] <- quantile(pk_extracted_pop[,i], probs=c(0.975))
  q2.5[i] <- quantile(pk_extracted_pop[,i], probs=c(0.025))}
q2.5
q97.5


#-------------------------------------------------------------------------------
# Plotting:
#-------------------------------------------------------------------------------
  

library(gridExtra)

# Create a list to store the plots
pk_pop_plots <- list()
pk_pop_big_plots <- list()


# Base plot loads theme & data.frame
base = ggplot(pk_extracted_pop) + theme_minimal()

# Parameter plots for 5 population parameters
for (i in 1:5){
  pk_pop_plots[[i]] <- base +
    geom_histogram(aes_string(x=pk_extracted_pop[,i]),
                   fill="#83a2f7", binwidth=function(x){2 * IQR(x) / (length(x)^(1/1.5))}) +
    geom_vline(xintercept = pk_pop_mean[i], size=0.7, col="darkorange") +
   # geom_vline(xintercept = pk_pop_mean[i], size=0.7) +
    geom_vline(xintercept = mean(pk_extracted_pop[,i]), col = "darkblue", linetype="dashed", size=0.7) +
    # geom_text(x=quantile(pk_extracted_pop[,i],0.8), y=3.3, label="Actual Mean") +
    # geom_text(x=quantile(pk_extracted_pop[,i],0.8), y=4.3, label="Sample Mean", col="darkblue") +
    xlab(paste("Median Pop", colnames(pk_extracted_pop)[i], "Values", sep=" "))
  
  pk_pop_big_plots[[i]] <- base +  geom_histogram(aes_string(x=pk_extracted_pop[,14+i]),
                         fill="lightgreen", binwidth=function(x){2 * IQR(x) / (length(x)^(1/1.5))}) +
    geom_histogram(aes_string(x=pk_extracted_pop[,9+i]),
                   fill="pink", binwidth=function(x){2 * IQR(x) / (length(x)^(1/1.5))}) +
    geom_histogram(aes_string(x=pk_extracted_pop[,i]),
                   fill="lightblue", binwidth=function(x){2 * IQR(x) / (length(x)^(1/1.5))}) +
    geom_vline(xintercept = pk_pop_mean[i], size=0.7) +   
    geom_vline(xintercept = mean(pk_extracted_pop[,9+i]), col = "red", linetype="dashed", size=0.7) +
    geom_vline(xintercept = mean(pk_extracted_pop[,14+i]), col = "darkgreen", linetype="dashed", size=0.7) +
    xlab(paste("2.5%, Median & 97.5% Pop", colnames(pk_extracted_pop)[i], "Values", sep=" ")) +
  geom_vline(xintercept = mean(pk_extracted_pop[,i]), col = "blue", linetype="dashed", size=0.7) 
}




#Plot of histograms for pop values first is just histogram of posterior medians
(medians <- do.call(grid.arrange,c(pk_pop_plots,nrow=1)))
# second is plot of posterior 2.5% and 97.5% quantiles, as well as medians:
(medians2 <- do.call(grid.arrange,c(pk_pop_big_plots,nrow=2)))

ggsave(paste("Median_PK_pop_plot_n",length(filesnames),".png", sep=""), medians)
ggsave(aste("Median_PK_95pctile_pop_plot_n",length(filesnames),".png", sep=""), medians2, width=15, height =9)

#-------------------------------------------------------------------------------

# extra visual checks can be completed for a random subset of 'a', plot full chains
a <- 2

(subset <- floor(runif(a, 0, n_datasets)))

#function to load whole outpu (i.e. all the samples from the chain, no just the median & key percentiles)
read.pk.load <- function(name){
  if(is.na(file.info(name)$size)==FALSE & file.info(name)$size >= 30000000){
    load(name)
    print(paste("sim", name, "did load", sep=" "))
    return(PKMVN_)}
  else{
    print(paste("sim", name, "did not load", sep=" "))
  }
}

simulations <- lapply(filenames[subset],read.pk.load)

#apply this as names to our list of stanfits 'simulations'
names(simulations) <- seeds[subset] #assuming correct order is maintained

#function to save extra diagnostic plots
extra_pk_plots <- function(seed){
  tr <- traceplot(simulations[[seed]], pars=c("thetaPop", "thetaInd","sig"))
  ggsave(paste("traceInd_",seed,".png",sep=""),tr, height=15, width=10)
}

lapply(seeds[subset], extra_pk_plots)

# if we want to select a specific seed/file we can do it like this:
# simulations[["301050"]] <- read.pk.load("MVN_pk_301050.stan")
# simulations[["318807"]] <- read.pk.load("MVN_pk_318807.stan")
# simulations[["919688"]] <- read.pk.load("MVN_pk_919688.stan")
# simulations[["142131"]] <- read.pk.load("MVN_pk_142131.stan")

posterior_plots <- function(seed){
  
  plots_list <- list()
  
  read_data <- readRDS(paste("simu_MVN_pk_", seed,".rdata", sep=""))
  data_conc = data.frame("t" = read_data$t, "conc" = read_data$conc)
  
  pos_data <- as.data.frame(summary(simulations[[seed]], pars="predict_conc")[[1]])
 
  #linear scale plot of data vs posterior median (95% posterior interval)
  for (i in 1:8){
    plots_list[[i]] <- ggplot(pos_data[(168+(288*(i-1))):(288+(288*(i-1))),], aes(x=1:121, y=`50%`)) +
      geom_ribbon(aes(x=1:121, min=`2.5%`,max=`97.5%`),fill="lightblue", lwd=2) + geom_line() +
      theme_light() + geom_point(aes(x=t-168,y=conc), data=data_conc[(13*(i-1)+1):(13*i),]) +
      ylab("ng/mL") + xlab("Time (hours") + ggtitle(label=paste("Patient", i, sep=" "))
    
  }
  
  (combo <- do.call(grid.arrange,c(plots_list,nrow=2)))
  ggsave(paste("posteriors_",seed,".png", sep=""),plot=combo, width=15, height=8)
  
}

# create these posterior plots for our subset
sapply(seeds[subset],posterior_plots)

#-------------------------------------------------------------------------------

# lets plot two PK profiles, one with with the real & one with posterior median of he simulated values
(sim_param <- as.numeric(colMeans(pk_extracted_pop)[1:5]))
#         Cl          Vc           Q          Vp          ka 
#5.4054847  61.8835608  12.3623718 111.4398121   0.9963475

#actual values
pk_pop_mean

# load in SAME pk function from PK sims, if not already in environment
two_cpt_oral_firstabs2 = function(dose,  t,  tLag,  tDose,  CL,   Vc,  Q,  Vp,   ka) {
  
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

# create the profiles, applying across time 1:100 for firstly the 'esimated' or posterior medians, and then the 'real' input values set with pk_pop_mean
profiles <- c(mapply(t=c(1:100),function(t){
  two_cpt_oral_firstabs2(10,t,5,0,sim_param[1],sim_param[2], 
                         sim_param[3], sim_param[4], sim_param[5])}),
  mapply(t=c(1:100),function(t){
    two_cpt_oral_firstabs2(10,t,0,5,pk_pop_mean[1],pk_pop_mean[2], 
                           pk_pop_mean[3], pk_pop_mean[4], pk_pop_mean[5])}))


# load this all into a data frame
compare_data <- data.frame(Parameters=as.factor(rep(c("Estimated", "'Actual'"), each=100)) , conc = 1000*profiles, t=rep(1:100,2)) 

# plot
ggplot(compare_data, aes(x=t, y=conc, col=Parameters, linetype=Parameters)) + theme_light() + 
  #geom_line(alpha=0.2,size=5)  +
  geom_line(size=0.8)  +
  ylab("Conc (ng/mL)") + xlab("Time (hours)") +
  scale_colour_manual(values=c("darkorange","darkblue"))

#-------------------------------------------------------------------------------

# Paper Figure 2 - 1 profile per simulated dataset's posterior medians, compared to 1 profile for 'real' underlying input values


# make all the profiles, across each dataset (i.e. each row of pk_extracted_pop)
Profiles <- c(mapply(i=1:n_datasets, function(i){
  mapply(t=c(1:100),function(t){
  two_cpt_oral_firstabs2(10,t,0,5, #tlag is 5h, this can be nice to plot a straight line, that then goes up, instead of immediate rise
                         pk_extracted_pop[i,1],
                         pk_extracted_pop[i,2], 
                         pk_extracted_pop[i,3], 
                         pk_extracted_pop[i,4], 
                         pk_extracted_pop[i,5])})}))

#combine into 1 big data frame, convert time into DAYS and minus 5 so dose time is t=0 
combo <- data.frame(Conc = Profiles,
                    ID = as.factor(rep(1:n_datasets, each=100)),
                    time = (rep(1:100, n_datasets)-5)/24)

# full profile
(PK_manyprof <- ggplot() + theme_light() +
  geom_line(data=combo, aes(x=time, y=1000*Conc, group=ID), lwd=0.6, 
            col="grey40", alpha=0.03) +
  ylab("Cipargamin Concentration (ng/mL)") + xlab("Time (days)") +

  geom_line(data=compare_data[101:200,], aes(x=(t-5)/24, y=conc), 
            col="black", lwd=1) +  scale_x_continuous(breaks=c(0:5)#, limits=c(0,5)
                                                      ) +
  theme(text=element_text(size=25)))

# absorption phase inset
(PK_manyprof_day1 <- ggplot() + theme_light() +
    geom_line(data=combo, aes(x=(time*24), y=1000*Conc, group=ID), lwd=0.5, 
              col="grey40", alpha=0.1) +
    ylab("Cipargamin Concentration (ng/mL)") + xlab("Time (hours)") +
    
    geom_line(data=compare_data[101:200,], aes(x=t-5, y=conc), 
              col="black", lwd=1) +  scale_x_continuous(breaks=c(0, 2,4,6), 
                                                        limits=c(0,6)
              ) +
    theme(text=element_text(size=20),
          panel.border = element_rect(colour = "black", fill=NA, size=3)))

# combine & save
(comp <- PK_manyprof + inset_element(PK_manyprof_day1, left = 0.3, bottom = 0.3, right = 1, top = 1))
ggsave(paste("PK_fig1_", n_datasets,"_sims.png", sep=""), comp, width = 20, height=10)


#-------------------------------------------------------------------------------
# now repeat with a 95% ribbon for each time point, if we switch to a quantile visualisation
# also option to add at a 25% - 75% (i.e. central 50%) ribbon, commented out

#scale combo time back up to hours
combo$time <- round(combo$time*24)+5
ymin <- c()
ymax <- c()
ymin2 <- c()
ymax2 <- c()
for (i in 1:100){
  ymin[i] <- quantile(combo$Conc[combo$time==i], probs=c(0.025))
  ymax[i] <-quantile(combo$Conc[combo$time==i], probs=c(0.975))
  ymin2[i] <- quantile(combo$Conc[combo$time==i], probs=c(0.25))
  ymax2[i] <-quantile(combo$Conc[combo$time==i], probs=c(0.75))
}

#scale up to match compare data units
ribbon <- data.frame(t = 1:100,
                     min = 1000*ymin,
                     min2 = 1000*ymin2,
                     max2 = 1000*ymax2,
                     max = 1000*ymax)


ggplot() + theme_light() +
  geom_ribbon(data = ribbon, 
              aes(ymin =min, 
                  ymax =max, 
                  x=t),
              fill="lightblue", col="lightblue") +
  #geom_ribbon(data = ribbon,  aes(ymin =min2,  ymax =max2,  x=t), fill="lightgreen", col="lightgreen") +
  geom_line(data=compare_data, aes(x=t, y=conc, col=Parameters, linetype=Parameters), lwd=1.5) +  scale_colour_manual(values=c("black","blue")) +
  ylab("Conc (ng/mL)") + xlab("Time (hours)") +
  theme(text = element_text(size=25))



