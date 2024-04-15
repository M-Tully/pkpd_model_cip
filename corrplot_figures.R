#import fits and do diagnostic plots on them
library(rstan)
library(rstan)
#library(bayesplot)
library(GGally)
library(ggpubr)
library(deeptime)

library(extrafont)
font_import()
loadfonts(device = "win")

#-----------PK ---------------------------

setwd("XXX")
stanfiles_PK <- list.files(pattern=".stan")
rfiles_PK <- list.files(pattern=".rdata")

id <- substring(stanfiles_PK[1], 8, 13)
load(stanfiles_PK[1])
extr_PK <- as.data.frame(extract(PKMVN_, pars=c("thetaPop"))[[1]])
colnames(extr_PK) <- c("Cl", "Vc", "Q", "Vp", "ka")
A = 1.3

real_vals <- data.frame(Cl = 5.5, Vc = 64.4, Q = 12.9, Vp = 107, ka = 0.919)

extr_PK <- rbind(extr_PK, real_vals)
extr_PK$type <- c(rep("Sim", 4500), "Real")

#pair_pk <- ggpairs(extr_PK, upper=list(continuous=wrap(ggally_cor, stars=F))) + 
#  theme_light() 

plots <- list()

for(i in  c("Cl", "Vc", "Q", "Vp", "ka")){
  for(j in  c("Cl", "Vc", "Q", "Vp", "ka")){
    (plots[[paste(i,j,sep="_x_")]] <-  ggplot()  + 
      geom_point(data=extr_PK[extr_PK$type=="Sim",], aes_string(x=i, y=j), alpha=0.25, pch=1, size=2.5) +
      geom_point(data=extr_PK[extr_PK$type=="Real",], aes_string(x=i, y=j), fill="red3", col="white", size=6, pch=21)+ theme_classic() +
       theme(text=element_text(size=22),
             axis.title.y=element_text(face="italic", size=26),
             axis.title.x=element_text(face="italic", size=26),
             plot.margin = unit(c(0,0,0,0), 'lines')) +
       labs(x="", y=""))
  }
}

(PK_pairs <- ggarrange2(plots[["Cl_x_Vc"]] + labs(y = bquote(italic(V[c]))),
                        plots[["Cl_x_Q"]] + labs(y = bquote(italic(Q))),
                        plots[["Cl_x_Vp"]] + labs(y = bquote(italic(V[p]))),
                        plots[["Cl_x_ka"]] + labs(y = bquote(italic(k[a])),
                                                  x = bquote(italic(Cl))),
           plots[["Vc_x_Q"]],plots[["Vc_x_Vp"]],plots[["Vc_x_ka"]]+ labs(x = bquote(italic(V[c]))),
           plots[["Q_x_Vp"]],plots[["Q_x_ka"]]+ labs(x = bquote(italic(Q))),
           plots[["Vp_x_ka"]]+ labs(x = bquote(italic(V[p]))),
           layout = matrix(c(1:4, 0, 5:7, 0, 0, 8, 9, 0, 0, 0, 10), nrow = 4, byrow = F),
           widths = c(1.015, 1,1,1),
           heights = c(1,1,1,1.015)))


ggsave("compile_pk.png", PK_pairs, width=45, height=45, units="cm")


#-----------PD ---------------------------

setwd("C:/Users/tully/OneDrive - The University of Melbourne/Documents/MBiostatPaper/pd/")
stanfiles_PD <- list.files(pattern=".stan")
rfiles_PD <- list.files(pattern=".rdata")

id <- substring(stanfiles_PD[1], 8, 13)
load(stanfiles_PD[1])
stan806034 <- PDMVN_
extr_PD <- as.data.frame(extract(stan806034, pars=c("thetaPop"))[[1]])
colnames(extr_PD) <- c("ipl", "iplMu", "iplSig", "PMF", "Emax", "EC_50", "gamma")

real_vals <- data.frame(ipl = 1800, iplMu = 2, iplSig = 3, PMF = 13,
                        Emax = 0.23, EC_50 = 15.1, gamma = 5)

extr_PD <- rbind(extr_PD, real_vals)
extr_PD$type <- c(rep("Sim", 1800), "Real")

#pair_pd <- ggpairs(extr_PD, upper=list(continuous=wrap(ggally_cor, stars=F))) + 
#  theme_light() 
extr_PD$ipl <- extr_PD$ipl/1000



plots <- list()

for(i in  c("ipl", "iplMu", "iplSig", "PMF", "Emax", "EC_50", "gamma")){
  for(j in  c("ipl", "iplMu", "iplSig", "PMF", "Emax", "EC_50", "gamma")){
    (plots[[paste(i,j,sep="_x_")]] <-  ggplot()  + 
       geom_point(data=extr_PD[extr_PD$type=="Sim",], aes_string(x=i, y=j), alpha=0.4, pch=1, size=2.5) +
       geom_point(data=extr_PD[extr_PD$type=="Real",], aes_string(x=i, y=j), fill="red3", col="white", size=6, pch=21)+ theme_classic() +
       theme(text=element_text(size=18),
             axis.title.y=element_text(face="italic", size=26),
             axis.title.x=element_text(face="italic", size=26),
             plot.margin = unit(c(0.1,0.1,0.1,0.1), 'lines'))+
     labs(x="", y=""))
  }
}

PD_pairs <- ggarrange2(plots[["ipl_x_iplMu"]] + labs(y=expression(mu["ipl"])),
                       plots[["ipl_x_iplSig"]] + labs(y = expression(sigma["ipl"])),
                       plots[["ipl_x_PMF"]]+ labs(y = "PMF"),
                       plots[["ipl_x_Emax"]]+ labs(y =  expression("E"[max])),
                       plots[["ipl_x_EC_50"]]+ labs(y = expression("EC"[50])),
                       plots[["ipl_x_gamma"]]+ labs(x = "ipl",
                                                    y = expression(gamma)),
                       plots[["iplMu_x_iplSig"]],
                       plots[["iplMu_x_PMF"]], 
                       plots[["iplMu_x_Emax"]],
                       plots[["iplMu_x_EC_50"]], 
                       plots[["iplMu_x_gamma"]]+ labs(x=expression(mu["ipl"])),
                       plots[["iplSig_x_PMF"]], 
                       plots[["iplSig_x_Emax"]], 
                       plots[["iplSig_x_EC_50"]], 
                       plots[["iplSig_x_gamma"]] + labs(x = expression(sigma["ipl"])),
                       plots[["PMF_x_Emax"]],
                       plots[["PMF_x_EC_50"]], 
                       plots[["PMF_x_gamma"]] + labs(x = "PMF"),
                       plots[["Emax_x_EC_50"]], 
                       plots[["Emax_x_gamma"]]+ labs(x =  expression("E"[max])),
                       plots[["EC_50_x_gamma"]] + labs(x =  expression("EC"[50])),
                       
                       
                       
                       layout = matrix(c(1:6, 0, 
                                         7:11, 0, 0, 
                                         12:15, 0, 0, 0, 
                                         16:18, 0, 0, 0, 0, 
                                         19, 20, 0, 0, 0, 0, 0 ,
                                         21), nrow = 6, byrow = F),
                       widths = c(1.005, 1,1,1,1,1),
                       heights = c(1, 1,1,1,1,1.005))


ggsave("compile_pd.png", PD_pairs, width=50, height=50, units="cm")
