These files can be used to simulate and fit a Bayesian heirarchical Pharmacokinetic-Pharmacodynamic (PK-PD) model

The framework has been customised to Cipargamin as per the Phase 2 clinical trial McCarthy et al (2021 10.1128/AAC.01423-20 PMCID: PMC7849011)
However it could be adjusted to different drugs or sampling schemes. Many of the functions are copies or adaptations of those presented in Dini et al (2018 10.1128/AAC.01068-18 PMCID: PMC6201091)

The file "PK_sims.R" simulates PK (drug-concentration) data sets and fits them using the Stan model "pk_hierarchy.stan"
The file "PK_sims_extract.R" compiles some basic diagnostic summary measures of the model fits and creates some visualisations

The file "PD_sims.R" simulates PD (parasite-concentration) data sets and fits them using the Stan model "pd_hierarchy.stan"
The file "PD_sims_extract.R" compiles some basic diagnostic summary measures of the model fits and creates some visualisations

Lastly "corrplot_figures.R" creates figures 3 & 5 in the manuscript, looking at a single dataset's posterior samples
