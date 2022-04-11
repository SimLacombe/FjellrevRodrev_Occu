##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##      Arctic Foxes in Varanger
##               --
##            Fit DCOMs
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##***************************
## INSERT DESCRIPTION HERE
##***************************
##
# install.packages("foreach")
# install.packages("stringr")
# install.packages("data.table")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("LaplacesDemon")
# install.packages("runjags")
# install.packages("rjags")
# install.packages("coda")
# install.packages("doParallel")
rm(list=ls())
##PARAMETERS TO SET MANUALY----
YEARMIN = 2006
YEARMAX = 2016
Nobs_min = 3 #min number of days sampled to consider week valid
Nweek_min = 4 #min number of weeks sampled a year to consider the year valid
CUTOFF = 35 #min number of photo a day to consider the observation valid

##GET LIBRARIES, PATHS AND FILENAMES-----
sapply(packages <- c('stringr', 'foreach', 'data.table', 'dplyr','tidyr', 'LaplacesDemon', 'runjags', 'rjags', 'coda', 'doParallel'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))
## GET DATA AND UTILITY FUNCTIONS ----
source("src/format_data_keepdates.R")
source("src/fit_models_utility_functions.R")

## MODELS ----
cat("##START PARALLEL RUNS ----------------\n")
###GLOBAL MODEL PARAMETERS

n.chains <- 4
adapt <- 1000
burnin <- 5000
sample <- ceiling(300)
thin <- 5
  
cat("###M----------------\n")
######################################################
#                   MODEL M
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
psi_covs <- c("int", "CLG", "FTG", "rodents_fall")
gam_covs <- c("int", "CLG", "FTG", "rodents_fall")
eps_covs <- c("int", "CLG", "FTG", "rodents_fall")
pi_covs  <- c("int")
tau_covs <- c("int")
rho_covs <- c("rodents_fall")
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
data_list <- list(psi_cov = covs[,psi_covs]%>%as.matrix(),
                  gam_cov = covs[,gam_covs]%>%as.matrix(),
                  eps_cov = covs[,eps_covs]%>%as.matrix(),
                  pi_cov = covs[,pi_covs]%>%as.matrix(),
                  tau_cov = covs[,tau_covs]%>%as.matrix(),
                  rho_cov = covs[,rho_covs]%>%as.matrix(),
                  year_cov = covs[,"year"]%>%as.matrix(),
                  ncov_psi = length(psi_covs), ncov_gam = length(gam_covs),
                  ncov_eps = length(eps_covs), ncov_pi = length(pi_covs),
                  ncov_tau = length(tau_covs), ncov_rho = length(pi_covs),
                  nspec = 2, nseason = T, nsite = M,
                  nsurvey = K, nout = 4, nyear = length(allyears),
                  y = ob_state,
                  bait=bait)


M <- run.jags(model = "src/dcom.R",
               monitor = c("a", "b", "d","f","g","h", "rho_bait", "yr_rho", "z"),
               data = data_list,
               n.chains = n.chains,
               inits = inits,
               adapt = adapt,
               burnin = burnin,
               sample = sample,
               thin = thin,
               summarise = TRUE,
               plots = FALSE,
               method = "parallel")
    
    
M_matrix <- as.matrix(as.mcmc.list(M), chains = TRUE)
    
saveRDS(M, "outputs/M.rds")
 