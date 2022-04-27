##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##      Arctic Foxes in Varanger
##               --
##      Calculate CPO scores
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##***************************
##***************************
##

get_gof_opened <- function(mm = NULL, data_list = NULL, e = .0001, Nsim = 1000){
  
  #locate each element in data_list to a single variable
  for(i in 1:length(data_list)){
    assign(names(data_list)[i], data_list[[i]])
  }
  
  chi2_l <- chi2.sim_l <- rep(0,Nsim)
   
  cat(paste0("-------------------------------------------------| ", Nsim, "\n"))
  O <- sample(1:nrow(mm),Nsim)
  for(i in 1:Nsim){
    if (!i%%(Nsim/50)) cat("+")
    o <- O[i]
    a <- matrix(mm[o,grep("a\\[", colnames(mm))], ncol = ncov_psi, nrow = nspec)
    b <- matrix(mm[o,grep("b\\[", colnames(mm))], ncol = ncov_gam, nrow = nspec)
    d <- matrix(mm[o,grep("d\\[", colnames(mm))], ncol = ncov_eps, nrow = nspec)
    g <- matrix(mm[o,grep("g\\[", colnames(mm))], ncol = ncov_pi, nrow = nspec)
    h <- matrix(mm[o,grep("h\\[", colnames(mm))], ncol = ncov_tau, nrow = nspec)
    
    z <- matrix(mm[o,grep("z", colnames(mm))], ncol = nseason, nrow = nsite)
      
    psinit  <- a %*% t(psi_cov)
    gam <-     b %*% t(gam_cov)
    eps <-     d %*% t(eps_cov)

    pi <- tau <-  gam_one <- eps_one <- matrix(0, nrow=nspec, ncol=nsite)
    if("pi_cov" %in% names(data_list)){
      pi <- g %*%t(pi_cov)
      tau <- h %*%t(tau_cov)
    }
    gam_one <- b %*% t(gam_cov) + pi
    eps_one <- d %*% t(eps_cov) + tau
    
    # This is the numerator of the softmax function for each of the 4 community
    # states a site can be in.
    fsm <- matrix(0, ncol = 4, nrow = nsite)
    fsm[, 1] <- 1 #---------------------------------------------------------|U
    fsm[, 2] <- exp( psinit[1, ] ) #----------------------------------------|A
    fsm[, 3] <- exp( psinit[2, ] ) #----------------------------------------|B
    fsm[, 4] <- exp( psinit[1, ] + psinit[2, ] ) #--------------------------|AB

    tpm <- array(0, dim = c(nsite, 4, 4))
    # U to ...
    tpm[ , 1, 1] <- 1 #-----------------------------------------------------|U
    tpm[ , 2, 1] <- exp( gam[1, ] ) #---------------------------------------|A
    tpm[ , 3, 1] <- exp(            gam[2, ] ) #----------------------------|B
    tpm[ , 4, 1] <- exp( gam[1, ] + gam[2, ] ) #----------------------------|AB
    # A to ...
    tpm[, 1, 2] <- exp( eps[1, ] ) #----------------------------------------|U
    tpm[, 2, 2] <- 1 #------------------------------------------------------|A
    tpm[, 3, 2] <- exp( eps[1, ] + gam_one[2, ] ) #----------------------|B
    tpm[, 4, 2] <- exp(            gam_one[2, ] ) #----------------------|AB
    # B to ...
    tpm[, 1, 3] <- exp(                   eps[2, ] ) #----------------------|U
    tpm[, 2, 3] <- exp( gam_one[1, ] + eps[2, ] ) #----------------------|A
    tpm[, 3, 3] <- 1 #------------------------------------------------------|B
    tpm[, 4, 3] <- exp( gam_one[1, ] ) #---------------------------------|A
    # AB to ..
    tpm[, 1, 4] <- exp( eps_one[1, ] + eps_one[2, ] ) #---------------|U
    tpm[, 2, 4] <- exp(                   eps_one[2, ] ) #---------------|A
    tpm[, 3, 4] <- exp( eps_one[1, ] ) #---------------------------------|B
    tpm[, 4, 4] <- 1 #------------------------------------------------------|AB
    
    #SIMULATED NUMBER OF DETECTIONS : y.sim
    #EXPECTED NUMBER OF DETECTIONS : Y.exp

    z.sim <- array(NA, dim = c(nsite,nseason))
    
    psi <- array(0,dim = c(nsite,nseason,4))
    
    for(m in 1:nsite){
      psi[m,1,] <- fsm[m,]/sum(fsm[m,])
      z.sim[m,1] <- rcat(1,fsm[m,])
      for(t in 2:nseason ){
        psi[m,t,] <- tpm[m,,]%*%psi[m,t-1,]
        psi[m,t,] <- psi[m,t,]/sum(psi[m,t,])
        z.sim[m,t] <- rcat(1,tpm[m, , z.sim[m,t-1]])
      }
    }
    z_freq <- apply(z, 1, function(x) {table(factor(x, levels = 1:4))/length(which(!is.na(x)))})
    z.sim_freq <- apply(z.sim, 1, function(x) {table(factor(x, levels = 1:4))/length(which(!is.na(x)))})
    z.exp_freq <- apply(psi, 1, function(x) {apply(x,2,function(x) sum(x, na.rm=T))/sum(x, na.rm=T)})
    
  #COMPUTE CHI²
    
    chi2 <- sum((z_freq - z.exp_freq)^2/(z.exp_freq+e))
    chi2.sim <- sum((z.sim_freq - z.exp_freq)^2/(z.exp_freq+e))
    chi2_l[i] <- chi2
    chi2.sim_l[i] <- chi2.sim

  }
  return(data.frame(chi2.obs = chi2_l, chi2.sim=chi2.sim_l))
}

################################################################################

get_gof_closed <- function(mm = NULL, data_list = NULL, e = .0001, Nsim = 1000){
  
  #locate each element in data_list to a single variable
  for(i in 1:length(data_list)){
    assign(names(data_list)[i], data_list[[i]])
  }
  
  detfreq <- apply(y, c(1,2), function(x) {table(factor(x, levels = 1:4))})
  
  chi2_l <- chi2.sim_l <- rep(0,Nsim)
  
  cat(paste0("-------------------------------------------------| ", Nsim, "\n"))
  O <- sample(1:nrow(mm),Nsim)
  for(i in 1:Nsim){
    if (!i%%(Nsim/50)) cat("+")
    o <- O[i]
    f <- matrix(mm[o,grep("f\\[", colnames(mm))], ncol = ncov_rho, nrow = nspec)
    yr_rho <- matrix(mm[o, grep("yr_rho", colnames(mm))], ncol = nyear, nrow = nspec)
    rho_bait <- matrix(mm[o, grep("rho_bait", colnames(mm))], ncol = 1, nrow = nspec)
    # estimated community state at time t and site j
    z <- matrix(mm[o,grep("z", colnames(mm))], ncol = nseason, nrow = nsite)
    
    rho <-     yr_rho[,year_cov] + f %*% t(rho_cov)
    
    rdm <- array(0, dim = c(nsite, 4, 4,2))
    ## Bait absent
    # TS = U
    rdm[, 1, 1, 1] <- 1 #----------------------------------------------------|OS = U
    rdm[, 2, 1, 1] <- 0 #----------------------------------------------------|OS = A
    rdm[, 3, 1, 1] <- 0 #----------------------------------------------------|OS = B
    rdm[, 4, 1, 1] <- 0 #----------------------------------------------------|OS = AB
    # TS = A
    rdm[, 1, 2, 1] <- 1 #----------------------------------------------------|OS = U
    rdm[, 2, 2, 1] <- exp( rho[1, ] ) #-------------------------------------|OS = A
    rdm[, 3, 2, 1] <- 0 #----------------------------------------------------|OS = B
    rdm[, 4, 2, 1] <- 0 #----------------------------------------------------|OS = AB
    # TS = B
    rdm[, 1, 3, 1] <- 1 #----------------------------------------------------|OS = U
    rdm[, 2, 3, 1] <- 0 #----------------------------------------------------|OS = A
    rdm[, 3, 3, 1] <- exp( rho[2, ] ) #-------------------------------------|OS = B
    rdm[, 4, 3, 1] <- 0 #----------------------------------------------------|OS = AB
    # TS = AB
    rdm[, 1, 4, 1] <- 1 #----------------------------------------------------|OS = U
    rdm[, 2, 4, 1] <- exp( rho[1, ] ) #-------------------------------------|OS = A
    rdm[, 3, 4, 1] <- exp(             rho[2, ] ) #-------------------------|OS = B
    rdm[, 4, 4, 1] <- exp( rho[1, ] + rho[2, ] ) #-------------------------|OS = AB
    
    ## Bait present
    # TS = U
    rdm[, 1, 1, 2] <- 1 #------------------------------------------------------|OS = U
    rdm[, 2, 1, 2] <- 0 #------------------------------------------------------|OS = A
    rdm[, 3, 1, 2] <- 0 #------------------------------------------------------|OS = B
    rdm[, 4, 1, 2] <- 0 #------------------------------------------------------|OS = AB
    # TS = A
    rdm[, 1, 2, 2] <- 1 #------------------------------------------------------|OS = U
    rdm[, 2, 2, 2] <- exp( rho[1, ] + rho_bait[1, 1]) #------------------------|OS = A
    rdm[, 3, 2, 2] <- 0 #------------------------------------------------------|OS = B
    rdm[, 4, 2, 2] <- 0 #------------------------------------------------------|OS = AB
    # TS = B
    rdm[, 1, 3, 2] <- 1 #------------------------------------------------------|OS = U
    rdm[, 2, 3, 2] <- 0 #------------------------------------------------------|OS = A
    rdm[, 3, 3, 2] <- exp( rho[2, ] + rho_bait[2, 1]) #------------------------|OS = B
    rdm[, 4, 3, 2] <- 0 #------------------------------------------------------|OS = AB
    # TS = AB
    rdm[, 1, 4, 2] <- 1 #---------------------------------------------------------|OS = U
    rdm[, 2, 4, 2] <- exp( rho[1, ] + rho_bait[1,1]) #----------------------------|OS = A
    rdm[, 3, 4, 2] <- exp(                             rho[2, ] + rho_bait[2,1]) #|OS = B
    rdm[, 4, 4, 2] <- exp( rho[1, ] + rho_bait[1,1] + rho[2, ] + rho_bait[2,1]) #-|OS = AB
    
    #SIMULATED NUMBER OF DETECTIONS : y.sim
    
    y.sim <- array(NA, dim = c(nsite,nseason, nsurvey))
    y.exp <- array(NA, dim = c(nsite,nseason, nsurvey, 4))
    for(m in 1:nsite){
      for(t in 1:nseason ){
        for(k in 1:nsurvey){
          if(!is.na(y[m,t,k])){
            bait_ <- bait[m,t,k]
            y.sim[m,t,k] <- rcat(1, rdm[m, ,z[m,t], bait_])
            y.exp[m,t,k,] <- rdm[m, ,z[m,t], bait_]/sum(rdm[m, ,z[m,t], bait_])
          }
        }
      }
    }

    detfreq_sim <- apply(y.sim, c(1,2), function(x) {table(factor(x, levels = 1:4))})
    detfreq_exp <- apply(y.exp, c(1,2), function(x) {apply(x,2,function(x) sum(x, na.rm=T))})
    
    #COMPUTE CHI²
    
    chi2 <- sum((detfreq - detfreq_exp)^2/(detfreq_exp+e))
    chi2.sim <- sum((detfreq_sim - detfreq_exp)^2/(detfreq_exp+e))
    chi2_l[i] <- chi2
    chi2.sim_l[i] <- chi2.sim
    
  }
  return(data.frame(chi2.obs = chi2_l, chi2.sim=chi2.sim_l))
}

