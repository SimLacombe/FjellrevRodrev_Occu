
sapply(c('stringr', 'foreach',
         'data.table', 'dplyr',
         'tidyr', 'LaplacesDemon',
         'runjags', 'rjags',
         'coda', 'rlist', 
         'mcmcplots', 'bayesplots',
         'ggplot2'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

rm(list=ls())

years_1 <- as.character(2006:2016)
years_2 <- as.character(2017:2021)

###LOAD MODELS ------
model.path <- "outputs/"
model.filenames <- list.files(path = model.path, pattern = ".rds")

M_1 <- readRDS(paste0(model.path, "/", model.filenames[1]))
M_2 <- readRDS(paste0(model.path, "/", model.filenames[2]))

M.mat_1 <- as.matrix(as.mcmc.list(M_1), chains = T)
M.mat_1 <- M.mat_1[, -grep("z", colnames(M.mat_1))]

M.mat_2 <- as.matrix(as.mcmc.list(M_2), chains = T)
M.mat_2 <- M.mat_2[, -grep("z", colnames(M.mat_2))]

###CHECK CONVERGENCE ------
# M.theta_1 <- as.mcmc.list(M_1)[,c(1:54)]
# gelman.diag(M.theta_1)
# denplot(M.theta_1, parms= c("h[1,1]", "h[2,1]", "rho_bait[1]", "rho_bait[2]"))
# 
# M.theta_2 <- as.mcmc.list(M_2)[,c(1:32)]
# gelman.diag(M.theta_2)
# denplot(M.theta_2, parms= c("h[1,1]", "h[2,1]", "rho_bait[1]", "rho_bait[2]"))

###GET CI FOR ALL PARAMS -----
M.sum_1 <- M.mat_1%>%
  as.data.frame()%>%
  pivot_longer(col = !CHAIN)%>%
  group_by(name)%>%
  summarise(mean = quantile(value,0.5),
            lower = quantile(value,0.025),
            upper= quantile(value,0.975))%>%
  mutate(model = "2006 - 2016")

M.sum_2 <- M.mat_2%>%
  as.data.frame()%>%
  pivot_longer(col = !CHAIN)%>%
  group_by(name)%>%
  summarise(mean = quantile(value,0.5),
            lower = quantile(value,0.025),
            upper= quantile(value,0.975))%>%
  mutate(model = "2017 - 2021")

M.sum <- rbind(M.sum_1, M.sum_2)

###GET LABELS FOR PLOTS -----
param <- c("psi", "gamma", "epsilon", "pi", "tau", "rho", "yr_rho", "rho_bait")
param_ch <- c("a\\[","b\\[","d\\[","g\\[","h\\[", "f\\[", "yr_rho", "rho_bait")

M.sum$param <- ""
for(i in 1:length(param)){
  M.sum[grep(param_ch[i], M.sum$name),]$param <- param[i]
}

M.sum$species <- "ArcticFox"
M.sum[grep("\\[1", M.sum$name),]$species <- "RedFox"

M.sum$cov <- ""
for(i in unique(M.sum$param)){
  if (i %in% c("psi", "gamma","epsilon", "pi", "tau")){
    cov <- c("int", "CLG", "FTG", "rodents")
    cov_ch <- c(",1\\]",",2\\]",",3\\]",",4\\]")
  }
  if (i %in% c("rho")){
    cov <- c("rodents")
    cov_ch <- c(",1\\]")
  }
  if (i %in% c("rho_bait")){
    cov <- c("bait")
    cov_ch <- c("")
  }
  if (i %in% c("yr_rho")){
    year_ch <- paste0(as.character(1:max(length(years_1),length(years_2))),"\\]")
    for (j in 1:length(years_1)){
      M.sum[rownames(M.sum)%in%as.character(grep(year_ch[j], M.sum$name))&
              M.sum$param == i&
              M.sum$model ==  "2006 - 2016",]$cov <- years_1[j]
    }
    for (j in 1:length(years_2)){
      M.sum[rownames(M.sum)%in%as.character(grep(year_ch[j], M.sum$name))&
              M.sum$param == i&
              M.sum$model == "2017 - 2021",]$cov <- years_2[j]
    }
  }
  if (!i %in% c("yr_rho"))
    {
    for(j in 1:length(cov)){
      M.sum[rownames(M.sum)%in%as.character(grep(cov_ch[j], M.sum$name))&M.sum$param == i,]$cov <- cov[j]
    }
  }
}

M.sum[M.sum$param %in% c("rho_bait", "yr_rho"),]$param <- "rho"

#########PLOT############

ggplot(M.sum[M.sum$species == "ArcticFox",])+
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_pointrange(aes(x = cov, y = mean, ymin=lower, ymax=upper, color = model), 
                  position = position_jitterdodge(), show.legend = F, size = .8)+
  geom_line(aes(x=0,y=0,color=model))+
  geom_point(aes(x=0,y=0,color=model))+
  scale_color_discrete(direction=-1)+
  facet_wrap(~param, scale = "free")+
  coord_flip()+
  xlab("")+
  ylab("")+
  theme(panel.background = element_rect(color="darkgrey",fill = "white"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color="darkgrey" ),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12, hjust = 0),
        strip.placement = "outside",
        axis.text.y = element_text(face = "bold"),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(face = "bold", size = 11),
        legend.title=element_text(face = "bold", size = 11))+
  guides(colour = guide_legend(title.position="top", title.hjust = 0))
