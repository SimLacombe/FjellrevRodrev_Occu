##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##      Arctic Foxes in Varanger
##               --
##    Prepare data for DCOM (Fedini 2018)
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##***************************
##Execute from fit_models.R
##***************************
##

# rm(list=ls())
# 
# library(devtools)
# install_github("vqv/ggbiplot")
# 
# sapply(packages <- c('stringr', 'foreach', 'data.table', 'dplyr','tidyr', 'LaplacesDemon', 'runjags', 'rjags', 'coda', 'doParallel', 'ggbiplot'),
#        function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))
# 
# ##PARAMETERS TO SET MANUALY----
# YEARMIN = 2006
# YEARMAX = 2021
# Nobs_min = 3 #min number of days sampled to consider week valid
# Nweek_min = 4 #min number of weeks sampled a year to consider the year valid
# CUTOFF = 35 #min number of photo a day to consider the observation valid

data.path <- "data/main_dat"
data.filenames <- list.files(path = data.path, pattern = ".rds")
cov.path <- "data/covs/camera_attributes_2019.txt"
cov.path2 <- "data/covs/cam-coords-vegprod-hum.txt"
rodents.path <- "data/covs/storskala_04-21.txt"

allyears <- as.character(2005:2021)
years <- as.character(YEARMIN:YEARMAX)
years_id <- which(allyears %in% years)

control_region <- c("komag", "nyborg", "stjernevann")
##LOAD DATA -----
cat("##LOAD DATA ----- \n")

dat.l <- foreach(i = data.filenames[years_id]) %do%{
  readRDS(file = paste0(data.path, "/", i))%>%as.data.table()
}
allyears <- years
names(dat.l) <- years
dat.df <- rbindlist(dat.l, fill = TRUE)%>%
  filter(site %in% control_region)%>%
  filter(vis==1)

dat.df$site.year <- paste0(dat.df$loc, ".", dat.df$year)
dat.df$site.year <- factor(dat.df$site.year, levels = unique(dat.df$site.year))

##FILTER DATA ----
cat("##FILTER DATA ----- \n")
  
###condense to daily observations 
daily_dat <- dat.df%>%
  dplyr::group_by(site.year, year, loc, julian2)%>%
  dplyr::summarise(RedFox  = as.numeric(any(RedFox>0)),
                  ArcticFox  = as.numeric(any(ArcticFox>0)),
                  Bait = as.numeric(any(bait_corr>0)),
                  newbait = as.numeric(any(newbait==1)),
                  Npic = n())%>%
  dplyr::filter(Npic>=CUTOFF)

daily_dat <- daily_dat%>%
  dplyr::group_by(site.year)%>%
  dplyr::mutate(week = (julian2 - julian2[1]) %/% 7+1)%>%
  dplyr::group_by(site.year, week)%>%
  dplyr::mutate(day = 1:n())

### filter based on number of observation per weeks and number of weeks

weeks_tokeep <- daily_dat %>% 
  dplyr::group_by(site.year, week) %>% 
  dplyr::summarise(Nobs = length(unique(julian2)))%>%
  dplyr::filter(Nobs>=Nobs_min)

daily_dat <- daily_dat%>%
  dplyr::inner_join(weeks_tokeep)
  
sites_tokeep <- daily_dat %>% 
  dplyr::group_by(site.year) %>% 
  dplyr::summarise(Nweek = length(unique(week)))%>%
  dplyr::filter(Nweek>=Nweek_min)

daily_dat <- daily_dat%>%
  dplyr::inner_join(sites_tokeep)

###Get the occupancy state (1 = no animal, 2 = RF only, 3 = AF only, 4 = both, NA = NA)
daily_dat$state <- paste0(daily_dat$RedFox,daily_dat$ArcticFox)
categories <- data.frame(state = c("00", "10", "01", "11"),
                         lvls = c(seq(1:4)))
daily_dat <- dplyr::left_join(daily_dat, categories, by = "state")

daily_dat$lvls <- factor(daily_dat$lvls, levels = unique(daily_dat$lvls))

##GET  PARAMETERS ----
cat("##GET PARAMETERS ----- \n")
###Number of site.year
site_infos <- daily_dat%>%
  dplyr::group_by(site.year,loc,year)%>%
  dplyr::summarize()
M <- nrow(site_infos)

###Number of preliminary seasons
T <- nweeks

### Number of days per week
K <- 7

##GET DATA ARRAY M*T*K----
cat("##DATA AS ARRAY----- \n")
#fill daily_dat with NA for julian days not sampled 
daily_dat_all <- data.frame(site.year = rep(site_infos$site.year, each = K*T),
                           week = rep(1:T, each = K, M),
                           day= rep(1:7,T*M))
daily_dat <- merge(daily_dat,daily_dat_all, by = c("site.year", "week", "day"),all=T)

###Get array M*T*K
ob_state <- array(daily_dat$lvls, dim=c(K,T,M))
ob_state <- aperm(ob_state, c(3,2,1))

###Get minimal occupancy state
min_occ <- daily_dat%>%
  dplyr::group_by(site.year, week)%>%
  dplyr::summarise(RedFox = as.numeric(any(RedFox>0, na.rm=T)),
            ArcticFox = as.numeric(any(ArcticFox>0, na.rm=T)))

min_occ$state <- paste0(min_occ$RedFox, min_occ$ArcticFox)
min_occ <- left_join(min_occ, categories, by = "state")
min_occ$lvls[is.na(min_occ$lvls)] <- NA #all species potentially present

init <- matrix(min_occ$lvls, nrow = M, ncol = T, byrow = T)

###Get observation covariates 

bait <- array(c(daily_dat$Bait), dim=c(K,T,M))
bait <- aperm(bait, c(3,2,1))
bait <- bait + 1 #1: nobait 2: bait
bait[is.na(bait)] <- 1

##GET  COVARIATES ----
cat("##GET COVARIATES ----- \n")
###get site infos
covs <- read.table(file = cov.path,
                   header = T, row.names = 1)
covs$region <- str_to_lower(covs$region)

covs <- covs[site_infos$loc,]%>%
  mutate(site.year = site_infos$site.year,
         loc = site_infos$loc,
         year = factor(site_infos$year, levels = allyears),
         int = 1)%>%
  dplyr::select(c("site.year", "loc", "year", "region", "int", "altitude", "distroad", "distforest", "distcoast"))

covs$region[covs$region%in%c("nyborg-vj","komag")] <- "varangerS"
covs$region[covs$region%in%c("stjernevann")] <- "varangerN"
covs$region[covs$region%in%c("gaissene")] <- "ifjord"

rownames(covs) <- 1:nrow(covs)

###Get proportion of productive habitat

prod <- read.table(file = cov.path2,
                   header = T, row.names = 2)
row.names(prod) <- str_to_lower(row.names(prod))

covs$propprod5 <- prod[site_infos$loc,]$propprod5
covs$propprod10 <- prod[site_infos$loc,]$propprod10

###Get coord. on the PCA space based on geographical features

covs_pca <- prcomp(covs[,c("altitude","distroad", "distforest", "distcoast", "propprod5")], scale. = TRUE)
#ggbiplot(covs_pca)

covs <- cbind(covs, covs_pca$x[,c(1,2)])
##PC1: coast to land gradient (CLG)
##PC2: Forest to toundra gradient (FTG)
setnames(covs,c("PC1", "PC2"), c("CLG", "FTG"))

###get lemming abundance

rodents_data <- read.table(file = rodents.path,
                           header = T)

old_region_names <- c("bekkarfjord", "nordkynn", "ifjordfjellet", "komagdalen", "stjernevann", "vestre_jakobselv")
new_region_names <- c("bekkarfjord", "nordkyn", "ifjord", "varangerS", "varangerN", "varangerS")

rodents_data <- rodents_data%>%
  dplyr::rowwise()%>%
  dplyr::mutate(region = new_region_names[which(old_region_names == region)])

covs <- covs%>%
  dplyr::rowwise()%>%
  dplyr::mutate(rodents_fall = sum(apply(as.matrix(rodents_data[rodents_data$season=="fall"&
                                                          rodents_data$year+1==year&
                                                          rodents_data$region==region,5:9]),2,mean)),
         rodents_spr = sum(apply(as.matrix(rodents_data[rodents_data$season=="spring"&
                                                        rodents_data$year==year&
                                                        rodents_data$region==region,5:9]),2,mean)))
covs$rodents_mean <- (covs$rodents_fall+covs$rodents_spr)/2



covs[,c("altitude", "distroad", "distforest", "distcoast", "propprod5", "propprod10",  "CLG", "FTG", "rodents_fall",
        "rodents_spr","rodents_mean")] <- scale(covs[,c("altitude","distroad", "distforest","distcoast",
                                                        "propprod5", "propprod10", "CLG", "FTG",
                                                        "rodents_fall","rodents_spr", "rodents_mean")])


covs$year <- as.numeric(covs$year)



