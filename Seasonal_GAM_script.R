
# Building a BBS model that accounts for the date of survey using  --------

library(bbsBayes)
library(tidyverse)
library(lubridate)
source("GAM_basis_function.R")

####################################
# 2: Prepare Data
####################################

# Stratify the data
strat = "bbs_cws"
bbs_strat <- stratify(by = strat)
species = "Wood Thrush"
model = "slope"
# Prepare data for JAGS.
# This includes subsetting based on species of interest,
# adding in zeros, and performing model-specific calculations
jags_data <- prepare_jags_data(strat_data = bbs_strat,
                               species_to_run = species,
                               min_max_route_years = 3,
                               model = model)

#need to add date column and convert to julian date

#model_to_file("gam",filename = "gam_model_w_season.txt")


jags_data$date <- paste(jags_data$r_year,jags_data$month, 
                        jags_data$day, sep = "-") %>%
  ymd() %>%
  as.Date()

earliest_date <- paste(2018,5, 
                        24, sep = "-") %>%
  ymd() %>%
  as.Date()
as.numeric(format(earliest_date, "%j"))


latest_date <- paste(2018,7, 
                       7, sep = "-") %>%
  ymd() %>%
  as.Date()
as.numeric(format(latest_date, "%j"))



jags_data$yday <- as.numeric(format(jags_data$date, "%j"))
jags_data$yday <- jags_data$yday-(min(jags_data$yday)-1)



gam_output <- gam.basis.func(orig.preds = jags_data$yday,
                             npredpoints = max(jags_data$yday),
                             nknots = 6,
                             random = F,
                             sm_name = "day",
                             standardize = "range",
                             even_gaps = F)



#model_to_file(model = model,filename = paste0(model,"add GAM.R"))



jags_GAM_data <- c(jags_data, gam_output[c(1,2,6)])

jags_GAM_data[["date"]] = NULL
jags_GAM_data[["yday"]] = NULL

####################################
# 3: Run JAGS model
####################################

#first run seasonal_GAM.R function
#defining new model
t1 = Sys.time()
mod <- run_model(jags_data = jags_GAM_data,
                 model_file_path = "slope_model_w_GAM.txt", ######new model with seasonal GAM
                 n_burnin = 10000,
                 n_thin = 10,
                 n_iter = 10000,
                 #n_adapt = 1000,
                 parallel = T,
                 #inits = new.inits,
                 parameters_to_save = c("n","beta","BETA","vis.sm_day","beta_day","strata"))
t2 = Sys.time()
t2-t1


df = data.frame(mod$summary)
df$node = row.names(df)
names(df)[3:7] <- c("lci","lqrt","med","uqrt","uci") 

dfsm = df[paste0("vis.sm_day[",1:gam_output$npredpoints_day,"]"),]
dfsm$day = 1:gam_output$npredpoints_day
dfsm$mean_exp = exp(dfsm$mean)
dfsm$lci_exp = exp(dfsm$lci)
dfsm$uci_exp = exp(dfsm$uci)


smp = ggplot(data = dfsm,aes(x = day, y = mean_exp))+
  geom_line()+
  geom_ribbon(aes(ymin = lci_exp,ymax = uci_exp),alpha = 0.2)
x11()
print(smp)













# GAM version -------------------------------------------------------------


model = "gam"
# Prepare data for JAGS.
# This includes subsetting based on species of interest,
# adding in zeros, and performing model-specific calculations
jags_data <- prepare_jags_data(strat_data = bbs_strat,
                               species_to_run = species,
                               min_max_route_years = 3,
                               model = model)

#need to add date column and convert to julian date

#model_to_file("gam",filename = "gam_model_w_season.txt")


jags_data$date <- paste(jags_data$r_year,jags_data$month, 
                        jags_data$day, sep = "-") %>%
  ymd() %>%
  as.Date()

earliest_date <- paste(2018,5, 
                       24, sep = "-") %>%
  ymd() %>%
  as.Date()
as.numeric(format(earliest_date, "%j"))


latest_date <- paste(2018,7, 
                     7, sep = "-") %>%
  ymd() %>%
  as.Date()
as.numeric(format(latest_date, "%j"))



jags_data$yday <- as.numeric(format(jags_data$date, "%j"))
jags_data$yday <- jags_data$yday-(min(jags_data$yday)-1)



gam_output <- gam.basis.func(orig.preds = jags_data$yday,
                             npredpoints = max(jags_data$yday),
                             nknots = 6,
                             random = F,
                             sm_name = "day",
                             standardize = "range",
                             even_gaps = F)



#model_to_file(model = model,filename = paste0(model,"add GAM.R"))



jags_GAM_data <- c(jags_data, gam_output[c(1,2,6)])

jags_GAM_data[["date"]] = NULL
jags_GAM_data[["yday"]] = NULL

####################################
# 3: Run JAGS model
####################################

#first run seasonal_GAM.R function
#defining new model
t1 = Sys.time()
mod <- run_model(jags_data = jags_GAM_data,
                 model_file_path = "gam_model_w_season.txt", ######new model with seasonal GAM
                 n_burnin = 5000,
                 n_thin = 10,
                 n_iter = 5000,
                 #n_adapt = 1000,
                 parallel = T,
                 #inits = new.inits,
                 parameters_to_save = c("n","beta.X","B.X","vis.sm_day","beta_day","strata"))
t2 = Sys.time()
t2-t1


df = data.frame(mod$summary)
df$node = row.names(df)
names(df)[3:7] <- c("lci","lqrt","med","uqrt","uci") 

dfsm = df[paste0("vis.sm_day[",1:gam_output$npredpoints_day,"]"),]
dfsm$day = 1:gam_output$npredpoints_day
dfsm$mean_exp = exp(dfsm$mean)
dfsm$lci_exp = exp(dfsm$lci)
dfsm$uci_exp = exp(dfsm$uci)


smp = ggplot(data = dfsm,aes(x = day, y = mean_exp))+
  geom_line()+
  geom_ribbon(aes(ymin = lci_exp,ymax = uci_exp),alpha = 0.2)
x11()
print(smp)


