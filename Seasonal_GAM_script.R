
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
                               min_max_route_years = 1,
                               model = model)

#need to add date column and convert to julian date

jags_data$date <- paste(jags_data$r_year,jags_data$month, 
                        jags_data$day, sep = "-") %>%
  ymd() %>%
  as.Date()

jags_data$yday <- as.numeric(format(jags_data$date, "%j"))
jags_data$yday <- jags_data$yday-(min(jags_data$yday)-1)



gam_output <- gam.basis.func(orig.preds = jags_data$yday,
                             npredpoints = max(jags_data$yday),
                             nknots = 6,
                             random = F,
                             sm_name = "day",
                             standardize = "range")


#model_to_file(model = model,filename = paste0(model,"add GAM.R"))



jags_GAM_data <- c(jags_data, gam_output[c(1,2,6)])

jags_GAM_data[["date"]] = NULL
jags_GAM_data[["yday"]] = NULL

####################################
# 3: Run JAGS model
####################################

#first run seasonal_GAM.R function
#defining new model

mod <- run_model(jags_data = jags_GAM_data,
                 model_file_path = "slope w GAM.R", ######new model with seasonal GAM
                 n_burnin = 20,
                 n_thin = 1,
                 n_iter=20,
                 #n_adapt = 1000,
                 parallel = T,
                 #inits = new.inits,
                 parameters_to_save = c("n","beta","BETA","gam.smooth"))
#new.inits = mod$mcmc.info$end.values
####################################
# 4: Model Analysis
####################################

# Calculate WAIC for the model (assuming "lambda" parameter was tracked)
#waic(jags_data = jags_data, jags_mod = mod)

# Generate and plot continental indices
cont_indices <- generate_cont_indices(mod)
pdf(paste0(species,"continental indices.pdf"))
c_plot <- plot_cont_indices(cont_indices)
print(c_plot)
dev.off()

# Generate and plot strata indices
strata_indices <- generate_strata_indices(mod)
s_plots <- plot_strata_indices(strata_indices)
pdf(paste0(species,"strata indices.pdf"))
for(j in 1:length(s_plots)){
print(s_plots[[j]])
}
dev.off()

# Generate continental trend
generate_cont_trend(indices = cont_indices)

# Generate strata-specific trends and plot on a map
s_trend <- generate_strata_trends(strata_indices)
map <- generate_map(s_trend, stratify_by = strat)
pdf(paste0(species,"trend map.pdf"))

print(map)
dev.off()

save.image("WOTH full image.RData")
save(list = c("mod","jags_data"),
     file = "WOTH model and data.RData")
