library(tidyr)
library(ggplot2)

# Read in all lorenz results
setwd('C:/Users/evsop/OneDrive/Documents/Master/4th Term/TFM/dynamical_fluid_forecast_TFM/Fluidic Pinball')
nb_trials = 500              # number of different trials/initial prediction points
nb_series = 2                # number of univariate time series
nb_time_steps = 50
pred_methods = c('NARX', 'NARMAX', 'SINDy')
noise_levels = c('n', '10', '20', '30')

# Calculate NRMSE distribution per coordinate
for (l in 1:length(pred_methods)) {
  for (k in 1:length(noise_levels)) {
    NRMSEs = matrix(0, nb_trials, nb_series)
    preds = readRDS(paste0('Results/', tolower(pred_methods[l]), '_preds_n', noise_levels[k],'.RData'))
    for(i in 1:nb_trials) {
      for(j in 1:nb_series) {
        obs_sd = sd(preds[[i]][[j]]$Y[1:nb_time_steps])
        NRMSE = sqrt(mean((preds[[i]][[j]]$yh[1:nb_time_steps] - preds[[i]][[j]]$Y[1:nb_time_steps])^2))/obs_sd
        NRMSEs[i, j] = NRMSE
      }
    }
    NRMSEs_df = data.frame(NRMSEs)
    colnames(NRMSEs_df) = c('Cl', 'Cd')
    
    # Calculation of mean, sd
    assign(paste0(tolower(pred_methods[l]), '_means_n', noise_levels[k]), apply(NRMSEs_df, 2, mean))
    assign(paste0(tolower(pred_methods[l]), '_vars_n', noise_levels[k]), apply(NRMSEs_df, 2, var))
    assign(paste0(tolower(pred_methods[l]), '_sd_n', noise_levels[k]), apply(NRMSEs_df, 2, sd))
    assign(paste0(tolower(pred_methods[l]), '_NRMSE_n', noise_levels[k]), NRMSEs_df)
  }
}

# Compare standard deviations
for (j in 1:nb_series) {
  df_joined = data.frame(c(narx_sd_nn[j], narx_sd_n10[j], narx_sd_n20[j], narx_sd_n30[j]),
                         c(narmax_sd_nn[j], narmax_sd_n10[j], narmax_sd_n20[j], narmax_sd_n30[j]),
                         c(sindy_sd_nn[j], sindy_sd_n10[j], sindy_sd_n20[j], sindy_sd_n30[j]))
  
  colnames(df_joined) = pred_methods
  
  NRMSEs_x_long = gather(df_joined, key="Method", value="NRMSE", 1:3)
  NRMSEs_x_long$Noise = c('0', '0.1', '0.2', '0.3')
  
  print(ggplot(NRMSEs_x_long, aes(x = factor(Noise, level = c('0', '0.1', '0.2', '0.3')), y = NRMSE, group = Method)) +
          theme_minimal() +
          theme(plot.title = element_text(size = 50,face="bold"), 
                axis.text=element_text(size=40),
                axis.title=element_text(size=45,face="bold"), 
                legend.text = element_text(size=40),
                legend.title = element_text(size = 45,face="bold")) +
          geom_line(aes(color=Method), size = 1.7) +
          geom_point(aes(color=Method), size = 1.7) +
          xlab('Noise Level') +
          ylab('SD of NRMSE')
  )
}

# Compare means
for (j in 1:nb_series) {
  df_joined = data.frame(c(narx_means_nn[j], narx_means_n10[j], narx_means_n20[j], narx_means_n30[j]),
                         c(narmax_means_nn[j], narmax_means_n10[j], narmax_means_n20[j], narmax_means_n30[j]),
                         c(sindy_means_nn[j], sindy_means_n10[j], sindy_means_n20[j], sindy_means_n30[j]))
  
  colnames(df_joined) = pred_methods
  
  NRMSEs_x_long = gather(df_joined, key="Method", value="NRMSE", 1:3)
  NRMSEs_x_long$Noise = c('0', '0.1', '0.2', '0.3')
  
  print(ggplot(NRMSEs_x_long, aes(x = factor(Noise, level = c('0', '0.1', '0.2', '0.3')), y = NRMSE, group = Method)) +
          theme_minimal() +
          theme(plot.title = element_text(size = 50,face="bold"), 
                axis.text=element_text(size=40),
                axis.title=element_text(size=45,face="bold"), 
                legend.text = element_text(size=40),
                legend.title = element_text(size = 45,face="bold")) +
          geom_line(aes(color=Method), size = 1.7) +
          geom_point(aes(color=Method), size = 1.7) +
          xlab('Noise Level') +
          ylab('Mean of NRMSE')
  )
}
