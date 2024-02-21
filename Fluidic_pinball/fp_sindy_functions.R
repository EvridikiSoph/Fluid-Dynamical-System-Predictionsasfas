###########################
### Necessary functions ###
###########################
# Sindy prediction
# Lorenz simulation
# Function that simulates a Lorenz model
# y: initial condition
# u: exogenous parameter
sindy_pred = function(x, x_lag, u, coeff, order) {
  
  col_names = c("Cl1", "Cd1", "Cl2", "Cd2", "u")
  x_pred = data.frame(t(x), t(x_lag), u)
  colnames(x_pred) = col_names

  # Add sinus of each term as feature
  for(j in 1:((nb_of_lags-1)*nb_series)) {
    x_pred = cbind(x_pred, sin(x_pred[, j]))
    col_name = paste0('sin(', colnames(x_pred)[j], ')')
    col_names = c(col_names, col_name)
  }
  
  # Add cosinus of each term as feature
  for(j in 1:((nb_of_lags-1)*nb_series)) {
    x_pred = cbind(x_pred, cos(x_pred[, j]))
    col_name = paste0('cos(', colnames(x_pred)[j], ')')
    col_names = c(col_names, col_name)
  }
  
  colnames(x_pred) = col_names
  
  # Creates data matrix with terms up to degree x
  features = sindyr::features(x_pred, polyorder = order)
  
  pred = c(features %*% coeff[, 1],
           features %*% coeff[, 2])
  
  return(pred)
}


# Integration function
# Function that uses Runge-Kutta to integrate the result of a given function
# for a given amount ot time steps ahead
# v: function the result of which needs to be integrated
# y: initial condition
# u: exogenous parameter of length n
# h: size of time step
# p: further parameters of function v
rk4u = function(v, y, y_lag, u, h, p, r) { 
  
  # RK4U   Runge-Kutta scheme of order 4 for control system
  k1 = v(y, y_lag, u, p, r)
  k2 = v(y + (h/2)*k1, y_lag, u, p, r)
  k3 = v(y + (h/2)*k2, y_lag, u, p, r)
  k4 = v(y + h*k3, y_lag, u, p, r)
  y = y + h*(k1 + 2*k2 + 2*k3 + k4)/6

  return(y)
}


# NRMSE statistic function
# Description: calculates the NRMSE btw real and predictedvalues for a certain amount of observations 
# (used for bootstraping)
# data: dataset with one column for real and one for predicted values
# indices: amount of rows to be used in calculation
# Output: NRMSE value
NRMSE_function <- function(data, indices) {
  
  d <- data[indices, ]
  
  obs_sd = sd(d[, 1])
  pred_errors = (d[, 1] - d[, 2])^2
  
  NRMSE = sqrt(mean(pred_errors))/obs_sd
  
  NRMSE
}
