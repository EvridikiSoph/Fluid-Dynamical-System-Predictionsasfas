# Choose optimal parameters
opt_n = data.frame(degree=2, lambda=0.5)
lags = c(1,18) # First and 18th lag are going to be used in the prediction
nb_of_lags = 3 # number of lags plus 1
max_lag = 19   # maximum lag plus 1

# Initiate lists to save predictions in
preds_nn = list()
preds_n10 = list()
preds_n20 = list()
preds_n30 = list()

# Create predictions for n different initial points
for (i in 1:nb_trials) {
  
  print(i)
  
  # Create training set
  if (initial_points[i] == 1) {
    train_set = as.matrix(xyzu_train)
  }
  
  if (initial_points[i] > 1) {
    train_set = as.matrix(rbind(xyzu_train, xyzu_pred[1:(initial_points[i]-1), ]))
  }
  
  # Add 18th lag to trainset
  train_set_new = matrix(0, nrow=nrow(train_set)-max_lag+1, ncol = 2*(ncol(train_set)-1)+1)
  
  k=0
  col_names=NULL
  
  for(j in 1:(nb_of_lags-1)) {
    train_set_new[, (k + j):(j + k + nb_series - 1)] = train_set[(max_lag - lags[j]):(nrow(train_set) - lags[j]), 1:2]
    col_names = c(col_names, paste0(colnames(train_set)[1:2],j))
    k = k + nb_series - 1
  }
  
  # Add u
  train_set_new[, ((nb_of_lags - 1)*nb_series + 1)] = train_set[(max_lag - 1):(nrow(train_set) - 1), 3]
  col_names = c(col_names, "u")
  
  colnames(train_set_new) = col_names
  
  # Add sinus of each term as feature
  for(j in 1:(nb_series*(nb_of_lags-1))) {
    train_set_new = cbind(train_set_new, sin(train_set_new[, j]))
    col_names = c(col_names, paste0("sin(", colnames(train_set_new)[j], ")"))
  }
  
  # Add cosinus of each term as feature
  for(j in 1:(nb_series*(nb_of_lags-1))) {
    train_set_new = cbind(train_set_new, cos(train_set_new[, j]))
    col_names = c(col_names, paste0("cos(", colnames(train_set_new)[j], ")"))
  }
  
  colnames(train_set_new) = col_names
  
  # Create prediction set
  pred_set = as.matrix(xyzu_pred[initial_points[i]:(initial_points[i]+99), ])
  x0 = pred_set[1, 1:2]                # initial point
  u = pred_set[, 3]
  steps = 100                          # number of time steps to take
  
  # Create training feature matrix
  feat_mat_train = sindyr::features(train_set_new, polyorder = opt_n$degree)
  
  # Train model
  sindy_model = sindy(train_set_new, dt = dt, Theta = feat_mat_train, lambda = opt_n$lambda)
  coeff = data.frame(sindy_model$B)


  ################
  ### No noise ###
  ################
  # Create prediction an real values
  xhat_nn = matrix(0, steps, nb_series)
  xhat_nn[1, ] = x0
  x_lag_nn = rbind(rep(0,2), train_set[(nrow(train_set)-(max_lag-3)):nrow(train_set), 1:nb_series])

  for (j in 1:(steps-1)) {
    xhat_nn[j+1, ] = rk4u(sindy_pred, xhat_nn[j, ], x_lag_nn[j+1, ], u[j], dt, coeff, opt_n$degree)
    x_lag_nn = rbind(x_lag_nn, xhat_nn[j+1, ])
  }
  
  # Creation of final dataframes
  pred_nn_aux = list()
  
  for(k in 1:nb_series) {
    pred_nn_aux[[k]] = data.frame(cbind(xhat_nn[, k], pred_set[, k]))
    colnames(pred_nn_aux[[k]]) = c('yh', 'Y')
  }
  
  # Smoothing the prediction
  for(j in 1:nb_series) {
    pred = data.frame(pred = xhat_nn[, j], time = seq(0.001, 0.1, 0.001))
    smoothed_pred = np.est(data = as.matrix(pred), h.seq = 0.03, newt = NULL, estimator = "LLP", kernel = "triweight")
    pred_nn_aux[[j]]$yh = smoothed_pred
  }
  
  preds_nn[[i]] = pred_nn_aux
  
  
  ##################
  ### 10% noise ###
  ##################
  # Add 10% of noise to initial point
  x0n10 = rep(0, length(x0))
  x_lag_n10 = rbind(rep(0,2), train_set[(nrow(train_set)-(max_lag-3)):nrow(train_set), 1:nb_series])
  
  for (j in 1:nb_series) {
    x0n10[j] = x0[j] + 0.1*sd(pred_set[, j])
  }
  
  # Create prediction an real values
  xhat_n10 = matrix(0, steps, nb_series)
  xhat_n10[1, ] = x0n10

  for (j in 1:(steps-1)) {
    xhat_n10[j+1, ] = rk4u(sindy_pred, xhat_n10[j, ], x_lag_n10[j+1, ], u[j], dt, coeff, opt_n$degree)
    x_lag_n10 = rbind(x_lag_n10, xhat_n10[j+1, ])
  }
  
  # Creation of final dataframes
  pred_n10_aux = list()
  
  for(k in 1:nb_series) {
    pred_n10_aux[[k]] = data.frame(cbind(xhat_n10[, k], pred_set[, k]))
    colnames(pred_n10_aux[[k]]) = c('yh', 'Y')
  }
  
  # Smoothing the prediction
  for(j in 1:nb_series) {
    pred = data.frame(pred = xhat_n10[, j], time = seq(0.001, 0.1, 0.001))
    smoothed_pred = np.est(data = as.matrix(pred), h.seq = 0.03, newt = NULL, estimator = "LLP", kernel = "triweight")
    pred_n10_aux[[j]]$pred_smooth = smoothed_pred
  }
  
  preds_n10[[i]] = pred_n10_aux
  
    
  #################
  ### 20% noise ###
  #################
  # Add 20% of noise to initial point
  x0n20 = rep(0, length(x0))
  x_lag_n20 = rbind(rep(0,2), train_set[(nrow(train_set)-(max_lag-3)):nrow(train_set), 1:nb_series])
  
  for (j in 1:nb_series) {
    x0n20[j] = x0[j] + 0.2*sd(pred_set[, j])
  }
  
  # Create prediction an real values
  xhat_n20 = matrix(0, steps, nb_series)
  xhat_n20[1, ] = x0n20

  for (j in 1:(steps-1)) {
    xhat_n20[j+1, ] = rk4u(sindy_pred, xhat_n20[j, ], x_lag_n20[j+1, ], u[j], dt, coeff, opt_n$degree)
    x_lag_n20 = rbind(x_lag_n20, xhat_n20[j+1, ])
  }
  
  # Creation of final dataframes
  pred_n20_aux = list()
  
  for(k in 1:nb_series) {
    pred_n20_aux[[k]] = data.frame(cbind(xhat_n20[, k], pred_set[, k]))
    colnames(pred_n20_aux[[k]]) = c('yh', 'Y')
  }
  
  # Smoothing the prediction
  for(j in 1:nb_series) {
    pred = data.frame(pred = xhat_n20[, j], time = seq(0.001, 0.1, 0.001))
    smoothed_pred = np.est(data = as.matrix(pred), h.seq = 0.03, newt = NULL, estimator = "LLP", kernel = "triweight")
    pred_n20_aux[[j]]$yh = smoothed_pred
  }
  
  preds_n20[[i]] = pred_n20_aux
  
    
  ################
  ### 30% noise ###
  ################
  # Add 30% of noise to initial point
  x0n30 = rep(0, length(x0))
  x_lag_n30 = rbind(rep(0,2), train_set[(nrow(train_set)-(max_lag-3)):nrow(train_set), 1:nb_series])
  
  for (j in 1:nb_series) {
    x0n30[j] = x0[j] + 0.3*sd(pred_set[, j])
  }
  
  # Create prediction an real values
  xhat_n30 = matrix(0, steps, nb_series)
  xhat_n30[1, ] = x0n30

  for (j in 1:(steps-1)) {
    xhat_n30[j+1, ] = rk4u(sindy_pred, xhat_n30[j, ], x_lag_n30[j+1, ], u[j], dt, coeff, opt_n$degree)
    x_lag_n30 = rbind(x_lag_n30, xhat_n30[j+1, ])
  }
  
  # Creation of final dataframes
  pred_n30_aux = list()
  
  for(k in 1:nb_series) {
    pred_n30_aux[[k]] = data.frame(cbind(xhat_n30[, k], pred_set[, k]))
    colnames(pred_n30_aux[[k]]) = c('yh', 'Y')
  }
  
  # Smoothing the prediction
  for(j in 1:nb_series) {
    pred = data.frame(pred = xhat_n30[, j], time = seq(0.001, 0.1, 0.001))
    smoothed_pred = np.est(data = as.matrix(pred), h.seq = 0.03, newt = NULL, estimator = "LLP", kernel = "triweight")
    pred_n30_aux[[j]]$yh = smoothed_pred
  }
  
  preds_n30[[i]] = pred_n30_aux
}


#######################################
### Save prediction and model lists ###
#######################################
saveRDS(preds_nn, file="Results/sindy_preds_nn.RData")
saveRDS(preds_n10, file="Results/sindy_preds_n10.RData")
saveRDS(preds_n20, file="Results/sindy_preds_n20.RData")
saveRDS(preds_n30, file="Results/sindy_preds_n30.RData")
