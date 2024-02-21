# Choose optimal parameters
opt_n = data.frame(ny=4, nu=4, ne=3, nl=4, rho_p1=1e-05, rho_p2=1e-09, rho_e1=1e-09, rho_e2=1e-09)

# Initiate lists to save predictions in
preds_nn = list()
preds_n10 = list()
preds_n20 = list()
preds_n30 = list()

# Create predictions for 100 different initial points
for (i in 1:nb_trials) {
  
  print(i)
  
  # Create training set
  if (initial_points[i] == 1) {
    train_set = as.matrix(xyzu_train)
  }
  
  if (initial_points[i] > 1) {
    train_set = as.matrix(rbind(xyzu_train, xyzu_pred[1:(initial_points[i]-1), ]))
  }
  
  # Create prediction set
  pred_set = as.matrix(xyzu_pred[initial_points[i]:(initial_points[i]+99), ])
  
  # Train model
  narmax_model = narmax(opt_n$ny, opt_n$nu, opt_n$ne, opt_n$nl)
  narmax_model = estimate.narmax(narmax_model, train_set[, 1:2], as.matrix(train_set[, 3]), c(opt_n$rho_p1, opt_n$rho_p2), c(opt_n$rho_e1, opt_n$rho_e2))
  print(narmax_model$terms)
  print(narmax_model$coefficients)
  print(narmax_model$nb_terms)
  
  ################
  ### No noise ###
  ################
  # Create predictions with no noise present
  preds_nn_aux = predict.narmax(narmax_model, pred_set[, 1:2], as.matrix(pred_set[, 3]), K=0)
  
  # Smoothing the prediction
  for(j in 1:nb_series) {
    pred = data.frame(pred = preds_nn_aux[[j]]$yh, time = seq(0.01, 1, 0.01))
    smoothed_pred = np.est(data = as.matrix(pred), h.seq = 0.08, newt = NULL, estimator = "LLP", kernel = "triweight")
    preds_nn_aux[[j]]$yh = smoothed_pred
  }
  
  preds_nn[[i]] = preds_nn_aux
  
  
  ##################
  ### 0.5% noise ###
  ##################
  # Add 0.5% of noise to initial condition
  pred_set_n10 = pred_set
  
  for(j in 1:nb_series) {
    for(k in 1:opt_n$ny) {
      pred_set_n10[k, j] = pred_set_n10[k, j] + 0.005*sd(pred_set[, j])
    }
  }
  
  # Create predictions
  preds_n10_aux = predict.narmax(narmax_model, pred_set_n10[, 1:2], as.matrix(pred_set_n10[, 3]), K=0)
  
  # Smoothing the prediction
  for(j in 1:nb_series) {
    pred = data.frame(pred = preds_n10_aux[[j]]$yh, time = seq(0.01, 1, 0.01))
    smoothed_pred = np.est(data = as.matrix(pred), h.seq = 0.08, newt = NULL, estimator = "LLP", kernel = "triweight")
    preds_n10_aux[[j]]$yh = smoothed_pred
  }
  
  preds_n10[[i]] = preds_n10_aux
  

  ################
  ### 1% noise ###
  ################
  # Add 1% of noise to initial condition
  pred_set_n20 = pred_set
  
  for(j in 1:nb_series) {
    for(k in 1:opt_n$ny) {
      pred_set_n20[k, j] = pred_set_n20[k, j] + 0.01*sd(pred_set[, j])
    }
  }
  
  # Create predictions
  preds_n20_aux = predict.narmax(narmax_model, pred_set_n20[, 1:2], as.matrix(pred_set_n20[, 3]), K=0)
  
  # Smoothing the prediction
  for(j in 1:nb_series) {
    pred = data.frame(pred = preds_n20_aux[[j]]$yh, time = seq(0.01, 1, 0.01))
    smoothed_pred = np.est(data = as.matrix(pred), h.seq = 0.08, newt = NULL, estimator = "LLP", kernel = "triweight")
    preds_n20_aux[[j]]$yh = smoothed_pred
  }
  
  preds_n20[[i]] = preds_n20_aux
  

  ################
  ### 2% noise ###
  ################
  # Add 2% of noise to initial condition
  pred_set_n30 = pred_set
  
  for(j in 1:nb_series) {
    for(k in 1:opt_n$ny) {
      pred_set_n30[k, j] = pred_set_n30[k, j] + 0.02*sd(pred_set[, j])
    }
  }
  
  # Create predictions
  preds_n30_aux = predict.narmax(narmax_model, pred_set_n30[, 1:2], as.matrix(pred_set_n30[, 3]), K=0)
  
  # Smoothing the prediction
  for(j in 1:nb_series) {
    pred = data.frame(pred = preds_n30_aux[[j]]$yh, time = seq(0.01, 1, 0.01))
    smoothed_pred = np.est(data = as.matrix(pred), h.seq = 0.08, newt = NULL, estimator = "LLP", kernel = "triweight")
    preds_n30_aux[[j]]$yh = smoothed_pred
  }
  
  preds_n30[[i]] = preds_n30_aux
}


#######################################
### Save prediction and model lists ###
#######################################
saveRDS(preds_nn, file="Results/narmax_preds_nn.RData")
saveRDS(preds_n10, file="Results/narmax_preds_n10.RData")
saveRDS(preds_n20, file="Results/narmax_preds_n20.RData")
saveRDS(preds_n30, file="Results/narmax_preds_n30.RData")
