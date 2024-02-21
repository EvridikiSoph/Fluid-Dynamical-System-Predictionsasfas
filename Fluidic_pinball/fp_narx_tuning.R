############################
### Tuning of parameters ###
############################
# Generate dataframe with all possible inputs
gs <- list(ny = c(3, 4),
           nu = c(3, 4),
           nl = c(3, 4),
           rho1 = c(0.000001, 0.00001, 0.0001, 0.001),
           rho2 = c(0.0000000001, 0.000000001, 0.00000001, 0.0000001)) %>%
  cross_df() # Convert to data frame grid

# Fit model with different parameters
MSEs = rep(0, length(gs))

# Calculate MAE between prediction and real value for different combinations
for (i in 1:nrow(gs)) {

  MSEs_aux = rep(0, 5)
  k = 0
  
  # Iterate through the five different segments for crossvalidation
  for (l in 1:5) {
    
    # Create training and test sets
    train_set = as.matrix(xyzu_train[1:(9200+k), ])
    val_set = as.matrix(xyzu_train[(9200+k+1):(9200+k+100), ])
    
    # Train model
    model_narx = narx(ny = gs[i, 1]$ny, nu = gs[i, 2]$nu, nl = gs[i, 3]$nl)
    
    tryCatch({
      model_narx = estimate.narx(model_narx, train_set[, 1:2], as.matrix(train_set[, 3]), c(gs[i, 4]$rho1, gs[i, 5]$rho2))
      print(model_narx$nb_terms)
      print(model_narx$terms)
    
      # Create predictions
      pred_narx = predict(model_narx, val_set[, 1:2], as.matrix(val_set[, 3]), K=0)
      MSE_sum = 0
      
      # Calculate MSE for each prediction and add results for all components
      for(j in 1:nb_series) {
        MSE = sqrt(mean((val_set[, j] - pred_narx[[j]]$yh)^2))/sd(val_set[, j])
        MSE_sum = MSE_sum + MSE
      }
      
      # Save to matrix
      MSEs_aux[l] = MSE_sum/nb_series
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

    # Take next step in crossvalidation
    k = k + 100
  }
  MSEs[i] = mean(MSEs_aux)
  print(MSEs[i])
}  


######################################
### Tuning of smoothing parameters ###
######################################
# Get optimal parameters
opt_n = data.frame(ny=3, nu=3, nl=3, rho1=1e-5, rho2=1e-9)

# Tuning smoothing parameters using PLR
MSEs_aux_n005 = rep(0, nb_series)
MSEs_aux_n01 = rep(0, nb_series)
MSEs_aux_n02 = rep(0, nb_series)

# Create training and test sets
train_set = as.matrix(xyzu_train[1:10097, ])
val_set = as.matrix(xyzu_train[10098:nrow(xyzu_train), ])

# Train model
narx_model = narx(opt_n$ny, opt_n$nu, opt_n$nl)
narx_model = estimate.narx(narx_model, train_set[, 1:nb_series], as.matrix(train_set[, (nb_series+1)]), c(opt_n$rho1, opt_n$rho2))


##################
### 0.5% noise ###
##################
# Add 0.5% of noise to initial condition
val_set_n005 = val_set

for(j in 1:nb_series) {
  for(f in 1:opt_n$ny) {
    val_set_n005[f, j] = val_set_n005[f, j] + 0.005*sd(val_set[, j])
  }
}

# Create predictions
pred_n005_aux = predict.narx(narx_model, val_set_n005[, 1:nb_series], as.matrix(val_set_n005[, (nb_series+1)]), K=0)

for(j in 1:nb_series) {
  
  pred = data.frame(pred = pred_n005_aux[[j]]$yh, time = seq(0.01, 3.03, 0.01))
  
  aux = np.cv(data = as.matrix(pred), h.seq = NULL, num.h = 1000, w = NULL, num.ln = 1,
              ln.0 = 0, step.ln = 2, estimator = "LLP", kernel = "triweight")
  
  MSEs_aux_n005[j] = aux$h.opt[2, ]
}


################
### 1% noise ###
################
# Add 1% of noise to initial condition
val_set_n01 = val_set

for(j in 1:nb_series) {
  for(f in 1:opt_n$ny) {
    val_set_n01[f, j] = val_set_n01[f, j] + 0.01*sd(val_set[, j])
  }
}

# Create predictions
pred_n01_aux = predict.narx(narx_model, val_set_n01[, 1:nb_series], as.matrix(val_set_n01[, (nb_series+1)]), K=0)

for(j in 1:nb_series) {
  
  pred = data.frame(pred = pred_n01_aux[[j]]$yh, time = seq(0.01, 3.03, 0.01))
  
  aux = np.cv(data = as.matrix(pred), h.seq = NULL, num.h = 1000, w = NULL, num.ln = 1,
              ln.0 = 0, step.ln = 2, estimator = "LLP", kernel = "triweight")
  
  MSEs_aux_n01[j] = aux$h.opt[2, ]
}


################
### 2% noise ###
################
# Add 2% of noise to initial condition
val_set_n02 = val_set

for(j in 1:nb_series) {
  for(f in 1:opt_n$ny) {
    val_set_n02[f, j] = val_set_n02[f, j] + 0.02*sd(val_set[, j])
  }
}

# Create predictions
pred_n02_aux = predict.narx(narx_model, val_set_n02[, 1:nb_series], as.matrix(val_set_n02[, (nb_series+1)]), K=0)

for(j in 1:nb_series) {
  
  pred = data.frame(pred = pred_n02_aux[[j]]$yh, time = seq(0.01, 3.03, 0.01))
  
  aux = np.cv(data = as.matrix(pred), h.seq = NULL, num.h = 1000, w = NULL, num.ln = 1,
              ln.0 = 0, step.ln = 2, estimator = "LLP", kernel = "triweight")
  print(aux$h.opt[2, ])

  MSEs_aux_n02[j] = aux$h.opt[2, ]
}
