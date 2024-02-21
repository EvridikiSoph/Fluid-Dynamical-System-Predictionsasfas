############################
### Tuning of parameters ###
############################
gs <- list(degree = c(2,3),
           lambda = seq(0.1, 0.5, 0.1),
           lag = seq(2, 22, 1))  %>% 
  cross_df() # Convert to data frame grid

nb_of_lags = 3

# Fit model with different parameters
MSEs = rep(0, length(gs))

for (i in 1:nrow(gs)) {
  
  MSEs_aux = rep(0, 5)
  k = 0
  lags = c(1,gs$lag[i])
  max_lag = gs$lag[i] + 1
  
  # Iterate through the five different segments for cross-validation
  for (l in 1:5) {
    
    # Create training and test sets
    train_set = as.matrix(xyzu_train[1:(9200+k), ])
    
    # Add lags to trainset
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
    
    val_set = as.matrix(xyzu_train[(9200+k+1):(9200+k+100), ])
    
    x0 = val_set[1, 1:2]                # initial point
    x_lag = rbind(rep(0,2), train_set[(nrow(train_set)-(max_lag-3)):nrow(train_set), 1:nb_series])
    u = val_set[, 3]
    steps = nrow(val_set)

    # Train model 
    feat_mat_train = sindyr::features(train_set_new, polyorder = gs$degree[i])
    snd = sindy(train_set_new, Theta = feat_mat_train, lambda = gs$lambda[i], dt = dt)
    coeff = data.frame(snd$B)

    # Create prediction an real values
    xhat = matrix(0, steps, nb_series)
    xhat[1, ] = x0

    tryCatch({
      for (j in 1:(steps-1)) {
        xhat[j+1, ] = rk4u(sindy_pred, xhat[j, ], x_lag[j+1, ], u[j], dt, coeff, gs$degree[i])
        x_lag = rbind(x_lag, xhat[j+1, ])
      }
      
      MSE_sum = 0
    
      for(j in 1:nb_series) {
        MSE = sqrt(mean((val_set[, j] - xhat[, j])^2))/sd(val_set[, j])
        MSE_sum = MSE_sum + MSE
      }
    
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
# Choose optimal parameters
opt_n = data.frame(degree=2, lambda=0.5)
lags = c(1,18) # First and 18th lag are going to be used in the prediction
nb_of_lags = 3 # number of lags plus 1
max_lag = 19   # maximum lag plus 1

#################################################################################
# Tuning smoothing parameters using PLR
MSEs_aux_n10 = rep(0, nb_series)
MSEs_aux_n20 = rep(0, nb_series)
MSEs_aux_n30 = rep(0, nb_series)

# Create training and test sets
train_set = as.matrix(xyzu_train[1:10097, ])
val_set = as.matrix(xyzu_train[10098:nrow(xyzu_train), ])

# Add lags to trainset
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

# Train model 
feat_mat_train = sindyr::features(train_set_new, polyorder = opt_n$degree)
snd = sindy(train_set_new, Theta = feat_mat_train, lambda = opt_n$lambda, dt = dt)
coeff = data.frame(snd$B)

x0 = val_set[1, 1:2]                # initial point
u = val_set[, 3]
steps = 100              # number of time steps to take


##################
### 10% noise ###
##################
# Add 10% of noise to initial point
x0n10 = rep(0, length(x0))
x_lag_n10 = rbind(rep(0,2), val_set[(nrow(val_set)-(max_lag-3)):nrow(val_set), 1:nb_series])

for (j in 1:nb_series) {
  x0n10[j] = x0[j] + 0.1*sd(val_set[, j])
}

# Create prediction an real values
xhat_n10 = matrix(0, steps, nb_series)
xhat_n10[1, ] = x0n10

for (j in 1:(steps-1)) {
  xhat_n10[j+1, ] = rk4u(sindy_pred, xhat_n10[j, ], x_lag_n10[j+1, ], u[j], dt, coeff, opt_n$degree)
  x_lag = rbind(x_lag_n10, xhat_n10[j+1, ])
}

# Creation of final dataframes
pred_n10_aux = list()

for(k in 1:nb_series) {
  pred_n10_aux[[k]] = data.frame(cbind(xhat_n10[, k], val_set[1:steps, k]))
  colnames(pred_n10_aux[[k]]) = c('pred', 'real')
}

for(j in 1:nb_series) {
  
  pred = data.frame(pred = pred_n10_aux[[j]]$pred, time = seq(0.01, 1, 0.01))
  
  aux = np.cv(data = as.matrix(pred), h.seq = NULL, num.h = 1000, w = NULL, num.ln = 1,
              ln.0 = 0, step.ln = 2, estimator = "LLP", kernel = "triweight")
  
  MSEs_aux_n10[j] = aux$h.opt[2, ]
}


################
### 20% noise ###
################
# Add 20% of noise to initial condition
x0n20 = rep(0, length(x0))
x_lag_n20 = rbind(rep(0,2), val_set[(nrow(val_set)-(max_lag-3)):nrow(val_set), 1:nb_series])

for (j in 1:nb_series) {
  x0n20[j] = x0[j] + 0.2*sd(val_set[, j])
}

# Create prediction an real values
xhat_n20 = matrix(0, steps, nb_series)
xhat_n20[1, ] = x0n20

for (j in 1:(steps-1)) {
  xhat_n20[j+1, ] = rk4u(sindy_pred, xhat_n20[j, ], x_lag_n20[j+1, ], u[j], dt, coeff, opt_n$degree)
  x_lag_n20 = rbind(x_lag, xhat_n20[j+1, ])
}

# Creation of final dataframes
pred_n20_aux = list()

for(k in 1:nb_series) {
  pred_n20_aux[[k]] = data.frame(cbind(xhat_n20[, k], val_set[1:steps, k]))
  colnames(pred_n20_aux[[k]]) = c('pred', 'real')
}

for(j in 1:nb_series) {
  
  pred = data.frame(pred = pred_n20_aux[[j]]$pred, time = seq(0.01, 1, 0.01))
  
  aux = np.cv(data = as.matrix(pred), h.seq = NULL, num.h = 1000, w = NULL, num.ln = 1,
              ln.0 = 0, step.ln = 2, estimator = "LLP", kernel = "triweight")
  
  MSEs_aux_n20[j] = aux$h.opt[2, ]
}


################
### 30% noise ###
################
# Add 30% of noise to initial condition
x0n30 = rep(0, length(x0))
x_lag_n30 = rbind(rep(0,2), val_set[(nrow(val_set)-(max_lag-3)):nrow(val_set), 1:nb_series])

for (j in 1:nb_series) {
  x0n30[j] = x0[j] + 0.3*sd(val_set[, j])
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
  pred_n30_aux[[k]] = data.frame(cbind(xhat_n30[, k], val_set[1:steps, k]))
  colnames(pred_n30_aux[[k]]) = c('pred', 'real')
}

for(j in 1:nb_series) {
  
  pred = data.frame(pred = pred_n30_aux[[j]]$pred, time = seq(0.01, 1, 0.01))
  
  aux = np.cv(data = as.matrix(pred), h.seq = NULL, num.h = 1000, w = NULL, num.ln = 1,
              ln.0 = 0, step.ln = 2, estimator = "LLP", kernel = "triweight")
  
  MSEs_aux_n30[j] = aux$h.opt[2, ]
}
