##############
### Models ###
##############

# ARX Model
# Description: Creates an autoregressive with exogenous inputs model
# ny: Number of autoregressive lags
# nu: Number of input lags
# Output: Object representing an ARX mode
arx = function (ny, nu) {
  model = list(
    ny = ny,
    nu = nu,
    maxLag = max(ny, nu) + 1,
    terms = NULL,
    coefficients = NULL, # Variable name required by base::
    nb_terms = NULL,
    call = match.call()
  )
  class(model) = 'arx'
  return(model)
}

# NARMAX Model
# Description: Creates a nonlinear autogressive with exogenous inputs model
# ny: Number of autoregressive lags
# nu: Number of input lags
# nl: Nonlinearity polynomial length
# Output: Object representing a NARX model
narx = function (ny, nu, nl) {
  model = list(
    ny = ny,
    nu = nu,
    nl = nl,
    maxLag = max(ny, nu) + 1,
    terms = NULL,
    coefficients = NULL,
    nb_terms = NULL,
    call = match.call()
  )
  class(model) = 'narx'
  return(model)
}


####################################
### Regression matrix generators ###
####################################

sizeGuard = function (Y, U, E = NULL) {
  if (!is.null(U) && nrow(Y) != nrow(U)) stop('Input-Output vector must have the same size')
  if (!is.null(E)) {
    if (nrow(Y) != nrow(E)) stop('Error-Output vector must have the same size')
  }
}

# Subset matrix
# Description: Subset a matrix like mat[rows, cols]. Return a matrix even when the
#              subset generates a lower dimensional structure, i.e., never loose rownames or colnames.
# mat: Target matrix
# rows: Rows subsetter
# cols: Cols subsetter
subsetMatrix = function (mat, rows, cols) {
  names = dimnames(mat)
  rows = if (is.null(rows)) 1:nrow(mat) else rows
  cols = if (is.null(cols)) 1:ncol(mat) else cols
  rownames = if (typeof(rows) == 'character') rows else names[[1]][rows]
  colnames = if (typeof(cols) == 'character') cols else names[[2]][cols]
  submat = matrix(
    mat[rows, cols],
    nrow = length(rows),
    ncol = length(cols),
    dimnames = list(rownames, colnames)
  )
  return(submat)
}

# General function creation
genRegMatrix = function (model, ...) UseMethod('genRegMatrix', model)

genRegMatrix.default = function (model, ...) {
  stop(sprintf('Unknown class %s for regression matrix generation', class(model)))
}

# Description: Generates a regression matrix for an ARMAX Model
# model: ARMAX Model
# Y: The target matrix
# U: The input vector
# Output: Object containing:
#         P: Regression matrix with all terms
#         Pp: Regression matrix with only process terms
#         Pnp: Regression matrix without process terms
genRegMatrix.arx = function (model, Y, U, E = NULL) {
  na = model$ny
  nb = model$nu
  p = model$maxLag
  N = nrow(Y)
  M = ncol(Y)
  
  # Check if sizes of Y, U and E are compatible
  sizeGuard(Y,U, E = NULL)
  
  # Initialize empty matrix that will contain potential termans of the model
  obj = list()
  #obj$P = matrix(0, nrow = N - p + 1, ncol = 2*M*na + nb)
  obj$P = matrix(0, nrow = N - p + 1, ncol = M*na + nb)
  
  colPhi = NULL         # Will contain column names
  
  #Y_der = matrix(0, N, M)
  
  #for(j in 1:M) {
  #  lk = glkerns(seq(0.01, N*0.01, 0.01), Y[, j], deriv = 1, n.out=N)
  #  Y_der[, j] = lk$est
  #}
  
  # Iterate through lags of Y and add them to the regression matrix
  # for all time series present 
  k = 0
  
  for(i in 1:na[1]) {
    obj$P[, (k + i):(i + k + M - 1)] = -Y[(p - i):(N - i),]
    for(j in 1:M) {
      colPhi = c(colPhi, paste0("-y(k-",i,")",j))
    }
    k = k + M - 1
  }
  
  # Iterate through derivatives of the lags of Y and add them to the regression matrix
  # for all time series present 
  #k = M*na
  
  #for(i in 1:na[1]) {
  #  obj$P[, (k + i):(i + k + M - 1)] = -Y_der[(p - i):(N - i), ]
  #  for(j in 1:M) {
  #    colPhi = c(colPhi, paste0("-dy(k-",i,")",j))
  #  }
  #  k = k + M - 1
  #}
  
  # Iterate through lags of U and add them to the regression matrix
  for(i in 1:nb[1]) {
    #obj$P[, (2*M*na + i)] = U[(p - i):(N - i)]
    obj$P[, (M*na + i)] = U[(p - i):(N - i)]
    colPhi = c(colPhi, paste0("u(k-",i,")"))
  }
  
  # Give each feature in the regression matirx its corresponding name
  rowPhi = paste0(rep("k=", N - p + 1), p:N)
  
  colnames(obj$P) = colPhi
  rownames(obj$P) = rowPhi
  
  # Create final complete regression matrix
  if (!is.null(model$terms)) {
    obj$P = subsetMatrix(obj$P, NULL, model$terms)
  }
  
  obj$Pp = obj$P
  return(obj)
}


# Description: Generates a regression matrix for an NARX Model
# model: NARX Model
# Y: The target matrix
# U: The input vector
# Output: Object containing:
#         P: Regression matrix with all terms
#         Pp: Regression matrix with only process terms
#         Pnp: Regression matrix without process terms
genRegMatrix.narx = function (model, Y, U, E = NULL) {
  ny = model$ny
  nu = model$nu
  nl = model$nl
  N = nrow(Y)
  M = ncol(Y)
  #n = 2*M*ny + nu
  n = M*ny + nu
  p = model$maxLag
  
  # Check if sizes of Y, U and E are compatible
  sizeGuard(Y, U, E = NULL)
  obj = list()
  
  # Create a list of all possible combinations of 1:nl, as a basis for the polynomial
  auxExp = list()
  candList = list()
  
  for (i in 1:nl[1]) {
    eval(parse(text = paste0('auxExp$x', i, '= 1:n')))
    cand = expand.grid(auxExp)
    candList[[i]] = unique(t(apply(cand, 1, sort)))
  }

  # Generate the regression matrix for armax
  P0 = genRegMatrix(arx(ny, nu), Y, U, E)$P
  P0[, 1:M*ny] = -P0[, 1:M*ny]
  #P0[, 1:2*M*ny] = -P0[, 1:2*M*ny]
  
  # Assign column names
  col_names = NULL
  
  for(j in 1:M) {
    col_names = c(col_names, paste0('y(k-', 1:ny[1], ')',j))
  }
  
  #for(j in 1:M) {
  #  col_names = c(col_names, paste0('dy(k-', 1:ny[1], ')',j))
  #}
  
  col_names = c(col_names, paste0('u(k-', 1:nu[1], ')'))
  
  # Add constant term to polynomial
  NP0 = nrow(P0)
  obj$P = NULL
  obj$P = cbind(rep(1, NP0), obj$P) # Constant
  colnames(obj$P) = "constant"
  obj$P = cbind(obj$P, P0)
  
  # Combine linear terms of arx feature matrix into nonlinear ones of order max nl
  if (nl >= 2) {
    for (i in 2:nl) {
      ncand = nrow(candList[[i]])
      for (j in 1:ncand[1]) {
        Pcand_a = subsetMatrix(P0, NULL, candList[[i]][j, ]) # P0[, candList[[i]][j, ]]
        names = colnames(Pcand_a)
        Pcand_b = matrix(apply(Pcand_a, 1, prod), ncol = 1)
        colnames(Pcand_b) = stringr::str_c(names, collapse = '')
        obj$P = cbind(obj$P, Pcand_b)
      }
    }
  }
  
  # Create final complete regression matrix
  if (!is.null(model$terms)) {
    obj$P = subsetMatrix(obj$P, NULL, model$terms)
  }
  
  obj$Pp = obj$P
  obj$Pnp = NULL

  return(obj)
}


# Title: Generate target vector
# Description: Generate target vector based on maximum lag. Works with any model
# model: Any model containing a $maxLag property
# Y: Original target matrix
genTarget = function (model, ...) UseMethod('genTarget')

genTarget.default = function (model, Y) {
  N = nrow(Y)
  M = ncol(Y)
  p = model$maxLag
  target = matrix(Y[p:N,], ncol = M)
  rownames(target) = paste0(rep("k=", N - p + 1), p:N)
  col_names = NULL
  for(j in 1:M) {
    col_names = c(col_names, paste0("y(k)",j))
  }
  colnames(target) = col_names
  return(target)
}

###########################
### Estimation function ###
###########################

# FROLS Algorithm
# Description: Forward Regression with Orthogonal Least Squares
# P: Matrix where each column is a potential feature
# Y: The target matrix
# rho: Stop criteria (run while (1 - sum(EER)) > rho)
# Output: List containing:
#         List of coefficient vector for each time series
#         List of selected features matrix for each time series
#         List of g vector for each time series
#         List of orthogonal matrix for each time series
#         List of upper-triangular matrix for each time series
#         List of error vector for each time series
frols = function (P, Y, rho) {
  MP = ncol(P)
  M = ncol(Y)
  NP = nrow(P)
  
  # Initialize lists to save the final results for each time series
  A_lst = list()
  gvec_lst = list()
  Qs_lst = list()
  ERRvec_lst = list()
  selectTerms = list()
  th_FROLS_lst = list()
  Psel_lst = list()
  
  for (j in 1:M) {
    # 1st step ----------------------------------------------------------------
    sig = Y[, j] %*% Y[, j]
    selectTerms = NULL
    ERRvec = NULL
    gvec = NULL
    
    Qs = P
    g = rep(0, MP)
    ERR = rep(0, MP)
    for (m in 1:MP) {
      g[m] = (Y[, j] %*% Qs[, m]) / (Qs[, m] %*% Qs[, m])
      ERR[m] = (g[m] ^ 2 * (Qs[, m] %*% Qs[, m])) / sig
    }
    
    l1 = which(ERR == max(ERR))[1]
    selectTerms = l1[1] # vector keeping all selected terms
    
    # init
    # A = diag(1,M)
    A = 1
    Qs = matrix(P[, l1], ncol = 1)
    gvec = g[l1]
    ERRvec = ERR[l1]
      
    # s-th step --------------------------------------------------------------
    for (s in 2:MP) {
      gm  = rep(0, MP)
      Qm  = matrix(0, NP, MP)
      ERR = rep(0, MP)
      A = cbind(rbind(A, 0), 0)
      
      for (m in (1:MP)[-selectTerms]) {
        sumQm = rep(0, NP)
        for (r in 1:(s - 1)){
          sumQm = sumQm + ((P[, m] %*% Qs[, r]) /  (Qs[, r] %*% Qs[, r])) %*% Qs[, r]
        }
        Qm[, m] = P[, m] - sumQm
        gm[m] = (Y[, j] %*% Qm[, m]) / (Qm[, m] %*% Qm[, m])
        ERR[m] = (gm[m] ^ 2 * (Qm[, m] %*% Qm[, m])) / sig
        #print(m)
        #print(colnames(P)[m])
        #print(ERR[m])
      }
      
      ls = which(ERR == max(ERR))
      selectTerms = cbind(selectTerms, ls) # vector keeping all selected terms
      
      Qs = cbind(Qs, Qm[, ls]) # keep set of orthogonal bases
      gvec = rbind(gvec, gm[ls])
      for (r in 1:(s - 1)){
        A[r, s] = (Qs[, r] %*% P[, ls]) / (Qs[, r] %*% Qs[, r])
      }
      A[s, s] = 1
      
      # Assing results to lists
      A_lst[[j]] = A
      gvec_lst[[j]] = gvec
      Qs_lst[[j]] = Qs
      ERRvec_lst[[j]] = ERRvec
      
      ERRvec = rbind(ERRvec, ERR[ls])
      ESR = 1 - sum(ERRvec)

      if (ESR <= rho[j]){
        M0 = s
        break
      }
    }
  
    th_FROLS = solve(A, gvec)
    th_FROLS_lst[[j]] = th_FROLS[, ]
    Psel = P[, selectTerms]
    Psel_lst[[j]] = Psel
  }
  
  return(list(
    th = th_FROLS_lst,
    Psel = Psel_lst,
    g = gvec_lst,
    W = Qs_lst,
    A = A_lst,
    ERR = ERRvec_lst
  ))
}


# Estimate NARX model
# Description: Uses the FROLS algorithm to estimate the models coefficients
# model: NARX model
# Y: training time series
# U: training input
# rho: FROLS threshold fo the selection of the model terms
# Output: model with defined coefficients and model terms
estimate.narx = function (model, Y, U, rho) {
  
  M = ncol(Y)
  P = genRegMatrix(model, Y, U, E)$P
  Target = genTarget(model, Y)
  
  resultNarx = frols(P, Target, rho)
  
  terms = NULL
  nb_terms = NULL
  
  for(j in 1:M) {
    terms =  c(terms, colnames(resultNarx$Psel[[j]]))
    nb_terms = c(nb_terms, length(colnames(resultNarx$Psel[[j]])))
  }
  
  model$terms = terms
  model$nb_terms = nb_terms
  model$coefficients = resultNarx$th

  return(model)
}


############################
### Prediction functions ###
############################

# Description: predictions using a NARX model, by choosing the prediction method depending on the value of K
# model: NARX model
# Y: time series to be predicted
# U: exogenous input
# K: number indicating the steps ahead to take
# Output: list containing predictions for each time series
predict.narx = function (model, Y, U, K = 1, ...) {
  cat('Running narx prediction ... \n')
  prediction = predict.default(model, Y, U, K)
  return(prediction)
}

# Prediction default function
predict.default = function (model, Y, U, K = 1, ...) {
  if (K < 0) stop('K must be greater or equal to zero')
  method = switch(
    as.character(K),
    "1" = oneStepAhead,
    "0" = freeRun,
    kStepAhead
  )
  return(method(model, Y, U, K))
}

# Description: predicts next time step of time series by using previous predictions
# model: model
# Y: time series to be predicted
# U: exogenous parameter input
# ignore K
# Output: list containing predictions for each time series
freeRun = function (model, Y, U, K) {
  
  p = model$maxLag
  E = matrix(0, nrow(Y), ncol(Y))
  N = nrow(Y)
  M = ncol(Y)
  yp_lst = list()
  type = "free-run"
  ySlice = Y[1:(p - 1), ]
  uSlice = U[1:(p - 1)]
  eSlice = E
  
  pb = progress::progress_bar$new(total = N-p+1)
  for (k in p:N) {
    pb$tick()
    
    auxY = rbind(ySlice[(k - p + 1):(k - 1), ], rep(0, M))
    auxU = as.matrix(c(uSlice[(k - p + 1):(k - 1)], 0))
    auxE = rbind(eSlice[(k - p + 1):(k - 1), ], rep(0, M))
    phiK = genRegMatrix(model, auxY, auxU, auxE)$P
    preds = NULL
    m = 0
    
    for (j in 1:M) {
      theta = as.matrix(model$coefficients[[j]])
      new_pred = (phiK[, (j + m):(j + m + length(theta) - 1)] %*% theta)[1]
      preds = cbind(preds, new_pred)
      m = m + length(theta) - 1
    }
    ySlice = rbind(ySlice, preds)
    uSlice[k] = U[k]
    eSlice[k, ] = rep(0, ncol(Y))
  }
  
  # Create prediction-real dataframes for each coordinate
  for (j in 1:M) {
    df = data.frame(time = 1:N,
                    Y = Y[1:N, j],
                    U = U[1:N],
                    yh = ySlice[1:N, j],
                    e = Y[1:N, j] - ySlice[1:N, j])
    yp_lst[[j]] = df
  }
  
  return(yp_lst)
}

kStepAhead = function (model, Y, U, K) {
  
  p = model$maxLag
  E = matrix(0, nrow(Y), ncol(Y))
  N = nrow(Y)
  M = ncol(Y)
  type = paste0(K,"-steps ahead")
  
  pb = progress::progress_bar$new(total = N-p-K+2)
  
  time = (p+K-1):N
  yh = NULL
  
  for (k in p:(N-K+1)) {
    pb$tick()
    
    ySlice = Y[(k - p + 1):(k - 1), ]
    uSlice = U[(k - p + 1):(k - 1)]
    eSlice = E[(k - p + 1):(k - 1), ]
    
    for (nsteps in 1:K) {
      auxY = rbind(ySlice[(nsteps):(p+nsteps-2), ], rep(0, M))
      auxU = as.matrix(c(uSlice[(nsteps):(p+nsteps-2)], 0))
      auxE = rbind(eSlice[(nsteps):(p+nsteps-2), ], rep(0, M))
      phiK = genRegMatrix(model, auxY, auxU, auxE)$P
      
      uSlice[p + nsteps - 1] = U[p + nsteps - 1]
      ySlice_aux = rep(0, M)
      eSlice_aux = rep(0, M)
      m = 0
      
      for(j in 1:M) {
        theta = as.matrix(model$coefficients[[j]])
        ySlice_aux[j] = (phiK[, (j + m):(j + m + length(theta) - 1)] %*% theta)[1]
        eSlice_aux[j] = 0
        m = m + length(theta) - 1
      }
      ySlice = rbind(ySlice, ySlice_aux)
      eSlice = rbind(eSlice, eSlice_aux)
    }
    yh = rbind(yh, ySlice[p + K - 1, ])
  }
  yp_lst = list()
  
  for(j in 1:M) {
    df = data.frame(time = (p+K-1):N,
                    Y = Y[(p+K-1):N, j],
                    U = U[(p+K-1):N],
                    yh = yh[, j],
                    E = Y[(p+K-1):N, j] - yh[, j])
    yp_lst[[j]] = df
  }
  
  return(yp_lst)
}

oneStepAhead = function (model, Y, U,...) {
  
  p = model$maxLag
  E = matrix(0, nrow(Y), ncol(Y))
  N = nrow(Y)
  M = ncol(Y)
  m = 0
  yp_lst = list()
  type = "one-step-ahead"
  
  if (class(model) %in% c("armax","arx")) nl = 1 # linear with X
  else if (class(model) %in% "narmax")    nl = 2 # nonlinear with x
  else                                    nl = 3 # undefined
  
  for(j in 1:M) {
    theta = as.matrix(model$coefficients[[j]])
    
    # If e[k] does not exist on model, return the prediction
    if (!any(grepl('e(', model$terms[(j + m):(j + m + 1)], fixed = TRUE))) {
      P = genRegMatrix(model, Y, U, E)$P
      yp = P[, (j + m):(j + m + 1)] %*% theta
      yp_lst[[j]] = yp[,]
    }
    else {
      ySlice = Y[1:(p-1), j]
      eSlice = matrix(0, N, M)
      
      pb = progress::progress_bar$new(total = N-p+1)
      for (k in p:N) {
        pb$tick()
        
        auxY = rbind(Y[(k - p + 1):(k - 1), ], rep(0, M))
        auxU = as.matrix(c(U[(k - p + 1):(k - 1)], 0))
        auxE = rbind(eSlice[(k - p + 1):(k - 1), ], rep(0, M))
        phiK = genRegMatrix(model, auxY, auxU, auxE)$P
        ySlice[k] = (phiK[, (j + m):(j + m + 1)] %*% theta)[1]
        eSlice[k, j] = Y[k, j] - ySlice[k]
      }
    }
    m = m + 1
  }
  
  if (any(grepl('e(', model$terms, fixed = TRUE))) {
    for(j in 1:M) {
      df = data.frame(time = p:N,
                      Y = Y[p:N, j],
                      U = U[p:N],
                      yh = ySlice[p:N],
                      e = Y[p:N, j] - ySlice[p:N])
      
      yp_lst[[j]] = df
    }
  }
  
  return(yp_lst)
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
