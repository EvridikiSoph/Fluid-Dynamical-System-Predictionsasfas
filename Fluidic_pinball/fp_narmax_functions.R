##############
### Models ###
##############

# ARMAX Model
# Description: Creates an autoregressive moving average with exogenous inputs model
# ny: Number of autoregressive lags
# nu: Number of input lags
# ne: Number of moving average lags
# Output: Object representing an ARMAX mode
armax = function (ny, nu, ne) {
  model = list(
    ny = ny,
    nu = nu,
    ne = ne,
    maxLag = max(ny, nu, ne) + 1,
    terms = NULL,
    coefficients = NULL, # Variable name required by base::
    nb_terms = NULL,
    call = match.call()
  )
  class(model) = 'armax'
  return(model)
}

# NARMAX Model
# Description: Creates a nonlinear autogressive moving average with exogenous inputs model
# ny: Number of autoregressive lags
# nu: Number of input lags
# ne: Number of moving average lags
# nl: Nonlinearity polynomial length
# Output: Object representing a NARMAX model
narmax = function (ny, nu, ne, nl) {
  model = list(
    ny = ny,
    nu = nu,
    ne = ne,
    nl = nl,
    maxLag = max(ny, nu, ne) + 1,
    terms = NULL,
    coefficients = NULL,
    nb_terms = NULL,
    call = match.call()
  )
  class(model) = 'narmax'
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
# E: The error matrix (can be NULL)
# Output: Object containing:
#         P: Regression matrix with all terms
#         Pp: Regression matrix with only process terms
#         Pnp: Regression matrix without process terms
genRegMatrix.armax = function (model, Y, U, E) {
  na = model$ny
  nb = model$nu
  nc = model$ne
  p = model$maxLag
  N = nrow(Y)
  M = ncol(Y)
  
  # Check if sizes of Y, U and E are compatible
  sizeGuard(Y,U, E)
  obj = list()
  
  # If there is no error matrix present, set lags to 0
  if (is.null(E)) nc = 0
  
  # Initialize empty matrix that will contain potential termans of the model
  obj = list()
  #obj$P = matrix(0, nrow = N - p + 1, ncol = 2*M*na + nb + M*nc)
  obj$P = matrix(0, nrow = N - p + 1, ncol = M*na + nb + M*nc)
  
  colPhi = NULL     # Will contain column names 
  
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
  
  # Iterate through the polynomial order and add all combinations of all lags up to that order to the regression matrix
  #k = 2*M*na + nb
  k = M*na + nb
  
  if (nc > 0) {
    for(i in 1:nc) {
      #obj$P[, (M*na + nb + k + i):(M*na + nb + i + k + M - 1)] = E[(p - i):(N - i),]
      obj$P[, (k + i):(i + k + M - 1)] = E[(p - i):(N - i),]
      for(j in 1:M) {
        colPhi = c(colPhi, paste0("e(k-",i,")",j))
      }
      k = k + M - 1
    }
  }
  
  # Give each feature in the regression matirx its corresponding name
  rowPhi = paste0(rep("k=", N - p + 1), p:N)
  
  colnames(obj$P) = colPhi
  rownames(obj$P) = rowPhi
  
  # Create final complete regression matrix
  if (!is.null(model$terms)) {
    obj$P = subsetMatrix(obj$P, NULL, model$terms)
  }
  
  # Create regression matrix with everything but the error terms
  errIndexes = grepl('e(', colnames(obj$P), fixed = TRUE)
  
  obj$Pp = matrix(
    obj$P[, !errIndexes],
    ncol = sum(!errIndexes),
    dimnames = list(rownames(obj$P), colnames(obj$P)[!errIndexes]))
  
  # Create regression matrix with only the error terms
  obj$Pnp = matrix(
    obj$P[, errIndexes],
    ncol = sum(errIndexes),
    dimnames = list(rownames(obj$P), colnames(obj$P)[errIndexes]))
  
  return(obj)
}

# Description: Generates a regression matrix for an NARMAX Model
# model: NARMAX Model
# Y: The target matrix
# U: The input vector
# E: The error matrix (can be NULL)
# Output: Object containing:
#         P: Regression matrix with all terms
#         Pp: Regression matrix with only process terms
#         Pnp: Regression matrix without process terms
genRegMatrix.narmax = function (model, Y, U, E) {
  ny = model$ny
  nu = model$nu
  ne = model$ne
  nl = model$nl
  N = nrow(Y)
  M = ncol(Y)
  #n = 2*M*ny + nu + ne
  n = M*ny + nu + ne
  p = model$maxLag
  
  # Check if sizes of Y, U and E are compatible
  sizeGuard(Y, U, E)
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
  P0 = genRegMatrix(armax(ny, nu, ne), Y, U, E)$P
  P0[, 1:2*M*ny] = -P0[, 1:2*M*ny]
  
  # Assign column names
  col_names = NULL
  
  for(j in 1:M) {
    col_names = c(col_names, paste0('y(k-', 1:ny, ')',j))
  }
  
  #for(j in 1:M) {
  #  col_names = c(col_names, paste0('dy(k-', 1:ny[1], ')',j))
  #}
  
  col_names = c(col_names, paste0('u(k-', 1:nu, ')'))
  
  for(j in 1:M) {
    col_names = c(col_names, paste0('e(k-', 1:ny, ')',j))
  }
  
  # Add constant term to polynomial
  NP0 = nrow(P0)
  obj$P = NULL
  obj$P = cbind(rep(1, NP0), obj$P) # Contant
  colnames(obj$P) = "constant"
  obj$P = cbind(obj$P, P0)
  
  # Combine linear terms of armax feature matrix into nonlinear ones of order max nl
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
  
  # Create regression matrix with everything but the error terms
  errIndexes = grepl('e(', colnames(obj$P), fixed = TRUE)
  obj$Pp = matrix(
    obj$P[, !errIndexes],
    ncol = sum(!errIndexes),
    dimnames = list(rownames(obj$P), colnames(obj$P)[!errIndexes]))
  
  # Create regression matrix with only the error terms
  if (any(errIndexes)){
    obj$Pnp = matrix(
      obj$P[, errIndexes],
      ncol = sum(errIndexes),
      dimnames = list(rownames(obj$P), colnames(obj$P)[errIndexes]))
  } else{
    obj$Pnp = NULL
  }
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
    
    l1 = which(ERR == max(ERR))
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
    
    #th_FROLS = solve(A, gvec)
    th_FROLS = backsolve(A, gvec, upper.tri = TRUE)
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

# Estimate NARMAX model
# Description: Uses the FROLS algorithm to estimate the models coefficients
# model: NARMAX model
# Y: training time series
# U: training input
# rho_p: FROLS threshold fo the selection of the process terms only
# rho_e: FROLS threshold for the selection of all terms (including errors)
# Output: model with defined coefficients and model terms
estimate.narmax = function (model, Y, U, rho_p, rho_n) {
  
  # Generate regression matrix with error = 0
  M = ncol(Y)
  E = matrix(0, nrow(Y), ncol(Y))
  Pp = genRegMatrix(model, Y, U, E)$Pp
  Target = genTarget(model, Y)
  
  # Apply FROLS
  resultNarx = frols(Pp, Target, rho_p)
  
  # Calculate errors
  for(j in 1:M) {
    E[, j] = c(rep(0, model$maxLag - 1), as.vector(Target[, j]) - resultNarx$W[[j]] %*% resultNarx$g[[j]])
  }
  
  # Concatenate selected features for all univariate time series into a new regression matrix
  Psel_conc = NULL
  
  for(j in 1:M) {
    Psel_conc = as.data.frame(cbind(Psel_conc, resultNarx$Psel[[j]]))
  }
  
  Psel_conc = as.matrix(Psel_conc[!duplicated(colnames(Psel_conc))])
  
  # Also add the error terms
  Reg = cbind(Psel_conc, genRegMatrix(model, Y, U, E)$Pnp)
  
  # Re-apply FROLS
  resultNarmax = frols(Reg, Target, rho_n)
  
  # Add names of the selected terms as model terms
  terms = NULL
  nb_terms = NULL
  
  for(j in 1:M) {
    terms =  c(terms, colnames(resultNarmax$Psel[[j]]))
    nb_terms = c(nb_terms, length(colnames(resultNarmax$Psel[[j]])))
  }
  
  model$terms = terms
  model$nb_terms = nb_terms

  # Repeat a least squares process to calculate final coefficients
  theta_lst = list()
  k = 0
  
  for(j in 1:M) {
    nth = length(resultNarmax$th[[j]])
    E[, j] = c(rep(0, model$maxLag - 1), as.vector(Target[, j]) - resultNarmax$W[[j]] %*% resultNarmax$g[[j]])
    
    iterELS = 10
    thNarmaxHat = matrix(0, nrow = nth, ncol = iterELS)
    dlt = rep(0, iterELS)
    theta = NULL
    
    P = genRegMatrix(model, Y, U, E)$P[, (j + k):(j + k + ncol(resultNarmax$A[[j]]) - 1)]
    
    theta = MASS::ginv(P) %*% Target[, j]
    
    E1 = E[, j]
    E[, j] = c(rep(0, model$maxLag - 1), Target[, j] - (P %*% theta)[, ])
    
    for (s in 1:iterELS) {
      
      dlt[s] = sqrt(sum((E[, j] - E1) ^ 2))
      thNarmaxHat[, s] = theta
    }
    theta_lst[[j]] = theta
    k = k + ncol(resultNarmax$A[[j]]) - 1
  }
  
  model$coefficients = theta_lst
  
  return(model)
}


############################
### Prediction functions ###
############################

# Description: predictions using a NARMAX model, by choosing the prediction method depending on the value of K
# model: NARMAX model
# Y: time series to be predicted
# U: exogenous input
# K: number indicating the steps ahead to take
# Output: list containing predictions for each time series
predict.narmax = function (model, Y, U, K = 1, ...) {
  cat('Running narmax prediction ... \n')
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
