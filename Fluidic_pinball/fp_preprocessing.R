library(ggplot2)
library(tidyr)
library("plot3D")
library(boot)
library(PLRModels)
library(sindyr)
library(tseries)
library(lokern)
library(purrr)
library(astsa)


####################
### Read in data ###
####################
Cl = read.table('C:/Users/evsop/OneDrive/Documents/Master/4th Term/TFM/dynamical_fluid_forecast_TFM/Fluidic Pinball/Input/ClValues.txt', header = FALSE)
lift = apply(Cl[, 2:4], 1, sum)

Cd = read.table('C:/Users/evsop/OneDrive/Documents/Master/4th Term/TFM/dynamical_fluid_forecast_TFM/Fluidic Pinball/Input/CdValues.txt', header = FALSE)
drag = apply(Cd[, 2:4], 1, sum)

ClCd = data.frame(Cl = lift, Cd = drag)

# Drop the first 101 lines
ClCd = ClCd[401:nrow(ClCd), ]

# Generate random initial points between 1000
nb_trials = 500              # number of different trials/initial prediction points
nb_series = 2                # number of univariate time series
dt = 0.01
initial_points = sort(sample(1:1000, nb_trials), decreasing = FALSE)
setwd('C:/Users/evsop/OneDrive/Documents/Master/4th Term/TFM/dynamical_fluid_forecast_TFM/Fluidic Pinball')
col_names = c("Cl", "Cd")

# Vector of rotation values in training data
u_list = c(1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 
           1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 
           1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 
           1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39)

u = NULL

for (i in u_list) {
  u_temp = rep(i, 300)
  u = c(u, u_temp)
}

ClCd_final = data.frame(ClCd, u)
col_names = c(col_names, "u")
colnames(ClCd_final) = col_names

xyzu_train = ClCd_final[1:10400, ]
xyzu_pred = ClCd_final[10401:nrow(ClCd_final), ]


#############
### Plots ###
#############
# Plot drag and lift coefficients
ts_xyz = ts(ClCd_final[1:1000, ], start = c(0,0.001),frequency = 1)
par(cex.axis=2,cex.lab=2.2, cex.main=3)
plot.ts(ts_xyz, main="")
