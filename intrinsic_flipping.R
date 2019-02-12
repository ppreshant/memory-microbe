# Code for memory-microbe
# Author: Prashant Kalvapalle
# Date started: 1 feb 2019

# initializing variables
r0 <- 1e10 # copies of plasmid (unflipped) in 1 ml culture
f0 <- 0
k1 <- 100 # sec-1
k2 <- 200 # sec-1

i0 <- 1e3 # number of integrase proteins
kcat <- 1 # sec-1
km <- 1e5 # copies
n <- 4 # effective hill coefficient - Integrase tetramerization and binding

# other initializations
t <- seq(0,.05,.0001) # time sequence

# calling libraries
library(deSolve)
# library(magrittr)
# library(ggplot2)

# differential equations
func_i <- function(t,y, params) # differential equation for intrinsic flipping - integrase independant
{with(as.list(c(y,params)),
      {f <- y[1]; r <- y[2]
      df <- k1*r - k2*f
      dr <- -df
      dy <- list(c(df,dr))}
)}

func_int <- function(t,y, params) # integrase dependant flipping differential eqn
{with(as.list(c(y,params)),
      {f <- y[1]; r <- y[2]
      df <- k1*r + kcat*i0*r^n/(km + r^n) - k2*f
      dr <- -df
      dy <- list(c(df,dr))}
)}

# Solving ODE for intrinsic flipping
out_i <- ode(c(f = f0, r = r0), times = t, func_i, c(k1,k2))
# View(out_i)
plot(out_i[,1], out_i[,2])