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

# calling libraries
library(deSolve)

# other initializations
t <- seq(0,10,.01)

# differential equations
func1 <- function(t,y, params)
{with(as.list(c(y,params)),
      {f <- y[1]; r <- y[2]
      df <- k1*r - k2*f
      dr <- -df
      dy <- list(df,dr)}
)}
  

# Solving ODE for intrinsic flipping
out1 <- ode(c(f = f0, r = r0), times = t, func1, c(k1,k2))
View(out1)