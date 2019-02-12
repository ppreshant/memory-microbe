# Code for memory-microbe
# Author: Prashant Kalvapalle
# Date started: 1 feb 2019

# initializing variables
r0 <- 1e10 # copies of plasmid (unflipped) in 1 ml culture
f0 <- 0
k1 <- 10 # hr-1
k2 <- 50 # hr-1

i0 <- 1e6 # number of integrase proteins
kcat <- 1e5 # hr-1
km <- 1e5 # copies
n <- 4 # effective hill coefficient - Integrase tetramerization and binding

# other initializations
t <- seq(0,.05,.001) # time sequence

# calling libraries
library(deSolve)
library(tidyverse)
library(ggplot2)
library(magrittr)

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
out_int <- ode(c(f = f0, r = r0), times = t, func_int, c(k1,k2,kcat,km,n,i0))

# polishing ODE output data
out_i1 <- out_i %>% as_tibble() %>% rename(flipped = f, unflipped = r) %>% gather('Orientation','# of plasmids', -time)
out_int1 <- out_int %>% as_tibble() %>% rename(flipped = f, unflipped = r) %>% gather('Orientation','# of plasmids', -time)
  
# plotting # of plasmids with time
plt <- out_i1 %>% ggplot() + aes(time,`# of plasmids`, color = Orientation) + geom_line() + geom_point() + 
       xlab('Time (hrs)') + ggtitle('Integrase independant flipping') +
       theme_classic() + scale_color_brewer(palette="Set1") + theme(plot.title = element_text(hjust = 0.5))  
# print(plt)

plt_int <- out_int1 %>% ggplot() + aes(time,`# of plasmids`, color = Orientation) + geom_line() + geom_point() + 
  xlab('Time (hrs)') + ggtitle('Baseline Integrase flipping') +
  theme_classic() + scale_color_brewer(palette="Set1") + theme(plot.title = element_text(hjust = 0.5))  

print(plt_int)