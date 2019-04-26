# Code for memory-microbe
# Author: Prashant Kalvapalle, Stadler lab/Rice University
# Date started: 1 feb 2019
# Description: Solving ODEs for integrase based flipping of DNA 

# initializing variables----
r0 <- 1e10 # copies of plasmid (unflipped) in 1 ml culture
f0 <- 0 # copies of flipped plasmid
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

# differential equations----
func_intrinsic_flip <- function(t,y, params) # differential equation for intrinsic flipping - integrase independent
{with(as.list(c(y,params)),
      {f <- y[1]; r <- y[2]
      df <- k1*r - k2*f
      dr <- -df
      dy <- list(c(df,dr))}
)}

func_integrase_flip <- function(t,y, params) # integrase dependant flipping differential eqn
{with(as.list(c(y,params)),
      {f <- y[1]; r <- y[2]
      df <- k1*r + kcat*i0*r^n/(km + r^n) - k2*f
      dr <- -df
      dy <- list(c(df,dr))}
)}

func_integrase_expression <- function(t,y, params) # integrase dependant flipping differential eqn
{with(as.list(c(y,params)),
      {IPTG <- y[1]; r <- y[2]
      dIntegrase <- k_leak + k_ind* IPTG -   k1*r + kcat*i0*r^n/(km + r^n) - k2*f
      dr <- -df
      dy <- list(c(df,dr))}
)}

# Solving ODE for intrinsic flipping
out_intrinsic_flip <- ode(c(f = f0, r = r0), times = t, func_intrinsic_flip, c(k1,k2))
out_integrase_flip <- ode(c(f = f0, r = r0), times = t, func_integrase_flip, c(k1,k2,kcat,km,n,i0))

# polishing ODE output data----
out_intrinsic_flip1 <- out_intrinsic_flip %>% as_tibble() %>% rename(flipped = f, unflipped = r) %>% gather('Orientation','# of plasmids', -time) %>% add_column(experiment = 'Intrinsic flip') # rename columns, gather into long array, add column to identify type
out_integrase_flip1 <- out_integrase_flip %>% as_tibble() %>% rename(flipped = f, unflipped = r) %>% gather('Orientation','# of plasmids', -time) %>% add_column(experiment = 'Baseline integrase activity')
timeseries <- full_join(out_intrinsic_flip1, out_integrase_flip1) # joining results from multiple ODES : integrase independant and dependant

# plotting # of plasmids with time----
# plots are faceted by experiment type : intrinsic flip vs basline integrase flip
plt <- timeseries %>% ggplot() + aes(time,`# of plasmids`, color = Orientation, facet = experiment) + geom_line() + geom_point() + 
       xlab('Time (hrs)') + ggtitle('Flipping with and without Integrase') + facet_wrap(~forcats::fct_inorder(experiment), scales = 'free_x') +  # fct_inorder ensures plotting in the order of joining of the experiments 
       theme_classic() + scale_color_brewer(palette="Set1") + theme(plot.title = element_text(hjust = 0.5))  # Setting theme, colour scale and need to centre the plot title 
print(plt)