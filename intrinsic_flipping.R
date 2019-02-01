# Code for memory-microbe
# Author: Prashant Kalvapalle
# Date started: 1 feb 2019

r = 1e10 # copies of plasmid (unflipped) in 1 ml culture
f = 0 
k1 = 100 # sec-1
k2 = 200 # sec-1

i0 = 1e3 # number of integrase proteins
kcat = 1 # sec-1
km = 1e5 # copies

# calling libraries
library(deSolve)

