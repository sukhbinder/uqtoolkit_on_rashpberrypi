#!/bin/bash -e

export PATH=../../../bin/:$PATH

######################################################

PCTYPE="LU"
DIM=2
ORDER=3

# Generate quadrature points
generate_quad -d$DIM -g$PCTYPE -x'full' -p7
# Evaluate a function
awk '{print $1*cos($2+0.1)+$1*$1}' qdpts.dat > ydata.dat
# Compute the PC response surface 
pce_resp -e -d$DIM -x$PCTYPE -o$ORDER 
