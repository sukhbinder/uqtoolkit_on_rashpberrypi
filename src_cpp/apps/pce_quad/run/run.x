#!/bin/bash -e

export PATH=../../../bin/:$PATH

######################################################

NSAM=100000
PCTYPE="HG"
SDIM=1
ORDER=3

# Generate random variables with a given PC
pce_rv -w'PC' -d1 -n$NSAM -x$PCTYPE -p$SDIM -o$ORDER -f'pccf.dat'

# Find a PC coefficients coresponding to the random variable given samples above
pce_quad -f'rvar.dat' -w-1 -x$PCTYPE -o$ORDER

