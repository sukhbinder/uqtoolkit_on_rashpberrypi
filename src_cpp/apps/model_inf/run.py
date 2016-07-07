#!/usr/bin/env python

import os
import shutil
import sys
import numpy as np
import math
from scipy import stats, mgrid, reshape, random
from scipy.stats.mstats import mquantiles
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab


sys.path.append(os.environ['HOME']+"/research/lib_dist/pytools")
from uqtools import *

from pylab import *
import cPickle as pick

rc('legend',loc='upper right', fontsize=32)
rc('lines', linewidth=4, color='r')
rc('axes',linewidth=3,grid=True,labelsize=32)
rc('xtick',labelsize=30)
rc('ytick',labelsize=30)


#############################################################
#############################################################
def model_example(modelPar):
    
    npar=modelPar.shape[0]
    mdim=modelPar.shape[1]
    
    nout=7
    
    # Create the design parameters
    designPar=np.array(range(nout)).reshape(-1,1)/float(nout-1)
    np.savetxt('designPar.dat',designPar)
    
    output=np.empty((npar,nout))
    for i in range(npar):
        print i+1, ": Running the model with parameter setting ", modelPar[i,:]
        for j in range(nout):
            aa=np.dot(modelPar[i,:]+modelPar[i,:]**2,pow(np.arange(1,mdim+1),-designPar[j]-1.))
            output[i,j]=aa *(sum(modelPar[i,:]))
    
    return output

#############################################################
nout=11
designPar=np.array(range(nout)).reshape(-1,1)/float(nout-1)
np.savetxt('xdata.dat',designPar)
data=designPar**3
np.savetxt('ydata.dat',data.T)


ppccf=np.zeros((2,nout))
for i in range(nout):
    ppccf[1,i]=designPar[i]

np.savetxt('ppccf.dat',ppccf)


cmd='../model_inf'
os.system(cmd)

mapparams=np.loadtxt('MAPparams.dat')

plot(designPar,data,'b*')
ngr=100
xgrid=np.array(range(ngr)).reshape(-1,1)/float(ngr-1)

plot(xgrid,mapparams*xgrid,'r-')
show()





