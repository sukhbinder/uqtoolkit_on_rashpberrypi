#!/usr/bin/env python


import os
import getopt
import shutil
import sys
import numpy as np
import math


#############################################################
#############################################################

def model(modelPar):
    return model_example(modelPar)

#############################################################
     
def model_example(modelPar):

    npar=modelPar.shape[0]
    mdim=modelPar.shape[1]
    
    nout=7
    
    # Create the design parameters #TODO maybe should be done outside
    designPar=np.array(range(nout)).reshape(-1,1) #/float(nout-1)
    np.savetxt('designPar.dat',designPar)
    
    output=np.empty((npar,nout))
    for i in range(npar):
        print i+1, ": Running the model with parameter setting ", modelPar[i,:]
        for j in range(nout):
            aa=np.dot(modelPar[i,:]+modelPar[i,:]**2,pow(np.arange(1,mdim+1),-designPar[j]-1.))
            output[i,j]=aa *(sum(modelPar[i,:]))

    return output

#############################################################

def main(argv):
    modelPar_file=argv[0]
    output_file=argv[1]
    #designPar_file=argv[2]
    modelPar=np.loadtxt(modelPar_file) # TODO what if 2d?

    output=model(modelPar)
    np.savetxt(output_file,output)

if __name__ == "__main__":
    main(sys.argv[1:])






