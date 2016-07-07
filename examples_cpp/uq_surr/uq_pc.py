#!/usr/bin/env python


import os
import getopt
import shutil
import sys
import numpy as np
import math
import random as rnd
from scipy import stats, mgrid, reshape, random
import cPickle as pick

from model import model

#############################################################
def usage():
    # TODO more informative usage()
    # Usage string
    usage_str='uq_pc.py -r <run_regime> -p <pdomain_file> -m <method> -o <ord> -s <sam_method> -n <nqd> -v <nval>'
    def_str='uq_pc.py -r online -p pdomain_3d.dat -m proj -o 3 -s quad -n 7 -v 0'
    print "Correct syntax is"
    print usage_str
    print "Default values are"
    print def_str
    print "For more detailed information on parameters, refer to the UQTk Manual"
    

#############################################################

def model_pc(modelParam, pdom, pcparams):
    """PC surrogate evaluator"""
    #print "Running the surrogate model with parameters ", mparam
    pctype='LEG_N'

    np.savetxt('mindex.dat',pcparams[0],fmt='%d')
    np.savetxt('pccf.dat',pcparams[1])

    modelParam_scaled=( 2.*modelParam-(pdom[:,1]+pdom[:,0]) ) / (pdom[:,1]-pdom[:,0])
    np.savetxt('xdata.dat',modelParam_scaled)
    cmd="pce_eval -x'PC_mi' -f'pccf.dat' -s"+pctype+" -r'mindex.dat' > fev.log"
    os.system(cmd)
    pcoutput=np.loadtxt('ydata.dat')
   
    return pcoutput

#######################################################################################
#######################################################################################
#######################################################################################

def main(argv):
    
    # Defaults
    run_regime="online"  # Running regime, "online", "offline_prep" or "offline_post"
    pdomain_file="pdomain_3d.dat" # Parameter domain file
    method="proj"              # Surrogate construction method, "proj", "bcs" or "lsq"
    ord=3                      # Order of PC surrogate
    sam_method="quad"          # Parameter sampling method, "quad" or "unif"
    nqd=7                      # Number of quadrature points per dim (if "quad" method), or
                               # number of uniformly random samples (if "unif" method)
    nval=0                     # Number of validation samples, uniform random
    
    
    # Flags for input checks
    rflag=False
    pflag=False
    sflag=False
    nflag=False
    mflag=False
    oflag=False
    vflag=False
    
    
    # Hardwired names
    input_train='ptrain.dat'
    input_val='pval.dat'
    pctype='LEG_N'
    
    try:
        opts, args = getopt.getopt(argv,"hr:p:m:o:s:n:v:",["regime=","pdom=","method=","ord=","sampl=","nqd=","nval="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-r", "--regime"):
            run_regime = arg
            rflag=True
        elif opt in ("-p", "--pdom"):
            pdomain_file = arg
            pflag=True
        elif opt in ("-m", "--method"):
            method = arg
            mflag=True
        elif opt in ("-o", "--ord"):
            ord = int(arg)
            oflag=True
        elif opt in ("-s", "--sampl"):
            sam_method = arg
            sflag=True
        elif opt in ("-n", "--nqd"):
            nqd = int(arg)
            nflag=True
        elif opt in ("-v", "--nval"):
            nval = int(arg)
            vflag=True
      
    # Load parameter domain file
    if (os.path.isfile(pdomain_file)):
        pdom=np.loadtxt(pdomain_file).reshape(-1,2)
    else:
        print "Error: The requested domain file %s does not exist. Exiting." % pdomain_file
        sys.exit()

    # Get the dimensions
    dim=pdom.shape[0]
    two=pdom.shape[1]

    # Sanity checks
    assert(two==2)
    for i in range(dim):
        if(pdom[i,0]>pdom[i,1]):
            print "Error: The domain file %s contains wrong bounds. Check the row number %d. Exiting." % (pdomain_file,i+1)
            sys.exit()

    # Print the inputs for reference
    print "Run regime                        ", run_regime
    print "Input parameter domain file       ", pdomain_file
    print "The number of input parameters    ", dim
    print "Surrogate construction method     ", method
    print " with order                       ", ord
    print "Sampling method                   ", sam_method
    print " with parameter                   ", nqd
    print "Number of validation points       ", nval

    # Sanity check to ensure projection uses quadrature points
    # TODO in principle, we can implement MC integration for the projection integral 
    if (run_regime=="offline_post" or run_regime=="online") and method=="proj" and sam_method!="quad":
        print "Projection requires quadrature sampling. Exiting."
        sys.exit(1)

    # (1) Generate sample points for online or offline_prep regimes 
    if run_regime=="online" or run_regime=="offline_prep":
        
        if sam_method=="quad":
            cmd="generate_quad -d"+str(dim)+"  -g'LU' -x'full' -p"+str(nqd)+" -s"+pdomain_file+" > gq.log; mv xqdpts.dat "+input_train
            print "Running "+cmd
            os.system(cmd)

        elif sam_method=="unif":
            print "Error: uniform sampling method is not implemented yet. Exiting"
            # TODO implement this
            sys.exit()

        else:
            print "Error: Sampling method is not recognized. Should be 'quad' or 'unif'. Exiting."
            sys.exit()


        ptrain_unsc=np.loadtxt(input_train)
        npt=ptrain_unsc.shape[0]
        print "Parameter samples for surrogate construction are in %s in a format %d x %d " % (input_train,npt,dim)

        # Generate points, if requested, for the validation of the surrogate
        if nval>0:
            pval=2.*np.random.rand(nval,dim)-1.
            pval_unsc=0.5*(pdom[:,1]-pdom[:,0])*pval+0.5*(pdom[:,1]+pdom[:,0])
            np.savetxt(input_val,pval_unsc)
            print "Parameter samples for surrogate validation are in %s in a format %d x %d " % (input_val,nval,dim)

       
        # Cleanup
        cmd="rm -rf indices.dat"
        os.system(cmd)

        # Exit if only sample preparation is required
        if run_regime=="offline_prep":
            print "Preparation of samples is done."
            sys.exit()
                
    ############################################################################
                
    # (2) Load sample points for online or offline_post regimes
    ptrain_unsc=np.loadtxt(input_train).reshape(-1,dim)
    if (nval>0):
        pval_unsc=np.loadtxt(input_val).reshape(-1,dim)  

    npt=ptrain_unsc.shape[0]
    print "Number of training points for surrogate construction : "+str(npt)

    ############################################################################

    # (3) Get model outputs

    # Run the model online or....
    if run_regime=="online":
        ytrain=model(ptrain_unsc)
        if (nval>0):
            yval=model(pval_unsc)

    # ...or read the results from offline simulations
    elif run_regime=="offline_post":
        ytrain=np.atleast_2d(np.loadtxt('ytrain.dat'))
        if (nval>0):
            yval=np.atleast_2d(np.loadtxt("yval.dat"))

    # Read the number of output observables or the number of values of deisgn parameters (e.g. location, time etc..) 
    nout=ytrain.shape[1]
    print "Number of output observables of the model is ", nout

    ############################################################################

    # (4) Obtain the PC surrogate using model simulations

    # Emplty arrays and lists to store results
    pccf_all=[]
    mindex_all=[]
    allsens=np.empty((nout,dim))
    allsens_sc=np.empty((nout,dim))
    ytrain_pc=np.empty((npt,nout))
    yval_pc=np.empty((nval,nout))

    # Generate PC multiindex
    cmd="gen_mi -x'TO' -p"+str(ord)+" -q"+str(dim)+" > gmi.log; mv mindex.dat mi.dat"
    print "Running "+cmd
    os.system(cmd)

    # Loop over all output observables/locations
    for i in range(nout):
        
        ################################
        
        # (4a) Build PC surrogate
        print "##################################################"
        print "Building surrogate for observable %d / %d" % (i+1,nout)
        np.savetxt('ydata.dat',ytrain[:,i])
        ptrain=2.*(ptrain_unsc-pdom[:,0])/(pdom[:,1]-pdom[:,0])-1.
        np.savetxt('xdata.dat',ptrain)

        if method=="proj":
            cmd="pce_resp -x"+pctype+" -o"+str(ord)+" -d"+str(dim)+" -e > pcr.log; mv PCcoeff_quad.dat PCcoeff.dat"
            print "Running "+cmd
            os.system(cmd)

        elif method=="bcs":
            # TODO Ideally this could be done by lin_reg
            print "BCS is yet to be implemented. Exiting. "
            sys.exit()
            cmd="pce_bcs -x"+pctype+"-m'mi.dat' -t1 -e1e-5 > pcbcs.log; mv mindex.dat mi.dat"
            print "Running "+cmd
            os.system(cmd)

        elif method=="lsq":
            
            print "least-squares is yet to be implemented. Exiting. "
            sys.exit()
            ##cmd="pce_infer -f'data_const' -n"+str(npt)+" -o"+str(ord)+" -d"+str(dim)+" > pcinf.log"
            cmd="lin_reg -p"+pctype+" -t'lsq' -m'mi.dat'"
            print "Running "+cmd
            os.system(cmd)

        # Get the PC coefficients and multiindex
        pccf=np.loadtxt('PCcoeff.dat')
        mindex=np.loadtxt('mi.dat')

        # Append the results
        pccf_all.append(pccf)
        mindex_all.append(mindex)

        ################################

        # (4b) Evaluate the PC surrogate at training and validation points
        print "Evaluating surrogate at %d training points" % (npt)
        ytrain_pc[:,i]=model_pc(ptrain_unsc,pdom,[mindex,pccf])
        err_training=np.linalg.norm(ytrain[:,i]-ytrain_pc[:,i])/np.linalg.norm(ytrain[:,i])
        print "Surrogate relative error at training points : ", err_training


        if (nval>0):
            print "Evaluating surrogate at %d validation points" % (nval)
            yval_pc[:,i]=model_pc(pval_unsc,pdom,[mindex,pccf])
            err_val=np.linalg.norm(yval[:,i]-yval_pc[:,i])/np.linalg.norm(yval[:,i])
            print "Surrogate relative error at validation points : ", err_val
            #np.savetxt('yval_pc.'+str(i+1)+'.dat',yval_pc)
        
        ################################
        
        # (4c) Compute sensitivities
        cmd="pce_sens -m'mi.dat' -f'PCcoeff.dat' -x"+pctype+" > pcsens.log"
        print "Running "+cmd
        os.system(cmd)
        allsens[i,:]=np.loadtxt('mainsens.dat')
        allsens_sc[i,:]=allsens[i,:]/sum(allsens[i,:])
        print "Sum of main sensitivities :",sum(allsens[i,:])

    ############################################################################

    # Store for quick visualization
    np.savetxt('allsens.dat',allsens)
    np.savetxt('allsens_sc.dat',allsens_sc)

    # Results container
    if(nval>0):
        results = {'training':(pdom,ptrain_unsc,ytrain,ytrain_pc),'validation':(pdom,pval_unsc,yval,yval_pc),'pcmi':(pccf_all,mindex_all),'sens':(allsens,allsens_sc),'err':(err_training,err_val)}
    else:
        results = {'training':(pdom,ptrain_unsc,ytrain,ytrain_pc),'pcmi':(pccf_all,mindex_all),'sens':(allsens,allsens_sc),'err':(err_training)}

    # Save results
    pick.dump(results,open('results.pk','wb'),-1)


    # Cleanup of unneeded leftovers
    del_cmd='rm -rf ydata_pc.dat ydata.dat xdata.dat pccf.dat varfrac.dat totsens.dat mainsens.dat jointsens.dat sp_mindex*dat mindex.dat PCcoeff.dat xwghts.dat wghts.dat qdpts.dat'
    os.system(del_cmd)


# search for TODO
# Enrich -h for both model.py and uq_pc.py
# add verbosity flag.

# Nobody comes with their own points! We generate them! The former can be separate example! Maybe unif is too much?
# run without any flags? maybe put some restriction here to ensure some inputs are given

# pce_resp does not use custom mindex
# clean lin_reg Add uniform sampling
# Maybe store all joint sensitivities?
# Make multi_genz example
# Sparse quadrature needed?
# Make more flexible error catches
# Plotting tools rely on lib_dist plottools, need to move them to uqtk
# check 1D case separately, check online/offline case

if __name__ == "__main__":
   main(sys.argv[1:])
