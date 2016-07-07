#!/usr/bin/env python
#=====================================================================================
#                     The UQ Toolkit (UQTk) version 2.1.1
#                    Copyright (2013) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#    Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
#    with Sandia Corporation, the U.S. Government retains certain rights in this software.
#
#    This file is part of The UQ Toolkit (UQTk)
#
#    UQTk is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    UQTk is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
#
#    Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
#    Sandia National Laboratories, Livermore, CA, USA
#====================================================================================


import os
import shutil
import sys
import math
import getopt

try:
    import numpy as np
except ImportError:
    print "\nline_infer.py requires numpy package -> Exit\n"
    quit()

import random as rnd
from scipy.stats.mstats import mquantiles
import matplotlib.pyplot as plt

import fileinput

try:
    import pyUQTk.inference.postproc as postp
except ImportError:
    print "\nline_infer.py requires the pyUQTk package. Make sure the UQTk root",
    print "directory is included in the PYTHONPATH environment variable.\n"
    quit()

from pylab import *

rc('legend',loc='upper right', fontsize=20)
rc('lines', linewidth=4, color='r')
rc('axes',linewidth=3,grid=True,labelsize=22)
rc('xtick',labelsize=16)
rc('ytick',labelsize=16)

###########################################################################
# define uqtkbin
if os.environ.get("UQTK_SRC") is None:
    print "Error: Need to set path to uqtk src as environment variable UQTK_SRC -> Abort"
    quit()
else:
    if ( not os.path.isdir(os.environ["UQTK_SRC"]) ):
        print "\"",os.environ["UQTK_SRC"],"\" is not a valid path -> Abort"
        quit()

uqtkbin=os.environ["UQTK_SRC"]+"/src_cpp/bin"
pcequad=uqtkbin+"/pce_quad"
###########################################################################

help_string = """
Usage:
  line_infer.py [-h] [--nd <number of points>] [--noplots] [--stats]
what
  Run a Bayesian line inference example; generate plots and statistics of the output
where
  -h    = print help info
  --nd  = number of data points to generate [defaults to 5]
  --noplots = do not generate chain and posterior plots [defaults to False]
  --stats = generate MCMC chain statistics [defaults to False]
""" 

#
# Process inputs 
#
try:
    # opts,v_names = getopt.getopt(sys.argv[1:],"hi:s:",["nb="])
    opts,extra_arguments = getopt.getopt(sys.argv[1:],"h",["nd=","noplots","stats"])
except getopt.GetoptError, err:
    print str(err)
    print help_string
    sys.exit(1)

# Default values
npt = 5    # Number of data points to infer line from
generate_plots = True # Generate chain and posterior pdf plots
generate_stats = False # By default, do not generate MCMC statistics

# samples_file_name=""
# n_burnin = 0
# stride = 1

print opts

for o,a in opts:
    if o == "-h":
        print help_string
        sys.exit(0)
    elif o == "--nd":
        npt = int(a)
    elif o == "--noplots":
        generate_plots = False
    elif o == "--stats":
        generate_stats = True
    else:
        assert False, "unhandled command line parsing option"

# error checking
if (npt < 2):
    assert False, "The number of data points needs to be >= 0"


# Other input settings 
model='lin'             # Type of the model: 'cos' or 'lin'
stdev= 0.1              # Added noise magnitude

#defparam=[-2, 4]        # Default parameters, used if only one is inferred, or to generate artificial data
defparam=[5, -1] 
# Do not change number of inferred parameters for this test problem, or it will break postprocessing
chdimlist=[1,2]               # Indices of chain parameters
chstart=np.array([2, -2])  # Starting value of the chain
chainfile='line_infer.chain.dat'      # Name of the chain file
input_data_file = "line_infer.inputdata.dat"    # Name of file with input data
n_burnin=10000              # Burnin samples
stride=5                   # Stride with which to read in MCMC samples file (reading only every stride samples)

pctype='HG'                # PC type
pcord=3                    # PC order of the chain samples 
bw=-1                    # Bandwidth for Rosenblatt (if negative, the program picks the optimal value) 

debug = 0               # Set to 1 or higher to get more output
dense_plots = True      # Set to True to get dense format for posterior plots in triangular format
## Prepare the xml file ##################################################
shutil.copyfile('line_infer.xml.templ','line_infer.xml')
for line in fileinput.input('line_infer.xml', inplace = 1):
    print line.replace('PCTYPE', pctype),
for line in fileinput.input('line_infer.xml', inplace = 1):
    print line.replace('PCORDER', str(pcord)),
for line in fileinput.input('line_infer.xml', inplace = 1):
    print line.replace('CHAINFILE', chainfile),
for line in fileinput.input('line_infer.xml', inplace = 1):
    print line.replace('INPUTDATAFILE', input_data_file),
for line in fileinput.input('line_infer.xml', inplace = 1):
    print line.replace('NOISE', str(stdev)),
    

## Synthetic data generation ###########################################
print "Generating synthetic data"
xin=np.random.uniform(-1.0,1.0, npt)
yin=np.zeros(npt)
for i in range(npt):
    if model=='cos':
        yin[i]=math.cos(xin[i]+math.exp(1))+np.random.normal(0,stdev)
    elif model=='lin':
        yin[i]=defparam[0]+defparam[1]*xin[i]+np.random.normal(0,stdev)
    else:
        print "Model unrecognized"
        sys.exit()

np.savetxt(input_data_file,np.transpose([xin,yin]))

## Run the inference code #####################################################
ii=0
for idim in chdimlist:
    for line in fileinput.input("line_infer.xml", inplace = 1):
        print line.replace('PAR_'+str(idim), str(chstart[ii])),
    ii=ii+1

print "Running the parameter inference"
os.system('./line_infer.x')


## Import data from MCMC file ###################################################
print "Loading in chain file",chainfile
all_samples, vnames = postp.extract_all_vars(chainfile,n_burnin,0,stride)
n_all_vars = len(vnames)

# Extract all MCMC chain variables in separate array
chn  = all_samples[:,0:1+n_all_vars]
nchn = chn.shape[0]

if generate_plots:
    ## Find posterior predictive #######################################################
    # ideally we should have a deterministic forward model code,
    # but for a line it is simple enough to compute here
    print "\nSampling chain to compute posterior predictive distribution\n"

    # First set parameters to the default parameters; then replace the inferred parameters
    # with samples from the chain.
    param = defparam

    # Set up x-locations to sample posterior predictive at and initialize arrays to keep
    # track of moments at those locations
    ngr=100 # number of locations in x
    xlin=np.linspace(-1.,1.,ngr)
    ypp=np.zeros((ngr))
    y2psh=np.zeros((ngr))

    npp = min(1000,nchn) # max number of posterior predictive samples wanted, regardless of chain length
    i_pp = 0 # counter of number of posterior predictive samples

    for ip in range(0,nchn,nchn/npp):
        ii=0
        for idim in chdimlist:
            # grab variable samples from posterior chain
            # (The first inferred parameter is in the second column of the matrix chn with samples)
            param[idim-1]=chn[ip,ii+1]
            ii=ii+1

        # run forward model for these sampled parameter values at the chosen x-locations
        for j in range(ngr):
            yval=param[0]+param[1]*xlin[j]
            # running sum for y
            ypp[j]=ypp[j]+yval
            # running sum for y^2
            y2psh[j]=y2psh[j]+yval**2

        i_pp += 1


    # pushed forward mean
    ypp=ypp/float(i_pp)
    # pushed forward ave y^2
    y2psh=y2psh/float(i_pp)
    # pushed forward std dev
    std_push=(y2psh-ypp**2)**0.5
    # Posterior predictive std dev (pushed forward variance + data noise variance)
    std_pp=(std_push**2 + stdev**2)**0.5


    # Plot pushed forward posterior and posterior predictive
    fig = plt.figure(figsize=(10,7))
    ax=fig.add_axes([0.10,0.15,0.85,0.75])
    plt.fill_between(xlin,ypp-std_pp,ypp+std_pp,color='lightgrey',label='Post predictive stdev')
    plt.fill_between(xlin,ypp-std_push,ypp+std_push,color='grey',label='Pushed forward stdev')
    plt.plot(xlin, ypp, linewidth=1, label='Mean prediction')
    plt.plot(xin,yin,'o', markersize=8, color='black',label='Data')
    ax.set_xlabel("x",fontsize=22)
    ax.set_ylabel("y",fontsize=22)
    plt.legend(loc='best')

    plt.savefig('line_infer.postpred.pdf')
    plt.clf()


    ## Plot chains ##################################################################
    print "\nPlotting chains ... \n"
    for i in range(n_all_vars):
        fig = plt.figure(figsize=(10,7))
        ax=fig.add_axes([0.10,0.15,0.85,0.75])

        plt.plot(chn[:,0],chn[:,i+1],color='black',linewidth=2)
        ax.set_xlabel("MCMC step",fontsize=22)
        ax.set_ylabel(vnames[i],fontsize=22)

        plt.savefig('line_infer.chn_'+vnames[i]+'.pdf')
        plt.clf()

    for i in range(n_all_vars):
        for j in range(i):
            fig = plt.figure(figsize=(10,7))
            ax=fig.add_axes([0.10,0.15,0.85,0.75])
        
            plt.plot(chn[:,j+1],chn[:,i+1],'ko',markeredgecolor='black',markersize=5)
            ax.set_xlabel(vnames[j],fontsize=22)
            ax.set_ylabel(vnames[i],fontsize=22)

            plt.savefig('line_infer.chn_'+vnames[j]+'_'+vnames[i]+'.pdf')
            plt.clf()


    ## Plot posterior PDF 'triangle' ################################################
    np_kde = 100
    postp.plot_all_posteriors(chn[:,1:],vnames,np_kde,"line_infer.posteriors",debug,dense_plots)

if generate_stats:
    postp.get_mcmc_stats(all_samples,vnames,"line_infer",debug)

# ## Find PC coefficients corresponding to the chain samples #######################
# samfile="chain_samples.dat"
# # save chain samples minus the step column 
# np.savetxt(samfile,chn[:,1:])
# print "Running KDE-Rosenblatt transformation to build input PCE"
# cmd = pcequad+' -o' + str(pcord)+' -f ' + samfile + ' -x ' + pctype + ' -w' + str(bw) + ' > pcequad.log'
# print "Running",cmd
# os.system(cmd)

# for res_file in [chainfile, samfile, 'quadpts.dat', 'quadpts_mapped.dat', 'PCcoeff.dat']:
#     shutil.copyfile(res_file, 'prob5_'+res_file)

# # Plots for demontration of Rosenblatt transformation
# fig = plt.figure(figsize=(8,6))
# ax=fig.add_axes([0.12,0.10,0.83,0.85])
# qdpts=np.loadtxt("quadpts.dat")
# qdpts_mapped=np.loadtxt("quadpts_mapped.dat")

# plt.plot(chn[:,1],chn[:,2],'*', color='blue', markersize=3)
# plt.plot(qdpts_mapped[:,0],qdpts_mapped[:,1],'o',color='red', markersize=9)
# ax.set_xlabel("param_a",fontsize=22)
# ax.set_ylabel("param_b",fontsize=22)

# plt.savefig('prob5_rosen.pdf')
# plt.clf()

# # Forward uncertainty propagation of PC intrusively through the line model a+bx
# print "Forward uncertainty propagation of PC"
# pcin=np.loadtxt("PCcoeff.dat")
# npc=pcin.shape[0]
# pcout=np.empty((ngr,npc))
# for j in range(ngr):
#     pcout[j,:] = pcin[:,0] + pcin[:,1] * xlin[j]

# fig = plt.figure(figsize=(10,7))
# ax=fig.add_axes([0.10,0.15,0.85,0.75])
# for i in range(npc):
#     plt.plot(xlin, pcout[:,i], linewidth=1, label='c_'+str(i))

# ax.set_xlabel("x",fontsize=22)
# ax.set_ylabel("PC coefficients",fontsize=22)
# plt.legend()

# plt.savefig('prob5_pccf_out.pdf')
# plt.clf()

print "END: Line Inference example done: check out line_infer.*.pdf files for results."
