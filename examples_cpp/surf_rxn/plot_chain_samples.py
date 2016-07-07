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
import sys
import shutil
import fileinput
import numpy as npy
import matplotlib.pyplot as plt
from scipy import stats, mgrid, c_, reshape, random, rot90
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from pylab import *
import file_utils

rc('legend',loc='upper left', fontsize=12)
rc('lines', linewidth=4, color='r')
rc('axes',linewidth=3,grid=True,labelsize=22)
rc('xtick',labelsize=20)
rc('ytick',labelsize=20)


chainfile="chain.dat"; # chain file
nskip=5000;            # skip first 'nskip' states
nthin=10;            # pick every 'nthin' state

all_samples, vnames = file_utils.extract_all_vars(chainfile,nskip,0,1)
n_all_vars = len(vnames)
n_cols = len(all_samples[0,:])
chain = all_samples[:,0:1+n_all_vars]


for i in range(n_all_vars):
    fig = plt.figure(figsize=(10,7))
    ax=fig.add_axes([0.10,0.15,0.85,0.75])

    plt.plot(chain[:,0],chain[:,i+1],color='black',linewidth=2)
    ax.set_xlabel("MCMC step",fontsize=22)
    ax.set_ylabel(vnames[i],fontsize=22)

    plt.savefig('chain_'+vnames[i]+'.eps')
    plt.clf()

for i in range(n_all_vars):
    for j in range(i):
        fig = plt.figure(figsize=(10,7))
        ax=fig.add_axes([0.10,0.15,0.85,0.75])
    
        plt.plot(chain[:,j+1],chain[:,i+1],'ko',markeredgecolor='black',markersize=5)
        ax.set_xlabel(vnames[j],fontsize=22)
        ax.set_ylabel(vnames[i],fontsize=22)

        plt.savefig('chain_'+vnames[j]+'_'+vnames[i]+'.eps')
        plt.clf()
