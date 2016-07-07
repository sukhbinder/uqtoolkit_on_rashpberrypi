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


from pyUQTk.plotting.surrogate import *

from pylab import *

import cPickle as pick

rc('legend',loc='best', fontsize=22)
rc('lines', linewidth=2, color='r')
rc('axes',linewidth=3,grid=True,labelsize=28)
rc('xtick',labelsize=20)
rc('ytick',labelsize=20)

    
#############################################################
#############################################################

plotid=(sys.argv[1])


results=pick.load(open('results.pk', 'rb'))

if(plotid=='sens'):
    sensdata=results['sens'][0]
    
    # TODO fix this via an argument
    #names_file='params'
    #text_file = open(names_file, "r")
    #par_labels = text_file.readlines()
    #print par_labels
    #text_file.close()
        
    #names_file='cases'
    #text_file = open(names_file, "r")
    #case_labels = text_file.readlines()
    #print case_labels
    #text_file.close()
    
    npar=sensdata.shape[1]
    ncases=sensdata.shape[0]
    pars=range(npar) #[0,1,2,5]
    cases=range(ncases)
    
    plot_sens(sensdata,pars,cases,vis="bar",reverse=False,par_labels=False,case_labels=False,showplot=True)

elif(plotid=='dm'):
    # TODO fix this via an argument
    icase=3
    val=False
    
    datas=[results['training'][2][:,icase]]
    models=[results['training'][3][:,icase]]
    labels=['training points']
    axes_labels=['Model','Polynomial Surrogate']

    if (val):
        datas.append(results['validation'][2][:,icase])
        models.append(results['validation'][3][:,icase])
        labels.append('validation points')
    
    plot_dm(datas,models,labels,axes_labels,showplot=True)

elif(plotid=='idm'):
    # TODO fix this via an argument
    icase=3
    trval='training' #'validation'
    
    data=results[trval][2][:,icase]
    model=results[trval][3][:,icase]
    label=trval+' points'
    
    plot_idm(data,model,sort='data',showplot=True)


