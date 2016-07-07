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
import numpy as np
import math
from scipy import stats, mgrid, reshape, random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import pyUQTk.utils as ut
from pylab import *

rc('legend',loc='upper left', fontsize=12)
rc('lines', linewidth=1, color='r')
rc('axes',linewidth=3,grid=True,labelsize=22)
rc('xtick',labelsize=20)
rc('ytick',labelsize=20)

#############################################################
#############################################################
def plot_micf(mindex,cfs=[],showplot=True):
    """Plots 2d or 3d multiindices"""
    
    custom_xlabel='Dim 1'
    custom_ylabel='Dim 2'
    custom_zlabel='Dim 3'

    npc=mindex.shape[0]
    ndim=mindex.shape[1]

    if (ndim==2):
        if cfs==[]:
            #plot(mindex[:,0],mindex[:,1],'bo',markersize=13)
            scatter(mindex[:,0],mindex[:,1],s=150, marker = 'o',cmap = cm.jet )
        else:
            scatter(mindex[:,0],mindex[:,1],s=150, c=cfs, marker = 'o',cmap = cm.jet )
    elif (ndim==3):
        ax = figure().add_subplot(111, projection='3d')
        if cfs==[]:
            ax.scatter(mindex[:,0],mindex[:,1],mindex[:,2],s=50)
        else:
            ax.scatter(mindex[:,0],mindex[:,1],mindex[:,2],c=cfs,s=50)

        ax.set_zlabel(custom_zlabel)
        ax.set_zlim([-.5, max(mindex[:,2])+.5])

    else:
        raise NameError("Multi-index should be 2d or 3d")

    xlabel(custom_xlabel)
    ylabel(custom_ylabel)
    xlim([-.5, max(mindex[:,0])+.5])
    ylim([-.5, max(mindex[:,1])+.5])

    if (showplot):
        show()
    else:
        savefig('micf.eps')
        clf()


#############################################################
def plot_idm(data,model,sort='none',showplot=True):
    """Plots data and model on the same axis"""
    
    axes_labels=['Run Id','Model / Surrogate']

    custom_xlabel=axes_labels[0]
    custom_ylabel=axes_labels[1]
    
    figure(figsize=(12,8))
    
    npts=data.shape[0]
    neach=1
    if (data.ndim>1):
        neach=data.shape[1]
    
    #neb=model.shape[1]-1# errbars not implemented yet
    
    if (sort=='model'):
        ind=argsort(model)
    elif (sort=='data'):
        ind=argsort(data)
    elif (sort=='none'):
        ind=range(npts)


    ddata=data.reshape(npts,neach)
    
    
    plot(range(1,npts+1),model[ind], 'bo', label='Surrogate')
    for j in range(neach):
        plot(range(1,npts+1),ddata[ind,j], 'ro',label='Model')
    
    
    xlabel(custom_xlabel)
    ylabel(custom_ylabel)
    #title('Data vs Model')
    legend()
    
    #xscale('log')
    #yscale('log')
    
    if (showplot):
        show()
    else:
        savefig('idm.eps')
        clf()


#############################################################
def plot_dm(datas,models,labels,axes_labels,showplot=True):
    """Plots data-vs-model and overays y=x"""
    
    custom_xlabel=axes_labels[0]
    custom_ylabel=axes_labels[1]

    figure(figsize=(10,10))
    ncase=len(datas)
    # Create colors list
    colors=ut.set_colors(ncase)
    yy=np.empty((0,1))
    for i in range(ncase):
        data=datas[i]
        model=models[i]
        
        npts=data.shape[0]
        neach=1
        if (data.ndim>1):
            neach=data.shape[1]

        #neb=model.shape[1]-1# errbars not implemented yet

        
    
        ddata=data.reshape(npts,neach)


        for j in range(neach):
            yy=np.append(yy,ddata[:,j])
            plot(ddata[:,j],model, 'o',color=colors[i],label=labels[i])

    delt=0.1*(yy.max()-yy.min())
    minmax=[yy.min()-delt, yy.max()+delt]
    plot(minmax,minmax,'r-',label='y=x')

    xlabel(custom_xlabel)
    ylabel(custom_ylabel)
    #title('Data vs Model')
    legend()

    #xscale('log')
    #yscale('log')

    if (showplot):
        show()
    else:
        savefig('dm.eps')
        clf()

#############################################################

def plot_sens(sensdata,pars,cases,vis="bar",reverse=False,par_labels=False,case_labels=False,showplot=True):
    """Plots sensitivity for multiple observables"""
    
    ncases=sensdata.shape[0]
    npar=sensdata.shape[1]
    
    wd=0.6
    xlbl=''
    ylbl='Sensitivity'
    legend_show=True #False
    assert set(pars) <= set(range(npar))
    assert set(cases) <= set(range(ncases))
    
    # Set up the figure
    # TODO need to scale figure size according to the expected amount of legends
    fig = plt.figure(figsize=(16,14))
    fig.add_axes([0.1,0.3,0.8,0.65])
    
    #########
    
    # Default parameter names
    if (~par_labels):
        par_labels = []
        for i in range(npar):
            par_labels.append(('par_'+str(i+1)))
    # Default case names
    if (~case_labels):
        case_labels=[]
        for i in range(ncases):
            case_labels.append(('case_'+str(i+1)))


    if(reverse):
        tmp=par_labels
        par_labels=case_labels
        case_labels=tmp
        tmp=pars
        pars=cases
        cases=tmp
        sensdata=sensdata.transpose()
    ##############################################################################

    npar_=len(pars)
    ncases_=len(cases)

    # Create colors list
    colors=ut.set_colors(npar_)

    case_labels_=[]
    for i in range(ncases_):
        case_labels_.append(case_labels[cases[i]])

    if (vis=="graph"):
        for i in range(npar_):
            plot(np.array(range(1,ncases_+1)),sensdata[cases,i], '-o',color=colors[pars[i]], label=par_labels[pars[i]])
    elif (vis=="bar"):
        curr=np.zeros((ncases_))
        for i in range(npar_):
            bar(np.array(range(1,ncases_+1)),sensdata[cases,i], width=wd,color=colors[pars[i]], bottom=curr, label=par_labels[pars[i]])
            curr=sensdata[cases,i]+curr
        xticks(np.array(range(1,ncases_+1))+wd/2.,case_labels_)
        xlim(1-wd/2.,ncases_+1.5*wd)
    #ylim(0.0,1.0)
    #elif (vis=="stack"):
    #    fnx = lambda : np.random.randint(5, 50, 10)
    #    y = np.row_stack((fnx(), fnx(), fnx()))
    #    x = np.arange(10)
    #    print y
    #    fig, ax = plt.subplots()
    #    plt.stackplot(x, y)

    xlabel(xlbl)
    ylabel(ylbl)
    if (legend_show):
        legend(bbox_to_anchor=(1.0, -0.05),fancybox=True, shadow=True,ncol=5,labelspacing=-0.1)



    if (showplot):
        show()
    else:
        savefig('sens.eps')
        clf()





