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
"""
Plotting scripts for samples, covariance matrices, and KL modes for
1D and 2D configurations. It takes 2-4 command line arguments depending
on the value of the first argument
"""
import os
import sys
import numpy as np
from   scipy import stats
import matplotlib
import matplotlib.pylab as plt
import matplotlib.tri as tri
#from   pylab import *

from pyutils import readfile, column, checkPtInside

sigma=5.0
Npl=5

rtype="anlcov"
if (len(sys.argv) > 1):
    rtype  = sys.argv[1]
    if rtype == "samples":
        clen  = sys.argv[2]
    if rtype == "anlcov":
        ctype = sys.argv[2]
        clen  = sys.argv[3]
    if rtype == "numcov":
        clen  = sys.argv[2]
        nreal = sys.argv[3]
    if rtype == "anlKLevec":
        ctype = sys.argv[2]
        clen  = sys.argv[3]
    if rtype == "numKLevec":
        clen  = sys.argv[2]
        nreal = sys.argv[3]
    if rtype == "xidata":
        clen  = sys.argv[2]
        nreal = sys.argv[3]
    if (rtype == "samples2D") | (rtype == "samples2Du"):
        clen  = sys.argv[2]
        nreal = sys.argv[3]
        nspl  = int(sys.argv[4])
    if ( rtype == "anlcov2D" ) | ( rtype == "anlcov2Du" ):
        ctype = sys.argv[2]
        clen  = sys.argv[3]
    if ( rtype == "numcov2D" ) | ( rtype == "numcov2Du" ) :
        clen  = sys.argv[2]
        nreal = sys.argv[3]
    if ( rtype == "anlKLevec2D" ) | ( rtype == "anlKLevec2Du" ):
        ctype = sys.argv[2]
        clen  = sys.argv[3]
    if ( rtype == "numKLevec2D" ) | (rtype == "numKLevec2Du" ) :
        clen  = sys.argv[2]
        nreal = sys.argv[3]

if rtype == "samples":
    fname = "klsampl_"+clen+"_512/samples_"+clen+"_512.dat"
    print "Processing file ",fname
    din,nliles=readfile(fname);
    # Plot samples
    fs1=18
    xp=column(din,0);
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_axes([0.15, 0.15, 0.75, 0.75]) 
    for i in range(1,Npl+1):
        plt.plot(xp,column(din,i))
    plt.xlabel('x',fontsize=fs1)
    plt.ylabel('f(x)',fontsize=fs1)
    ax.set_ylim([-4*sigma,4*sigma])
    ax.set_yticks([-4*sigma,-2*sigma,0,2*sigma,4*sigma])
    plt.savefig("rf1D_"+clen+".eps")

if rtype == "anlcov":
    fname = "klcov_"+ctype+"_"+clen+"/cov_"+clen+"_"+ctype+"_anl.dat"
    print "Processing file ",fname
    cov,nlines=readfile(fname);
    vmax = np.array(cov).max()
    fig=plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.05, 0.05, 0.85, 0.85]) 
    plt.imshow(cov, interpolation='nearest', vmin=-0.5*vmax, vmax=vmax,
              cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    plt.title("$c_l=$"+clen)
    plt.savefig("cov_"+ctype+"_"+clen+"_anl.eps")

if rtype == "numcov":
    fname = "klsampl_"+clen+"_"+nreal+"/cov_"+clen+"_"+nreal+".dat"
    print "Processing file ",fname
    cov,nlines=readfile(fname);
    vmax = np.array(cov).max()
    fig=plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.05, 0.05, 0.85, 0.85]) 
    plt.imshow(cov, interpolation='nearest', vmin=-0.5*vmax, vmax=vmax,
              cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    plt.title(r"$c_l="+clen+", N_{\Theta}="+nreal+"$")
    plt.savefig("cov_"+clen+"_"+nreal+"_num.eps")

if rtype == "anlKLevec":
    # plot KL modes
    fname = "klcov_"+ctype+"_"+clen+"/KLmodes_"+clen+"_"+ctype+"_anl.dat"
    din,nliles=readfile(fname);
    #parameters
    lw=2
    fs=16
    Nkl=4
    xp=column(din,0);
    fig = plt.figure(figsize=(4,6))
    ax=fig.add_axes([0.17, 0.15, 0.75, 0.75]) 
    pleg=[]
    for i in range(1,Nkl+1):
    	y=column(din,i);
        if ( clen == "0.10" ):
            if ( i == 1 ):
                y=(-1)*np.array(y)
        if ( clen == "0.20" ):
            if ( i == 2 ) | (i == 3) :
                y=(-1)*np.array(y)
    	pleg.append(plt.plot(xp,y,linewidth=lw))
    plt.xlabel("x",fontsize=fs)
    plt.ylabel(r"$\sqrt{\lambda_k}f_k$",fontsize=fs)
    ax.set_ylim([-4,4])
    #ax.set_yticks([1.e-4,1.e-2,1.,1.e2])
    #ax.set_xlim([0,40])
    #ax.set_xticks([0,10,20,30,40])
    leg=plt.legend( (pleg[0][0], pleg[1][0], pleg[2][0], pleg[3][0]),
                    (r"$f_1$", r"$f_2$", r"$f_3$", r"$f_4$"),'lower right' )
    plt.title("$c_l=$"+clen)
    plt.savefig("KLmodes_"+ctype+"_"+clen+"_anl.eps")

if rtype == "numKLevec":
    # plot KL modes
    fname = "klsampl_"+clen+"_"+nreal+"/KLmodes_"+clen+"_"+nreal+".dat"
    din,nliles=readfile(fname);
    #parameters
    lw=2
    fs=16
    Nkl=4
    xp=column(din,0);
    fig = plt.figure(figsize=(4,6))
    ax=fig.add_axes([0.17, 0.15, 0.75, 0.75]) 
    pleg=[]
    for i in range(1,Nkl+1):
    	y=column(din,i);
        if ( clen == "0.20" ):
            if ( nreal == "131072" ) | ( nreal == "8192" ):
                if ( i == 1 ):
                    y=(-1)*np.array(y)
            if ( nreal == "512" ):
                if ( i == 2 ):
                    y=(-1)*np.array(y)
        if ( clen == "0.05" ):
            if ( nreal == "8192" ):
                if ( i == 2 ):
                    y=(-1)*np.array(y)
            if ( nreal == "131072" ):
                if ( i == 2 ) | ( i == 3 ):
                    y=(-1)*np.array(y)
        if ( clen == "0.10" ):
            if ( nreal == "512" ):
                if ( i == 1 ):
                    y=(-1)*np.array(y)
                if ( i == 3 ):
                    y=(-1)*np.array(y)
    	pleg.append(plt.plot(xp,y,linewidth=lw))
    plt.xlabel("x",fontsize=fs)
    plt.ylabel(r"$\sqrt{\lambda_k}f_k$",fontsize=fs)
    ax.set_ylim([-4,4])
    #ax.set_yticks([1.e-4,1.e-2,1.,1.e2])
    #ax.set_xlim([0,40])
    #ax.set_xticks([0,10,20,30,40])
    leg=plt.legend( (pleg[0][0], pleg[1][0], pleg[2][0], pleg[3][0]),
                    (r"$f_1$", r"$f_2$", r"$f_3$", r"$f_4$"),'lower right' )
    #plt.title(r"$c_l="+clen+", N_{\Theta}="+nreal+"$")
    plt.savefig("KLmodes_"+clen+"_"+nreal+".eps")

if rtype == "xidata":
    fname = "klsampl_"+clen+"_"+nreal+"/xidata_"+clen+"_"+nreal+".dat"
    din,nliles=readfile(fname);
    #parameters
    lw=2
    fs=16
    fig = plt.figure(figsize=(6,4))
    ax=fig.add_axes([0.15, 0.15, 0.75, 0.75]) 
    pleg=[]
    Nxi=4
    for i in range(Nxi):
        x=np.linspace(np.array(column(din,i)).min(),np.array(column(din,i)).max(),100)
        kerns=stats.kde.gaussian_kde(column(din,i))
        pleg.append(plt.plot(x,kerns(x),linewidth=lw))
    #pleg.append(plt.plot(x,np.exp(-x**2/(2*sigma*sigma))/(np.sqrt(2*np.pi)*sigma),linewidth=lw))
    plt.xlabel(r"$\xi(\theta)$",fontsize=fs)
    plt.ylabel(r"$PDF(\xi)$",fontsize=fs)
    ax.set_xlim([-4,4])
    ax.set_ylim([0.0,0.45])
    ax.set_yticks([0.0, 0.10, 0.20, 0.30, 0.40])
    ax.set_xticks([-4,-2,0,2,4])
    leg=plt.legend( (pleg[0][0], pleg[1][0], pleg[2][0], pleg[3][0]),
                    (r"$\xi_1$", r"$\xi_2$", r"$\xi_3$", r"$\xi_4$"),'upper right' )
    plt.title(r"$c_l="+clen+", N_{\Theta}="+nreal+"$")
    plt.savefig("xidata_"+clen+"_"+nreal+".eps")


if rtype == "samples2D":
    #parameters
    lw=2
    fs=16
    # plot KL modes
    fname = "klsampl2D_"+clen+"_"+nreal+"/samples2D_"+clen+"_"+nreal+".dat"
    din,nlines=readfile(fname);
    x,nx=readfile("klsampl2D_"+clen+"_"+nreal+"/xgrid.dat");
    y,ny=readfile("klsampl2D_"+clen+"_"+nreal+"/ygrid.dat");
    x=column(x,0);
    y=column(y,0);
    X,Y=np.meshgrid(x,y)
    for i in range(nspl):
        z=column(din,i);
        Z=[[z[k*ny+j] for j in range(ny)] for k in range(nx)]
        fig = plt.figure(figsize=(4,4))
        ax=fig.add_axes([0.15, 0.15, 0.75, 0.75]) 
        plt.contourf(X,Y,Z,101)
        plt.xlabel("x",fontsize=fs)
        plt.ylabel("y",fontsize=fs)
        ax.set_xlim([0.0,1.0])
        ax.set_ylim([0.0,1.0])
        ax.set_aspect('equal')
        plt.savefig("samples2D_"+clen+"_"+nreal+"_s"+str(i+1)+".eps")

if rtype == "anlcov2D":
    fname = "klcov2D_"+ctype+"_"+clen+"/cov2D_"+clen+"_"+ctype+"_anl.dat"
    print "Processing file ",fname
    cov,nlines=readfile(fname);
    vmax = np.array(cov).max()
    fig=plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.05, 0.05, 0.85, 0.85]) 
    plt.imshow(cov, interpolation='nearest', vmin=-0.5*vmax, vmax=vmax,
              cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    plt.title("$c_l=$"+clen)
    plt.savefig("cov2D_"+ctype+"_"+clen+"_anl.eps")

if rtype == "numcov2D":
    fname = "klsampl2D_"+clen+"_"+nreal+"/cov2D_"+clen+"_"+nreal+".dat"
    print "Processing file ",fname
    cov,nlines=readfile(fname);
    vmax = np.array(cov).max()
    fig=plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.05, 0.05, 0.85, 0.85]) 
    plt.imshow(cov, interpolation='nearest', vmin=-0.5*vmax, vmax=vmax,
              cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    plt.title(r"$c_l="+clen+", N_{\Theta}="+nreal+"$")
    plt.savefig("cov2D_"+clen+"_"+nreal+"_num.eps")

if rtype == "anlKLevec2D":
    #parameters
    lw=2
    fs=16
    # plot KL modes
    Nkl=6
    fname = "klcov2D_"+ctype+"_"+clen+"/KLmodes2D_"+clen+"_"+ctype+"_anl.dat"
    din,nlines=readfile(fname);
    x=[din[i][0] for i in range(65)];
    y=[din[65*i][1] for i in range(65)];
    X,Y=np.meshgrid(x,y)
    for i in range(2,Nkl+2):
        z=column(din,i);
        Z=[[z[k*65+j] for j in range(65)] for k in range(65)]
        fig = plt.figure(figsize=(4,4))
        ax=fig.add_axes([0.15, 0.15, 0.75, 0.75]) 
        plt.contourf(X,Y,Z,101)
        plt.xlabel("x",fontsize=fs)
        plt.ylabel("y",fontsize=fs)
        ax.set_xlim([0.0,1.0])
        ax.set_ylim([0.0,1.0])
        ax.set_aspect('equal')
        plt.title("$c_l=$"+clen+", Mode "+str(i-1))
        plt.savefig("KLmodes2D_"+ctype+"_"+clen+"_m"+str(i-1)+".eps")

if rtype == "numKLevec2D":
    #parameters
    lw=2
    fs=16
    # plot KL modes
    Nkl=8
    fname = "klsampl2D_"+clen+"_"+nreal+"/KLmodes2D_"+clen+"_"+nreal+".dat"
    din,nlines=readfile(fname);
    x=[din[i][0] for i in range(65)];
    y=[din[65*i][1] for i in range(65)];
    X,Y=np.meshgrid(x,y)
    for i in range(2,Nkl+2):
        z=column(din,i);
        Z=[[z[k*65+j] for j in range(65)] for k in range(65)]
        fig = plt.figure(figsize=(4,4))
        ax=fig.add_axes([0.15, 0.15, 0.75, 0.75]) 
        plt.contourf(X,Y,Z,101)
        plt.xlabel("x",fontsize=fs)
        plt.ylabel("y",fontsize=fs)
        ax.set_xlim([0.0,1.0])
        ax.set_ylim([0.0,1.0])
        ax.set_aspect('equal')
        #plt.title("$c_l=$"+clen+", Mode "+str(i-1))
        plt.savefig("KLmodes2D_"+clen+"_"+nreal+"_m"+str(i-1)+".eps")

if rtype == "samples2Du":
    #parameters
    lw=2
    fs=16
    # plot KL modes
    fname = "klsampl2Du_"+clen+"_"+nreal+"/samples2Du_"+clen+"_"+nreal+".dat"
    din,nlines=readfile(fname);
    xy,nxy=readfile("data/cali_grid.dat");
    x = np.array(column(xy,0))
    y = np.array(column(xy,1))
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    triGrid = tri.Triangulation(x,y)
    # Mask off unwanted triangles.
    xyc,nl=readfile('data/cali.dat')
    xc = np.array(column(xyc,0))
    yc = np.array(column(xyc,1))
    xmid = x[triGrid.triangles].mean(axis=1)
    ymid = y[triGrid.triangles].mean(axis=1)
    mask = checkPtInside(xc,yc,xmid,ymid);
    #mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
    triGrid.set_mask(mask)
    for i in range(nspl):
        fig = plt.figure(figsize=(4,4))
        ax=fig.add_axes([0.10, 0.10, 0.8, 0.8]) 
        ax.set_aspect('equal')
        plt.tricontourf(triGrid,column(din,i),101)
        ax.set_xlabel("lon",fontsize=fs)
        ax.set_ylabel("lat",fontsize=fs)
        ax.set_xlim([x.min(),x.max()])
        ax.set_ylim([y.min(),y.max()])
        ax.set_xticks([])
        ax.set_yticks([])
        #plt.colorbar()
        plt.title("$c_l=$"+clen)
        plt.savefig("cali_samples2D_"+clen+"_"+nreal+"_s"+str(i+1)+".eps")

if rtype == "anlcov2Du":
    fname = "klcov2Du_"+ctype+"_"+clen+"/cov2Du_"+ctype+"_"+clen+"_anl.dat"
    print "Processing file ",fname
    cov,nlines=readfile(fname);
    vmax = np.array(cov).max()
    fig=plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.05, 0.05, 0.85, 0.85]) 
    plt.imshow(cov, interpolation='nearest', vmin=-vmax, vmax=vmax,
              cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    plt.title("$c_l=$"+clen)
    plt.savefig("cov2Du_"+ctype+"_"+clen+"_anl.eps")

if rtype == "numcov2Du":
    fname = "klsampl2Du_"+clen+"_"+nreal+"/cov2Du_"+clen+"_"+nreal+".dat"
    print "Processing file ",fname
    cov,nlines=readfile(fname);
    vmax = np.array(cov).max()
    fig=plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.05, 0.05, 0.85, 0.85]) 
    plt.imshow(cov, interpolation='nearest', vmin=-vmax, vmax=vmax,
              cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    plt.title(r"$c_l="+clen+", N_{\Theta}="+nreal+"$")
    plt.savefig("cov2Du_"+clen+"_"+nreal+"_num.eps")

if rtype == "anlKLevec2Du":
    #parameters
    lw=2
    fs=16
    # plot KL modes
    Nkl=12
    fname = "klcov2Du_"+ctype+"_"+clen+"/KLmodes2Du_"+ctype+"_"+clen+"_anl.dat"
    din,nlines=readfile(fname);
    x = np.array(column(din,0))
    y = np.array(column(din,1))
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    triGrid = tri.Triangulation(x,y)
    # Mask off unwanted triangles.
    xyc,nl=readfile('data/cali.dat')
    xc = np.array(column(xyc,0))
    yc = np.array(column(xyc,1))
    xmid = x[triGrid.triangles].mean(axis=1)
    ymid = y[triGrid.triangles].mean(axis=1)
    mask = checkPtInside(xc,yc,xmid,ymid);
    #mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
    triGrid.set_mask(mask)
    for i in range(2,Nkl+2):
        fig = plt.figure(figsize=(4,4))
        ax=fig.add_axes([0.10, 0.10, 0.8, 0.8]) 
        ax.set_aspect('equal')
        plt.tricontourf(triGrid,column(din,i),101)
        ax.set_xlabel("lon",fontsize=fs)
        ax.set_ylabel("lat",fontsize=fs)
        ax.set_xlim([x.min(),x.max()])
        ax.set_ylim([y.min(),y.max()])
        ax.set_xticks([])
        ax.set_yticks([])
        #plt.colorbar()
        plt.title("$c_l=$"+clen+", Mode "+str(i-1))
        plt.savefig("cali_KLmodes2D_"+ctype+"_"+clen+"_m"+str(i-1)+".eps")


if rtype == "numKLevec2Du":
    #parameters
    lw=2
    fs=16
    # plot KL modes
    Nkl=12
    fname = "klsampl2Du_"+clen+"_"+nreal+"/KLmodes2Du_"+clen+"_"+nreal+".dat"
    din,nlines=readfile(fname);
    x = np.array(column(din,0))
    y = np.array(column(din,1))
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    triGrid = tri.Triangulation(x,y)
    # Mask off unwanted triangles.
    xyc,nl=readfile('data/cali.dat')
    xc = np.array(column(xyc,0))
    yc = np.array(column(xyc,1))
    xmid = x[triGrid.triangles].mean(axis=1)
    ymid = y[triGrid.triangles].mean(axis=1)
    mask = checkPtInside(xc,yc,xmid,ymid);
    #mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
    triGrid.set_mask(mask)
    for i in range(2,Nkl+2):
        fig = plt.figure(figsize=(4,4))
        ax=fig.add_axes([0.10, 0.10, 0.8, 0.8]) 
        ax.set_aspect('equal')
        plt.tricontourf(triGrid,column(din,i),101)
        ax.set_xlabel("lon",fontsize=fs)
        ax.set_ylabel("lat",fontsize=fs)
        ax.set_xlim([x.min(),x.max()])
        ax.set_ylim([y.min(),y.max()])
        ax.set_xticks([])
        ax.set_yticks([])
        #plt.colorbar()
        plt.title("$c_l=$"+clen+", Mode "+str(i-1))
        plt.savefig("cali_KLmodes2D_"+clen+"_"+nreal+"_m"+str(i-1)+".eps")





