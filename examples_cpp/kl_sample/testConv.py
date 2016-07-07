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
import numpy as np
from   scipy import stats
import matplotlib
import matplotlib.pylab as plt
from   pylab import *

def readfile(filename):
    d0  = []
    nlines = 0
    for line in file(filename):
        line = line.rstrip('\n')
        line_list = [float(x) for x in line.split()]
        d0.append(line_list)
        nlines = nlines+1
    return d0,nlines

def column(matrix, i):
    return [row[i] for row in matrix]

clen="0.05"
dir="cl_"+clen
K1modes=[]
for i in range(9,18):
	KLmodes,nl=readfile(dir+"/KLmodes_"+clen+"_"+str(2**i)+".dat")
	K1=np.array(column(KLmodes,1));
	if K1.max()<0:
		K1=K1*(-1);
	K1modes.append(K1);

KLmodes,nl=readfile(dir+"/KLmodes_"+clen+"_SqExp.dat")
K1=np.array(column(KLmodes,1));
if K1.max()<0:
	K1=K1*(-1);

errors1=[]
for K1m in K1modes:
	ksum=0.0;
	nsum=0.0
	for j in range(len(K1m)):
		ksum=ksum+(K1m[j]-K1[j])**2;
		nsum=nsum+K1[j]**2;
	errors1.append(np.sqrt(ksum/nsum))


clen="0.10"
dir="cl_"+clen
K1modes=[]
for i in range(9,18):
	KLmodes,nl=readfile(dir+"/KLmodes_"+clen+"_"+str(2**i)+".dat")
	K1=np.array(column(KLmodes,1));
	if K1.max()<0:
		K1=K1*(-1);
	K1modes.append(K1);

KLmodes,nl=readfile(dir+"/KLmodes_"+clen+"_SqExp.dat")
K1=np.array(column(KLmodes,1));
if K1.max()<0:
	K1=K1*(-1);

errors2=[]
for K1m in K1modes:
	ksum=0.0;
	nsum=0.0
	for j in range(len(K1m)):
		ksum=ksum+(K1m[j]-K1[j])**2;
		nsum=nsum+K1[j]**2;
	errors2.append(np.sqrt(ksum/nsum))


clen="0.20"
dir="cl_"+clen
K1modes=[]
for i in range(9,18):
	KLmodes,nl=readfile(dir+"/KLmodes_"+clen+"_"+str(2**i)+".dat")
	K1=np.array(column(KLmodes,1));
	if K1.max()<0:
		K1=K1*(-1);
	K1modes.append(K1);

KLmodes,nl=readfile(dir+"/KLmodes_"+clen+"_SqExp.dat")
K1=np.array(column(KLmodes,1));
if K1.max()<0:
	K1=K1*(-1);

errors3=[]
for K1m in K1modes:
	ksum=0.0;
	nsum=0.0
	for j in range(len(K1m)):
		ksum=ksum+(K1m[j]-K1[j])**2;
		nsum=nsum+K1[j]**2;
	errors3.append(np.sqrt(ksum/nsum))


lw=2
fs=16
fig = plt.figure(figsize=(6,4))
ax=fig.add_axes([0.15, 0.15, 0.75, 0.75]) 
pleg=[]
pleg.append(plt.plot(errors1,linewidth=lw))
pleg.append(plt.plot(errors2,linewidth=lw))
pleg.append(plt.plot(errors3,linewidth=lw))

plt.xlabel("No samples",fontsize=fs)
plt.ylabel(r"$L_2(\sqrt{\lambda_1}f_1)$",fontsize=fs)
ax.set_xlim([-0.5,8.5])
ax.set_xticks([0,1,2,3,4,5,6,7,8])
ax.set_xticklabels(["$2^9$","$2^{10}$","$2^{11}$","$2^{12}$","$2^{13}$","$2^{14}$","$2^{15}$","$2^{16}$","$2^{17}$"])
ax.set_yscale("log")

leg=plt.legend( (pleg[0][0], pleg[1][0], pleg[2][0]),
                    (r"$0.05$", r"$0.10$", r"$0.20$"),'upper right' )
plt.savefig("L2error_KL1.eps")
