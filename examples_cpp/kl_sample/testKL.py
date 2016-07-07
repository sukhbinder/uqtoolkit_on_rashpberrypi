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

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math
from scipy import linalg

from pyutils import readfile, column, checkPtInside

Npts=129
Nspl=4096
Nk=21
cov_anl,nl=readfile('klcov_SqExp_0.05/cov_0.05_SqExp_anl.dat')
cov_num,nl=readfile('save/klsampl_0.05_'+str(Nspl)+'/cov_0.05_'+str(Nspl)+'.dat')

cov_anl=np.matrix(cov_anl)
cov_num=np.matrix(cov_num)

x = np.arange(0.0,1.0000001,1.0/(Npts-1))

# compute weights
w=np.array(np.zeros(Npts))
w[0]     =0.5*(x[1]     -x[0]    )
w[Npts-1]=0.5*(x[Npts-1]-x[Npts-2])
for i in range(1,Npts-1):
    w[i] = 0.5*(x[i+1]-x[i-1]);

# no change
ws=np.matrix(np.diag(np.sqrt(w)));
amt_anl=ws * cov_anl * ws ;
amt_num=ws * cov_num * ws ;
eva_anl,eve_anl=linalg.eigh(amt_anl, b=None, lower=True, eigvals_only=False);
eva_num,eve_num=linalg.eigh(amt_num, b=None, lower=True, eigvals_only=False);

# 1-diagonal
for i in range(Npts):
    cov_num[i,i] = cov_anl[i,i];

amt_num=ws * cov_num * ws ;
eva_num1,eve_num1=linalg.eigh(amt_num, b=None, lower=True, eigvals_only=False);

# 2-diagonal
for k in range(1,Nk):
    for i in range(Npts-k):
        cov_num[i,i+k] = cov_anl[i,i+k];
        cov_num[i+k,i] = cov_anl[i+k,i];

amt_num=ws * cov_num * ws ;
eva_numk,eve_numk=linalg.eigh(amt_num, b=None, lower=True, eigvals_only=False);


plt.figure()
for i in range(1,5):
    z_a=eve_anl[:,Npts-i]
    z_n=eve_num[:,Npts-i]
    z_n1=eve_num1[:,Npts-i]
    z_nk=eve_numk[:,Npts-i]
    if i == 1:
        if ( z_a[Npts/2]*z_n[Npts/2] < 0):
            z_n = z_n*(-1);
        if ( z_a[Npts/2]*z_n1[Npts/2] < 0):
            z_n1 = z_n1*(-1);
        if ( z_a[Npts/2]*z_nk[Npts/2] < 0):
            z_nk = z_nk*(-1);
    elif (i>1):
        ichk=30
        if ( z_a[ichk]*z_n[ichk] < 0):
            z_n  = z_n*(-1);
        if ( z_a[ichk]*z_n1[ichk] < 0):
            z_n1 = z_n1*(-1);
        if ( z_a[ichk]*z_nk[ichk] < 0):
            z_nk = z_nk*(-1);
    plt.subplot(int("22"+str(i)))
    plt.plot(x,z_a,color='r')
    plt.plot(x,z_n,color='b')
    plt.plot(x,z_n1,color='k')
    plt.plot(x,z_nk,color='g')
    plt.title('Mode '+str(i))

plt.savefig('num2anl_'+str(Nk)+'_0.05.eps')

