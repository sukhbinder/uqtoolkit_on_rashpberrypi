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
#!/usr/bin/env python

"""
Load the locations of the grid points, created triangulation, and plots
the unstructured grid to file
"""

import matplotlib.pyplot as plt
import matplotlib.tri    as tri
from   matplotlib.lines  import Line2D                        
import numpy as np
import math

from pyutils import readfile, column, checkPtInside

def trArea(x,y):
    return (abs(x[0]*(y[1]-y[2])+x[1]*(y[2]-y[0])+x[2]*(y[0]-y[1]))/2.0) 

xyc,nl=readfile('data/cali.dat')
xc = np.array(column(xyc,0))
yc = np.array(column(xyc,1))

xy,nl=readfile('data/cali_cvt_4096.dat')
x = np.array(column(xy,0))
y = np.array(column(xy,1))

# Create the Triangulation; no triangles so Delaunay triangulation created.
triang = tri.Triangulation(x, y)

# Mask off unwanted triangles.
xmid = x[triang.triangles].mean(axis=1)
ymid = y[triang.triangles].mean(axis=1)
mask = checkPtInside(xc,yc,xmid,ymid);
triang.set_mask(mask)

# get valid traingles
forOut=triang.get_masked_triangles()

fig = plt.figure(figsize=[8,8])                                
ax = fig.add_axes([0.1,0.1,0.85,0.85])                             
ax.set_aspect('equal')
plt.plot(x,y,".",markersize=2)

ax.set_xlim([x.min(),x.max()])
ax.set_ylim([y.min(),y.max()])
ax.set_xlabel("lon",fontsize=18)
ax.set_ylabel("lat",fontsize=18)

xlo=-120.1; xhi=-119
ylo=37.4; yhi=38.4
ax.add_line(Line2D([xlo,-118],[ylo,38.55],color='r'));
ax.add_line(Line2D([xhi,-114.15],[ylo,38.55],color='r'));
ax.add_line(Line2D([xlo,-118],[yhi,41.9913],color='r'));
ax.add_line(Line2D([xhi,-114.15],[yhi,41.9913],color='r'));

ax.add_line(Line2D([xlo,xhi],[ylo,ylo],color='b',linewidth=2));
ax.add_line(Line2D([xhi,xhi],[ylo,yhi],color='b',linewidth=2));
ax.add_line(Line2D([xhi,xlo],[yhi,yhi],color='b',linewidth=2));
ax.add_line(Line2D([xlo,xlo],[yhi,ylo],color='b',linewidth=2));

ax1 = fig.add_axes([0.65,0.635,0.3,0.3])                             
ax1.set_aspect('equal')
ax1.set_xlim([xlo,xhi])
ax1.set_ylim([ylo,yhi])
ax1.set_xticks([])
ax1.set_yticks([])

for k in range(len(forOut)):
    x1=x[forOut[k,0]]
    x2=x[forOut[k,1]]
    x3=x[forOut[k,2]]
    y1=y[forOut[k,0]]
    y2=y[forOut[k,1]]
    y3=y[forOut[k,2]]
    drawtr=True
    if ( x1 < xlo ) | ( x1 > xhi ): drawtr=False
    if ( x2 < xlo ) | ( x2 > xhi ): drawtr=False
    if ( x3 < xlo ) | ( x3 > xhi ): drawtr=False
    if ( y1 < ylo ) | ( y1 > yhi ): drawtr=False
    if ( y2 < ylo ) | ( y2 > yhi ): drawtr=False
    if ( y3 < ylo ) | ( y3 > yhi ): drawtr=False
    if drawtr:
        aaa=ax1.add_line(Line2D([x1,x2],[y1,y2]));
        aaa=ax1.add_line(Line2D([x2,x3],[y2,y3]));
        aaa=ax1.add_line(Line2D([x3,x1],[y3,y1]));

plt.savefig("cali_trian.eps")
plt.show()

