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

import numpy as np


def set_colors(npar):
    """ Sets a list of different colors of requested length"""
    colors = []
    pp=1+npar/6
    for i in range(npar):
        c=1-(float) (i/6)/pp
        b=np.empty((3))
        for jj in range(3):
            b[jj]=c*int(i%3==jj)
        a=int(i%6)/3
        colors.append(((1-a)*b[2]+a*(c-b[2]),(1-a)*b[1]+a*(c-b[1]),(1-a)*b[0]+a*(c-b[0])))
    
    return colors

