/* =====================================================================================
                      The UQ Toolkit (UQTk) version 2.1.1
                     Copyright (2013) Sandia Corporation
                      http://www.sandia.gov/UQToolkit/

     Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
     with Sandia Corporation, the U.S. Government retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
 
     You should have received a copy of the GNU Lesser General Public License
     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
 
     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
#ifndef ROSENBLATT_H
#define ROSENBLATT_H

/// \note There is a version ov invRos() with an automatic bandwidth selection. However, this rule of thumb is not always reliable.
/// For the time being, that code stays in the general toolbox probability.cpp 
/// It is recommended to pick bandwidth off-line according to standard deviation per dimension.

/// \brief Generates new xi samples by inverse Rosenblatt using the ones we already have(xi) and uniform samples (given the bandwidth vector)
void invRos(Array1D<double>& unif, Array2D<double>& xi, Array1D<double>& newXi,Array1D<double>& sig);

/// \brief Generates new xi samples by inverse Rosenblatt using the ones we already have(xi) and uniform samples (given the same bandwidth for all dimensions)
void invRos(Array1D<double>& unif, Array2D<double>& xi, Array1D<double>& newXi,double bw);


/// \brief Generates new xi samples by inverse Rosenblatt using the ones we already have(xi) and uniform samples (picking the optimal bandwidth)
void invRos(Array1D<double>& unif, Array2D<double>& xi, Array1D<double>& newXi);

/// \brief Calculates 'rule of thumb' optimal KDE bandwidth for the Multi-D data
void get_opt_KDEbdwth(const Array2D<double>& data,Array1D<double>& bdwth);

/// \brief An outdated version with a homemade bisection algorithm that seems slightly faster
void invRos_old(Array1D<double>& unif, Array2D<double>& xi, Array1D<double>& newXi, Array1D<double>& sig);

#endif // ROSENBLATT_H
