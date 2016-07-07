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
#ifndef MINMAX_H
#define MINMAX_H

/** \file minmax.h
 * Tools to find min/max values of arrays.
 */

/// \brief Returns the maximum value of a 1d double array
double maxVal(const Array1D<double> &vector) ;
/// \brief Returns the maximum value of a 1d int array
int    maxVal(const Array1D<int>    &vector) ;
/// \brief Returns the maximum value of a 2d double array
double maxVal(const Array2D<double> &vector) ;
/// \brief Returns the maximum value of a 2d int array
int    maxVal(const Array2D<int>    &vector) ;

/// \brief Returns the minimum value of a 1d double array
double minVal(const Array1D<double> &vector) ;
/// \brief Returns the minimum value of a 1d int array
int    minVal(const Array1D<int>    &vector) ;
/// \brief Returns the minimum value of a 2d double array
double minVal(const Array2D<double> &vector) ;
/// \brief Returns the minimum value of a 2d int array
int    minVal(const Array2D<int>    &vector) ;

/// \brief Returns the index of the maximal value of a 1d double array
int maxIndex(const Array1D<double>& vector);
/// \brief Returns the index of the maximal value of a 1d int array
int maxIndex(const Array1D<int>& vector);
/// \brief Returns the index of the minimal value of a 1d double array
int minIndex(const Array1D<double>& vector);
/// \brief Returns the index of the minimal value of a 1d int array
int minIndex(const Array1D<int>& vector);


/// \brief Returns the column number of the maximal element in the irow-th row of a 2d double array
int maxIndexR_2D(const Array2D<double>& vector, const int irow);

/// \brief Returns the column number of the minimal element in the irow-th row of a 2d double array
int minIndexR_2D(const Array2D<double>& vector, const int irow);

/// \brief Returns the row number of the maximal element in the icol-th column of a 2d double array
int maxIndexC_2D(const Array2D<double>& vector, const int icol);

/// \brief Returns the row number of the minimal element in the icol-th column of a 2d double array
int minIndexC_2D(const Array2D<double>& vector, const int icol);

#endif // MINMAX_H