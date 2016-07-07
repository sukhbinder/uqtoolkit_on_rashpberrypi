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
#ifndef MULTIINDEX_H
#define MULTIINDEX_H

#include "Array1D.h"
#include "Array2D.h"


/** \file multiindex.h
 * Functions that deal with integer multiindices.
 * \todo Multiindex could be a separate class and a part of core UQTk.
 */


/// \brief Computes the number of PC basis terms 
/// for Total-Order truncation with a given dimensionality and order
/// \note The formula is (ndim+norder)!/(ndim!norder!)  
int computeNPCTerms(int ndim,int norder);

/// \brief Computes the multiindex file of a PC basis 
/// for Total-Order truncation with a given dimensionality and order
/// Also, returns the number of terms.
int computeMultiIndex(int ndim,int norder, Array2D<int>& mi);

/// \brief Computes the number of PC basis terms 
/// for HDMR truncation with a given dimensionality and maxorders array 
/// that contains maximal orders per interaction dimensionalities.
int computeNPCTermsHDMR(int ndim, Array1D<int>& maxorders);

/// \brief Computes  the multiindex file of a PC basis 
/// for HDMR truncation with a given dimensionality and maxorders array 
/// that contains maximal orders per interaction dimensionalities.
int computeMultiIndexHDMR(int ndim, Array1D<int>& maxorders,Array2D<int>& mindex);

/// \brief Decode a multiindex from a sparse format to a regular format
/// \note For encoding and for more details on the format, see encodeMindex function of PCSet class
/// \sa PCSet.h
void decodeMindex(Array1D< Array2D<int> >& sp_mindex, int ndim, Array2D<int>& mindex);


/// \brief Given a multiindex it computes a new multiindex where only 'admissible' bases are added
/// \note A new basis is admissible, if by subtracting one order from any of the dimensions with 
/// non-zero order, one never leaves the set of old multiindices
void upOrder(Array2D<int>& mindex,Array2D<int>& new_mindex);

/// \brief A boolean check to see if a new basis term is admissible or not
bool is_admis(Array1D<int>& mindex_try,Array2D<int>& mindex);

/// \brief Given a multiindex array, it returns the orders of each basis term
/// \note Essentially, this function performs sums of each rows
void getOrders(Array2D<int>& mindex,Array1D<int>& orders);

int get_invmindex(Array1D<int> mi);

int get_invmindex_ord(Array1D<int> mi);


#endif // MULTIINDEX_H
