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
#ifndef ARRAYTOOLS_H
#define ARRAYTOOLS_H

#include <stdlib.h>
#include "Array1D.h"
#include "Array2D.h"
#include "Array3D.h"

/** \file arraytools.h
 * Tools to manipulate array objects. Some tools mimick MATLAB functionalities.
 * \todo Some functions are not optimal in terms of array access. They could be slower than MATLAB counterparts.
 */


/// \brief Store a given 1d array in a 2d array with a single second dimension
template <typename T> void array1Dto2D(Array1D<T>& arr_1d,Array2D<T>& arr);

/// \brief Store a given 2d array with a single second dimension in a 1d array
template <typename T> void array2Dto1D(Array2D<T>& arr_2d,Array1D<T>& arr);

/// \brief Paste two 1d arrays of same size into a single 2d array with second dimension equal to two
template <typename T> void paste(Array1D<T>& arr1,Array1D<T>& arr2,Array2D<T>& arr);

/// \brief Merge 2d double arrays
/// \todo Proper size-checks needed
void merge(Array2D<double>& x, Array2D<double>& y, Array2D<double>& xy);
/// \brief Merge 1d double arrays
/// \todo Proper size-checks needed
void merge(Array1D<double>& x, Array1D<double>& y, Array1D<double>& xy);
/// \brief Merge 1d int arrays
/// \todo Proper size-checks needed
void merge(Array1D<int>&    x, Array1D<int>&    y, Array1D<int>&    xy);

/// \brief Append array y to array x in place (double format)
void append(Array1D<double>& x, Array1D<double>& y);
/// \brief Append array y to array x in place (int format)
void append(Array1D<int>&    x, Array1D<int>&    y);

/// \brief Copy array data_in to an array data_out
/// \todo A simple data_out=data_in also should work
void copy(Array2D<double>& data_in, Array2D<double>& data_out);

/// \brief Transpose a 2d double array x and return the result in xt
void transpose(Array2D<double>& x, Array2D<double>& xt);

/// \brief Unfold/flatten a 2d array into a 1d array (double format)
/// \note does the same as flatten
void unfold_2dto1d(Array2D<double>& x2, Array1D<double>& x1);

/// \brief Unfold/flatten a 2d array into a 1d array (double format)
/// \note does the same as unfold_2dto1d
void flatten(Array2D<double>& arr_2, Array1D<double>& arr_1);

/// \brief Fold a 1d array into a 2d array (double format)
/// \note The dimension of the 1d array needs to be equal to
///  the product of the dimensions of the 2d array
void fold_1dto2d(Array1D<double>& x1, Array2D<double>& x2);

/// \brief Swap i-th and j-th elements of the array arr
void swap(Array1D<double>& arr,int i,int j);

/// \brief Access element \f$j+i\times ny\f$ from 1D array 'arr_1'
double access(int nx, int ny, Array1D<double>& arr_1, int i, int j);

/// \brief Creates a 2d array 'arr' \f$n\times 1\f$ from 1D array 'arr_1d'
/// \todo this duplicates the function array1Dto2D above
void moveArray1Dto2D(Array1D<double>& arr_1d,Array2D<double>& arr);
/// \brief Creates a 2d array 'arr' \f$n\times 1\f$ from 1D array 'arr_1d'
/// \todo this duplicates the function array1Dto2D above
void moveArray1Dto2D(Array1D<int>&    arr_1d,Array2D<int>&    arr);

/// \brief Retrieves row 'k' from 2D array 'arr2d' and returns it in 1D array 'arr1d'
void getRow(Array2D<double>& arr2d, int k, Array1D<double>& arr1d);
/// \brief Retrieves column 'k' from 2D array 'arr2d' and returns it in 1D array 'arr1d'
void getCol(Array2D<double>& arr2d, int k, Array1D<double>& arr1d);
/// \brief Retrieves row 'k' from 2D array 'arr2d' and returns it in 1D array 'arr1d'
void getRow(Array2D<int>& arr2d, int k, Array1D<int>& arr1d);
/// \brief Retrieves column 'k' from 2D array 'arr2d' and returns it in 1D array 'arr1d'
void getCol(Array2D<int>& arr2d, int k, Array1D<int>& arr1d);

/// \brief Adds 'val' to all elements of 1D array arr1d
void addVal(Array1D<double>& arr1d, double val) ;
/// \brief Adds 'val' to all elements of 1D array arr1d
void addVal(Array1D<int>&    arr1d, int    val) ;
/// \brief Adds 'val' to all elements of 2D array arr2d
void addVal(Array2D<double>& arr2d, double val) ;
/// \brief Adds 'val' to all elements of 2D array arr2d
void addVal(Array2D<int>&    arr2d, int    val) ;

/// \brief Extracts from 'vector', elements corresponding to indices 'ind' and returns them in 'subvector'
void subVector(Array1D<double> &vector, Array1D<int> &ind, Array1D<double>& subvector);
/// \brief Extracts from 'vector', elements corresponding to indices 'ind' and returns them in 'subvector'
void subVector(Array1D<double> &vector, Array1D<int> &ind, Array1D<double>& subvector);
/// \brief Extracts from 'matrix' rows corresponding to indices 'ind' and returns them in 'submatrix'
void subMatrix_row(Array2D<double> &matrix, Array1D<int> &ind, Array2D<double>& submatrix);
/// \brief Extracts from 'matrix' rows corresponding to indices 'ind' and returns them in 'submatrix'
void subMatrix_row(Array2D<int> &matrix, Array1D<int> &ind, Array2D<int>& submatrix);
/// \brief Extracts from 'matrix' columns corresponding to indices 'ind' and returns them in 'submatrix'
void subMatrix_col(Array2D<double> &matrix, Array1D<int> &ind, Array2D<double>& submatrix);
/// \brief Extracts from 'matrix' columns corresponding to indices 'ind' and returns them in 'submatrix'
void subMatrix_col(Array2D<int> &matrix, Array1D<int> &ind, Array2D<int>& submatrix);


/// \brief Returns maximum value in 'vector' and its location
double maxVal(const Array1D<double>& vector, int *indx) ;
/// \brief Returns maximum value in 'vector' and its location
int maxVal(const Array1D<int>& vector, int *indx) ;

/// \brief Returns \f$ C=A\backslash B\f$ ( C=Elements of A that are not in B); C is sorted in ascending order
void setdiff(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C) ;

/// \brief Returns \f$ C=A\backslash B\f$ ( C=Elements of A that are not in B); C is sorted in ascending order
/// Assumes A is sorted and uses a faster algorithm than setdiff
/// \todo In future, this should sort A too and replace setdiff 
/// \note B is sorted on output as well
void setdiff_s(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C) ;

/// \brief Sorts array in ascending order
void shell_sort(Array1D<int>& array);

/// \brief Sorts array in ascending order
void shell_sort(Array1D<double>& array);

/// \brief Sorts the array in ascending order according to a certain column
void shell_sort_col(Array2D<double>& array,int col,Array1D<int>& newInd, Array1D<int>& oldInd);

void shell_sort_all(Array2D<double>& array,Array1D<int>& newInd, Array1D<int>& oldInd);

void shell_sort_incr(Array2D<double>& array,int col, Array1D<int>& newInd, Array1D<int>& oldInd);


/// \brief Finds common entries in 1D arrays 'A' and 'B' and returns them in 'C', sorted in ascending order. It also 
/// returns the original locations of these entries in 1D arrays 'iA' and 'iB', respectively
/// \note Currently, duplicated entries in either 'A' and 'B' will be duplicated in 'C'
void intersect(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C, Array1D<int> &iA,Array1D<int> &iB) ;
/// \brief Find common entries in 1D arrays 'A' and 'B' and return them in 'C', sorted in ascending order
/// \note Currently, duplicated entries in either 'A' and 'B' will be duplicated in 'C'
void intersect(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C) ;

/// \brief return list of indices corresponding to elements of 1D array theta that are: larger ( type="gt" ), 
/// larger or equal ( type="ge" ), smaller ( type="lt" ), smaller or equal ( type="le" ) than lambda
void find(Array1D<int>    &theta, int    lambda, string type, Array1D<int> &indx) ;
/// \brief return list of indices corresponding to elements of 1D array theta that are: larger ( type="gt" ), 
/// larger or equal ( type="ge" ), smaller ( type="lt" ), smaller or equal ( type="le" ) than lambda
void find(Array1D<double> &theta, double lambda, string type, Array1D<int> &indx) ;


/// \brief Returns \f$y=\alpha Ax\f$, where 'A' is a \f$\left[n\times m\right]\f$ 2D array, 'x' is 
/// 1D array of size \f$m\f$ and 'alpha' is a scalar. The 1D array 'y' has \f$n\f$ elements
void prodAlphaMatVec (Array2D<double>& A, Array1D<double>& x, double alpha, Array1D<double>& y) ;
/// \brief Returns \f$y=\alpha A^Tx\f$, where 'A' is a \f$\left[m\times n\right]\f$ 2D array, 'x' is 
/// 1D array of size \f$m\f$ and 'alpha' is a scalar. The 1D array 'y' has \f$n\f$ elements
void prodAlphaMatTVec(Array2D<double>& A, Array1D<double>& x, double alpha, Array1D<double>& y) ;
/// \brief Returns \f$C=\alpha AB\f$, where 'A' and 'B' are \f$\left[m\times n\right]\f$ 2D arrays
/// and 'alpha' is a scalar. The 2D array 'C' has \f$m\times m\f$ elements
void prodAlphaMatMat(Array2D<double>& A, Array2D<double>& B, double alpha, Array2D<double>& C);
/// \brief Returns \f$C=\alpha A^TB\f$, where 'A' and 'B' are \f$\left[m\times n\right]\f$ 2D arrays
/// and 'alpha' is a scalar. The 2D array 'C' has \f$m\times m\f$ elements
void prodAlphaMatTMat(Array2D<double>& A, Array2D<double>& B, double alpha, Array2D<double>& C) ;
/// \brief Returns \f$y_i=\alpha x_i^ip\f$, where 'x' and 'y' are 1D arrays with \f$n\f$ elements
void addVecAlphaVecPow(Array1D<double>& x, double alpha, Array1D<double>& y, int ip) ;


/// \brief Deletes row 'irow' from 2D array 'A'
void delRow(Array2D<double>& A, int irow) ;
/// \brief Deletes column 'icol' from 2D array 'A'
void delCol(Array2D<double>& A, int icol) ;
/// \brief Deletes element 'icol' from 1D array 'A'
void delCol(Array1D<double>& x, int icol) ;
/// \brief Deletes column 'icol' from 2D array 'A'
void delCol(Array2D<int>& A, int icol) ;
/// \brief Deletes element 'icol' from 1D array 'A'
void delCol(Array1D<int>& x, int icol) ;

/// \brief Padds 2D array 'A' with the row 'x'
/// \note the number of elements in 'x' should be the same as the number of rows of 'A'
void paddMatRow(Array2D<double>& A, Array1D<double>& x) ;
/// \brief Padds 2D array 'A' with the column 'x'
/// \note the number of elements in 'x' should be the same as the number of rows in 'A'
void paddMatCol(Array2D<double>& A, Array1D<double>& x) ;
/// \brief Padds 2D array 'A' with the row 'x'
/// \note the number of elements in 'x' should be the same as the number of rows of 'A'
void paddMatRow(Array2D<int>& A, Array1D<int>& x) ;
/// \brief Padds 2D array 'A' with the column 'x'
/// \note the number of elements in 'x' should be the same as the number of rows in 'A'
void paddMatCol(Array2D<int>& A, Array1D<int>& x) ;
/// \brief Padds square 2D array 'A' \f$\left[n\times n\right]\f$ with the elements of 'x' and 'scal' as follows:
/// \f$A_{n,i}=A_{i,n}=x_i\f$ and \f$A_{n,n}=scal\f$
void paddMatColScal(Array2D<double>& A, Array1D<double>& x, double scal) ;

/// \brief Checks if two 1d int arrays are equal
bool is_equal(Array1D<int>& a, Array1D<int>& b);

/// \brief Checks if vec matches with any of the rows of array
/// Returns the row number, or -1 if vec is not equal to any of the rows of array
int vecIsInArray(Array1D<int>& vec, Array2D<int>& array);

/// \brief Select the k-th smallest element of an array arr
double select_kth(int k, Array1D<double>& arr);


//---------------------------------------------------------------------------------------
#endif // ARRAYTOOLS_H
