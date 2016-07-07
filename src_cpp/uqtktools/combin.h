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
#ifndef COMBIN_H
#define COMBIN_H

#include "Array2D.h"

/** \file combin.h
 * Tools to evaluate combinatorial quantities.
 * \note Some functions are likely not optimal and could have been computed more efficiently.
 */

/// \brief Calculates n-choose-k
int choose(int n,int k);

/// \brief Calculates the factorial of a number
int factorial(int number);

/// \brief Computes all possible k-combinations of the first n non-negative integers 
/// and returns them in fullInd
void chooseComb(int n, int k,Array2D<int>& fullInd);


/// \brief Computes a random permutation of the first n non-negative integers
/// and returns is in perm 
void get_perm(int n, int* perm,int seed);

/// \brief Computes a random permutation of the first n non-negative integers
/// and returns is in perm 
/// \note n is the size of the array argument perm
void get_perm(Array1D<int>& perm, int seed);

 
/// \brief Compute the incomplete Gamma function with parameter a at point x
/// \note This is a slightly modified version of a code distributed by John Burkardt
/// \note see http://people.sc.fsu.edu/~jburkardt/cpp_src/asa147/asa147.html
/// \note see comments under the function definition
double gammai(const double a, const double x);

/// \brief Compute the Beta function at the point pair (z,w) 
double beta(const double z, const double w);

/// \brief Compute the incomplete Beta function with parameters a and b at point x
/// \note This is a slightly modified version of a code distributed by John Burkardt
/// \note see http://people.sc.fsu.edu/~jburkardt/cpp_src/asa063/asa063.html
/// \note see comments under the function definition
double betai(const double p, const double q, const double x);

/// \brief Computes the digamma, or psi, function, i.e. derivative of the logarithm of gamma function
/// \note This is a slightly modified version of a code distributed by John Burkardt
/// \note see http://people.sc.fsu.edu/~jburkardt/cpp_src/asa103/asa103.cpp
double digama ( double x );



//---------------------------------------------------------------------------------------
#endif // COMBIN_H
