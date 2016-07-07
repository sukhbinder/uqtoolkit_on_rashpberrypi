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
#ifndef PROBABILITY_H
#define PROBABILITY_H

#include "Array1D.h"
#include "Array2D.h"

#define DSFMT_DO_NOT_USE_OLD_NAMES
#include "dsfmt_add.h"



/** \file probability.h
 * Probability and random number generation- related tools.
 * \todo There shuold be a RNG class as a part of core UQTk - most of these functions will fit there.
 */

/// \brief An implementation of error function using incomplete gamma function
double erff(const double x);

/// \brief Inverse error function, input rescaled to [-1,1]
/// \note Cephes Math Library Release 2.8:  June, 2000.
/// Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
// The website (http://www.boutell.com/lsm/lsmbyid.cgi/000626) states copying policy=freely distributable as of July 2012
//  \note modified by Sandia UQTk group to scale the input to [-1,1]
double inverf(double y0);

/// \brief Inverse of the CDF of the normal random variable, uses inverf
double invnormcdf(double y);

/// \brief Normal random variable CDF
double normcdf(double y);

/// \brief Complementary function for normcdf
double normcdfc(double y);

/// \brief Generates a vector of i.i.d. uniform(0,1) random variable samples of size ns*nd, given integer seed
void generate_uniform(double* rvar,int ns, int nd, int zSeed);

/// \brief Generates a matrix of i.i.d. uniform(0,1) random variable samples, given integer seed
void generate_uniform(Array2D<double>& rvar,int zSeed);

/// \brief Generates a vector of i.i.d. uniform(0,1) random variable samples of size ns*nd, given pointer to the state of current random number generator
void generate_uniform(double *rvar, int ns, int nd, dsfmt_t *rnstate);

/// \brief Generates a matrix of i.i.d. uniform(0,1) random variable samples, given pointer to the state of current random number generator
void generate_uniform(Array2D<double> &rvar, dsfmt_t *rnstate);

/// \brief Generates a vector of i.i.d. uniform(0,1) random variable LHS samples of size ns*nd, given integer seed
void generate_uniform_lhs(double* rvar,int ns, int nd, int zSeed);

/// \brief Generates a matrix of i.i.d. uniform(0,1) random variable LHS samples, given integer seed
void generate_uniform_lhs(Array2D<double>& rvar,int zSeed);

/// \brief Generates a vector of i.i.d. uniform(0,1) random variable LHS samples of size ns*nd, given pointer to the state of current random number generator
void generate_uniform_lhs(double *rvar, int ns, int nd, dsfmt_t *rnstate);

/// \brief Generates a matrix of i.i.d. uniform(0,1) random variable LHS samples, given pointer to the state of current random number generator
void generate_uniform_lhs(Array2D<double> &rvar, dsfmt_t *rnstate);

/// \brief Generates a matrix of i.i.d. normal(0,1) random variable samples 
void generate_normal(Array2D<double>& rvar,int zSeed);

/// \brief Generates a matrix of i.i.d. normal(0,1) random variable LHS samples
/// \todo LHS generation is far from optimal, it is quite slow
void generate_normal_lhs(Array2D<double>& rvar,int zSeed);

/// \brief Returns the median of a data array
double get_median(const Array1D<double>& data);
/// \brief Returns the mean of a data array
double get_mean(const Array1D<double>& data);
/// \brief Returns the std of a data array
double get_std(const Array1D<double>& data);


void ihsU(Array2D<double> &rndnos, int dfac, dsfmt_t *rnstate) ;

void ihsU(int ns, int np, double *rndnos, int dfac, dsfmt_t *rnstate) ;

void ihsP(int ns, int np, int *rpos, int dfac, dsfmt_t *rnstate);

/// \todo move and merge this into arraytools.h tools
void shell_sort (int *a, int n) ;

//---------------------------------------------------------------------------------------
#endif // PROBABILITY_H
