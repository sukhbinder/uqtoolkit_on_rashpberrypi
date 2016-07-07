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
#ifndef PCMAPS_H
#define PCMAPS_H

/** \file pcmaps.h
 * Suite of functions to help map one kind of a PC variable to another.
 * \todo Perhaps use more robust tools, like dcdflib.
 */

/// \brief Implements a map y=f(x), where f is a function mapping one PC domain(pcIn with parameters in1,in2) 
/// to another (pcOut with parameters out1,out2) 
/// \note It also incorporates a 'TG' truncated-gaussian variable
/// \note Given cdf.dat it maps an arbitrary distribution to PC variables as well
double PCtoPC(double x, string pcIn, double in1, double in2, string pcOut, double out1, double out2);

/// \brief Bisection method modified to invert PCtoPC maps
double rtbis_mod(double func(double,string,double,double,string,double,double), const double x1, const double x2, const double xacc,double x, string pcIn, double in1, double in2, string pcOut, double out1, double out2);

/// \brief Implements an entriwise map yy=f(xx), where f is a function mapping one PC domain(pcIn with parameters in1,in2) 
/// to another (pcOut with parameters out1,out2) 
/// \note It also incorporates a 'TG' truncated-gaussian variable
/// \note Given cdf.dat it maps an arbitrary distribution to PC variables as well
void PCtoPC(Array2D<double>& xx, string pcIn, double in1, double in2, Array2D<double>& yy, string pcOut, double out1, double out2);

/// \brief Auxiliary linear interpolation function
void linint( Array2D<double> &xydata, const double x, double &y ) ;
/// \brief Auxiliary linear interpolation function
void linint( Array2D<double> &xydata, const double x, double &y, int col ) ;

//---------------------------------------------------------------------------------------
#endif // PCMAPS_H
