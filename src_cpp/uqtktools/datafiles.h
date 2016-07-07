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
#ifndef DATAFILES_H
#define DATAFILES_H

#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "Array1D.h"
#include "Array2D.h"
#include "Array3D.h"

/** \file datafiles.h
 * Read/write capabilities from/to matrix or vector form arrays/files
 * \todo Perhaps it makes sense to move these functions to array class
 */


/// \brief Read a datafile from filename in a matrix form 
/// and store it in the double-array data
/// \note The array data needs to have the correct sizes
void read_datafile   (Array2D<double>& data, const char* filename);

/// \brief Read a datafile from filename in a matrix form 
/// and store it in the int-array data
/// \note The array data needs to have the correct sizes
void read_datafile   (Array2D<int>   & data, const char* filename);

/// \brief Read a datafile from filename in a matrix form 
/// and store it in the double-array data
/// \note The array data is resized to match the file contents
void read_datafileVS (Array2D<double>& data, const char* filename);

/// \brief Read a datafile from filename in a matrix form 
/// and store it in the int-array data
/// \note The array data is resized to match the file contents
void read_datafileVS (Array2D<int>   & data, const char* filename);

/// \brief Read a datafile from filename in a vector form 
/// and store it in the double-array data
/// \note The array data needs to have the correct size
void read_datafile_1d(Array1D<double>& data, const char* filename);

/// \brief Read a datafile from filename in a vector form 
/// and store it in the int-array data
/// \note The array data needs to have the correct size
void read_datafile_1d(Array1D<int>   & data, const char* filename);

/// \brief Read a datafile from filename in a vector form 
/// and store it in the double-array data
/// \note The array data is resized to match the file contents
void read_datafileVS(Array1D<double>& data, const char* filename);

/// \brief Read a datafile from filename in a vector form 
/// and store it in the int-array data
/// \note The array data is resized to match the file contents
void read_datafileVS(Array1D<int>& data, const char* filename);

/// \brief Write the contents of a 2d double array data
/// to a file filename in a matrix form 
void write_datafile   (const Array2D<double>& data, const char* filename);

/// \brief Write the contents of a 2d int array data
/// to a file filename in a matrix form 
void write_datafile   (const Array2D<int>   & data, const char* filename);

/// \brief Write the contents of a 1d double array data
/// to a file filename in a vector form 
void write_datafile_1d(const Array1D<double>& data, const char* filename);

/// \brief Write the contents of a 1d int array data
/// to a file filename in a vector form 
void write_datafile_1d(const Array1D<int>   & data, const char* filename);


#endif // DATAFILES_H
