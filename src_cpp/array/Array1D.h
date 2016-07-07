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
/// \file Array1D.h
#ifndef ARRAY1D_H_SEEN
#define ARRAY1D_H_SEEN

#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>

#include "error_handlers.h"

using namespace std;

/// \class  Array1D
/// \brief  Stores data of any type T in a 1D array
///
/// This class also provides a Fortran-like access operator ()
/// as well as a function to access the data in the array through a pointer that
/// can be passed to F77 or C routines.
/// \author Bert Debusschere <bjdebus@sandia.gov>
/// \date Apr 2005 - Nov 2007
/// \note Inpired by Helgi Adalsteinsson's Array class implementation
/// \todo double check copy constructor
template <typename T>
class Array1D {
  public:
    /// \brief Default constructor, which does not allocate any memory
    Array1D(): xsize_(0) {};

    /// \brief Constructor that allocates the memory
    Array1D(const size_t& nx): xsize_(nx) {
      data_.resize(xsize_);
    }

    /// \brief Constructor that allocates and initializes the data to a value t
    Array1D(const size_t& nx, const T& t): xsize_(nx) {
      data_.resize(xsize_, t);
    }

    /// \brief Assignment operator copies the data structure by value
    Array1D& operator=(const Array1D &obj) {
      xsize_ = obj.xsize_;
      data_ = obj.data_;
      return *this;
    }

    /// \brief Copy constructor
    Array1D(const Array1D &obj): xsize_(obj.xsize_), data_(obj.data_) {};

    /// \brief Destructor that frees up the memory
    ~Array1D() {data_.clear();}

    /// \brief Function to clear the memory
    void Clear() {
      xsize_ = 0;
      data_.clear();
    }

    /// \brief Returns size in the x-direction
    size_t XSize() const {return xsize_;}

    /// \brief Returns length (i.e. size in the x-direction)
    size_t Length() const {return xsize_;}

    /// \brief Resizes the array
    void Resize(const size_t& nx) {
      xsize_ = nx;
      data_.resize(xsize_);
    }

    /// \brief Resizes the array and sets ALL entries to the specified value
    /// \warning All original data will get lost if this function is used!
    /// \todo Write an implementation that is more closely follows the resize
    /// command in the vector class, which keeps the original elements and only
    /// initializes the new elements.
    void Resize(const size_t& nx, const T& t) {
      data_.clear();
      xsize_ = nx;
      data_.resize(xsize_, t);
    }

    /// \brief Set all values in the array to the given value
    void SetValue(const T& t){
      for(size_t i=0; i < data_.size(); i++){
        data_[i] = t;
      }
    }

    /// \brief Add element to the end of the vector
    void PushBack(const T& t){
       xsize_ += 1;
       data_.push_back(t);
    }

    /// \brief Return a pointer to the first element of the data in the
    /// vector so we can use it access the data in array format (e.g. for
    /// passing it to a Fortran program).
    T* GetArrayPointer() {
      return &(data_[0]);
    }

    /// \brief Return a const point to the first element of the data in the
    /// vector so we can use it access the data in array format (e.g. for
    /// passing it to a Fortran program).
    const T* GetConstArrayPointer() const {
      return &(data_[0]);
    }

    /// \brief Fortran-like () operator to access values in the 1D data array
    ///
    /// If "my_data" is an object of type Array1D, then its array values can
    /// be accessed as my_data(ix), where ix is the index in the
    /// x dimension.
    T& operator()(size_t ix) {return data_[ix];}

    /// \brief Fortran-like () const operator to access values in the 1D data array
    ///
    /// If "my_data" is an object of type Array1D, then its array values can
    /// be accessed as my_data(ix), where ix is the index in the
    /// x dimension.
    const T& operator()(size_t ix) const {return data_[ix];}

    /// \brief Insert a given array to the position ix
    /// \note ix=0 means insert at the beginning, ix=xsize_ means insert at the end
    void insert(Array1D<T>& insarr,size_t ix){
        if (ix<0 || ix>xsize_)
            throw Tantrum("Array1D:insert():: insert index out of bounds.");
        size_t addsize=insarr.Length();
        xsize_+=addsize;
        T* ptr=insarr.GetArrayPointer();
        data_.insert(data_.begin()+ix,ptr,ptr+addsize);
    }
    
    /// \brief Insert a given value to the position ix
    /// \note ix=0 means insert at the beginning, ix=xsize_ means insert at the end
    void insert(const T& insval,size_t ix){
        if (ix<0 || ix>xsize_)
            throw Tantrum("Array1D:insert():: insert index out of bounds.");
        xsize_+=1;
        data_.insert(data_.begin()+ix,insval);
    }
    
    /// \brief Erase the value from the position ix
    void erase(size_t ix){
        if (ix<0 || ix>=xsize_)
            throw Tantrum("Array1D:erase():: erase index out of bounds.");
        xsize_-=1;
        data_.erase(data_.begin()+ix);
    }
    

    /// \brief Dump contents of the array to a file in binary format
    void DumpBinary(FILE* f_out) const {
      fwrite(&xsize_,sizeof(xsize_),1,f_out);
      fwrite(this->GetConstArrayPointer(),sizeof(T),xsize_,f_out);
    }
    
   
    /// \brief Read contents of the array from a file in binary format
    void ReadBinary(FILE* f_in){
      fread(&xsize_,sizeof(xsize_),1,f_in);
      data_.resize(xsize_);
      fread(this->GetArrayPointer(),sizeof(T),xsize_,f_in);
    }

  private:

    /// \brief Copy constructor, which is made private so it would not be used inadvertently
    /// (until we define a proper copy constructor)
    ///Array1D(const Array1D &obj) {};

    /// \brief Number of elements in the x-dimension
    size_t xsize_;

    /// \brief Data in the array with size = xsize_
    ///
    /// The data is stored with the fastest running index in the x-dimension
    /// The index runs from 0 to "size-1"
    vector<T> data_;
};

#endif /* ARRAY1D_H_SEEN */
