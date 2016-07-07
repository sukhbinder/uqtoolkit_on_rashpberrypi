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
/// \file Array2D.h
#ifndef ARRAY2D_H_SEEN
#define ARRAY2D_H_SEEN

#include <stddef.h>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include "Array1D.h"

using namespace std;

/// \class  Array2D
/// \brief  Stores data of any type T in a 2D array
///
/// This class also provides a Fortran-like access operator ()
/// as well as a function to access the data in the array through a pointer that
/// can be passed to F77 or C routines.
/// \author  Bert Debusschere <bjdebus@sandia.gov>
/// \date  Jan 2005
/// \note  Inpired by Helgi Adalsteinsson's Array class implementation
/// \todo  Define copy constructor
template <typename T>
class Array2D {
  public:
    /// \brief Default constructor, which does not allocate any memory
    Array2D(): xsize_(0), ysize_(0) {};

    /// \brief Constructor that allocates the memory
    Array2D(const size_t& nx, const size_t& ny):  xsize_(nx), ysize_(ny) {
    data_.resize(xsize_*ysize_);
    }

    /// \brief Constructor that allocates and initializes the data to a constant t
    Array2D(const size_t& nx, const size_t& ny, const T& t):  xsize_(nx), ysize_(ny) {
    data_.resize(xsize_*ysize_ , t);
    }

    /// \brief Copy constructor
    Array2D(const Array2D &obj): xsize_(obj.xsize_), ysize_(obj.ysize_), data_(obj.data_) {};
    
      /// \brief Destructor that frees up the memory
    ~Array2D() {data_.clear();}

    /// \brief Function to clear the memory
    void Clear() {
      xsize_ = 0;
      ysize_ = 0;
      data_.clear();
    }

    /// \brief Returns size in the x-direction
    size_t XSize() const {return xsize_;}
    /// \brief Returns size in the y-direction
    size_t YSize() const {return ysize_;}

    /// \brief Resizes the array
    /// \warning In its current implementation, most of the original data
    /// will get lost if the xsize changes as this changes the indexing for all entries.
    /// \todo Write a better implementation that preserves the original data by
    /// copying it to a temporary array and putting the elements back where they were before.
    /// This would bring this resize() command more closely in line with vector::resize()
    /// function in the original vector class.
    void Resize(const size_t& nx, const size_t& ny) {
      xsize_ = nx;
      ysize_ = ny;
      data_.resize(xsize_*ysize_);
    }

    /// \brief Resizes the array and sets ALL entries to the specified value
    /// \warning All original data will get lost if this function is used!
    /// \todo Write an implementation that is more closely follows the resize
    /// command in the vector class, which keeps the original elements and only
    /// initializes the new elements.
    void Resize(const size_t& nx, const size_t& ny, const T& t) {
      data_.clear();
      xsize_ = nx;
      ysize_ = ny;
      data_.resize(xsize_*ysize_, t);
    }

    /// \brief Set all values in the array to the given value
    void SetValue(const T& t){
      for(size_t i=0; i < data_.size(); i++){
        data_[i] = t;
      }
    }

    /// \brief Return a pointer to the first element of the data in the
    /// vector so we can use it access the data in array format (e.g. for
    /// passing it to a Fortran program).
    T* GetArrayPointer() {
      return &(data_[0]);
    }

    /// \brief Return a cont point to the first element of the data in the
    /// vector so we can use it access the data in array format (e.g. for
    /// passing it to a Fortran program).
    const T* GetConstArrayPointer() const {
      return &(data_[0]);
    }

    /// \brief Fortran-like () operator to access values in the 2D data array
    ///
    /// If "my_data" is an object of type Array2D, then its array values can
    /// be accessed as my_data(ix,iy), where ix, and iy are the indices in the
    /// x, and y dimensions respectively.
    T& operator()(size_t ix,size_t iy) {return data_[ix+xsize_*iy];}

    /// \brief Fortran-like () const operator to access values in the 2D data array
    ///
    /// If "my_data" is an object of type Array2D, then its array values can
    /// be accessed as my_data(ix,iy), where ix, and iy are the indices in the
    /// x, and y dimensions respectively.
    const T& operator()(size_t ix,size_t iy) const {return data_[ix+xsize_*iy];}

    /// \brief Insert array insarr as a row into position ix
    void insertRow(Array1D<T>& insarr,size_t ix){
        if (ix<0 || ix>xsize_)
            throw Tantrum("Array2D:insertRow():: insert index out of bounds.");
        if ( insarr.Length() != ysize_ )
            throw Tantrum("Array2D:insertRow():: insert row size does not match.");
        
	vector<T> data_old;
	data_old=data_;
	
	xsize_ += 1;
	data_.resize(xsize_*ysize_);

	for(size_t iy=0;iy<ysize_;iy++){	
	  for(size_t i=0; i < ix; i++)
	    data_[i+xsize_*iy] = data_old[i+(xsize_-1)*iy];
	  data_[ix+xsize_*iy]=insarr(iy);
	  for(size_t i=ix+1; i < xsize_; i++)
	    data_[i+xsize_*iy] = data_old[i-1+(xsize_-1)*iy];
	}
        
    }
    
    /// \brief Insert a 2d-array insarr into a row position ix
    void insertRow(Array2D<T>& insarr,size_t ix){
        if (ix<0 || ix>xsize_)
            throw Tantrum("Array2D:insertRow():: insert index out of bounds.");
        if ( insarr.YSize() != ysize_ )
            throw Tantrum("Array2D:insertRow():: insert row size does not match.");
        

	vector<T> data_old;
	data_old=data_;

	size_t insx=insarr.XSize();

	xsize_ += insx;
	data_.resize(xsize_*ysize_);

	for(size_t iy=0;iy<ysize_;iy++){	
	  for(size_t i=0; i < ix; i++)
	    data_[i+xsize_*iy] = data_old[i+(xsize_-insx)*iy];
	  for(size_t i=ix; i < ix+insx; i++)
	    data_[i+xsize_*iy]=insarr(i-ix,iy);
	  for(size_t i=ix+insx; i < xsize_; i++)
	    data_[i+xsize_*iy] = data_old[i-insx+(xsize_-insx)*iy];
	}

        
    }

   /// \brief Erase the row ix
    void eraseRow(size_t ix){
        if (ix<0 || ix>=xsize_)
            throw Tantrum("Array2D:eraseRow():: erase index out of bounds.");
      
	vector<T> data_old;
	data_old=data_;  

        xsize_-=1;
	data_.resize(xsize_*ysize_);

	for(size_t iy=0;iy<ysize_;iy++){	
	  for(size_t i=0; i < ix; i++)
	    data_[i+xsize_*iy] = data_old[i+(xsize_+1)*iy];
	  for(size_t i=ix; i < xsize_; i++)
	    data_[i+xsize_*iy] = data_old[i+1+(xsize_+1)*iy];
	}        

	if (xsize_==0)
	  printf("eraseRow(): WARNING: the xsize is zeroed!");

    }


    /// \brief Insert array insarr as a column into position iy
    void insertCol(Array1D<T>& insarr,size_t iy){
        if (iy<0 || iy>ysize_)
            throw Tantrum("Array2D:insertCol():: insert index out of bounds.");
        if ( insarr.Length() != xsize_ )
            throw Tantrum("Array2D:insertCol():: insert column size does not match.");
        
        
        T* ptr=insarr.GetArrayPointer();
        data_.insert(data_.begin()+xsize_*iy,ptr,ptr+xsize_);
        
        ysize_+=1;
        
    }

    /// \brief Insert a 2d-array insarr into a column position iy
    void insertCol(Array2D<T>& insarr,size_t iy){
        if (iy<0 || iy>ysize_)
            throw Tantrum("Array2D:insertCol():: insert index out of bounds.");
        if ( insarr.XSize() != xsize_ )
            throw Tantrum("Array2D:insertRow():: insert column size does not match.");
        
        size_t insy=insarr.YSize();
        
        T* ptr=insarr.GetArrayPointer();
        data_.insert(data_.begin()+xsize_*iy,ptr,ptr+xsize_*insy);
        
        
        ysize_+=insy;
        
    }
    
    /// \brief Erase the column iy
    void eraseCol(size_t iy){
        if (iy<0 || iy>=ysize_)
            throw Tantrum("Array2D:eraseCol():: erase index out of bounds.");
        
        data_.erase(data_.begin()+xsize_*iy,data_.begin()+xsize_*(iy+1));
        
        ysize_-=1;
        
	if (ysize_==0)
	  printf("eraseCol(): WARNING: the ysize is zeroed!");


    }


    
    /// \brief Dump contents of the array to a file in binary format
    void DumpBinary(FILE* f_out) const {
      fwrite(&xsize_,sizeof(xsize_),1,f_out);
      fwrite(&ysize_,sizeof(ysize_),1,f_out);
      fwrite(this->GetConstArrayPointer(),sizeof(T),xsize_*ysize_,f_out);
    }
    
   
    /// \brief Read contents of the array from a file in binary format
    void ReadBinary(FILE* f_in){
      fread(&xsize_,sizeof(xsize_),1,f_in);
      fread(&ysize_,sizeof(ysize_),1,f_in);
      data_.resize(xsize_*ysize_);
      fread(this->GetArrayPointer(),sizeof(T),xsize_*ysize_,f_in);
    }
    
   
  private:

    /// \brief Copy constructor, which is made private so it would not be used inadvertently
    /// (until we define a proper copy constructor)
    //Array2D(const Array2D &obj) {};

    /// \brief Number of elements in the x-dimension
    size_t xsize_;
    /// \brief Number of elements in the y-dimension
    size_t ysize_;

    /// \brief Data in the array with size = xsize_ * ysize_
    ///
    /// The data is stored with the fastest running index in the x-dimension
    /// and the slowest running index in the y-dimension. The indices in every dimension
    /// run from 0 to their respective "size-1"
    vector<T> data_;
};

#endif /* ARRAY2D_H_SEEN */
