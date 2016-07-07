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
#include "Array1D.h"
#include "Array2D.h"

double maxVal(const Array1D<double>& vector)
{

  double maxVal_ = vector(0);
  for(int i=1;i< (int) vector.XSize();i++)
    if (vector(i) > maxVal_) maxVal_ = vector(i);
    
  return maxVal_;

}

int maxVal(const Array1D<int>& vector)
{

  int maxVal_ = vector(0);
  for(int i=1;i< (int) vector.XSize();i++)
    if (vector(i) > maxVal_) maxVal_ = vector(i);
    
  return maxVal_;

}

double maxVal(const Array2D<double>& vector)
{

  double maxVal_ = vector(0,0);
  for(int j=0;j< (int) vector.YSize();j++)
    for(int i=0;i< (int) vector.XSize();i++)
      if (vector(i,j) > maxVal_) maxVal_ = vector(i,j);

  return maxVal_;

}

int maxVal(const Array2D<int>& vector)
{

  int maxVal_ = vector(0,0);
  for(int j=0;j< (int) vector.YSize();j++)
    for(int i=0;i< (int) vector.XSize();i++)
      if (vector(i,j) > maxVal_) maxVal_ = vector(i,j);

  return maxVal_;

}

double minVal(const Array1D<double>& vector)
{
  double minVal_ = vector(0);
  for(int i=1;i<(int) vector.XSize();i++){
    if (vector(i) < minVal_){
      minVal_ = vector(i);
    }
  }
  return minVal_;
}

int minVal(const Array1D<int>& vector)
{
  int minVal_ = vector(0);
  for(int i=1;i<(int) vector.XSize();i++){
    if (vector(i) < minVal_){
      minVal_ = vector(i);
    }
  }
  return minVal_;
}

double minVal(const Array2D<double>& vector)
{

  double minVal_ = vector(0,0);
  for(int j=0;j< (int) vector.YSize();j++)
    for(int i=0;i< (int) vector.XSize();i++)
      if (vector(i,j) < minVal_) minVal_ = vector(i,j);

  return minVal_;

}

int minVal(const Array2D<int>& vector)
{

  int minVal_ = vector(0,0);
  for(int j=0;j< (int) vector.YSize();j++)
    for(int i=0;i< (int) vector.XSize();i++)
      if (vector(i,j) < minVal_) minVal_ = vector(i,j);

  return minVal_;

}

int maxIndex(const Array1D<double>& vector)
{
  int maxInd_ = 0;
  for(int i=1;i< (int) vector.XSize();i++){
    if (vector(i) > vector(maxInd_)){
      maxInd_ = i;
    }
  }
  return maxInd_;
}

int minIndex(const Array1D<double>& vector)
{
  int minInd_ = 0;
  for(int i=1;i< (int) vector.XSize();i++){
    if (vector(i) < vector(minInd_)){
      minInd_ = i;
    }
  }
  return minInd_;
}

int maxIndex(const Array1D<int>& vector)
{
  int maxInd_ = 0;
  for(int i=1;i< (int) vector.XSize();i++){
    if (vector(i) > vector(maxInd_)){
      maxInd_ = i;
    }
  }
  return maxInd_;
}

int minIndex(const Array1D<int>& vector)
{
  int minInd_ = 0;
  for(int i=1;i< (int) vector.XSize();i++){
    if (vector(i) < vector(minInd_)){
      minInd_ = i;
    }
  }
  return minInd_;
}


int maxIndexR_2D(const Array2D<double> a2d, const int irow)
{
  if ( ( irow < 0 ) ||( irow >= (int) a2d.XSize() ) ) {
    printf("Error in maxIndexR_2D() : illegal row index %d\n",irow) ;
    exit(1) ;
  }
 
  int maxInd_ = 0;
  for( int j = 1; j < (int) a2d.YSize(); j++){
    if (a2d(irow,j) > a2d(irow,maxInd_)){
      maxInd_ = j;
    }
  }

  return ( maxInd_ ) ;
}


int minIndexR_2D(const Array2D<double> a2d, const int irow)
{
  if ( ( irow < 0 ) ||( irow >= (int) a2d.XSize() ) ) {
    printf("Error in minIndexR_2D() : illegal row index %d\n",irow) ;
    exit(1) ;
  }
 
  int minInd_ = 0;
  for( int j = 1; j < (int) a2d.YSize(); j++){
    if (a2d(irow,j) < a2d(irow,minInd_)){
      minInd_ = j;
    }
  }

  return ( minInd_ ) ;
}



int maxIndexC_2D(const Array2D<double>& a2d, const int icol)
{
  if ( ( icol < 0 ) ||( icol >= (int) a2d.YSize() ) ) {
    printf("Error in maxIndexC_2D() : illegal column index %d\n",icol) ;
    exit(1) ;
  }
 
  int maxInd_ = 0;
  for( int i = 1; i < (int) a2d.XSize(); i++){
    if (a2d(i,icol) > a2d(maxInd_,icol)){
      maxInd_ = i;
    }
  }

  return ( maxInd_ ) ;
}



int minIndexC_2D(const Array2D<double>& a2d, const int icol)
{
  if ( ( icol < 0 ) ||( icol >= (int) a2d.YSize() ) ) {
    printf("Error in minIndexC_2D() : illegal column index %d\n",icol) ;
    exit(1) ;
  }
 
  int minInd_ = 0;
  for( int i = 1; i < (int) a2d.XSize(); i++){
    if (a2d(i,icol) < a2d(minInd_,icol)){
      minInd_ = i;
    }
  }

  return ( minInd_ ) ;
}
