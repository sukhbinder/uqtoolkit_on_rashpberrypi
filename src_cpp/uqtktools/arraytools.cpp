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
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "assert.h"
#include <sstream>
#include <fstream>

#include "arraytools.h"
#include "ftndefs.h"
#include "blas.h"

using namespace std;


template <typename T>
void array1Dto2D(Array1D<T>& arr_1d,Array2D<T>& arr)
{
  int nd=arr_1d.XSize();
  arr.Resize(nd,1);

  for(int i=0;i<nd;i++)
    arr(i,0)=arr_1d(i);
       
  return;
}
template void array1Dto2D(Array1D<double>& arr_1d,Array2D<double>& arr);
template void array1Dto2D(Array1D<int>& arr_1d,Array2D<int>& arr);


template <typename T>
void array2Dto1D(Array2D<T>& arr_2d,Array1D<T>& arr)
{
  int nd=arr_2d.XSize();
  int one=arr_2d.YSize();
  if (one != 1) {printf("array2Dto1D: the second dimension must be one\n"); exit(1);}
  
  arr.Resize(nd);
  
  for(int i=0;i<nd;i++)
    arr(i)=arr_2d(i,0);


  return;
}
template void array2Dto1D(Array2D<double>& arr_2d,Array1D<double>& arr);
template void array2Dto1D(Array2D<int>& arr_2d,Array1D<int>& arr);


template <typename T>
void paste(Array1D<T>& arr1,Array1D<T>& arr2,Array2D<T>& arr)
{
  int n=arr1.XSize();
  int m=arr2.XSize();
  if (n != m) {printf("paste: the arrays must have the same size\n"); exit(1);}
    
  arr.Resize(n,2);
  for(int i=0;i<n;i++){
    arr(i,0)=arr1(i);
    arr(i,1)=arr2(i);
  }
  return;
}
template void paste(Array1D<double>& arr1,Array1D<double>& arr2,Array2D<double>& arr);
template void paste(Array1D<int>& arr1,Array1D<int>& arr2,Array2D<int>& arr);


void merge(Array2D<double>& x, Array2D<double>& y, Array2D<double>& xy)
{
  int nsample=xy.XSize();
  int ndim=xy.YSize();//=x(y).YSize(), too.
  int nsample1=x.XSize();
  int nsample2=y.XSize();
  if (nsample1+nsample2 != nsample) {printf("merge: dimension error\n"); exit(1);}
  
  for (int i=0;i<nsample1;i++){
    for(int idim=0;idim<ndim;idim++){
    xy(i,idim)=x(i,idim);
    }
  }
  for (int i=nsample1;i<nsample;i++){
    for(int idim=0;idim<ndim;idim++){
   xy(i,idim)=y(i-nsample1,idim);
  }
  }

  return;
}

void merge(Array1D<double>& x, Array1D<double>& y, Array1D<double>& xy)
{
  
  
  int ns1=x.XSize();
  int ns2=y.XSize();

  int ns12=ns1+ns2;

  xy.Resize(ns12,0.e0);

  for (int i=0;i<ns1;i++){
      xy(i)=x(i);
  }
  for (int i=ns1;i<ns12;i++){
      xy(i)=y(i-ns1);
  }

  return;
}

void merge(Array1D<int>& x, Array1D<int>& y, Array1D<int>& xy)
{
  
  
  int ns1=x.XSize();
  int ns2=y.XSize();

  int ns12=ns1+ns2;

  xy.Resize(ns12,0);

  for (int i=0;i<ns1;i++){
      xy(i)=x(i);
  }
  for (int i=ns1;i<ns12;i++){
      xy(i)=y(i-ns1);
  }

  return;
}


void append(Array1D<double>& x, Array1D<double>& y)
{
  
  
  int ns1 = x.XSize();
  int ns2 = y.XSize();
  int ns12= ns1+ns2;

  x.Resize(ns12);
  for (int i=ns1;i<ns12;i++) x(i)=y(i-ns1) ;
  return;
}

void append(Array1D<int>& x, Array1D<int>& y)
{
  
  
  int ns1 = x.XSize();
  int ns2 = y.XSize();
  int ns12= ns1+ns2;

  x.Resize(ns12);
  for (int i=ns1;i<ns12;i++) x(i)=y(i-ns1) ;
  return;
}

void copy(Array2D<double>& data_in, Array2D<double>& data_out)
{
int nsample=data_in.XSize();
 int ndim=data_in.YSize();

 for (int i=0;i<nsample;i++){
    for(int idim=0;idim<ndim;idim++){
      data_out(i,idim)=data_in(i,idim);
    }
    }
  return;
}


void transpose(Array2D<double>& x, Array2D<double>& xt)
{
  int nx=x.XSize();
  int ny=x.YSize();
  
xt.Resize(ny,nx,0.e0);

  for (int ix=0;ix<nx;ix++){
    for (int iy=0;iy<ny;iy++){
      xt(iy,ix)=x(ix,iy);
    }
  }

  return;
}


void unfold_2dto1d(Array2D<double>& x2, Array1D<double>& x1)
{
  int nx=x2.XSize();
  int ny=x2.YSize();
  int nxy=nx*ny;
  
  x1.Resize(nxy,0.e0);
  
  for(int i=0;i<nxy;i++){
   int i2= (i%ny);
   int i1= (i-i2)/ny;

    x1(i)=x2(i1,i2);
  
  }

  return;
}

void flatten(Array2D<double>& arr_2, Array1D<double>& arr_1)
{
    int nx=arr_2.XSize();
    int ny=arr_2.YSize();
    
    int nxy=nx*ny;
    
    arr_1.Resize(nxy,0.e0);
    
    for(int i=0;i<nx;i++){
	for(int j=0;j<ny;j++){
	    arr_1(j+i*ny)=arr_2(i,j);
	}
    }

    return;
}



void fold_1dto2d(Array1D<double>& x1, Array2D<double>& x2)
{
  int nx=x2.XSize();
  int ny=x2.YSize();
  int nxy=nx*ny;

  if (nxy != (int) x1.XSize()) {printf("fold_1dto2d: dimension error\n"); exit(1);}
   
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      x2(i,j)=x1(i*ny+j);

   
      }  
  }

  return;
}



void swap(Array1D<double>& arr,int i,int j)
{
  double h=arr(i);
  arr(i)=arr(j);
  arr(j)=h;
  return;
}

double access(int nx, int ny,Array1D<double>& arr_1, int i, int j)
{
    assert(nx*ny==(int) arr_1.XSize());

    

    return arr_1(j+i*ny);    
}


void moveArray1Dto2D(Array1D<double>& arr_1d,Array2D<double>& arr)
{
  int nd=arr_1d.XSize();
  arr.Resize(nd,1,0.e0);

  for(int i=0;i<nd;i++)
    arr(i,0)=arr_1d(i);
       
  

 return;
}

void moveArray1Dto2D(Array1D<int>& arr_1d,Array2D<int>& arr)
{
  int nd=arr_1d.XSize();
  arr.Resize(nd,1,0);

  for(int i=0;i<nd;i++)
    arr(i,0)=arr_1d(i);

 return;
}

void getRow(Array2D<double>& arr2d, int k, Array1D<double>& arr1d)
{
  arr1d.Clear();
  for(int i=0;i<(int)arr2d.YSize();i++)
    arr1d.PushBack(arr2d(k,i));

  return;
}

void getCol(Array2D<double>& arr2d, int k, Array1D<double>& arr1d)
{
  arr1d.Clear();
  for(int i=0;i<(int)arr2d.XSize();i++)
    arr1d.PushBack(arr2d(i,k));

  return;
}

void getRow(Array2D<int>& arr2d, int k, Array1D<int>& arr1d)
{
    arr1d.Clear();
    for(int i=0;i<(int)arr2d.YSize();i++)
        arr1d.PushBack(arr2d(k,i));
    
    return;
}

void getCol(Array2D<int>& arr2d, int k, Array1D<int>& arr1d)
{
    arr1d.Clear();
    for(int i=0;i<(int)arr2d.XSize();i++)
        arr1d.PushBack(arr2d(i,k));
    
    return;
}


void addVal(Array1D<double>& arr1d, double val)
{
  for(int i=0; i< (int)arr1d.XSize(); i++) arr1d(i) += val;
  return;
}

void addVal(Array1D<int>& arr1d, int val)
{
  for(int i=0; i< (int)arr1d.XSize(); i++) arr1d(i) += val;
  return;
}

void addVal(Array2D<double>& arr2d, double val)
{
  for(int j=0; j< (int)arr2d.YSize(); j++) 
    for(int i=0; i< (int)arr2d.XSize(); i++) 
      arr2d(i,j) += val;
  return;
}

void addVal(Array2D<int>& arr2d, int val)
{
  for(int j=0; j< (int)arr2d.YSize(); j++) 
    for(int i=0; i< (int)arr2d.XSize(); i++) 
      arr2d(i,j) += val;
  return;
}



void subVector(Array1D<double> &vector, Array1D<int> &ind, Array1D<double>& subvector)
{
    int n=vector.XSize();
    int k=ind.XSize();
    
    
    subvector.Resize(k,0.e0);
    for (int ik=0;ik<k;ik++){
        if (ind(ik)<0 || ind(ik)>=n){
            printf("subVector()::Index ind(%d)=%d is not allowed. Exiting.\n",ik,ind(ik));
            exit(1);
        }
        subvector(ik)=vector(ind(ik));
        
    }
    
    
    
    return;
}


void subVector(Array1D<int> &vector, Array1D<int> &ind, Array1D<int>& subvector)
{
    int n=vector.XSize();
    int k=ind.XSize();
    
    
    subvector.Resize(k,0);
    for (int ik=0;ik<k;ik++){
        if (ind(ik)<0 || ind(ik)>=n){
            printf("subVector()::Index ind(%d)=%d is not allowed. Exiting.\n",ik,ind(ik));
            exit(1);
        }
        subvector(ik)=vector(ind(ik));
        
    }
    
    
    
    return;
}

void subMatrix_row(Array2D<double> &matrix, Array1D<int> &ind, Array2D<double>& submatrix)
{
    int n=matrix.XSize();
    int m=matrix.YSize();
    int k=ind.XSize();
    
    
    submatrix.Resize(k,m,0.e0);
    for (int ik=0;ik<k;ik++){
        if (ind(ik)<0 || ind(ik)>=n){
            printf("subMatrix()::Index ind(%d)=%d is not allowed. Exiting.\n",ik,ind(ik));
            exit(1);
        }
        for (int im=0;im<m;im++){
            submatrix(ik,im)=matrix(ind(ik),im);
        }
    }
    
    
    
    return;
}


void subMatrix_row(Array2D<int> &matrix, Array1D<int> &ind, Array2D<int>& submatrix)
{
    int n=matrix.XSize();
    int m=matrix.YSize();
    int k=ind.XSize();
    
    
    submatrix.Resize(k,m,0);
    for (int ik=0;ik<k;ik++){
        if (ind(ik)<0 || ind(ik)>=n){
            printf("subMatrix()::Index ind(%d)=%d is not allowed. Exiting.\n",ik,ind(ik));
            exit(1);
        }
        for (int im=0;im<m;im++){
            submatrix(ik,im)=matrix(ind(ik),im);
        }
    }
    
    
    
    return;
}

void subMatrix_col(Array2D<double> &matrix, Array1D<int> &ind, Array2D<double>& submatrix)
{
    int n=matrix.XSize();
    int m=matrix.YSize();
    int k=ind.XSize();
    
    
    submatrix.Resize(n,k,0.e0);
    for (int ik=0;ik<k;ik++){
        if (ind(ik)<0 || ind(ik)>=m){
            printf("subMatrix()::Index ind(%d)=%d is not allowed. Exiting.\n",ik,ind(ik));
            exit(1);
        }
        for (int in=0;in<n;in++){
            submatrix(in,ik)=matrix(in,ind(ik));
        }
    }
    
    
    return;
}

void subMatrix_col(Array2D<int> &matrix, Array1D<int> &ind, Array2D<int>& submatrix)
{
    int n=matrix.XSize();
    int m=matrix.YSize();
    int k=ind.XSize();
    
    
    submatrix.Resize(n,k,0);
    for (int ik=0;ik<k;ik++){
        if (ind(ik)<0 || ind(ik)>=m){
            printf("subMatrix()::Index ind(%d)=%d is not allowed. Exiting.\n",ik,ind(ik));
            exit(1);
        }
        for (int in=0;in<n;in++){
            submatrix(in,ik)=matrix(in,ind(ik));
        }
    }
    
    
    return;
}



double maxVal(const Array1D<double>& vector, int *indx)
{
    
    double maxVal_ = vector(0);
    (*indx) = 0 ;
    for(int i=1;i< (int) vector.XSize();i++)
        if (vector(i) > maxVal_) 
        {
            maxVal_ = vector(i);
            (*indx) = i ;
        }
    return maxVal_;
    
}

int maxVal(const Array1D<int>& vector, int *indx)
{
    
    int maxVal_ = vector(0);
    (*indx) = 0 ;
    for(int i=1;i< (int) vector.XSize();i++)
        if (vector(i) > maxVal_) 
        {
            maxVal_ = vector(i);
            (*indx) = i ;
        }
    return maxVal_;
    
}

void setdiff(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C)
{
    C.Clear() ;
    bool fnd;
    for ( int i = 0; i < (int) A.XSize() ; i++ )
    {
        fnd = false;
        for ( int j = 0; j < (int) B.XSize() ; j++ )
            if ( A(i) == B(j) ) fnd = true ;
        if ( !fnd) C.PushBack(A(i));
    }

    /* order C in ascending order */
    shell_sort(C);
    return ;
    
}

void setdiff_s(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C)
{
  shell_sort(B);

  C.Clear() ;
  int j=0;

    for ( int i = 0; i < (int) A.XSize() ; i++ )
    {

      while(A(i)>B(j)){
	j++;
      }
      if ( A(i) < B(j) ) 
	C.PushBack(A(i));

    }

    return ;
    
}


void shell_sort(Array1D<int>& array)
{
     int flag = 1, length = array.XSize(), i;
     int temp;
     int d=length;
     while( flag || (d>1))      // boolean flag (true when not equal to 0)
     {
          flag = 0;           // reset flag to 0 to check for future swaps
         d = (d+1) / 2;
         for (i = 0; i < (length - d); i++)
        {
               if (array(i + d) < array(i))
              {
                      temp = array(i + d);      // swap items at positions i+d and d
                      array(i + d)= array(i);
                      array(i) = temp;
                      flag = 1;                  // indicate that a swap has occurred
              }
         }
     }
     return;
}

void shell_sort(Array1D<double>& array)
{
     int flag = 1, length = array.XSize(), i;
     double temp;
     int d=length;
     while( flag || (d>1))      // boolean flag (true when not equal to 0)
     {
          flag = 0;           // reset flag to 0 to check for future swaps
         d = (d+1) / 2;
         for (i = 0; i < (length - d); i++)
        {
               if (array(i + d) < array(i))
              {
                      temp = array(i + d);      // swap items at positions i+d and d
                      array(i + d)= array(i);
                      array(i) = temp;
                      flag = 1;                  // indicate that a swap has occurred
              }
         }
     }
     return;
}

void shell_sort_col(Array2D<double>& array,int col, Array1D<int>& newInd, Array1D<int>& oldInd)
{

     int flag = 1, length = array.XSize(), ncol=array.YSize(),i,j;
     double temp;
     int d=length, tmp;

    
     newInd.Resize(length,0);

    


     while( flag || (d>1))      // boolean flag (true when not equal to 0)
     {
          flag = 0;           // reset flag to 0 to check for future swaps
         d = (d+1) / 2;
         for (i = 0; i < (length - d); i++)
        {
               if (array(i + d,col) < array(i,col))
              {
		for (j=0;j<ncol;j++){
                      temp = array(i + d,j);      // swap items at positions i+d and d
                      array(i + d,j)= array(i,j);
                      array(i,j) = temp;
		}
		newInd(oldInd(i+d))=i;
		newInd(oldInd(i))=i+d;

		tmp=oldInd(i);
		oldInd(i)=oldInd(i+d);
		oldInd(i+d)=tmp;
                      flag = 1;                  // indicate that a swap has occurred
              }
         }
     }
     return;
}

void shell_sort_all(Array2D<double>& array,Array1D<int>& newInd, Array1D<int>& oldInd)
{

     int flag = 1, length = array.XSize(), ncol=array.YSize(),i,j;
     double temp;
     int d=length, tmp;

    
     newInd.Resize(length,0);

    


     while( flag || (d>1))      // boolean flag (true when not equal to 0)
     {
          flag = 0;           // reset flag to 0 to check for future swaps
         d = (d+1) / 2;
         for (i = 0; i < (length - d); i++)
        {
	  bool swflag=false;
	  for(int col=0;col<ncol;col++){
	    //if ( fabs(array(i + d,col) - array(i,col) ) > 1e-10 ){
	    if (fabs(array(i + d,col) - array(i,col) ) > 1e-10 && array(i + d,col) < array(i,col)){
	      swflag=true; break;
	    }
	    	    else if (fabs(array(i + d,col) - array(i,col) ) > 1e-10 && array(i + d,col) > array(i,col)){
		      swflag=false; break;
	    }
	    
	    else {}
	  }

               if (swflag)
              {
		for (j=0;j<ncol;j++){
                      temp = array(i + d,j);      // swap items at positions i+d and d
                      array(i + d,j)= array(i,j);
                      array(i,j) = temp;
		}
		newInd(oldInd(i+d))=i;
		newInd(oldInd(i))=i+d;

		tmp=oldInd(i);
		oldInd(i)=oldInd(i+d);
		oldInd(i+d)=tmp;
                      flag = 1;                  // indicate that a swap has occurred
              }
         }
     }
     return;
  
}


void shell_sort_incr(Array2D<double>& array,int col, Array1D<int>& newInd, Array1D<int>& oldInd)
{

     int flag = 1, length = array.XSize(), ncol=array.YSize(),i,j;
     double temp;
     int d=length, tmp;
     //int dim=array.YSize();
    
     newInd.Resize(length,0);
    
     while( flag || (d>1))      // boolean flag (true when not equal to 0)
     {
          flag = 0;           // reset flag to 0 to check for future swaps
         d = (d+1) / 2;
         for (i = 0; i < (length - d); i++)
        {
          
          bool pr_flag=true;
          if (col>0){
            if (array(i + d,col-1) > array(i,col-1)) pr_flag=false;
          }


               if (array(i + d,col) < array(i,col) && pr_flag)
              {
		for (j=0;j<ncol;j++){
                      temp = array(i + d,j);      // swap items at positions i+d and d
                      array(i + d,j)= array(i,j);
                      array(i,j) = temp;
		}
		newInd(oldInd(i+d))=i;
		newInd(oldInd(i))=i+d;

		tmp=oldInd(i);
		oldInd(i)=oldInd(i+d);
		oldInd(i+d)=tmp;
                      flag = 1;                  // indicate that a swap has occurred
              }
         }
     }
     return;
}




void intersect(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C,Array1D<int> &iA,Array1D<int> &iB)
{
    C.Clear() ;
    iA.Clear() ;
    iB.Clear() ;
    for ( int i = 0; i < (int) A.XSize() ; i++ )
        for ( int j = 0; j < (int) B.XSize() ; j++ )
            if ( A(i) == B(j) )
            {
                C.PushBack(A(i));
                iA.PushBack(i);
                iB.PushBack(j);
            }
    
    /* order C in ascending order */
    bool chgOrd=true;
    while (chgOrd)
    {
        chgOrd=false;
        for ( int i = 0; i < (int) C.XSize()-1 ; i++ )
        {
            if (C(i)>C(i+1))
            {
                chgOrd=true;
                int itmp ;
                itmp = C(i);  C(i)  = C(i+1);  C(i+1) =itmp ;
                itmp = iA(i); iA(i) = iA(i+1); iA(i+1)=itmp ;
                itmp = iB(i); iB(i) = iB(i+1); iB(i+1)=itmp ;
            }
        }
    }
    return ;
    
}

void intersect(Array1D<int> &A, Array1D<int> &B, Array1D<int> &C)
{
    C.Clear() ;
    for ( int i = 0; i < (int) A.XSize() ; i++ )
        for ( int j = 0; j < (int) B.XSize() ; j++ )
            if ( A(i) == B(j) )
                C.PushBack(A(i));
    
    /* order C in ascending order */
    bool chgOrd=true;
    while (chgOrd)
    {
        chgOrd=false;
        for ( int i = 0; i < (int) C.XSize()-1 ; i++ )
        {
            if (C(i)>C(i+1))
            {
                chgOrd=true;
                int itmp ;
                itmp = C(i);  C(i)  = C(i+1);  C(i+1) =itmp ;
            }
        }
    }
    return ;
    
}


void find(Array1D<double> &theta, double lambda, string type, Array1D<int> &indx)
{
    indx.Clear();
    if ( type == "gt" )
    {
        for ( int i = 0; i<(int) theta.XSize(); i++)
            if ( theta(i) > lambda ) indx.PushBack(i) ;
        return ;
    }
    if ( type == "ge" )
    {
        for ( int i = 0; i<(int) theta.XSize(); i++)
            if ( theta(i) >= lambda ) indx.PushBack(i) ;
        return ;
    }
    if ( type == "lt" )
    {
        for ( int i = 0; i<(int) theta.XSize(); i++)
            if ( theta(i) < lambda ) indx.PushBack(i) ;
        return ;
    }
    if ( type == "le" )
    {
        for ( int i = 0; i<(int) theta.XSize(); i++)
            if ( theta(i) <= lambda ) indx.PushBack(i) ;
        return ;
    }
    return ;
}

void find(Array1D<int> &theta, int lambda, string type, Array1D<int> &indx)
{
    indx.Clear();
    if ( type == "gt" )
    {
        for ( int i = 0; i<(int) theta.XSize(); i++)
            if ( theta(i) > lambda ) indx.PushBack(i) ;
        return ;
    }
    if ( type == "ge" )
    {
        for ( int i = 0; i<(int) theta.XSize(); i++)
            if ( theta(i) >= lambda ) indx.PushBack(i) ;
        return ;
    }
    if ( type == "lt" )
    {
        for ( int i = 0; i<(int) theta.XSize(); i++)
            if ( theta(i) < lambda ) indx.PushBack(i) ;
        return ;
    }
    if ( type == "le" )
    {
        for ( int i = 0; i<(int) theta.XSize(); i++)
            if ( theta(i) <= lambda ) indx.PushBack(i) ;
        return ;
    }
    if ( type == "eq" )
    {
        for ( int i = 0; i<(int) theta.XSize(); i++)
            if ( theta(i) == lambda ) indx.PushBack(i) ;
        return ;
    }
    return ;
}





void prodAlphaMatVec(Array2D<double>& A, Array1D<double>& x, double alpha, Array1D<double>& y)
{
    
    int n=A.XSize();
    int m=A.YSize();
    
    if ( m != (int) x.XSize() )
    {
        printf("prodAlphaMatVec() : Error : no. of columns in A and size of x are not the same : %d %d\n",
               m,(int)x.XSize());
        exit(1);
    }
    
    y.Resize(n,0.e0);
    
    char trans='n';
    double beta=0.e0;
    int xinc=1;
    int yinc=1;
    FTN_NAME(dgemv)(&trans, &n, &m, &alpha, A.GetArrayPointer(), &n, x.GetArrayPointer(), &xinc,  &beta, y.GetArrayPointer(), &yinc );
    

   /* 
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<m;j++)
            y(i) += A(i,j)*x(j) ;
        y(i) *= alpha ;
    }
    */
    
    
    return;
    
}

void prodAlphaMatTVec(Array2D<double>& A, Array1D<double>& x, double alpha, Array1D<double>& y)
{

    int n=A.YSize();
    int m=A.XSize();
    
    if ( m != (int) x.XSize() )
    {
        printf("prodAlphaMatTVec() : Error : no. of columns in A^T and size of x are not the same : %d %d",m,(int)x.XSize());
        exit(1);
    }


    y.Resize(n,0.e0);

    char trans='t';
    double beta=0.e0;
    int xinc=1;
    int yinc=1;
    FTN_NAME(dgemv)(&trans, &m, &n, &alpha, A.GetArrayPointer(), &m, x.GetArrayPointer(), &xinc,  &beta, y.GetArrayPointer(), &yinc );

      /*    
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<m;j++)
            y(i) += A(j,i)*x(j) ;
        y(i) *= alpha ;
    }
      */
    return;
    
}

void prodAlphaMatTMat(Array2D<double>& A, Array2D<double>& B, double alpha, Array2D<double>& C)
{
    
    int n=A.XSize();
    int m=A.YSize();
    int k=B.YSize();
    
    if ( n != (int) B.XSize() ) 
      {
	printf("prodAlphaMatTMat() : Error : no. of columns in A^T and no. of rows in B are not the same : %d %d",n,(int)B.XSize());
        exit(1);
      }
    
    C.Resize(m,k,0.e0);
    
    char transa='t';
    char transb='n';

    double beta=0.e0;

   
     FTN_NAME(dgemm)(&transa, &transb, &m, &k, &n, &alpha, A.GetArrayPointer(), &n, B.GetArrayPointer(), &n,  &beta, C.GetArrayPointer(), &m );

     /*
    for (int j=0;j<m;j++)
      for (int i=0;i<m;i++) {
        for ( int k=0; k<n; k++)
          C(i,j) += A(k,i)*B(k,j) ; // reversed indices for A
        C(i,j) *= alpha ;
      }
     */
    return;
    
}

void prodAlphaMatMat(Array2D<double>& A, Array2D<double>& B, double alpha, Array2D<double>& C)
{
    
    int n=A.YSize();
    int m=A.XSize();
    int k=B.YSize();
    
    if ( n != (int) B.XSize() )
    {
        printf("prodAlphaMatMat() : Error : no. of columns in A and no. of rows in B are not the same : %d %d",n,(int)B.XSize());
        exit(1);
    }
    
    C.Resize(m,k,0.e0);
    
    char transa='n';
    char transb='n';
    
    double beta=0.e0;
    
    
    FTN_NAME(dgemm)(&transa, &transb, &m, &k, &n, &alpha, A.GetArrayPointer(), &m, B.GetArrayPointer(), &n,  &beta, C.GetArrayPointer(), &m );
    
    /*
     for (int j=0;j<m;j++)
     for (int i=0;i<m;i++) {
     for ( int k=0; k<n; k++)
     C(i,j) += A(k,i)*B(k,j) ; // reversed indices for A
     C(i,j) *= alpha ;
     }
     */
    return;
    
}



void addVecAlphaVecPow(Array1D<double>& x, double alpha, Array1D<double>& y, int ip)
{
    
    if ( (int) x.XSize() != (int) y.XSize() )
    {
        printf("addVecAlphaVecPow() : Error : Array dimensions are not the same : %d %d",(int)x.XSize(),(int)y.XSize());
        exit(1);
    }
    
    for (int i = 0; i < (int) x.XSize(); i++)
        x(i) += alpha*pow(y(i),ip) ;
    
    return ;
    
}

void delRow(Array2D<double>& A, int irow)
{
    
    int n = A.XSize() ;
    int m = A.YSize() ;
    
    if ( n <= 1 || m == 0 ) return ;
    
    Array2D<double> B(n-1,m) ;
    for ( int i = 0; i < irow; i++ )
        for ( int j = 0; j < m; j++)   
            B(i,j) = A(i,j) ;
    for ( int i = irow+1; i < n; i++ )
        for ( int j = 0; j < m; j++)   
            B(i-1,j) = A(i,j) ;
    
    A.Resize(n-1,m);
    for ( int i = 0; i < n-1; i++ )
        for ( int j = 0; j < m; j++)   
            A(i,j) = B(i,j) ;
    
    return ;
    
}


void delCol(Array2D<double>& A, int icol)
{
    
    int n = A.XSize() ;
    int m = A.YSize() ;
    
    if ( n == 0 || m <= 1 ) return ;
    
    Array2D<double> B(n,m-1) ;
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < icol; j++)   
            B(i,j) = A(i,j) ;
    for ( int i = 0; i < n; i++ )
        for ( int j = icol+1; j < m; j++)   
            B(i,j-1) = A(i,j) ;
    
    A.Resize(n,m-1);
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < m-1; j++)   
            A(i,j) = B(i,j) ;
    
    return ;
    
}

void delCol(Array1D<double>& x, int icol)
{
    
    int n = x.XSize() ;
    
    if ( n == 0 ) return ;
    
    Array1D<double> y(n-1) ;
    for ( int i = 0; i<icol; i++)   y(i ) = x(i) ;
    for ( int i = icol+1; i<n; i++) y(i-1) = x(i) ;
    
    x.Resize(n-1);
    for ( int i = 0; i<n-1; i++) x(i) = y(i) ;
    
    return ;
    
}



void delCol(Array2D<int>& A, int icol)
{
    
    int n = A.XSize() ;
    int m = A.YSize() ;
    
    if ( n == 0 || m <= 1 ) return ;
    
    Array2D<int> B(n,m-1) ;
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < icol; j++)   
            B(i,j) = A(i,j) ;
    for ( int i = 0; i < n; i++ )
        for ( int j = icol+1; j < m; j++)   
            B(i,j-1) = A(i,j) ;
    
    A.Resize(n,m-1);
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < m-1; j++)   
            A(i,j) = B(i,j) ;
    
    return ;
    
}

void delCol(Array1D<int>& x, int icol)
{
    
    int n = x.XSize() ;
    
    if ( n == 0 ) return ;
    
    Array1D<int> y(n-1) ;
    for ( int i = 0; i<icol; i++)   y(i ) = x(i) ;
    for ( int i = icol+1; i<n; i++) y(i-1) = x(i) ;
    
    x.Resize(n-1);
    for ( int i = 0; i<n-1; i++) x(i) = y(i) ;
    
    return ;
    
}




void paddMatRow(Array2D<double>& A, Array1D<double>& x)
{
    
    int n = A.XSize() ;
    int m = A.YSize() ;
    
    assert( m == (int) x.XSize() ) ;
    
    Array2D<double> B ;
    B=A;
    A.Resize(n+1,m);
    for ( int i = 0; i < m; i++ )
    {
        for ( int j = 0; j < n; j++)   
            A(j,i) = B(j,i) ;
        A(n,i) = x(i) ;
    }
    
    return ;
    
}



void paddMatCol(Array2D<double>& A, Array1D<double>& x)
{
    
    int n = A.XSize() ;
    int m = A.YSize() ;
    
    assert( n == (int) x.XSize() ) ;
    
    Array2D<double> B ;
    B=A;
    A.Resize(n,m+1);
    for ( int i = 0; i < n; i++ )
    {
        for ( int j = 0; j < m; j++)   
            A(i,j) = B(i,j) ;
        A(i,m) = x(i) ;
    }
    
    return ;
    
}

void paddMatRow(Array2D<int>& A, Array1D<int>& x)
{
    
    int n = A.XSize() ;
    int m = A.YSize() ;
    
    assert( m == (int) x.XSize() ) ;
    
    Array2D<int> B ;
    B=A;
    A.Resize(n+1,m);
    for ( int i = 0; i < m; i++ )
    {
        for ( int j = 0; j < n; j++)   
            A(j,i) = B(j,i) ;
        A(n,i) = x(i) ;
    }
    
    return ;
    
}



void paddMatCol(Array2D<int>& A, Array1D<int>& x)
{
    
    int n = A.XSize() ;
    int m = A.YSize() ;
    
    assert( n == (int) x.XSize() ) ;
    
    Array2D<int> B ;
    B=A;
    A.Resize(n,m+1);
    for ( int i = 0; i < n; i++ )
    {
        for ( int j = 0; j < m; j++)   
            A(i,j) = B(i,j) ;
        A(i,m) = x(i) ;
    }
    
    return ;
    
}


void paddMatColScal(Array2D<double>& A, Array1D<double>& x, double scal)
{
    
    int n = A.XSize() ;
    int m = A.YSize() ;
    
    assert( n==m && n == (int) x.XSize() ) ;
    
    Array2D<double> B ;
    B=A;
    A.Resize(n+1,m+1);
    for ( int i = 0; i < n; i++ )
    {
        for ( int j = 0; j < m; j++)   
            A(i,j) = B(i,j) ;
        A(i,m) = A(n,i)=x(i) ;
    }
    A(n,m) = scal ;
    
    
    return ;
    
}



bool is_equal(Array1D<int>& a, Array1D<int>& b){
    
    int n=a.XSize();
    int m=b.XSize();
    if (n!=m)
        return false;
    
    for (int i=0;i<n;i++){
        
        if (a(i)!=b(i))
            return false;
        
    }
    
    return true;
}



int vecIsInArray(Array1D<int>& vec, Array2D<int>& array)
{
  int dd=array.XSize();
int dim=vec.XSize();
 if (dim != (int) array.YSize()) {printf("vecIsInArray: dimension error\n"); exit(1);}
  

  for(int i=0;i<dd;i++){
int fl=0;
for(int id=0;id<dim;id++){
    if (vec(id) == array(i,id))
      fl+=1;
    else break;
}
if (fl==dim)
return i;
  }
return -1;
}


double select_kth(int k, Array1D<double>& arr)
{
  int i,ir,j,l,mid;
  double a;

  int n=arr.XSize();
  l=0;
  ir=n-1;
  for(;;){
    if (ir<=l+1){
      if (ir==l+1 && arr(ir)<arr(l))
	swap(arr,l,ir);
      return arr(k);
    }
    else {
      mid=(l+ir) >> 1;
      swap(arr,mid,l+1);
      if (arr(l)>arr(ir))
	swap(arr,l,ir);
      if (arr(l+1)>arr(ir))
	swap(arr,l+1,ir);
      if (arr(l)>arr(l+1))
	swap(arr,l,l+1);
      i=l+1;
      j=ir;
      a=arr(l+1);
      for(;;){
	do i++; while (arr(i)<a);
	do j--; while (arr(j)>a);
	if (j<i) break;
	swap(arr,i,j);
      }
      arr(l+1)=arr(j);
      arr(j)=a;
      if (j>=k) ir=j-1;
      if (j<=k) l=i;
    }
  }
}
