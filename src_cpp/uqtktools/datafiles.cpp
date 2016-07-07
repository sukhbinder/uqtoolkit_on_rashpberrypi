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
#include "datafiles.h"

using namespace std;

/*
  Note: the *VS() functions make two passes: the first pass figures the no. or rows and columns, 
  then the data array is appropriately resized, and the filename is read during second pass
*/


// Read a datafile from filename in a matrix form and store it in the double-array data
void read_datafile(Array2D<double>& data, const char* filename)
{
  int nx=data.XSize();
  int ny=data.YSize();  

  if (nx==0 || ny==0){
    printf("read_datafile()    : the requested data array is empty\n") ;
    exit(1) ;
  }
  ifstream in(filename);
  
  if(!in){
    printf("read_datafile()    : the requested file %s does not exist\n",filename) ;
    exit(1) ;
  }
  
  string theLine="";
  int ix=0;

  while(in.good()){
    getline(in,theLine);

    if (theLine=="") break;
    if ( theLine.compare(0, 1, "#") == 0 ) continue ;

    istringstream s(theLine);
    int iy=0;
    double tmp;
    while(s>>tmp){
      data(ix,iy)=tmp;
      iy++;
    }
    if (iy!=ny) {printf("Error at line %d while reading %s; number of columns should be %d\n", ix+1, filename,ny); exit(1);}
    ix++;
  }
  if (ix!=nx) {printf("Error while reading %s; number of rows should be %d\n", filename,nx); exit(1);}

  in.close();

  return;
}


//  Read a datafile from filename in a matrix form and store it in the int-array data
void read_datafile(Array2D<int>& data, const char* filename)
{
  int nx=data.XSize();
  int ny=data.YSize();  

  if (nx==0 || ny==0){
    printf("read_datafile()    : the requested data array is empty\n") ;
    exit(1) ;
  }

 ifstream in(filename);
  
  if(!in){
    printf("read_datafile()    : the requested file %s does not exist\n",filename) ;
    exit(1) ;
  }
  
  string theLine="";

  int ix=0;

  
  while(in.good()){
    getline(in,theLine);
    
    if (theLine=="") break;
    if ( theLine.compare(0, 1, "#") == 0 ) continue ;

    istringstream s(theLine);
    int iy=0;
    int itmp;
    while(s>>itmp){
      data(ix,iy)=itmp;
      iy++;
    }
    if (iy!=ny) {printf("Error at line %d while reading %s; number of columns should be %d\n", ix+1, filename,ny); exit(1);}
    ix++;
  }
  if (ix!=nx) {printf("Error while reading %s; number of rows should be %d\n", filename,nx); exit(1);}


 return;
}


// Read a datafile from filename in a vector form and store it in the double-array data
void read_datafileVS(Array2D<double>& data, const char* filename)
{

  ifstream in(filename);
  
  if(!in){
    printf("read_datafileVS()    : the requested file %s does not exist\n",filename) ;
    exit(1) ;
  }
  
  string theLine="";


  // figure out number of lines and columns
  int nx, ny, ix = 0 ;
  while(in.good()){
    getline(in,theLine);
    
    if ( theLine == "" ) break;
    if ( theLine.compare(0,1,"#") == 0 ) continue ;

    istringstream s(theLine);
    int    iy = 0 ;
    double tmp    ;
    while( s >> tmp ) iy++ ;

    if ( ( ix > 0 ) && ( iy != ny ) )
    {
      printf("read_datafileVS() : Error at line %d !!!\n",ix+1) ;
      printf("                    no. of columns should be %d instead of %d\n",ny,iy) ;
      exit(1) ;
    }
    
    ny = iy ;

    ix++ ;

  }

  nx = ix ;

#ifdef VERBOSE
  printf("File \"%s\" contains %d rows and %d columns \n",filename,nx,ny) ;
#endif
  // Resize, goto beginning, and read again

  if ( ( (int) data.XSize() != nx ) || ( (int) data.YSize() != ny ) )
    data.Resize(nx,ny) ;

  //in.close() ;
  //in.open(filename);
  in.clear() ;
  in.seekg(0, ios::beg ) ;
  ix = 0 ;
  while( in.good() ){

    getline(in,theLine);

    if ( theLine == "" ) break;
    if ( theLine.compare(0,1,"#") == 0 ) continue ;

    istringstream s(theLine);
    int    iy = 0 ;
    double tmp ;
    while( s >> tmp ) {
      data(ix,iy)=tmp; 
      iy++;
    }
    if ( iy != ny ) {
      printf("read_datafileVS() : Error in file \"%s\" \n",filename);
      printf("                    -> at line %d while reading %s; number of columns should be %d\n", 
              ix+1, filename, ny); 
      exit(1) ;
    }
    ix++;
  }
  if ( ix != nx ) {
    printf("read_datafileVS() : Error while reading \"%s\" -> number of rows should be %d\n", filename,nx) ; 
    exit(1) ;
  }

  return ;

}

//  Read a datafile from filename in a vector form and store it in the int-array data
void read_datafileVS(Array2D<int>& data, const char* filename)
{

   ifstream in(filename);
  
  if(!in){
    printf("read_datafileVS()    : the requested file %s does not exist\n",filename) ;
    exit(1) ;
  }
  
  string theLine="";

  // figure out number of lines and columns
  int nx, ny, ix = 0 ;
  while(in.good()){
    getline(in,theLine);
    
    if ( theLine == "" ) break;
    if ( theLine.compare(0,1,"#") == 0 ) continue ;

    istringstream s(theLine);
    int iy = 0 ;
    int itmp   ;
    while( s >> itmp ) iy++ ;

    if ( ( ix > 0 ) && ( iy != ny ) )
    {
      printf("read_datafileVS() : Error at line %d !!!\n",ix+1) ;
      printf("                    no. of columns should be %d instead of %d\n",ny,iy) ;
      exit(1) ;
    }
    
    ny = iy ;

    ix++ ;

  }

  nx = ix ;

#ifdef VERBOSE
  printf("File \"%s\" contains %d rows and %d columns \n",filename,nx,ny) ;
#endif

  /* Resize, goto beginning, and read again */
  if ( ( (int) data.XSize() != nx ) || ( (int) data.YSize() != ny ) )
    data.Resize(nx,ny) ;

  //in.close() ;
  //in.open(filename);
  in.clear() ;
  in.seekg(0, ios::beg ) ;
  ix = 0 ;
  while( in.good() ){

    getline(in,theLine);

    if ( theLine == "" ) break;
    if ( theLine.compare(0,1,"#") == 0 ) continue ;

    istringstream s(theLine);
    int iy = 0 ;
    int itmp   ;
    while( s >> itmp ) {
      data(ix,iy)=itmp; 
      iy++;
    }
    if ( iy != ny ) {
      printf("read_datafileVS() : Error in file \"%s\" \n",filename);
      printf("                    -> at line %d while reading %s; number of columns should be %d\n", 
              ix+1, filename, ny); 
      exit(1) ;
    }
    ix++;
  }
  if ( ix != nx ) {
    printf("read_datafileVS() : Error while reading \"%s\" -> number of rows should be %d\n", filename,nx) ; 
    exit(1) ;
  }

  return ;

}

// Read a datafile from filename in a vector form and store it in the double-array data
void read_datafile_1d(Array1D<double>& data, const char* filename)
{
  int nx=data.XSize();

  if (nx==0){
    printf("read_datafile_1d()    : the requested data array is empty\n") ;
    exit(1) ;
  }

  int ny=1;

 ifstream in(filename);
  
  if(!in){
    printf("read_datafileVS()    : the requested file %s does not exist\n",filename) ;
    exit(1) ;
  }
  
  string theLine="";
  int ix=0;
  

  while(in.good()){
    getline(in,theLine);
    
    if (theLine=="") break;

    istringstream s(theLine);
    int iy=0;
    double tmp;
    //  while(s>>tmp){
    s>>tmp;
    data(ix)=tmp;
    iy++;
      // }
    if (s>>tmp) {printf("Error at line %d while reading %s; number of columns should be %d\n", ix+1, filename,ny); exit(1);}
    ix++;
  }
  if (ix!=nx) {printf("Error while reading %s; number of rows should be %d\n", filename,nx); exit(1);}

  return ;

}

// Read a datafile from filename in a vector form and store it in the int-array data
void read_datafile_1d(Array1D<int>& data, const char* filename)
{
  int nx=data.XSize();

  if (nx==0){
    printf("read_datafile_1d()    : the requested data array is empty\n") ;
    exit(1) ;
  }

  int ny=1;

  ifstream in(filename);
  
  if(!in){
    printf("read_datafileVS()    : the requested file %s does not exist\n",filename) ;
    exit(1) ;
  }
  
  string theLine="";
  int ix=0;
  
  
 while(in.good()){
    getline(in,theLine);
    
    if (theLine=="") break;

    istringstream s(theLine);
    int iy=0;
    int itmp;
    //  while(s>>itmp){
    s>>itmp;
    data(ix)=itmp;
    iy++;
      // }
    if (s>>itmp) {printf("Error at line %d while reading %s; number of columns should be %d\n", ix+1, filename,ny); exit(1);}
    ix++;
  }
  if (ix!=nx) {printf("Error while reading %s; number of rows should be %d\n", filename,nx); exit(1);}
  
  return;
}

//Read a datafile from filename in a vector form and store it in the double-array data
void read_datafileVS(Array1D<double>& data, const char* filename)
{
  data.Clear();

 ifstream in(filename);
  
  if(!in){
    printf("read_datafileVS()    : the requested file %s does not exist\n",filename) ;
    exit(1) ;
  }
  
  string theLine="";
  int ix=0;
  

  while(in.good()){
    getline(in,theLine);
    
    if (theLine=="") break;

    istringstream s(theLine);
    double tmp;
    //  while(s>>tmp){
    s>>tmp;
    data.PushBack(tmp);
      // }
    if (s>>tmp) {printf("Error at line %d while reading %s; number of columns should be 1\n", ix+1, filename); exit(1);}
    ix++;
  }

  return ;

}

// Read a datafile from filename in a vector form and store it in the int-array data
void read_datafileVS(Array1D<int>& data, const char* filename)
{
  data.Clear();
  
  ifstream in(filename);
  
  if(!in){
    printf("read_datafileVS()    : the requested file %s does not exist\n",filename) ;
    exit(1) ;
  }
  
  string theLine="";
  int ix=0;
  
  
 while(in.good()){
    getline(in,theLine);
    
    if (theLine=="") break;

    istringstream s(theLine);
    int itmp;
    //  while(s>>itmp){
    s>>itmp;
    data.PushBack(itmp);
      // }
    if (s>>itmp) {printf("Error at line %d while reading %s; number of columns should be 1\n", ix+1, filename); exit(1);}
    ix++;
  }
  
  return;
}

// Write the contents of a 2d double array data to a file filename in a matrix form
void write_datafile(const Array2D<double>& data, const char* filename)
{
  int nx=data.XSize();
  int ny=data.YSize();

  
  FILE* f_out;
  if(!(f_out = fopen(filename,"w"))){ 
    printf("write_datafile: could not open file '%s'\n",filename); 
    exit(1); 
  }

 for(int ix = 0 ; ix < nx ; ix++){
    for(int iy = 0 ; iy < ny ; iy++){
      fprintf(f_out, "%24.16lg ", data(ix,iy));
    }
    fprintf(f_out, "\n");
 }

 if(fclose(f_out)){ 
   printf("write_datafile: could not close file '%s'\n",filename); 
   exit(1); 
 }

#ifdef VERBOSE
 printf("Data written to '%s' in a matrix form [%d X %d]\n", filename,nx,ny);
#endif

 return ;

}

// Write the contents of a 2d int array data to a file filename in a matrix form 
void write_datafile(const Array2D<int>& data, const char* filename)
{
  int nx=data.XSize();
  int ny=data.YSize();

  
  FILE* f_out;
  if(!(f_out = fopen(filename,"w"))){ 
    printf("write_datafile: could not open file '%s'\n",filename); 
    exit(1); 
  }

 for(int ix = 0 ; ix < nx ; ix++){
    for(int iy = 0 ; iy < ny ; iy++){
      fprintf(f_out, "%d ", data(ix,iy));
    }
    fprintf(f_out, "\n");
 }

 if(fclose(f_out)){ 
   printf("write_datafile: could not close file '%s'\n",filename); 
   exit(1); 
 }

#ifdef VERBOSE
 printf("Data written to '%s' in a matrix form [%d X %d]\n", filename,nx,ny);
#endif

 return;
}

// Write the contents of a 1d double array data to a file filename in a vector form 
void write_datafile_1d(const Array1D<double>& data, const char* filename)
{

  int nx=data.XSize();

  FILE* f_out;
  if(!(f_out = fopen(filename,"w"))){ 
    printf("write_datafile: could not open file '%s'\n",filename); 
    exit(1); 
  }

  for(int ix = 0 ; ix < nx ; ix++) fprintf(f_out, "%24.16lg\n", data(ix));

  if(fclose(f_out)){ 
    printf("write_datafile: could not close file '%s'\n",filename); 
    exit(1); 
  }

#ifdef VERBOSE
  printf("Data written to '%s' in a matrix form [%d X 1]\n", filename,nx);
#endif

 return;
}

// Write the contents of a 1d int array data to a file filename in a vector form 
void write_datafile_1d(const Array1D<int>& data, const char* filename)
{
  int nx=data.XSize();

  FILE* f_out;
  if(!(f_out = fopen(filename,"w"))){ 
    printf("write_datafile: could not open file '%s'\n",filename); 
    exit(1); 
  }

 for(int ix = 0 ; ix < nx ; ix++) fprintf(f_out, "%d\n", data(ix));

 if(fclose(f_out)){ 
   printf("write_datafile: could not close file '%s'\n",filename); 
   exit(1); 
 }

#ifdef VERBOSE
 printf("Data written to '%s' in a matrix form [%d X 1]\n", filename,nx);
#endif 

 return;
}


