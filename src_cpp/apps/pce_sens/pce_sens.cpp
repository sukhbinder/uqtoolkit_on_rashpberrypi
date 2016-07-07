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
/* pce_sens.cpp */

#include "PCSet.h"
#include "uqtktools.h"
#include <unistd.h>

using namespace std;



/// default PC type
#define CHAOS "LU"
/// default multiindex file
#define MINDEX_FILE "mindex.dat" 
/// default coefficient file
#define COEF_FILE "PCcoeff.dat"



/// Displays information about this program
int usage(){
  printf("This program parses the information contained in given pair of multiindex-coefficients\n");
  printf("usage: pce_sens [-h] [-m<mindex_file>] [-f<coef_file>] [-x<which_chaos>]\n");
  printf(" -h               : print out this help message \n");
  printf(" -x <which_chaos> : define the PC type (default=%s) \n",CHAOS);
  printf(" -m <mindex_file> : define multiindex filename (default=%s) \n",MINDEX_FILE);
  printf(" -f <coef_file>   : define the coefficient filename (default=%s) \n",COEF_FILE);
  printf("================================================================================\n");
  printf("Input  : None \n");
  printf("Output : sp_mindex.X.dat - sparse format of multiindices, has 2*X columns, \n");
  printf("                         - e.g. for X=3, each row has a form [i j k a_i a_j a_k] \n");
  printf("                         - and corresponds to a PC term with order a_i in the i-th dimension etc. \n");
  printf("       : varfrac.dat     - a column file of variance fractions corresponding to each PC term \n");
  printf("                         - the size is Nx1, where N is the number of PC terms \n");
  printf("       : mainsens.dat    - a columm file of main sensitivities. The size is dx1, where d is the dimensionality\n");
  printf("       : totsens.dat     - a columm file of total sensitivities. The size is dx1\n");
  printf("       : jointsens.dat   - a matrix file of joint sensitivities. The size is dxd\n");
  printf("--------------------------------------------------------------------------------\n");
  printf("Comments  : Sparse format is useful in high-dimensional problems \n");
  printf("Complexity: Linear in the number of PC terms \n");
  printf("================================================================================\n");
  
  exit(0);
  return 0;
}


/// Main program: parses the information contained in given multiindices and corresponding coefficients
int main(int argc, char *argv[])
{

  /// Set the defaults and parse the input arguments
  int c;
  
  char* mindex_file=(char *)MINDEX_FILE;
  char* coef_file=(char *)COEF_FILE;
  char* which_chaos=(char *)CHAOS;
  
  while ((c=getopt(argc,(char **)argv,"hm:f:x:"))!=-1){
     switch (c) {
     case 'h':
       usage();
       break;
     case 'm':
       mindex_file =  optarg;	
       break;
     case 'f':
       coef_file =  optarg;	
       break;
     case 'x':
       which_chaos =  optarg;
       break;
     default :
       break;
     }
  }

  /// Print out input information
  fprintf(stdout,"---------------------------------\n") ;
  fprintf(stdout,"pce_sens() parameters : \n") ;
  fprintf(stdout,"mindex_file = %s \n",mindex_file);
  fprintf(stdout,"coef_file = %s \n",coef_file);
  fprintf(stdout,"which_chaos = %s \n",which_chaos);
  fprintf(stdout,"---------------------------------\n") ;


  /// Read the multiindex and coefficients' files
  Array2D<int> mindex; 
  read_datafileVS(mindex,mindex_file);
  int npc=mindex.XSize();
  int ndim=mindex.YSize();
  Array1D<double> coef(npc,0.e0);
  read_datafile_1d(coef,coef_file); 
  

  /// Declare PC in NISP formulation with no quadrature
  string which_chaos_str(which_chaos);
  PCSet PCModel("NISPnoq",mindex,which_chaos_str,0.0,1.0);

  /// Encode the multiindex in a sparse format and print to files
  Array1D< Array2D<int> > sp_mindex;
  PCModel.EncodeMindex(sp_mindex);    
    
  char filename[25];
  for (int i=1;i<(int)sp_mindex.XSize();i++){
    int nn=sprintf(filename,"sp_mindex.%d.dat",i);
    write_datafile(sp_mindex(i),filename);
  }

  
  /// Compute mean and variance of PC and variance fractions for each term
  Array1D<double> varfrac;
  double mean, var;
  
  mean=PCModel.ComputeMean(coef);
  var=PCModel.ComputeVarFrac(coef,varfrac);
  cout << "Mean = " << mean << endl;
  cout << "Var  = " << var << endl;
  write_datafile_1d(varfrac,"varfrac.dat");
  
  /// Compute main sensitivities
  Array1D<double> mainsens;
  PCModel.ComputeMainSens(coef,mainsens);
  write_datafile_1d(mainsens,"mainsens.dat");

  /// Compute total sensitivities
  Array1D<double> totsens;
  PCModel.ComputeTotSens(coef,totsens);
  write_datafile_1d(totsens,"totsens.dat");

  /// Compute joint sensitivities
  Array2D<double> jointsens;
  PCModel.ComputeJointSens(coef,jointsens);
  for(int id=0;id<ndim;id++)
    jointsens(id,id)=mainsens(id);
  write_datafile(jointsens,"jointsens.dat");
  
  return 0;
}


