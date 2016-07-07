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
/* pce_resp.cpp */

#include "PCSet.h"
#include "uqtktools.h"
#include <unistd.h>

using namespace std;


/// default data dimensionality
#define DIM 2	
/// default PC type
#define CHAOS "LEG" 
/// default PC order
#define ORD 3		
/// default quadrature type
#define QUADTYPE  "CC"   
/// default quadrature sparseness
#define FSTYPE "sparse" 
/// default quadrature parameter   
#define PARAM 3	
///default alpha parameter for PC	
#define ALPHA 0.0
/// default beta parameter for PC
#define BETA 1.0 



/******************************************************************************/
/// Displays information about this program
int usage(){
  printf("usage: pce_resp [-h] [-e] [-d<ndim>] [-x<which_chaos>] [-o<nord>] [-f<fsType>] [-p<param>] [-a<alpha>] [-b<beta>]\n");
  printf(" -h               : print out this help message \n");
  printf(" -e               : flag if we want to evaluate PC representation at the quadrature points and get the error \n");
  printf(" -x <which_chaos> : define the PC type; only LEG is implemented (default=%s) \n",CHAOS);
  printf(" -d <ndim>        : define the data dimensionality (default=%d) \n",DIM);
  printf(" -o <nord>        : define the PC order (default=%d) \n",ORD);
  printf(" -f <fsType>      : define the quadrature sparseness  (default=%s) \n",FSTYPE);
  printf(" -p <param>       : define the  (default=%d) \n",PARAM);
  printf(" -a <alpha>       : define the alpha parameter of the quadrature (default=%lg) \n",ALPHA);
  printf(" -b <beta>        : define the beta parameter of the quadrature (default=%lg) \n",BETA);
  printf("================================================================================\n");
  printf("Input  : qdpts.dat  - quadrature points on [-1,1]^d\n");
  printf("       : wghts.dat  - quadrature weights\n");
  printf("       : ydata.dat  - function evaluations corresponding to quadrature points\n");
  printf("Output : PCcoeff_quad.dat - PC coefficients\n");
  printf("       : ydata_pc.dat     - PC response at the quadrature points\n");
  printf("--------------------------------------------------------------------------------\n");
  printf("Todo   : Add an option to use the default quadrature for a given PC             \n");
  printf("================================================================================\n");
  exit(0);
  return 0;
}


///  Main program: gets PC coefficients of a response curve given function evaluations at a grid
int main (int argc, char *argv[]) {

  /// Set the defaults and parse the input arguments
  int ndim= DIM;
  int nord= ORD;
  char* which_chaos=(char *)CHAOS;
  double alpha=ALPHA;
  double beta =BETA;

  bool eflag=false;
  bool aflag=false;
  bool bflag=false;

  int c;

  while ((c=getopt(argc,(char **)argv,"hex:d:o:a:b:"))!=-1){
    switch (c) {
    case 'h':
      usage();
      break;
    case 'e':
      eflag=true;
      break;
    case 'x':
      which_chaos =  optarg;
      break;
    case 'd':
      ndim =  strtol(optarg, (char **)NULL,0);	
      break;
    case 'o':
      nord =  strtol(optarg, (char **)NULL,0);	
      break;
  case 'a':
      aflag=true;
      alpha = strtod(optarg, (char **)NULL);	
      break;
    case 'b':
      bflag=true;
      beta = strtod(optarg, (char **)NULL);	
      break;
    default :
      break;
    }
  }
  
  /// Print the input information
  fprintf(stdout,"---------------------------------\n") ;
  fprintf(stdout,"pce_resp() parameters : \n") ;
  fprintf(stdout,"which_chaos = %s \n",which_chaos);  
  fprintf(stdout,"ndim = %d \n",ndim);
  fprintf(stdout,"nord = %d \n",nord);
  if ( aflag )
    fprintf(stdout,"alpha     = %lg \n",alpha);
  if ( bflag )
    fprintf(stdout,"beta     = %lg \n",beta);
  if ( eflag )
    fprintf(stdout,"will compute the L_2 error, too.\n");
  fprintf(stdout,"---------------------------------\n") ;
 

  /// Read the ydata evaluated at quadrature points
  Array2D<double> ydata_2d;
  read_datafileVS(ydata_2d,"ydata.dat");
  Array1D<double> ydata;
  for(int i=0;i<(int) ydata_2d.XSize();i++)
    ydata.PushBack(ydata_2d(i,0));
  
  /// Read the quadrature rule
  Array2D<double> qdpts;
  Array1D<double> wghts;
  Array2D<int> indices;
  
  read_datafileVS(qdpts,"qdpts.dat");
  wghts.Resize(qdpts.XSize());
  indices.Resize(qdpts.XSize(),qdpts.YSize());
  read_datafile_1d(wghts,"wghts.dat");
  // read_datafile_1d(indices,"indices.dat");

  /// Set the chaos
  string which_chaos_str(which_chaos);
  /// Declare PC model with no quadrature
  PCSet currPCModel("NISPnoq",nord,ndim,which_chaos_str,alpha,beta);
  
 
  /// Set the quadrature rule
  Quad spRule;
  spRule.SetRule(qdpts,wghts,indices);
  currPCModel.SetQuadRule(spRule);
  int nup=currPCModel.GetNumberPCTerms()-1;
  int totquad=currPCModel.GetNQuadPoints();

  /// Get the multiindex for postprocessing
  Array2D<int> mindex;
  currPCModel.GetMultiIndex(mindex);
  write_datafile(mindex,"mindex.dat");

  /// Get PC coefficients by Galerkin Projection (c_k=<ydata(iq) * psi_k(iq)>)
  Array1D<double> c_k(nup+1,0.e0); 
  currPCModel.GalerkProjection(ydata,c_k);
  write_datafile_1d(c_k,"PCcoeff_quad.dat");

  /// If requested, compute the L2 error at quadrature points
  /// \note Note that the error is computed with the same quadrature rule, and may be inaccurate
  /// \todo Perhaps only compute simpler, l2 norm of the error without involving weights
  if (eflag){
    // Evaluate the PC expansion at quadrature points
    Array1D<double> pcxvalues(totquad,0.e0);
    currPCModel.EvalPCAtCustPoints(pcxvalues,qdpts,c_k);
  // Write-out
    write_datafile_1d(pcxvalues,"ydata_pc.dat");
    
    double sum1=0.0, sum2=0.0;
    for(int it=0;it<totquad;it++){
      sum1+=(  wghts(it)*( pcxvalues(it)-ydata(it) )*( pcxvalues(it)-ydata(it) ) );
      sum2+= ( wghts(it)* ydata(it)*ydata(it) );
    }
    printf("pce_resp() : Relative L2 error^2=%lg/%lg=%lg\n",sum1,sum2, sum1/sum2);

  } 
  

  


  return 0;
}

