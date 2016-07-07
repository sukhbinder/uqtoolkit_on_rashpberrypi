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
/* gen_mi.cpp */

#include "uqtktools.h"
#include <unistd.h>

using namespace std;


/// default multiindex type
#define MI_TYPE "TO" 
/// default first param
#define ORD 1
/// default second param		
#define DIM 3 
/// default parameter filename
#define PARAM_FILE "mi_param.dat" 


/******************************************************************************/
/// \brief Displays information about this program
int usage(){
  printf("This program to generate multiindex files given rules.\n");
  printf("usage: gen_mi [-h]  [-x<mi_type>] [-p<nord>] [-q<ndim>] [-f<param_file>]\n");
  printf(" -h               : print out this help message \n");
  printf(" -x <mi_type> : define the multiindex type (default=%s) \n",MI_TYPE);
  printf(" -p <nord>        : define the first parameter (default=%d) \n",ORD);
  printf(" -q <ndim>        : define the second parameter (default=%d) \n",DIM);
  printf(" -f <param_file>     : define the parameter filename for multiindex (default=%s) \n",PARAM_FILE);
  printf("================================================================================\n");
  printf("Input  : None \n");
  printf("Output : File 'mindex.dat'\n");
  printf("--------------------------------------------------------------------------------\n");
  printf("================================================================================\n");
  exit(0);
  return 0;
}


/// Main program: Generates multiindex of requested type with given parameters

int main(int argc, char *argv[])
{

  
  /// Set the default values
  int nord = ORD;
  int ndim=DIM;
  char* param_file=(char *)PARAM_FILE;
  char* mi_type=(char *)(char *)(char *)(char *)(char *)(char *)(char *)(char *)(char *)MI_TYPE;
    
  bool pflag = false;
  bool qflag = false;
  bool fflag = false;

  /// Read the user input
  int c;

  while ((c=getopt(argc,(char **)argv,"hx:p:q:f:"))!=-1){
     switch (c) {
     case 'h':
       usage();
       break;
     case 'x':
       mi_type =  optarg;
       break;
     case 'p':
       nord =  strtol(optarg, (char **)NULL,0);	
       pflag=true;
       break;
     case 'q':
       ndim =  strtol(optarg, (char **)NULL,0);
       qflag=true;
       break;
     case 'f':
       param_file =  optarg;
       fflag=true;
       break;
     default :
       break;
     }
  }

/*----------------------------------------------------------------------------*/ 
  /// Print the input information on screen 
  fprintf(stdout,"mi_type = %s \n",mi_type);
fprintf(stdout,"ndim = %d \n",ndim);
  if (pflag)  fprintf(stdout,"nord = %d \n",nord);
  if (fflag)  fprintf(stdout,"param_file = %s \n",param_file);

    /*----------------------------------------------------------------------------*/

   if(fflag && (pflag)){
     printf("gen_mi(): Can not specify both parameter file and order. Exiting.\n");
     exit(1);
   }


  
   /// Cast multiindex type as a string
   string mi_type_str(mi_type);

   int npc;
   Array2D<int> mindex;

   /// Total order
   if (mi_type_str=="TO"){
     npc=computeMultiIndex(ndim,nord,mindex);
   }

   /// HDMR ordering
   else if(mi_type_str=="HDMR"){
     Array1D<int> maxorders;
     Array2D<int> maxorders2d;
     read_datafileVS(maxorders2d,param_file);
     getCol(maxorders2d, 0, maxorders);
     npc=computeMultiIndexHDMR(ndim, maxorders, mindex);
   }
   
   else {
     printf("gen_mi():: Multiindex type %s is not recognized. \n", mi_type);
     exit(1);
   }
   
   
   write_datafile(mindex, "mindex.dat");
   cout << "Generated multiindex of size " << npc << " and stored in mindex.dat" << endl;
   
   return 0;
}


