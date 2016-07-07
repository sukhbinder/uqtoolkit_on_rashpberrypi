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
#include <cstdio>
#include <stddef.h>
#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include "assert.h"

#include <getopt.h>
#include "Array1D.h"
#include "Array2D.h"

#include "kldecompuni.h"
#include "uqtktools.h"
#include "lapack.h"
#include "error_handlers.h"

using namespace std;

#include "kl_sample.h"

#define NPTS  129
#define NKL   64
#define NSPL  128
#define CLEN  0.05
#define SIG   5.0

#define COVTYPE "SqExp"
#define XFILE   "xgrid.dat"


/// \brief Displays information about this program
int usage(){

  printf("usage: cor_kl [-h]  [-x<grid_file>] [-c<cov_type>] [-n<npts>] [-e<nkl>] [-p<nspl>] [-s<sigma>] [-l<cor_len>]\n");
  printf(" -h                 : print out this help message \n");
  printf(" -x <xfile>       : define the file from which the time grid is being read (default=%s) \n",XFILE);
  printf(" -c <cov_type>    : define wether using analytical covariance (type: %s) or compute it from samples \n",COVTYPE);
  printf(" -n <npts>        : define the number of grid points (default=%d) \n",NPTS);
  printf(" -e <nkl>         : define the number of KL modes retained (default=%d) \n",NKL);
  printf(" -p <nspl>        : define the number of samples (default=%d) \n",NSPL);
  printf(" -s <sigma>       : define standard deviation (default=%e) \n",SIG);
  printf(" -l <clen>        : define correlation length (default=%e) \n",CLEN);
  printf("================================================================================================================\n");
  printf("Input:: Nothing hard-coded.\n");
  printf("Output:: \n");
  printf("        - cov_out.dat:  covariance matrix\n");
  printf("        - eig.dat:      eigenvalues\n");
  printf("        - KLmodes.dat:  eigenmodes scaled with sqrt(eig)\n");
  printf("        - rel_diag.dat: pointwise covariance fraction explained by the finite KL expansion\n");
  printf("        - mean.dat:     pointwise mean based on samples\n");
  printf("        - xi_data.dat:  samples of eigenvalues\n");
  printf("================================================================================================================\n");
  exit(0);
  return 0;
}

/*!
   \brief Karhunen-Loeve decomposition of a UNIVARIATE process
   given covariance type or samples drawn from multivariate normal distribution
*/
int main(int argc, char *argv[])
{

  /* Set the default values */
  int    npts = NPTS ;
  int    nkl  = NKL  ;
  int    nspl = NSPL ;

  bool  xflag = false;    
  char* xfile = (char *) XFILE;
  Array1D<double> xgrid;

  bool   cflag = false ;
  double clen  = CLEN  ;
  double sigma = SIG   ;
  char* cov_type = (char *) COVTYPE;
  Array2D<double> cov ;

  /* Read the user input */
  int c;

  while ((c=getopt(argc,(char **)argv,"hx:c:n:e:p:s:l:"))!=-1){
    switch (c) {
    case 'h':
      usage();
      break;
    case 'x':
      xflag = true;
      xfile = optarg;
      break;
    case 'c':
      cflag=true;
      cov_type = optarg;
      break;
    case 'n':
      npts = strtol(optarg, (char **)NULL,0);	
      break;
    case 'e':
      nkl = strtol(optarg, (char **)NULL,0);	
      break;
    case 'p':
      nspl = strtol(optarg, (char **)NULL,0);  
      break;
    case 's':
      sigma = strtod(optarg, (char **)NULL);	
      break;
   case 'l':
      clen = strtod(optarg, (char **)NULL); 
      break;
    }
  }

  if ( nkl > npts )
    throw Tantrum("kl_sample.cpp::main(): cannot request more KL modes than number of grid points");

  /* Print the input information on screen */
  cout << " - Number of grid points: " << npts << endl<<flush;
  cout << " - Correlation length:    " << clen << endl<<flush;
  cout << " - Standard deviation:    " << sigma<< endl<<flush;
  cout << " - Number of KL modes:    " << nkl  << endl<<flush;
  if ( cflag )
    cout<<" - Will generate covariance of type "<<cov_type<<endl<<flush;
  else    
    cout << " - Will generate covariance from "<<nspl<<" samples" << endl<<flush;
  if ( xflag )
    cout<<" - Will read grid from file: "<<xfile<<endl<<flush;
  else    
    cout << " - Will create equally spaced grid with "<<npts<<" points in [0,1]" << endl<<flush;
 
  /* read grid from file or generate uniform grid on [0,1] */
  if ( xflag )
    /* read grid from file */
    read_datafile_1d(xgrid,xfile);
  else {
    genGrid(xgrid,npts);
    write_datafile_1d(xgrid,xfile);
  }

  Array2D<double> ySamples;

  cov.Resize(npts,npts,0.e0);	
  if ( cflag ) {
    for ( int i = 0; i < npts; i++)
      for ( int j = 0; j < npts; j++)
        cov(i,j)=genCovAnl(xgrid(i),xgrid(j),clen,sigma,cov_type);
  }
  else {

    for ( int i = 0; i < npts; i++)
      for ( int j = 0; j < npts; j++)
        cov(i,j)=genCovAnl(xgrid(i),xgrid(j),clen,sigma,string("SqExp"));
    for ( int i = 0; i < npts; i++)
      cov(i,i) += 1.e-13;

    /* Generate samples */
    char *lu = (char *) "L";
    int info ;
    FTN_NAME(dpotrf)( lu, &npts, cov.GetArrayPointer(), &npts, &info );

    /* Check the success in Cholesky factorization */
    if ( info != 0 ) {

      cout<<"Error in Cholesky factorization, info=" << info << endl << flush ;;
      exit(1);

    } /* done if Cholesky factorization fails */

    dsfmt_t  rnstate ;
    int rseed=20120828;
    dsfmt_init_gen_rand(&rnstate, (uint32_t) rseed );

    Array1D<double> randSamples(npts,0.0);
    ySamples.Resize(npts,nspl,0.0);
    for ( int j = 0; j < nspl; j++) {
      for (int i = 0 ; i < npts ; i++ )
        randSamples(i) = dsfmt_genrand_nrv(&rnstate);
      for ( int i = 0; i < npts; i++ ) {
        ySamples(i,j)=0.0;
        for ( int k = 0; k < i+1; k++) 
          ySamples(i,j) += (cov.GetArrayPointer())[i+k*npts]*randSamples(k);  
      } 
    }
    write_datafile( ySamples, "samples.dat" );

    /* Compute samples mean */
    Array1D<double> mean(npts,0.e0);
    for ( int i = 0 ; i < npts ; i++ )
      mean(i) = getMean( ySamples, string("L"), i );
    write_datafile_1d( mean, "mean.dat" );

    /* Compute covariance matrix */

    /* 1. subtract the mean from samples */
    for ( int j = 0 ; j < nspl ; j++)
      for ( int i = 0 ; i < npts ; i++)
        ySamples(i,j) -= mean(i) ;

    /* 2. compute the upper triangular part */
    for ( int i = 0; i < npts; i++ ) {
      for ( int j = i; j < npts; j++ ) {

        double dsum=0.0;
        for(int k = 0; k < nspl; k++ )
          dsum += ySamples(i,k)*ySamples(j,k);
    
        cov(i,j) = dsum/( (double) nspl );
      }
    }

    /* 3. transpose to fill out the lower triangle */
    for ( int i = 0; i < npts; i++ ) 
      for ( int j = 0; j < i; j++ ) 
        cov(i,j) = cov(j,i) ;

  }
  write_datafile( cov, "cov.dat" );
 
  /*  Performing KL decomposition */
  cout << " --> Starting KL decomposition " << endl;

  KLDecompUni decomposer(xgrid);  
  int n_eig = decomposer.decompose(cov,nkl);

  if(n_eig <  nkl){ 
    printf("There are only %d  eigenvalues available (requested %d) \n",n_eig, nkl); 
    nkl = n_eig; 
  }

  const Array1D<double>& eigs    = decomposer.eigenvalues();
  const Array2D<double>& KLmodes = decomposer.KLmodes();

  cout << " --> KL decomposition done " << endl << flush;

  cout << "      - Obtained " << n_eig << " eigenvalues:" << endl;
  //for ( int i = 0; i < n_eig; i++) 
  //  cout << "        " << i << " : " << eigs(i) << endl;
  write_datafile_1d(eigs,"eig.dat");

  cout << "      - Computing relative variances" << endl;
  Array1D<double> rel (npts,0.0);
  for ( int i = 0; i < npts; i++){
    for ( int j = 0; j < n_eig; j++){
      rel(i) += eigs(j)*KLmodes(i,j)*KLmodes(i,j);
    }
    rel(i) /= cov(i,i);
  }
  write_datafile_1d(rel,"relVar.dat");

  Array2D<double> scaledKLmodes(npts,n_eig+1,0.0);
  for ( int i = 0; i < npts; i++ ){
    scaledKLmodes(i,0) = xgrid(i);
    for ( int j = 0; j < n_eig; j++ )
      scaledKLmodes(i,j+1) = KLmodes(i,j)*sqrt(eigs(j));
  }
  write_datafile(scaledKLmodes,"KLmodes.dat");

  if ( !cflag) {

    /* Project realizations onto KL modes */
    cout << "      - Project realizations onto KL modes " << endl;
    Array2D<double> xi(nkl, nspl, 0.e0);
    decomposer.KLproject( ySamples, xi );

    Array2D<double> xit(nspl,nkl,0.e0);
    transpose(xi,xit);
    write_datafile(xit,"xi_data.dat");

  }

  cout << " --> KL example done " << endl;

  return ( 0 ) ;

}


double genCovAnl(const double x, const double y, const double clen, const double sigma, const string covtype) {

  double cov;

  if ( covtype == "SqExp" )
    cov = exp( - (x-y) * (x-y) / ( clen * clen ) ) * sigma * sigma ;
  else if ( covtype == "Exp" )
    cov = exp( - fabs(x-y) / clen ) * sigma * sigma ;
  else 
    throw Tantrum("kl_sample.cpp::genCovAnl(): covariance type is not recognized!");

  return ( cov );

}

void genGrid(Array1D<double> &xgrid, const int npts) {

  if (npts<=0)
    throw Tantrum("kl_sample::genGrid() : number of grid points needs to be greater than 0") ;

  xgrid.Resize(npts,0.0);

  for ( int i=0; i<npts; i++ ) xgrid(i) = (double) i / ((double) npts - 1.0);

  for ( int i=0; i<npts; i++ ) xgrid(i) = scomp(xgrid(i), 1.01, 1.0);
  
  return;

}

double getMean(Array2D<double> &y, const string idir, const int ij) {

  int nx = (int) y.XSize();
  int ny = (int) y.YSize();

  double dmean, dtmp = 0.0 ;
 
  if ( idir == "L" ) {

    int mst=ny%6;

    for ( int j = 0; j<mst; j++ ) 
      dtmp += y(ij,j);      

    if ( ny >= 6 ) {
      for ( int j = mst; j < ny; j += 6)
        dtmp += (y(ij,j)+y(ij,j+1)+y(ij,j+2)+y(ij,j+3)+y(ij,j+4)+y(ij,j+5));
    }

    dmean = dtmp / ((double) ny );

  } 
  else if ( idir == "C" ) {

    int mst=nx%6;

    for ( int i = 0; i<mst; i++ ) 
      dtmp += y(i,ij);      

    if ( nx >= 6 ) {
      for ( int i = mst; i < nx; i += 6)
        dtmp += (y(i,ij)+y(i+1,ij)+y(i+2,ij)+y(i+3,ij)+y(i+4,ij)+y(i+5,ij));
    }

    dmean = dtmp / ( (double) nx );
  }

  return ( dmean ) ;

}

double dcomp(double xi, double a, double b, double L) {
  double dtmp = pow((b+1.0)/(b-1.0),(xi-a)/(1.0-a));
  return ( 
           L * ((2.0*a+b)*dtmp+2.0*a-b)
	      /((2.0*a+1.0)*(1.0+dtmp)) 
         );

}

double scomp(double xi, double b, double L) {
  double dtmp = pow((b+1.0)/(b-1.0),1.0-xi);
  return ( 
	  L * (b+1.0-(b-1.0)*dtmp)
	     /(dtmp+1.0) 
         );

}
