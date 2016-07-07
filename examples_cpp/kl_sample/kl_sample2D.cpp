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

#include "kl_sample2D.h"

#define NX    65
#define NY    65
#define NKL   64
#define NSPL  128
#define CLEN  0.05
#define SIG   5.0

#define COVTYPE "SqExp"
#define XFILE   "xgrid.dat"
#define YFILE   "ygrid.dat"


/// \brief Displays information about this program
int usage(){

  printf("usage: cor_kl [-h]  [-x<grid_file>] [-c<cov_type>] [-n<nx>] [-m<ny>] [-e<nkl>] [-p<nspl>] [-s<sigma>] [-l<cor_len>]\n");
  printf(" -h                 : print out this help message \n");
  printf(" -x <xfile>       : define the file from which the time grid is being read (default=%s) \n",XFILE);
  printf(" -c <cov_type>    : define wether using analytical covariance (type: %s) or compute it from samples \n",COVTYPE);
  printf(" -n <nx>          : define the number of x-grid points (default=%d) \n",NX);
  printf(" -m <ny>          : define the number of y-grid points (default=%d) \n",NY);
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
  int    nx   = NX   ;
  int    ny   = NY   ;
  int    nkl  = NKL  ;
  int    nspl = NSPL ;

  bool  xflag = false;    
  char* xfile = (char *) XFILE;
  char* yfile = (char *) YFILE;
  Array1D<double> xgrid, ygrid;

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
      break;
    case 'c':
      cflag=true;
      cov_type = optarg;
      break;
    case 'n':
      nx   = strtol(optarg, (char **)NULL,0);	
      break;
    case 'm':
      ny   = strtol(optarg, (char **)NULL,0);	
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

  int nxy=nx*ny;
  if ( nkl > nxy )
    throw Tantrum("kl_sample.cpp::main(): cannot request more KL modes than number of grid points");

  /* Print the input information on screen */
  cout << " - Number of grid points: " << nx <<"x"<<ny << endl<<flush;
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
    cout << " - Will create equally spaced grid with "<<nx<<"x"<<ny<<" points in [0,1]^2" << endl<<flush;
 
  /* read grid from file or generate uniform grid on [0,1] */
  if ( xflag ) {
    /* read grid from file */
    read_datafile_1d(xgrid,xfile);
    read_datafile_1d(ygrid,yfile);
  }
  else {
    genGrid(xgrid,ygrid,nx,ny);
    write_datafile_1d(xgrid,xfile);
    write_datafile_1d(ygrid,yfile);
  }
  Array2D<double> xygrid(nxy,2);
  for (int k=0; k < nxy; k++ ) {
    int i = k%nx;
    int j = k/nx;
    xygrid(k,0) = xgrid(i);
    xygrid(k,1) = ygrid(j);
  }

  Array2D<double> ySamples;

  cov.Resize(nxy,nxy,0.e0);	
  if ( cflag ) {
    for ( int i = 0; i < nxy; i++)
      for ( int j = 0; j < nxy; j++)
        cov(i,j)=genCovAnl(xygrid,i,j,clen,sigma,cov_type);
  }
  // else {
  //   read_datafile(cov,"cov.dat");
  // }
  else {

    double dfac=1.0e-11;
    bool   tryagain = true;

    while ( ( tryagain ) && ( dfac < 1.e-6 ) ) {
      tryagain = false ;
      for ( int i = 0; i < nxy; i++)
        for ( int j = 0; j < nxy; j++)
          cov(i,j)=genCovAnl(xygrid,i,j,clen,sigma,string("SqExp"));
      for ( int i = 0; i < nxy; i++)
        cov(i,i) += dfac;
      //write_datafile( cov, "cov.dat" );

      /* Generate samples */
      char *lu = (char *) "L";
      int info ;
      FTN_NAME(dpotrf)( lu, &nxy, cov.GetArrayPointer(), &nxy, &info );

      /* Check the success in Cholesky factorization */
      if ( info != 0 ) {

        cout<<"Error in Cholesky factorization, info=" << info << endl << flush ;;
        dfac = dfac * 10.0 ;
        tryagain = true ;
        cout<<"  will try again with diagonal factor" << dfac << endl << flush ;;

      } /* done if Cholesky factorization fails */
    }
    dsfmt_t  rnstate ;
    int rseed=20120828;
    dsfmt_init_gen_rand(&rnstate, (uint32_t) rseed );

    Array1D<double> randSamples(nxy,0.0);
    ySamples.Resize(nxy,nspl,0.0);
    for ( int j = 0; j < nspl; j++) {
      for (int i = 0 ; i < nxy ; i++ )
        randSamples(i) = dsfmt_genrand_nrv(&rnstate);
      for ( int i = 0; i < nxy; i++ ) {
        ySamples(i,j)=0.0;
        for ( int k = 0; k < i+1; k++) 
          ySamples(i,j) += (cov.GetArrayPointer())[i+k*nxy]*randSamples(k);  
      } 
    }
    write_datafile( ySamples, "samples.dat" );

    /* Compute samples mean */
    Array1D<double> mean(nxy,0.e0);
    for ( int i = 0 ; i < nxy ; i++ )
      mean(i) = getMean( ySamples, string("L"), i );
    write_datafile_1d( mean, "mean.dat" );

    /* Compute covariance matrix */

    /* 1. subtract the mean from samples */
    for ( int j = 0 ; j < nspl ; j++)
      for ( int i = 0 ; i < nxy ; i++)
        ySamples(i,j) -= mean(i) ;

    /* 2. compute the upper triangular part */
    for ( int i = 0; i < nxy; i++ ) {
      for ( int j = i; j < nxy; j++ ) {

        double dsum=0.0;
        for(int k = 0; k < nspl; k++ )
          dsum += ySamples(i,k)*ySamples(j,k);
    
        cov(i,j) = dsum/( (double) nspl );
      }
    }

    /* 3. transpose to fill out the lower triangle */
    for ( int i = 0; i < nxy; i++ ) 
      for ( int j = 0; j < i; j++ ) 
        cov(i,j) = cov(j,i) ;

  }
  write_datafile( cov, "cov.dat" );
 
  /*  Performing KL decomposition */
  cout << " --> Starting KL decomposition " << endl;

  /* Create 1D equivalent grid*/
  Array1D<double> xg1d(nxy,0.0);
  getGrid1dEquiv(xgrid,ygrid,xg1d);

  KLDecompUni decomposer(xg1d);  
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
  Array1D<double> rel (nxy,0.0);
  for ( int i = 0; i < nxy; i++){
    for ( int j = 0; j < n_eig; j++){
      rel(i) += eigs(j)*KLmodes(i,j)*KLmodes(i,j);
    }
    rel(i) /= cov(i,i);
  }
  write_datafile_1d(rel,"relVar.dat");

  Array2D<double> scaledKLmodes(nxy,n_eig+2,0.0);
  for ( int k = 0; k < nxy; k++ ){
    int i = k%nx;
    int j = k/nx;
    scaledKLmodes(k,0) = xgrid(i);
    scaledKLmodes(k,1) = ygrid(j);
    for ( int l = 0; l < n_eig; l++ )
      scaledKLmodes(k,l+2) = KLmodes(k,l)*sqrt(eigs(l));
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


double genCovAnl( Array2D<double> &xy, const int i, const int j, const double clen, const double sigma, const string covtype) {

  double cov;

  double dxy = sqrt((xy(i,0)-xy(j,0))*(xy(i,0)-xy(j,0))+(xy(i,1)-xy(j,1))*(xy(i,1)-xy(j,1)));
  if ( covtype == "SqExp" )
    cov = exp( - dxy*dxy / ( clen * clen ) ) * sigma * sigma ;
  else if ( covtype == "Exp" )
    cov = exp( - dxy / clen ) * sigma * sigma ;
  else 
    throw Tantrum("kl_sample.cpp::genCovAnl(): covariance type is not recognized!");

  return ( cov );

}

void genGrid(Array1D<double> &xgrid,Array1D<double> &ygrid, const int nx, const int ny){

  if ( (nx<=0) || (ny<=0) )
    throw Tantrum("kl_sample::genGrid() : number of grid points needs to be greater than 0") ;

  xgrid.Resize(nx,0.0);
  for ( int i=0; i<nx; i++ ) xgrid(i) = (double) i / ((double) nx - 1.0);
  ygrid.Resize(ny,0.0);
  for ( int i=0; i<ny; i++ ) ygrid(i) = (double) i / ((double) ny - 1.0);

  for ( int i=0; i<nx; i++ ) xgrid(i) = dcomp(xgrid(i), 0.5, 1.1, xgrid(nx-1));
  for ( int i=0; i<ny; i++ ) ygrid(i) = dcomp(ygrid(i), 0.5, 1.1, ygrid(ny-1));

  return;

}

void getGrid1dEquiv(Array1D<double> &xgrid, Array1D<double> &ygrid, Array1D<double> &xg1d) {

  int nx = (int) xgrid.XSize();
  int ny = (int) ygrid.XSize();
  assert( nx*ny == (int) xg1d.XSize() );

  Array1D<double> w(nx*ny,0.0) ;

  double dxm,dxp,dym,dyp;
  for ( int k = 0; k < nx*ny; k++){
    int i=k%nx;
    int j=k/nx;
    if ( i==0 )
      dxm=0.0;
    else
      dxm=xgrid(i)-xgrid(i-1);
    if ( i==nx-1 )
      dxp=0.0;
    else
      dxp=xgrid(i+1)-xgrid(i);
    if ( j==0 )
      dym=0.0;
    else
      dym=ygrid(j)-ygrid(j-1);
    if ( j==ny-1 )
      dyp=0.0;
    else
      dyp=ygrid(j+1)-ygrid(j);
    w(k)=0.25*(dxm*(dym+dyp)+dxp*(dym+dyp));
  }

  xg1d(1)=w(0)*2.0+xg1d(0);
  for ( int k = 2; k < nx*ny; k++)
    xg1d(k)=2.0*w(k-1)+xg1d(k-2);

  write_datafile_1d(xg1d,"xg1d.dat");

  return ;

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
