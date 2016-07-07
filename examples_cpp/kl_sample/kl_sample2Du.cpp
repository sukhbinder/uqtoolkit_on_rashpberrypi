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

#include <string>
#include "dsfmt_add.h"

#include "kldecompuni.h"
#include "uqtktools.h"
#include "lapack.h"
#include "error_handlers.h"

using namespace std;

#define NKL   64
#define NSPL  128
#define CLEN  0.05
#define SIG   5.0

#define COVTYPE "SqExp"
#define XFILE   "data/cali_grid.dat"
#define TFILE   "data/cali_tria.dat"


/// \brief Displays information about this program
int usage(){

  printf("usage: cor_kl [-h]  [-x<grid_file>] [-c<cov_type>] [-n<nx>] [-m<ny>] [-e<nkl>] [-p<nspl>] [-s<sigma>] [-l<cor_len>]\n");
  printf(" -h                 : print out this help message \n");
  printf(" -x <xfile>       : define the file from which the time grid is being read (default=%s) \n",XFILE);
  printf(" -c <cov_type>    : define wether using analytical covariance (type: %s) or compute it from samples \n",COVTYPE);
  printf(" -e <nkl>         : define the number of KL modes retained (default=%d) \n",NKL);
  printf(" -p <nspl>        : define the number of samples (default=%d) \n",NSPL);
  printf(" -s <sigma>       : define standard deviation (default=%e) \n",SIG);
  printf(" -l <clen>        : define correlation length (default=%e) \n",CLEN);
  printf("================================================================================================\n");
  printf("Input:: Nothing hard-coded.\n");
  printf("Output:: \n");
  printf("        - cov_out.dat:  covariance matrix\n");
  printf("        - eig.dat:      eigenvalues\n");
  printf("        - KLmodes.dat:  eigenmodes scaled with sqrt(eig)\n");
  printf("        - rel_diag.dat: pointwise covariance fraction explained by the finite KL expansion\n");
  printf("        - mean.dat:     pointwise mean based on samples\n");
  printf("        - xi_data.dat:  samples of eigenvalues\n");
  printf("================================================================================================\n");
  exit(0);
  return 0;
}

double trArea(double x1,double x2,double x3, double y1,double y2,double y3) ;
double genCovAnl( Array2D<double> &xy, const int i, const int j, const double clen, 
                  const double sigma, const string covtype) ;
void getWeights(Array2D<double> &xgrid, Array2D<int> &tgrid, Array1D<double> &w) ;
double getMean(Array2D<double> &y, const string idir, const int ij) ;


/*!
   \brief Karhunen-Loeve decomposition of a UNIVARIATE process
   given covariance type or samples drawn from multivariate normal distribution
*/
int main(int argc, char *argv[])
{

  /* Set the default values */
  int    nkl  = NKL  ;
  int    nspl = NSPL ;

  char* xfile = (char *) XFILE;
  char* tfile = (char *) TFILE;
  Array2D<double> xgrid;
  Array2D<int>    tgrid;

  bool   cflag = false ;
  double clen  = CLEN  ;
  double sigma = SIG   ;
  char* cov_type = (char *) COVTYPE;
  Array2D<double> cov ;

  /* Read the user input */
  int c;

  while ((c=getopt(argc,(char **)argv,"hc:e:p:s:l:"))!=-1){
    switch (c) {
    case 'h':
      usage();
      break;
    case 'c':
      cflag=true;
      cov_type = optarg;
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

  /* Print the input information on screen */
  cout << " - Correlation length:    " << clen << endl<<flush;
  cout << " - Standard deviation:    " << sigma<< endl<<flush;
  cout << " - Number of KL modes:    " << nkl  << endl<<flush;
  if ( cflag )
    cout<<" - Will generate covariance of type "<<cov_type<<endl<<flush;
  else    
    cout << " - Will generate covariance from "<<nspl<<" samples" << endl<<flush;
 
  /* read grid from file */
  read_datafileVS(xgrid,xfile);
  read_datafileVS(tgrid,tfile);
  cout<<" - No. of grid points:    " <<xgrid.XSize() << endl << flush ;
  cout<<" - No. of triangles  :    " <<tgrid.XSize() << endl << flush ;

  int nxy = (int) xgrid.XSize() ;
  if ( nkl > nxy )
    throw Tantrum("kl_sample2Du::main(): cannot request more KL modes than number of grid points");

  Array2D<double> ySamples;

  cov.Resize(nxy,nxy,0.e0);	
  if ( cflag ) {
    for ( int i = 0; i < nxy; i++)
      for ( int j = 0; j < nxy; j++)
        cov(i,j)=genCovAnl(xgrid,i,j,clen,sigma,cov_type);
  }
  else {

    double dfac=1.0e-12;
    bool   tryagain = true;

    while ( ( tryagain ) && ( dfac < 1.e-6 ) ) {

      tryagain = false ;
      for ( int i = 0; i < nxy; i++)
        for ( int j = 0; j < nxy; j++)
          cov(i,j)=genCovAnl(xgrid,i,j,clen,sigma,string("SqExp"));
      for ( int i = 0; i < nxy; i++) cov(i,i) += dfac;
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
  Array1D<double> weights(nxy,0.0);
  getWeights(xgrid,tgrid,weights);

  KLDecompUni decomposer;  
  decomposer.SetWeights(weights);
  int n_eig = decomposer.decompose(cov,nkl);

  if(n_eig <  nkl){ 
    printf("There are only %d  eigenvalues available (requested %d) \n",n_eig, nkl); 
    nkl = n_eig; 
  }

  const Array1D<double>& eigs    = decomposer.eigenvalues();
  const Array2D<double>& KLmodes = decomposer.KLmodes();

  cout << " --> KL decomposition done " << endl << flush;

  cout << "      - Obtained " << n_eig << " eigenvalues:" << endl;
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
    scaledKLmodes(k,0) = xgrid(k,0);
    scaledKLmodes(k,1) = xgrid(k,1);
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

  double dxy2 = (xy(i,0)-xy(j,0))*(xy(i,0)-xy(j,0))+(xy(i,1)-xy(j,1))*(xy(i,1)-xy(j,1));
  if ( covtype == "SqExp" )
    cov = exp( - dxy2 / ( clen * clen ) ) * sigma * sigma ;
  else if ( covtype == "Exp" )
    cov = exp( - sqrt(dxy2) / clen ) * sigma * sigma ;
  else 
    throw Tantrum("kl_sample.cpp::genCovAnl(): covariance type is not recognized!");

  return ( cov );

}

void getWeights(Array2D<double> &xgrid, Array2D<int> &tgrid, Array1D<double> &w) {

  int nxy = xgrid.XSize() ;
  int nt  = tgrid.XSize() ;
 
  w.Resize(nxy,0.0) ;
  for ( int i = 0; i<nxy; i++)
    for ( int j = 0; j<nt; j++)
      if (( tgrid(j,0) == i ) || ( tgrid(j,1) == i ) || ( tgrid(j,2) == i ))
	w(i) += trArea(xgrid(tgrid(j,0),0),xgrid(tgrid(j,1),0),xgrid(tgrid(j,2),0),
		       xgrid(tgrid(j,0),1),xgrid(tgrid(j,1),1),xgrid(tgrid(j,2),1))/3.0;

  return ;

}

double trArea(double x1,double x2,double x3, double y1,double y2,double y3) {
  return (fabs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2.0) ;
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
