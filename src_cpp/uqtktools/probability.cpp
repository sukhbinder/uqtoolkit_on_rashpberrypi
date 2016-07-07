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
#include <math.h>
#include <iostream>
#include <float.h>
#include <limits.h>

#include "Array1D.h"
#include "Array2D.h"
#include "probability.h"
#include "combin.h"
#include "arraytools.h"
#include "gen_defs.h"




using namespace std;




double erff(const double x)
{
	return x < 0.0 ? -gammai(0.5,x*x) : gammai(0.5,x*x);
}


double inverf(double y0)
{
  double result_;
  
    double expm2;
    double s2pi;
    double x;
    double y;
    double z;
    double y2;
    double x0;
    double x1;
    int code;
    double p0;
    double q0;
    double p1;
    double q1;
    double p2;
    double q2;


    y0=0.5*(y0+1);
    
    expm2 = 0.13533528323661269189;
    s2pi = 2.50662827463100050242;
    if( y0<=0 or y0 >=1)
    {
      cout << "Error in inverting erf, the argument should be stricly between -1 and 1 " << y0 << endl;
      return 0;
      exit(1);
    }
            

    code = 1;
    y = y0;
    if( y>1.0-expm2 )
    {
        y = 1.0-y;
        code = 0;
    }
    if( y>expm2 )
    {
        y = y-0.5;
        y2 = y*y;
        p0 = -59.9633501014107895267;
        p0 = 98.0010754185999661536+y2*p0;
        p0 = -56.6762857469070293439+y2*p0;
        p0 = 13.9312609387279679503+y2*p0;
        p0 = -1.23916583867381258016+y2*p0;
        q0 = 1;
        q0 = 1.95448858338141759834+y2*q0;
        q0 = 4.67627912898881538453+y2*q0;
        q0 = 86.3602421390890590575+y2*q0;
        q0 = -225.462687854119370527+y2*q0;
        q0 = 200.260212380060660359+y2*q0;
        q0 = -82.0372256168333339912+y2*q0;
        q0 = 15.9056225126211695515+y2*q0;
        q0 = -1.18331621121330003142+y2*q0;
        x = y+y*y2*p0/q0;
        x = x*s2pi;
        result_ = x;
        return result_/sqrt(2.);
    }
    x = sqrt(-2.0*log(y));
    x0 = x-log(x)/x;
    z = 1.0/x;
    if( x<8.0 )
    {
        p1 = 4.05544892305962419923;
        p1 = 31.5251094599893866154+z*p1;
        p1 = 57.1628192246421288162+z*p1;
        p1 = 44.0805073893200834700+z*p1;
        p1 = 14.6849561928858024014+z*p1;
        p1 = 2.18663306850790267539+z*p1;
        p1 = -1.40256079171354495875*0.1+z*p1;
        p1 = -3.50424626827848203418*0.01+z*p1;
        p1 = -8.57456785154685413611*0.0001+z*p1;
        q1 = 1;
        q1 = 15.7799883256466749731+z*q1;
        q1 = 45.3907635128879210584+z*q1;
        q1 = 41.3172038254672030440+z*q1;
        q1 = 15.0425385692907503408+z*q1;
        q1 = 2.50464946208309415979+z*q1;
        q1 = -1.42182922854787788574*0.1+z*q1;
        q1 = -3.80806407691578277194*0.01+z*q1;
        q1 = -9.33259480895457427372*0.0001+z*q1;
        x1 = z*p1/q1;
    }
    else
    {
        p2 = 3.23774891776946035970;
        p2 = 6.91522889068984211695+z*p2;
        p2 = 3.93881025292474443415+z*p2;
        p2 = 1.33303460815807542389+z*p2;
        p2 = 2.01485389549179081538*0.1+z*p2;
        p2 = 1.23716634817820021358*0.01+z*p2;
        p2 = 3.01581553508235416007*0.0001+z*p2;
        p2 = 2.65806974686737550832*0.000001+z*p2;
        p2 = 6.23974539184983293730*0.000000001+z*p2;
        q2 = 1;
        q2 = 6.02427039364742014255+z*q2;
        q2 = 3.67983563856160859403+z*q2;
        q2 = 1.37702099489081330271+z*q2;
        q2 = 2.16236993594496635890*0.1+z*q2;
        q2 = 1.34204006088543189037*0.01+z*q2;
        q2 = 3.28014464682127739104*0.0001+z*q2;
        q2 = 2.89247864745380683936*0.000001+z*q2;
        q2 = 6.79019408009981274425*0.000000001+z*q2;
        x1 = z*p2/q2;
    }
    x = x0-x1;
    if( code!=0 )
    {
        x = -x;
    }
    result_ = x;
    
    return result_/sqrt(2.);
}


double invnormcdf(double y)
{
  return sqrt(2.)*inverf(2.*y-1.);
}

double normcdf(double y)
{
  return 0.5*(1.+erf(y/sqrt(2.)));
}

double normcdfc(double y)
{
  return 1.-normcdf(y);
}

void generate_uniform(double* rvar,int ns, int nd, int zSeed)
{
 

  dsfmt_gv_init_gen_rand(zSeed );
  
  
  for(int is = 0 ; is < ns*nd ; is++)
    rvar[is]=dsfmt_gv_genrand_urv();

      
 
  return;
}

void generate_uniform(Array2D<double>& rvar,int zSeed)
{
  int nsample = (int) rvar.XSize();
  int ndim    = (int) rvar.YSize();


  generate_uniform(rvar.GetArrayPointer(), nsample, ndim, zSeed);

  return;
}

void generate_uniform(double *rvar, int ns, int nd, dsfmt_t *rnstate)
{
  int i ;
  // Need to check allocation?
  for ( i=0; i<ns*nd; i++ )
      rvar[i] = dsfmt_genrand_open_open( rnstate ) ;

  return ;

}

void generate_uniform(Array2D<double> &rvar, dsfmt_t *rnstate)
{
  int nsample = (int) rvar.XSize() ;
  int ndim    = (int) rvar.YSize() ;

  generate_uniform(rvar.GetArrayPointer(),nsample, ndim, rnstate) ;

  return ;

}

void generate_uniform_lhs(double *rvar,int nsample, int ndim, int zSeed)
{
 
  int    *perm  = (int    *) malloc(nsample*   sizeof(int   )) ;


  dsfmt_gv_init_gen_rand(zSeed );
  
  int ii=0;
  for(int id=0;id<ndim;id++){
    
    int seed=(int) (dsfmt_gv_genrand_urv()*INT_MAX);
    get_perm(nsample,perm, seed);
    
    for(int is=0;is<nsample;is++){
      double urv=dsfmt_gv_genrand_urv();
      rvar[ii]=(urv+perm[is])/nsample;
      ii++;
    }
  }

  free(perm);

  return;
}

void generate_uniform_lhs(Array2D<double>& rvar,int zSeed)
{
  int ndim=rvar.YSize();
  int nsample=rvar.XSize();

  generate_uniform_lhs(rvar.GetArrayPointer(),nsample, ndim,zSeed);

  


  return;
}

void generate_uniform_lhs(double *rvar,int nsample, int ndim, dsfmt_t *rnstate)
{
 
  int    *perm  = (int    *) malloc(nsample*   sizeof(int   )) ;


  
  int ii=0;
  for(int id=0;id<ndim;id++){
    
    int seed = dsfmt_genrand_uint32(rnstate);
    get_perm(nsample,perm,seed);
    
    for(int is=0;is<nsample;is++){
      double urv=dsfmt_genrand_open_open( rnstate );;
      rvar[ii]=(urv+perm[is])/nsample;
      ii++;
    }
  }

  free(perm);

  return;
}

void generate_uniform_lhs(Array2D<double>& rvar,dsfmt_t *rnstate)
{
  int ndim=rvar.YSize();
  int nsample=rvar.XSize();

  generate_uniform_lhs(rvar.GetArrayPointer(),nsample, ndim,rnstate);

  


  return;
}


void generate_normal(Array2D<double>& rvar,int zSeed)
{
  int ndim=rvar.YSize();
  int nsample=rvar.XSize();

    dsfmt_gv_init_gen_rand( zSeed );


   
   for(int is = 0 ; is < nsample ; is++){
      for(int id = 0 ; id < ndim ; id++){
	rvar(is,id)=dsfmt_gv_genrand_nrv();
      }
 }
  return;
}

void generate_normal_lhs(Array2D<double>& rvar,int zSeed)
{
  int ndim=rvar.YSize();
  int nsample=rvar.XSize();

  generate_uniform_lhs(rvar,zSeed);
  for(int is = 0 ; is < nsample ; is++){
    for(int id = 0 ; id < ndim ; id++){
      rvar(is,id)=invnormcdf(rvar(is,id));
    }
  }


  return;
}




double get_median(const Array1D<double>& data)
{

  

  int ndata=data.XSize();
  int k=(int) ndata/2;
 
  Array1D<double> data_copy(ndata,0.e0);
  for (int i=0;i<ndata;i++){data_copy(i)=data(i);}
  double median=select_kth(k,data_copy);
 


  

  return median;

}

double get_mean(const Array1D<double>& data)
{
double mean=0.0;
  

  int ndata=data.XSize();

  for(int i=0;i<ndata;i++){
    mean=mean+data(i);
  
  }
  mean=mean/ndata;
 

  

  return mean;

}

double get_std(const Array1D<double>& data)
{
  double mean=get_mean(data); 
  double mean2=0.0;

  int ndata=data.XSize(); //for (int j=1;j<ndata;j++){printf("**stdd********%lg\n",data(j));}

  for(int i=0;i<ndata;i++){
    mean2=mean2+data(i)*data(i);
  }
 
  mean2=mean2/ndata;



  double std=sqrt(mean2-mean*mean);
  if (std==0) {printf("Heads Up! Standard Deviation=0\n");}

  return std;

}


void ihsU(Array2D<double> &rndnos, int dfac, dsfmt_t *rnstate) {

  int np = (int) rndnos.XSize() ;
  int ns = (int) rndnos.YSize() ;
  double *rndnosPNT = rndnos.GetArrayPointer() ;
  ihsU(ns, np, rndnosPNT, dfac, rnstate) ;

}
void ihsU(int ns, int np, double *rndnos, int dfac, dsfmt_t *rnstate) {

  int *ipos = (int *) malloc( ns*np*sizeof(int)) ;
  ihsP(ns, np, ipos, dfac, rnstate) ;

  for ( int j = 0; j < ns; j++)
    for ( int i = 0; i < np; i++)
      rndnos[j*np+i] = (((double) ipos[j*np+i])+dsfmt_genrand_open_open( rnstate ))
	              *2.0/((double) ns)-1.0 ;

  free(ipos) ;

  return ;

}

void ihsP(int ns, int np, int *rpos, int dfac, dsfmt_t *rnstate)
{

  /* Allocate space */
  int *binAvail = (int *) malloc( ns        *np*sizeof(int)) ;
  int *listTry  = (int *) malloc((ns-1)*dfac*np*sizeof(int)) ;
  int *tmpSort  = (int *) malloc( ns           *sizeof(int)) ;

  double dopt = ((double) ns) / pow((double) ns,1.0/np);
  for ( int i=0; i<ns*np; i++ ) rpos[i] = 0 ;

  for ( int i=0; i<np; i++ ) 
    rpos[(ns-1)*np+i] = dsfmt_genrand_uint32(rnstate) % ns ;

  for ( int j = 0; j < ns; j++ )
    for ( int i = 0; i < np; i++ ) 
      binAvail[j*np+i] = j ;

  for ( int i = 0; i < np; i++ ) 
    binAvail[rpos[(ns-1)*np+i]*np+i] = ns ;

  for ( int i = 0; i < np; i++ ) {
    for ( int j = 0; j < ns; j++ )
      tmpSort[j] = binAvail[j*np+i] ;
    shell_sort (tmpSort, ns);
    for ( int j = 0; j < ns; j++ )
      binAvail[j*np+i] = tmpSort[j] ;
  }

  for ( int iter = ns-2; iter > 0; iter--) {

    /* Generate (n-i)*dfac samples from list of available bins */
    for ( int i = 0; i < np; i++ )
      for ( int j = 0; j < iter*dfac; j++) {
	int idx = dsfmt_genrand_uint32(rnstate) % (iter+1) ;
	listTry[j*np+i] = binAvail[idx*np+i] ;
      }

    /* Look for sample with minimum distance closest to optimum */
    double minall = 1.0e30;
    int    sidx   = -1;
    for ( int j = 0; j < iter*dfac; j++) {
      double minspl = 1.0e30;
      for (int j1 = iter + 1; j < ns; j++) {
	double dist = 0.0 ;
	for ( int i = 0; i < np; i++ )
	  dist += pow(listTry[j*np+i]-rpos[j1*np+i],2);
        minspl = MIN(minspl,sqrt(dist)) ;
      }

      if (fabs(minspl - dopt) < minall) {
        minall = fabs(minspl - dopt); 
        sidx = j;
      }

    } /* done looking for minimum */

    /* Chose point with optimum distance */
    for ( int i = 0; i < np; i++ )
      rpos[iter*np+i] = listTry[sidx*np+i] ;

    /* Mark corresponding bins */
    for ( int i = 0; i < np; i++ )
      for ( int j = 0; j < iter; j++) 
        if (binAvail[j*np+i] == rpos[iter*np+i])
	  binAvail[j*np+i] = ns;

    /* Sort to move marked bins to the end */
    for ( int i = 0; i < np; i++ ) {
      for ( int j = 0; j < ns; j++ )
        tmpSort[j] = binAvail[j*np+i] ;
      shell_sort (tmpSort, ns);
      for ( int j = 0; j < ns; j++ )
        binAvail[j*np+i] = tmpSort[j] ;
    }

  }

  /* Last point - only one option */
  for ( int i = 0; i < np; i++ ) rpos[i] = binAvail[i] ;

  free(binAvail) ;
  free(listTry ) ;
  free(tmpSort ) ;

  return ;

}

/* 
  Returns a random permutation of 0..n-1
*/
void rperm(int n, int *a, dsfmt_t *rnstate) 
{
  int k;

  for (k = 0; k < n; k++) a[k] = k;

  for ( k = n-1; k > 0; k-- ) 
  {
    int j = dsfmt_genrand_uint32(rnstate) % (k+1) ;
    int temp = a[j];
    a[j] = a[k];
    a[k] = temp;
  }
    return ;
}

void shell_sort (int *a, int n) {

  int j ;
  for (int h = n/2; h>0; h = h/2) {

    for ( int i = h; i < n; i++) {

      int k = a[i];
      for ( j = i; j >= h && k < a[j - h]; j -= h)
	a[j] = a[j - h];
      a[j] = k;

    }

  }

  return ;

}
