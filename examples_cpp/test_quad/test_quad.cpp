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
#include "math.h"
#include "string.h"
#include "uqtktools.h"
using namespace std ;

#define SWAP(t,a,b) ( (t) = (a), (a) = (b), (b) = (t) )

string convertInt(int number);
void outXW( const char * froot, const int kind, Array1D<double> x, Array1D<double> w ) ;

int main() {

  int npts=9 ;
  Array1D<double> x(npts),w(npts);
  Array1D<double> tmp1,tmp2;
  double a=0.1, b=0.8;

  /* Check kinds*/
  for ( int k=1; k<=6; k++ ) {
#ifndef OLDQ
    //cout<<"Testing kind = "<<k<<endl;
    gq ( k, a, b, x, w ) ;
    outXW((char *)"newq",k,x,w) ;
#else
    gaussqC(k, x.XSize(), a, b, 0, tmp1, tmp2, x, w) ;
    outXW((char *)"oldq",k,x,w) ;
#endif
  }

  /* Check Vandermonde */
  if (npts==1){
    x(0)=0.0;
    w(0)=2.0;
  }
  else {
    for (int i=0;i<npts;i++)
      x(i) =  -1.0+2.0*double(i)/ double(npts-1);
    
    Array1D<double> q(npts,0.e0);
    for (int i=0; i<npts; i=i+2) q(i)=2.0/(i+1.0);
#ifndef OLDQ
    vandermonde_gq(x,w,q);
    outXW((char *)"newq",10,x,w) ;
#else
    vander(x,w,q);
    outXW((char *)"oldq",10,x,w) ;
#endif
  }

  /* Check Log-normal (as a generic recursive) */
  Array1D<double> al(npts,0.0), be(npts,0.0);
  double logmu = 0.1;
  double sigma = 0.4;

  double es2=exp(sigma*sigma);
  for(int i=0; i<npts;i++){
      al(i)=exp(logmu)*pow(es2,1.0*i-0.5)*((es2+1.0)*pow(es2,i)-1.e0);
      be(i)=exp(2.0*logmu)*pow(es2,3.0*i-2.0)*(pow(es2,i)-1.0);
  }

#ifndef OLDQ
  gq_gen(al, be, 1.0, x, w) ;
  double stmp ;
  for(int i=0; i<npts/2;i++) {
    SWAP(stmp,x(i),x(npts-1-i));
    SWAP(stmp,w(i),w(npts-1-i));
  }
  outXW((char *)"newq",11,x,w) ;
#else
  gaucof(al, be, 1.0, x, w) ;
  outXW((char *)"oldq",11,x,w) ;
#endif

  return (0) ;

}

void outXW( const char * froot, const int kind, Array1D<double> x, Array1D<double> w ) {

  string fname(froot);
  fname.append(convertInt(kind).c_str());
  fname.append(".dat");

  FILE *fout=fopen(fname.c_str(),"w");
  for ( int i=0; i<(int) x.XSize(); i++ )
    fprintf(fout,"%20.12e %20.12e\n",x(i),w(i));
  fclose(fout);

  return ;

}

string convertInt(int number)
{
   stringstream ss;
   ss << number;
   return ss.str();
}
