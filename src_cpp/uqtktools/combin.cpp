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
#include "Array1D.h"
#include "Array2D.h"
#include "gen_defs.h"
#include "probability.h"
#include "combin.h"
#include <math.h>
#include <float.h>
#include "error_handlers.h"


int choose(int n,int k)
{

//make sure input is nonnegative
  k=MIN(k,n-k);/////

int num=1, den=1;

    for(int i=0; i<k;i++){
      num = num*(n-i);
      den = den*(1+i);
    } 
  
  int nup = num/den;

if (k<0)
return 0;
else
  return nup;


}


int factorial(int number) {
	int temp;

	if(number <= 1) return 1;

	temp = number * factorial(number - 1);
	return temp;
}


void chooseComb(int n, int k,Array2D<int>& fullInd)
{
 
  
  int n_k=choose(n,k);
  fullInd.Resize(n_k,k,0);

  //  fullInd(0,0)=0;
  for(int ik=0;ik<k;ik++)
    fullInd(0,ik)=ik;
  fullInd(0,k-1)=k-2;

  int iii=0;
  int j=k-1;
  while(j>=0){
    if(fullInd(iii,j)<n-k+j){
      //printf("j=%d\n",j);
        fullInd(iii,j)++;
        // printf("********%d\n", fullInd(0,1));
        for(int jj=j+1;jj<k;jj++)
          fullInd(iii,jj)=fullInd(iii,j)+jj-j;
        
        j=k-1;
        iii++;
        if (iii==(int)fullInd.XSize()) {j=-1; break;}
        //        fullInd(iii,0)=iii;
        for(int ik=0;ik<k;ik++)
          fullInd(iii,ik)=fullInd(iii-1,ik);
    }
    else j--;

    
  }



  return;
}


void get_perm(Array1D<int>& perm,int seed)
{

    int nn=perm.XSize();
   
    
    get_perm(nn,perm.GetArrayPointer(),seed);


    return;
}

void get_perm(int nn, int* perm,int seed)
{

  dsfmt_gv_init_gen_rand(seed );
   
    int j,t;
    for(int is=0;is<nn;is++)
	perm[is]=is;

     for(int is=0;is<nn;is++){
	 j=is+(int) (dsfmt_gv_genrand_urv()*(nn-is));
	 t=perm[is];
	 perm[is]=perm[j];
	 perm[j]=t;
     }     


    return;
}




double gammai ( const double p, const double x )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMDS computes the incomplete Gamma integral.
//
//  Discussion:
//
//    The parameters must be positive.  An infinite series is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by Chi Leung Lau.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Chi Leung Lau,
//    Algorithm AS 147:
//    A Simple Series for the Incomplete Gamma Integral,
//    Applied Statistics,
//    Volume 29, Number 1, 1980, pages 113-114.
//
//  Parameters:
//
//    Input, double X, P, the arguments of the incomplete
//    Gamma integral.  X and P must be greater than 0.
//
//    Output, int *IFAULT, error flag.
//    0, no errors.
//    1, X <= 0 or P <= 0.
//    2, underflow during the computation.
//
//    Output, double GAMMDS, the value of the incomplete
//    Gamma integral.
/////////////////////////////////////////////////////////////////////////////////////////
// The code is taken from http://people.sc.fsu.edu/~jburkardt/cpp_src/asa147/asa147.html
// UQTk group has modified the code in the following way:
// - Renamed this function from gammds to gammai for backward compatibility
// - Instead of integer fault indicator, we throw Tantrum and exit when parameters or arguments are outside bounds
// - Log of complete gamma function is computed by lgamma instead of being computed by another function from the same code-suite
// - Instead of an underflow error-message, implemented an approximation (i.e. returning 1) for large arguments when the output value is near 1
{
  double a;
  double arg;
  double c;
  double e = 1.0E-09;
  double f;
  //int ifault2;
  double uflo = 1.0E-37;
  double value;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    throw Tantrum("gammai() error:: the argument of the incomplete gamma function is non-positive. Exit.");
    value = 0.0;
    return value;
  }

  if ( p <= 0.0 ) 
  {
    throw Tantrum("gammai() error:: the parameter of the incomplete gamma function is non-positive. Exit.");
    value = 0.0;
    return value;
  }
//
//
  arg = p * log ( x ) - lgamma ( p + 1.0 ) - x;

  

  f = exp ( arg );

  if ( f < uflo )
  {
    // This is added by UQTk group as an aproximation in the case of large arguments
    value = 1.0;
    return value;
  }

//
//  Series begins.
//
  c = 1.0;
  value = 1.0;
  a = p;

  for ( ; ; )
  {
    a = a + 1.0;
    c = c * x / a;
    value = value + c;

    if ( c <= e * value )
    {
      break;
    }
  }

  value = value * f;

  return value;
}

double beta(const double z, const double w)
{
  return exp(lgamma(z)+lgamma(w)-lgamma(z+w));
}




double betai ( const double p, const double q, const double x )
//****************************************************************************
//
//  Purpose:
//
//    BETAIN computes the incomplete Beta function ratio.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    KL Majumder, GP Bhattacharjee,
//    Algorithm AS 63:
//    The incomplete Beta Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 409-411.
//
//  Parameters:
//
//    Input, double X, the argument, between 0 and 1.
//
//    Input, double P, Q, the parameters, which
//    must be positive.
//
//    Input, double BETA, the logarithm of the complete
//    beta function.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    nonzero, an error occurred.
//
//    Output, double BETAIN, the value of the incomplete
//    Beta function ratio.
/////////////////////////////////////////////////////////////////////////////////////////
// The code is taken from http://people.sc.fsu.edu/~jburkardt/cpp_src/asa063/asa063.html
// UQTk group has modified the code in the following way:
// - Renamed this function from betain to betai for backward compatibility
// - Instead of integer fault indicator, we throw Tantrum and exit when parameters or arguments are outside bounds
// - Log of complete beta function is computed instead of being given as an argument
//
{
  double acu = 0.1E-14;
  double ai;
  //double betain;
  double cx;
  bool indx;
  int ns;
  double pp;
  double psq;
  double qq;
  double rx;
  double temp;
  double term;
  double value;
  double xx;




  value = x;
//
//  Check the input arguments.
//
  if ( p <= 0.0 || q <= 0.0 )
  {
    throw Tantrum("betai() error:: parameters of the incomplete beta function are negative. Exit.");

    return 0.0;
  }

  if ( x < 0.0 || 1.0 < x )
  {
    throw Tantrum("betai() error:: argument of the incomplete beta function are outside bounds [0,1]. Exit.");
    
    return 0.0;
  }
//
//  Special cases.
//
  if ( x == 0.0 || x == 1.0 )
  {
    return value;
  }

  // Added by Sandia UQTk group
  double lbeta=lgamma(p)+lgamma(q)-lgamma(p+q);

//
//  Change tail if necessary and determine S.
//
  psq = p + q;
  cx = 1.0 - x;

  if ( p < psq * x )
  {
    xx = cx;
    cx = x;
    pp = q;
    qq = p;
    indx = true;
  }
  else
  {
    xx = x;
    pp = p;
    qq = q;
    indx = false;
  }

  term = 1.0;
  ai = 1.0;
  value = 1.0;
  ns = ( int ) ( qq + cx * psq );
//
//  Use the Soper reduction formula.
//
  rx = xx / cx;
  temp = qq - ai;
  if ( ns == 0 )
  {
    rx = xx;
  }

  for ( ; ; )
  {
    term = term * temp * rx / ( pp + ai );
    value = value + term;;
    temp = fabs ( term );

    if ( temp <= acu && temp <= acu * value )
    {
      value = value * exp ( pp * log ( xx ) 
      + ( qq - 1.0 ) * log ( cx ) - lbeta ) / pp;

      if ( indx )
      {
        value = 1.0 - value;
      }
      break;
    }

    ai = ai + 1.0;
    ns = ns - 1;

    if ( 0 <= ns )
    {
      temp = qq - ai;
      if ( ns == 0 )
      {
        rx = xx;
      }
    }
    else
    {
      temp = psq;
      psq = psq + 1.0;
    }
  }

  return value;
}

double digama ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    DIGAMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2013
//
//  Author:
//
//    Original FORTRAN77 version by Jose Bernardo.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jose Bernardo,
//    Algorithm AS 103:
//    Psi ( Digamma ) Function,
//    Applied Statistics,
//    Volume 25, Number 3, 1976, pages 315-317.
//
//  Parameters:
//
//    Input, double X, the argument of the digamma function.
//    0 < X.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    1, X <= 0.
//
//    Output, double DIGAMA, the value of the digamma function at X.
//
/////////////////////////////////////////////////////////////////////////////////////////
// The code is taken from http://people.sc.fsu.edu/~jburkardt/cpp_src/asa103/asa103.cpp
// UQTk group has modified the code in the following way:
// - Instead of integer fault indicator, we throw Tantrum and exit when parameters or arguments are outside bounds
//
{
  double euler_mascheroni = 0.57721566490153286060;
  double r;
  double value;
  double x2;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    value = 0.0;
    throw Tantrum("digama() error:: argument of the digamma function is non-positive. Exit.");
    return value;
  }
//
//  Initialize.
//
  x2 = x;
  value = 0.0;
//
//  Use approximation for small argument.
//
  if ( x2 <= 0.00001 )
  {
    value = - euler_mascheroni - 1.0 / x2;
    return value;
  }
//
//  Reduce to DIGAMA(X + N).
//
  while ( x2 < 8.5 )
  {
    value = value - 1.0 / x2;
    x2 = x2 + 1.0;
  }
//
//  Use Stirling's (actually de Moivre's) expansion.
//
  r = 1.0 / x2;
  value = value + log ( x2 ) - 0.5 * r;
  r = r * r;
  value = value 
    - r * ( 1.0 / 12.0
    - r * ( 1.0 / 120.0 
    - r *   1.0 / 252.0 ) );

  return value;
}
//****************************************************************************
