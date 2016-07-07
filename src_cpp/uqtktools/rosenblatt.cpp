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
#include "Array1D.h"
#include "Array2D.h"
#include "uqtktools.h"

using namespace std;



void invRos(Array1D<double>& unif, Array2D<double>& xi, Array1D<double>& newXi, Array1D<double>& sig)
{
  // Accuracy, i.e. the stopping criterion in the bisection algorithm
  double xacc=1e-8;
  // The initial range of the bisection algorithm [-xmax,xmax]
  double xmax=1e+8;
	
	
  // Dimensionality of the transformation
  int ndim = unif.XSize();
  int ns	 = xi.YSize();
	
  // Output container
  newXi.Resize(ndim,0.e0);
	
  // Dimension check
  if (ndim != (int) xi.XSize()  || ndim != (int) sig.XSize())
    {printf("invRos: dimension error\n"); exit(1);}

	
  // Work arrays
  Array1D<double> kern(ns,1.e0);
  Array1D<double> numer(ndim,0.e0);
  Array1D<double> denom(ndim,0.e0);
	
  for(int id=0; id < ndim; id++){
    ////  including this is equivalent to dropping | conditionals in Rosenblatt transformation
    //  kern.Resize(ndim,1.e0); 
		
    // Get the corresponding uniform r.v. sample
    double Un=unif(id);
		
    // Starting point of bisection
    double xx=0.0; 
		
    // Cure against extremes
    if (Un==0.0)
      Un+=DBL_EPSILON;
    if (Un==1.0)
      Un-=DBL_EPSILON;
		
    // Build denominator
    for(int is=0; is < ns; is++)
      denom(id) += kern(is);
		
    // Bisection iteration counter
    int ii=0;
		
    // Bisection looks in the interval [x1,x2]
    double x1=-xmax;
    double x2=xmax;
		
    // Bisection loop
    do{
      numer(id)=0;
      for(int is=0; is < ns; is++)
        numer(id) += (kern(is)/denom(id)) *  ( 0.5+0.5*erf( (xx-xi(id,is))/(sqrt(2)*sig(id)) ) ) ;// (xx > xi(id,is));
			
      if (Un < numer(id)){ x1=x1; x2=xx;} 
      else {x1=xx; x2=x2;}
      xx=(x1+x2)/2.0;
      ii++;
    } while (fabs(x2-x1) >  xacc);
    // End of bisection loop
		
		
    // Just-in-case warning
    if (ii>1000)
      printf( "Warning: Inverted CDF in really many (%d) iterations\n", ii);
		
    // Select the solution as a new sample
    newXi(id)=xx;
		
    // Update the kernel for the next dimension
    for(int is=0; is < ns; is++)
      kern(is) = kern(is)*exp(-(newXi(id)-xi(id,is))*(newXi(id)-xi(id,is))/(2*sig(id)*sig(id)));           // (bw*sqrt(2*PI)) cancels;
		
  }
	
  return;
}





void invRos(Array1D<double>& unif,  Array2D<double>& xi, Array1D<double>& newXi, double bw)
{
  if (bw<=0)
    {printf("invRos: bandwidth needs to be positive"); exit(1);}

  int ndim = unif.XSize();
 
  // dimension check
  if (ndim != (int) xi.XSize())
    {printf("invRos: dimension error"); exit(1);}
  
  Array1D<double> sig(ndim,bw);
  
  invRos(unif,xi,newXi,sig);
  
  return;
}




void invRos(Array1D<double>& unif, Array2D<double>& xi, Array1D<double>& newXi)
{
  int nkl = unif.XSize();
  const int nr_proj = xi.YSize();

  // dimension check
  if (nkl != (int) newXi.XSize() || nkl > (int) xi.XSize())
    {printf("invRos: dimension error"); exit(1);}
  
  Array1D<double> sig(nkl,0.e0);
  
  
  Array2D<double> xi_t(nr_proj,nkl,0.e0);
  transpose(xi,xi_t);
  
  get_opt_KDEbdwth(xi_t,sig);
  
  invRos(unif,xi,newXi,sig);
  
  return;
}

void get_opt_KDEbdwth(const Array2D<double>& data,Array1D<double>& bdwth)
{
  int ndata=data.XSize();
  int ndim=data.YSize();
  if (ndim != (int) bdwth.XSize()) {printf("KDEBdwth dimension error\n"); exit(-1);}

  Array1D<double> flag(ndim,1.);
  

  Array1D<double> datamin(ndim,1000.0);
  Array1D<double> datamax(ndim,-1000.0);

  for(int idim=0;idim<ndim;idim++){
  int numBorder=0;
    for(int idata=0;idata<ndata;idata++){ 
      if (data(idata,idim)<datamin(idim)) {datamin(idim)=data(idata,idim);}
      if (data(idata,idim)>datamax(idim)) {datamax(idim)=data(idata,idim);}
    }
      double nearBorder=(datamax(idim)-datamin(idim))/20.;
     for(int idata=0;idata<ndata;idata++){ 
       if ( data(idata,idim)-datamin(idim)<nearBorder || datamax(idim)-data(idata,idim)<nearBorder ) {numBorder++;}
     }
     // printf("numBorder=%d\n",numBorder);
       if ( numBorder > ndata/20.) {flag(idim)=0.5;}
  } 
  
  Array1D<double> stdd(ndim,0.e0);
  Array1D<double> data_1d(ndata,0.e0);

  for(int idim=0;idim<ndim;idim++){
    for(int idata=0;idata<ndata;idata++){ 
      data_1d(idata)=data(idata,idim);

    }
    stdd(idim)=get_std(data_1d);
    //printf("**stdd********%lg\n",stdd(idim));
 
    bdwth(idim)=flag(idim)*pow(4./(ndim+2),1./(ndim+4.))*stdd(idim)*pow(ndata,-1./(ndim+4.));
  }

  return;
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


void invRos_old(Array1D<double>& unif, Array2D<double>& xi, Array1D<double>& newXi, Array1D<double>& sig)
{
 
  int nkl = unif.XSize();
  const int nr_proj = xi.YSize();
  
  // dimension check
  if (nkl != (int) newXi.XSize() || nkl > (int) xi.XSize() || nkl != (int) sig.XSize())
    {printf("invRos: dimension error"); exit(1);}
  
  Array1D<double> kern(nr_proj,1.e0);
  Array1D<double> numer(nkl,0.e0);
  Array1D<double> denom(nkl,0.e0);
  
  for(int ikl=0; ikl < nkl; ikl++){
    //  kern.Resize(nr_proj,1.e0); //including this is equivalent to dropping | conditionals in Rosenblatt transformation
    
    double Un=unif(ikl);
    double xx=xi(ikl,0); //start somewhere 'random'
    double step=2.;
    if (Un==0.0)
      Un+=.001/nr_proj;
    if (Un==1.0)
      Un-=.001/nr_proj;
    
    int fl;
    int oldfl;
    
    denom(ikl)=0;
    for(int ir=0; ir < nr_proj; ir++)
      denom(ikl) += kern(ir);
    
    int ii=0;
    do{
      numer(ikl)=0;
      for(int ir=0; ir < nr_proj; ir++)
	numer(ikl) += (kern(ir)/denom(ikl)) *  ( 0.5+0.5*erf( (xx-xi(ikl,ir))/(sqrt(2)*sig(ikl)) ) ) ;// (xx > xi(ikl,ir));
      
      if (Un < numer(ikl)){ xx=xx-step/2;fl=-1;} 
      else {xx=xx+step/2;fl=1;}
      if (oldfl != fl && ii>0) {step=step/2;}
      oldfl=fl;
      ii++;
      // cout << ii <<":  "<< step << " " << fl << endl;
    } while (step >  1./(2.*nr_proj));
    
    newXi(ikl) = xx;
    
    for(int ir=0; ir < nr_proj; ir++)
      kern(ir) = kern(ir)*exp(-(newXi(ikl)-xi(ikl,ir))*(newXi(ikl)-xi(ikl,ir))/(2*sig(ikl)*sig(ikl)));           // (bw*sqrt(2*PI)) cancels;
    
  }
     
  return;
}
