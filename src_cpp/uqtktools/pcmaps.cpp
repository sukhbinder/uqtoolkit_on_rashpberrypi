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
#include "error_handlers.h"
#include "probability.h"
#include "datafiles.h"
#include "pcmaps.h"
#include "combin.h"
#include <assert.h>
#include <math.h>
#include <float.h>


#define EPS 1e-16
#define MXEPS 1e+100

using namespace std;


double PCtoPC(double x, string pcIn, double in1, double in2, string pcOut, double out1, double out2)
{
double y;

///////// Set domain check for x!!!////////
 if ( pcIn=="HG" || pcIn=="pdf" ){
   x=x;
 }
 else if (pcIn=="LU" or pcIn=="JB"){
   if (fabs(x)>1.0) 
     throw Tantrum("pcmaps.cpp:: PCtoPC:: input x is out of the range!");
 }
 else if (pcIn=="SW" or pcIn=="LG"){
   if (x<0.0) 
     throw Tantrum("pcmaps.cpp:: PCtoPC:: input x is out of the range!");
 }
 else if (pcIn=="TG"){
   if (x>in1)
     throw Tantrum("pcmaps.cpp:: PCtoPC:: input x is out of the range!");
 }
 else if (pcIn=="RB"){
   if (x>in2/(1.-in1) or x<0)
     throw Tantrum("pcmaps.cpp:: PCtoPC:: input x is out of the range!");
 }
 else 
 {
   printf("pcIn = %s\n",pcIn.c_str()) ;
   throw Tantrum("pcmaps.cpp::Input PC is not recognized!\n");
 }
///////// Set input parameter check!!!////////
 if (pcIn=="LU" or pcIn=="HG" or pcIn=="TG" or pcIn=="pdf"){
   x=x;
 }
 else if (pcIn=="JB"){
   if (in1<-1. or in2>1.) 
     throw Tantrum("pcmaps.cpp:: PCtoPC:: input parameter is out of the range!");
 }
 else if (pcIn=="SW"){
   if (in2<=0.)
     throw Tantrum("pcmaps.cpp:: PCtoPC:: input parameter is out of the range!");
 }
 else if (pcIn=="LG"){
   if (in1<-1.)
     throw Tantrum("pcmaps.cpp:: PCtoPC:: input parameter is out of the range!");
 }
 else if (pcIn=="RB"){
   if (in2<=0. or in1<=0)
	throw Tantrum("pcmaps.cpp:: PCtoPC:: input parameter is out of the range!");
 }
 else 
   throw Tantrum("pcmaps.cpp::Input PC is not recognized!\n");

///////// Set output parameter check!!!////////
 if (pcOut=="LU" or pcOut=="HG" or pcOut=="TG" or pcOut=="pdf"){
   x=x;
 }
 else if (pcOut=="JB"){
   if (out1<-1. or out2>1.) 
     throw Tantrum("pcmaps.cpp:: PCtoPC:: output parameter is out of the range!");
 }
 else if (pcOut=="SW"){
   if (out2<=0.)
     throw Tantrum("pcmaps.cpp:: PCtoPC:: output parameter is out of the range!");
 }
 else if (pcOut=="LG") {
   if (out1<0.)
     throw Tantrum("pcmaps.cpp:: PCtoPC:: output parameter is out of the range!");
 }
 else if (pcOut=="RB"){
   if (out2<=0. or out1<=0)
     throw Tantrum("pcmaps.cpp:: PCtoPC:: output parameter is out of the range!");
 }
 else 
 {
   throw Tantrum("pcmaps.cpp::Output PC is not recognized!\n");
 }
 ////////////////////////////////////////////
 if (pcIn=="LU" && pcOut=="LU"){
   y=x;
 }
 else if (pcIn=="HG" && pcOut=="HG"){
   y=x;
 }
 else if (pcIn=="SW" && pcOut=="SW"){
   y=pow(x,out2/in2)*exp(out1-in1*out2/in2);
 }
 else if (pcIn=="JB" && pcOut=="JB"){
   if(fabs(x)==1) y=x;
   else
     y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"JB",out1,out2);
 }
 else if (pcIn=="LG" && pcOut=="LG"){
   if (x==0) y=0.;
   else y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"LG",out1,0);
 }
 else if (pcIn=="TG" && pcOut=="TG"){
   if (x==in1) y=out1;
   else if (in1==out1) y=x;
   else y=PCtoPC(PCtoPC(x,"TG",in1,0,"LU",0,0),"LU",0,0,"TG",out1,0);
 }
 else if (pcIn=="RB" && pcOut=="RB"){
   y=PCtoPC(PCtoPC(x,"RB",in1,in2,"LU",0,0),"LU",0,0,"RB",out1,out2);
 }
 ////////////////////////////////////////////////////////////////////
else if (pcIn=="LU" && pcOut=="HG"){
  //if(x==1) x-=2*EPS;
  //if(x==-1) x+=2*EPS;
  if(fabs(x)==1.) 
    throw Tantrum("pcmaps.cpp::LU->HG: the value at the domain boundary (would map to infinity)!\n");	
  y=sqrt(2.0)*inverf(x);
}
 else if (pcIn=="HG" && pcOut=="LU"){
   y=erf(x/sqrt(2.0));
 }
 
 else if (pcIn=="HG" && pcOut=="SW"){
   y=exp(out1+out2*x);
 }
 else if (pcIn=="SW" && pcOut=="HG"){
   y=(log(x)-in1)/in2;
 }
 
 else if (pcIn=="LU" && pcOut=="SW"){
   y=PCtoPC(PCtoPC(x,"LU",0,0,"HG",0,0),"HG",0,0,"SW",out1,out2);
 }
 else if (pcIn=="SW" && pcOut=="LU"){
   y=PCtoPC(PCtoPC(x,"SW",in1,in2,"HG",0,0),"HG",0,0,"LU",0,0);
 }
 
 else if (pcIn=="JB" && pcOut=="LU"){
   y=2.*betai(in1+1,in2+1,(x+1.)/2.)-1.;
 }
 
 else if (pcIn=="LU" && pcOut=="JB"){
   if(fabs(x)==1) y=x;
   else  
     y=rtbis_mod(PCtoPC,-1.,1.,EPS,x,"JB",out1,out2,"LU",in1,in2);
 }
 
 else if (pcIn=="JB" && pcOut=="HG"){
   y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"HG",0,0);
 }
 
 else if (pcIn=="HG" && pcOut=="JB"){
  y=PCtoPC(PCtoPC(x,"HG",0,0,"LU",0,0),"LU",0,0,"JB",out1,out2);
 }
 
 else if (pcIn=="JB" && pcOut=="SW"){
   y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"SW",out1,out2);
 }
 
 else if (pcIn=="SW" && pcOut=="JB"){
   y=PCtoPC(PCtoPC(x,"SW",in1,in2,"LU",0,0),"LU",0,0,"JB",out1,out2);
 }
 
 else if (pcIn=="LG" && pcOut=="LU"){
   y=2.*gammai(in1+1,x)-1.;
 }
 else if (pcIn=="LU" && pcOut=="LG"){
  y=rtbis_mod(PCtoPC,0.,MXEPS,EPS,x,"LG",out1,out2,"LU",in1,in2);
 }
 
 else if (pcIn=="LG" && pcOut=="HG"){
   y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"HG",0,0);
 }
 else if (pcIn=="HG" && pcOut=="LG"){
   y=PCtoPC(PCtoPC(x,"HG",0,0,"LU",0,0),"LU",0,0,"LG",out1,0);
 }
 
 else if (pcIn=="LG" && pcOut=="SW"){
   y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"SW",out1,out2);
 }
 else if (pcIn=="SW" && pcOut=="LG"){
  y=PCtoPC(PCtoPC(x,"SW",in1,in2,"LU",0,0),"LU",0,0,"LG",out1,0);
 }
 
 else if (pcIn=="LG" && pcOut=="JB"){
   y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"JB",out1,out2);
 }
 else if (pcIn=="JB" && pcOut=="LG"){
   y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"LG",out1,0);
 }
 else if (pcIn=="TG" && pcOut=="LU"){
   y=normcdf(x)/normcdf(in1)*2.-1.;
 }
 else if (pcIn=="LU" && pcOut=="TG"){
   //if(x==1) x-=2*EPS;
   //if(x==-1) x+=2*EPS;
   if(fabs(x)==1.) 
     throw Tantrum("pcmaps.cpp::LU->TG: the value at the domain boundary (would map to infinity)!\n");	
   y=invnormcdf(0.5*(x+1.)*normcdf(out1));
 }
 else if (pcIn=="TG" && pcOut=="HG"){
   y=PCtoPC(PCtoPC(x,"TG",in1,0,"LU",0,0),"LU",0,0,"HG",0,0);
 }
 else if (pcIn=="HG" && pcOut=="TG"){
   y=PCtoPC(PCtoPC(x,"HG",0,0,"LU",0,0),"LU",0,0,"TG",out1,0);
 }
 else if (pcIn=="TG" && pcOut=="LG"){
   y=PCtoPC(PCtoPC(x,"TG",in1,0,"LU",0,0),"LU",0,0,"LG",out1,0);
 }
 else if (pcIn=="LG" && pcOut=="TG"){
  y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"TG",out1,0);
 }
 else if (pcIn=="TG" && pcOut=="SW"){
   y=PCtoPC(PCtoPC(x,"TG",in1,0,"LU",0,0),"LU",0,0,"SW",out1,out2);
 }
 else if (pcIn=="SW" && pcOut=="TG"){
   y=PCtoPC(PCtoPC(x,"SW",in1,in2,"LU",0,0),"LU",0,0,"TG",out1,0);
 }
 else if (pcIn=="TG" && pcOut=="JB"){
   y=PCtoPC(PCtoPC(x,"TG",in1,0,"LU",0,0),"LU",0,0,"JB",out1,out2);
 }
 else if (pcIn=="JB" && pcOut=="TG"){
   y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"TG",out1,0);
 }
 /////////////////////////////////////////Roe-Baker
 else if (pcIn=="TG" && pcOut=="RB"){
   y=out2/(1.-PCtoPC(x,"TG",in1,0,"TG",out1,0));
 }
 else if (pcIn=="RB" && pcOut=="TG"){
   y=PCtoPC(1.-in2/x,"TG",in1,0,"TG",out1,0);
 }
 else if (pcIn=="RB" && pcOut=="LU"){
   y=PCtoPC(PCtoPC(x,"RB",in1,in2,"TG",in1,0),"TG",in1,0,"LU",0,0);
 }
 else if (pcIn=="LU" && pcOut=="RB"){
   if(fabs(x)==1.) 
     throw Tantrum("pcmaps.cpp::LU->RB: the value at the domain boundary (would map to infinity)!\n");	
   y=PCtoPC(PCtoPC(x,"LU",0,0,"TG",out1,0),"TG",out1,0,"RB",out1,out2);
 }
 else if (pcIn=="RB" && pcOut=="HG"){
   y=PCtoPC(PCtoPC(x,"RB",in1,in2,"LU",0,0),"LU",0,0,"HG",0,0);
 }
 else if (pcIn=="HG" && pcOut=="RB"){
   y=PCtoPC(PCtoPC(x,"HG",0,0,"LU",0,0),"LU",0,0,"RB",out1,out2);
 }
 else if (pcIn=="RB" && pcOut=="LG"){
   y=PCtoPC(PCtoPC(x,"RB",in1,in2,"LU",0,0),"LU",0,0,"LG",out1,0);
 }
 else if (pcIn=="LG" && pcOut=="RB"){
   y=PCtoPC(PCtoPC(x,"LG",in1,0,"LU",0,0),"LU",0,0,"RB",out1,out2);
 }
 else if (pcIn=="RB" && pcOut=="SW"){
   y=PCtoPC(PCtoPC(x,"RB",in1,in2,"LU",0,0),"LU",0,0,"SW",out1,out2);
 }
 else if (pcIn=="SW" && pcOut=="RB"){
   y=PCtoPC(PCtoPC(x,"SW",in1,in2,"LU",0,0),"LU",0,0,"RB",out1,out2);
 }
 else if (pcIn=="RB" && pcOut=="JB"){
   y=PCtoPC(PCtoPC(x,"RB",in1,in2,"LU",0,0),"LU",0,0,"JB",out1,out2);
 }
 else if (pcIn=="JB" && pcOut=="RB"){
   y=PCtoPC(PCtoPC(x,"JB",in1,in2,"LU",0,0),"LU",0,0,"RB",out1,out2);
 }
 else if (pcIn=="LU" && pcOut=="pdf"){
   Array2D<double> cdf ;
   read_datafileVS(cdf,"cdf.dat");
   linint( cdf, 0.5*(x+1) , y, 1 ) ;
 }
 else if (pcIn=="pdf" && pcOut=="LU"){
   Array2D<double> cdf ;
   read_datafileVS(cdf,"cdf.dat");
   linint( cdf, x , y, 0 ) ;
   y = y*2.0-1.0 ;
 }
 else if (pcIn=="pdf" && pcOut=="HG")
   y=PCtoPC(PCtoPC(x,"pdf",0,0,"LU",0,0),"LU",0,0,"HG",0,0);
 else if (pcIn=="HG" && pcOut=="pdf")
   y=PCtoPC(PCtoPC(x,"HG",0,0,"LU",0,0),"LU",0,0,"pdf",0,0);
 else if (pcIn=="pdf" && pcOut=="pdf"){
   y=x ;
 }
 ////////////////////////////////////////////
 else 
   throw Tantrum("pcmaps.cpp::Input-Output PC pair is not recognized!\n");
 

 
 return y;
}

double rtbis_mod(double func(double,string,double,double,string,double,double), const double x1, const double x2, const double xacc,double x, string pcIn, double in1, double in2, string pcOut, double out1, double out2)
{
	const int JMAX=400;
	int j;
	double dx,f,fmid,xmid,rtb;

	f=func(x1,pcIn,in1,in2,pcOut,out1,out2)-x;
	fmid=func(x2,pcIn,in1,in2,pcOut,out1,out2)-x;
	if (f*fmid > 0.0) {printf("Root must be bracketed for bisection in rtbis"); exit(1);}// made the inequality strict
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=0;j<JMAX;j++) {
                xmid=rtb+(dx *= 0.5);
		fmid=func(xmid,pcIn,in1,in2,pcOut,out1,out2)-x;
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	printf("Too many bisections in rtbis"); exit(1);
	return 0.0;
}

void PCtoPC(Array2D<double>& xx, string pcIn, double in1, double in2, Array2D<double>& yy, string pcOut, double out1, double out2)
{

  int n=xx.XSize();
  int m=xx.YSize();

  yy.Resize(n,m);

  for(int i=0;i<n;i++)
    for(int j=0;j<m;j++)
      yy(i,j)=PCtoPC(xx(i,j),pcIn,in1,in2,pcOut,out1,out2);

  return;
}

void linint( Array2D<double> &xydata, const double x, double &y )
{
  
  int n = xydata.XSize() ;

  if ( x < xydata(0,0) ) 
    y = 0.0 ;
  else if ( x > xydata(n-1,0) ) 
    y = 0.0 ;

  int i = 0 ;
  while ( ( x > xydata(i+1,0) ) && ( i < n-2 ) ) i++ ;

  y = xydata(i,1) + ( x - xydata(i,0) )
     *( xydata(i+1,1) - xydata(i,1) ) / ( xydata(i+1,0) - xydata(i,0) ) ;

  return ;

}

void linint( Array2D<double> &xydata, const double x, double &y, int col )
{
  
  assert(col==0 || col==1) ;

  int n = xydata.XSize() ;

  if ( x < xydata(0,col) ) 
    y = 0.0 ;
  else if ( x > xydata(n-1,col) ) 
    y = 0.0 ;

  int i = 0 ;
  while ( ( x > xydata(i+1,col) ) && ( i < n-2 ) ) i++ ;

  y = xydata(i,1-col) + ( x - xydata(i,col) )
     *( xydata(i+1,1-col) - xydata(i,1-col) ) / ( xydata(i+1,col) - xydata(i,col) ) ;

  return ;

}
