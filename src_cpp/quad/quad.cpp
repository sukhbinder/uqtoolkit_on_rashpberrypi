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
#include <cfloat>
#include "quad.h"

#include "error_handlers.h"
#include "combin.h"
#include "multiindex.h"
#include "gq.h"
#include "datafiles.h"
#include "arraytools.h"


Quad::Quad(Array1D<string>& grid_types, char *fs_type, Array1D<int>& param,Array1D<double>& alphas, Array1D<double>& bettas)//Default
{
  //  indexing_=false;  
  // Initialize the dimensionality, grid types and boundaries
  this->ndim_=grid_types.Length();
  aa_.Resize(this->ndim_,-1.e0);
  bb_.Resize(this->ndim_,1.e0);
  this->growth_rules_.Resize(this->ndim_,0);

  this->grid_types_=grid_types; // char* to string?
  this->fs_type_=fs_type;
  
  // Initialize parameters
  this->alphas_=alphas;
  this->betas_=bettas;

  // Set the quadrature level
  maxlevel_=param(0);
  param_=param;

  this->init();
 

}

Quad::Quad(char *grid_type, char *fs_type,int ndim, int param,double alpha, double betta)//Default
{
  
  // Initialize the dimensionality, grid types and boundaries
  this->ndim_=ndim;
  aa_.Resize(ndim,-1.e0);
  bb_.Resize(ndim,1.e0);
  this->grid_type_=grid_type;
  this->fs_type_=fs_type;
  
  // Initialize parameters
  this->alpha_=alpha;
  this->beta_=betta;

  this->alphas_.Resize(this->ndim_,alpha);
  this->betas_.Resize(this->ndim_,betta);
  this->growth_rules_.Resize(this->ndim_,0);
  this->grid_types_.Resize(this->ndim_,grid_type);
  this->param_.Resize(this->ndim_,param);
  // Set the quadrature level
  maxlevel_=param;
  
  this->init();
 

}

void Quad::init()
{
  
  // Require proper indexing for full quadrature
  if (this->fs_type_=="full"){
    //indexing_=true;
    if(this->maxlevel_==0)
      throw Tantrum("Quad::Quad(): 'full' does not make sense with parameter 0.");
  }
  
  // Require (or not) indexing for sparse quadrature, as well as set the proper growth rules
  else if (this->fs_type_=="sparse"){
    if(this->ndim_==1)
      throw Tantrum("Quad::Quad(): 'sparse' does not make sense in 1D, use 'full' instead.");
    
    for (int i=0;i<this->ndim_;i++){
      if (this->grid_types_(i)=="CC" or this->grid_types_(i)=="NC"){
	this->growth_rules_(i)=0;
	//indexing_=true;
    }
      else if (this->grid_types_(i)=="LU"){
	this->growth_rules_(i)=0;
	//      indexing_=false;
    } 
    else if (this->grid_types_(i)=="NCO" or this->grid_types_(i)=="CCO"){
	this->growth_rules_(i)=1;
	//indexing_=true;
    }
    else if (this->grid_types_(i)=="HG"){
	this->growth_rules_(i)=1;
	//      indexing_=false;
    }
    else if (this->grid_types_(i)=="JB"){
      	this->growth_rules_(i)=0;
	//      indexing_=false;
    } 
    else if (this->grid_types_(i)=="GLG"){
      	this->growth_rules_(i)=0;
      //indexing_=false;
    } 
    else if (this->grid_types_(i)=="SW"){
	this->growth_rules_(i)=0;
	//      indexing_=false;
    } 
    else if (this->grid_types_(i)=="pdf"){
	this->growth_rules_(i)=0;
      // indexing_=false;
    } 
    else {
      throw Tantrum("Quad::Quad(): Grid type unrecognized! Options are 'CC','CCO','NC','NCO','LU', 'HG', 'JB', 'GLG', 'SW' or 'pdf'");
    }
	  
    }
  }
  
  else 
    throw Tantrum("Quad::Quad(): Either 'full' or 'sparse' should be specified");

  return;

}

void Quad::SetDomain(Array1D<double>& aa, Array1D<double>& bb)
{
  // Dimensionality check
  if ( (int) aa.XSize() != ndim_ or (int) bb.XSize() != ndim_ )
    throw Tantrum("Quad::SetDomain(): Dimension error!");

  aa_=aa;
  bb_=bb;

  return;
}


void Quad::SetDomain(Array1D<double>& aa)
{
  // Dimensionality check
  if ( (int) aa.XSize() != ndim_ )
    throw Tantrum("Quad::SetDomain(): Dimension error!");
  
  aa_=aa;
  
  return;
}

void Quad::SetRule(Array2D<double>& q, Array1D<double>& w)
{
  // Set a rule
  this->SetQdpts(q);
  this->SetWghts(w);

  return;
}

// void Quad::SetRule(Array2D<double>& q, Array1D<double>& w, Array2D<int>& ind)
// {
//   // Set a rule
//   this->SetQdpts(q);
//   this->SetWghts(w);
//   this->SetIndices(ind);

//   return;
// }

// void Quad::SetRule(Array2D<double>& q, Array1D<double>& w, Array2D<int>& ind, int param)
// {
//   // Set a rule and a level parameter
//   this->SetQdpts(q);
//   this->SetWghts(w);
//   this->SetIndices(ind);

//   this->SetLevel(param);

//   return;
// }

// void Quad::GetRule(Array2D<double>& q, Array1D<double>& w, Array2D<int>& ind)
// {
//   // Get the rule
//   this->GetQdpts(q);
//   this->GetWghts(w);
//   this->GetIndices(ind);

//   return;
// }



void Quad::GetRule(Array2D<double>& q, Array1D<double>& w)
{
  // Get the rule
  this->GetQdpts(q);
  this->GetWghts(w);

  return;
}



void Quad::MultiplyTwoRules(QuadRule *rule1,QuadRule *rule2,QuadRule *rule_prod)
{ 

  // Get the sizes
  int n1=rule1->qdpts.XSize();
  int n2=rule2->qdpts.XSize();
  
  int d1=rule1->qdpts.YSize();
  int d2=rule2->qdpts.YSize();

  // Compute the sizes of the product rule
  int n12=n1*n2;
  int d12=d1+d2;

  // Resize the product rule containers
  rule_prod->qdpts.Resize(n12,d12,0.e0);
  rule_prod->wghts.Resize(n12,0.e0);
  //rule_prod->indices.Resize(n12,d12,0);

  
  for(int in2=0;in2<n2;in2++){
    for(int in1=0;in1<n1;in1++){
      
      // Quadrature points and indices are 'mesh'-ed
      for(int id1=0;id1<d1;id1++){
        rule_prod->qdpts(in1+in2*n1,id1)=rule1->qdpts(in1,id1);
      }
      for(int id2=0;id2<d2;id2++){
        rule_prod->qdpts(in1+in2*n1,id2+d1)=rule2->qdpts(in2,id2);
      }
      
      // Weights are multiplied
      rule_prod->wghts(in1+in2*n1)=rule1->wghts(in1)*rule2->wghts(in2);

     
	

    }
  }
  
  



  return;
}


void Quad::MultiplyManyRules(int nrules, QuadRule *rules, QuadRule *rule_prod)
{

  // Working rules
  QuadRule rule_1d;
  QuadRule rule_cur;
  QuadRule rule_prod_h;
  
    
  
  for(int i=0;i<nrules;i++){

      rule_1d=rules[i];
   
      // Recursively multiply all the rules
      if(i==0)
        rule_prod_h=rule_1d;
      else
        this->MultiplyTwoRules(&rule_cur,&rule_1d,&rule_prod_h);
      
      rule_cur=rule_prod_h;
      
  }

  
  rule_prod->qdpts=rule_prod_h.qdpts;
  rule_prod->wghts=rule_prod_h.wghts;
  
  return;
}

void Quad::SubtractTwoRules(QuadRule *rule1,QuadRule *rule2,QuadRule *rule_sum)
{

  for(int i=0;i<(int) rule2->wghts.XSize();i++)
    rule2->wghts(i) *= -1.;

  this->AddTwoRules(rule1,rule2,rule_sum);

  return;
}

void Quad::AddTwoRules(QuadRule *rule1,QuadRule *rule2,QuadRule *rule_sum)
{
  // Get the sizes
  int n1=rule1->qdpts.XSize();
  int n2=rule2->qdpts.XSize();

  int d1=rule1->qdpts.YSize();
  int d2=rule2->qdpts.YSize();

  // Sanity check
  if(d1!=d2){
    printf("quad.cpp:only rules of same dimensionality can be added to each other! %d %d\n",d1,d2);
    exit(1);
  }
  int d=d1;
  // The size of the full sum
  int n12=n1+n2;
  

  // Work arrays
  Array2D<double> t_qdpts(n12,d,0.e0);
  Array1D<double> t_wghts(n12,0.e0);


  for(int in1=0;in1<n1;in1++){
    for(int id1=0;id1<d;id1++){
    t_qdpts(in1,id1)=rule1->qdpts(in1,id1);
    }
    t_wghts(in1)=rule1->wghts(in1);
  }
  
  //int n2d=0;
  
  // Bookkeeping if indexing is true, do not have duplicate points
  for(int in2=0;in2<n2;in2++){ 
      for(int id1=0;id1<d;id1++){
        t_qdpts(n1+in2,id1)=rule2->qdpts(in2,id1);
      }
      t_wghts(n1+in2)=rule2->wghts(in2);
  }
  
  // Resize the sum rule container
  rule_sum->qdpts.Resize(n1+n2,d,0.e0);
  rule_sum->wghts.Resize(n1+n2,0.e0);
  
  // MAYBE COMPRESS HERE?
 
  // Set the values of the sum rule
  for(int i=0;i<n1+n2;i++){
    for(int id1=0;id1<d;id1++){
      rule_sum->qdpts(i,id1)=t_qdpts(i,id1);
    }
    rule_sum->wghts(i)=t_wghts(i);
  }

  return;
}

/**********************************************************************************/
/**********************************************************************************/
void Quad::create1DRule(string gridtype,Array1D<double>& qdpts,Array1D<double>& wghts,int ngr, double a, double b)
{
  if (gridtype=="CC"){
    this->create1DRule_CC(qdpts,wghts,ngr,a,b);
  }
  else if (gridtype=="NC"){
    this->create1DRule_NC(qdpts,wghts,ngr,a,b);
  }
  else if (gridtype=="NCO"){
    this->create1DRule_NCO(qdpts,wghts,ngr,a,b);
  }
  else if (gridtype=="CCO"){
    this->create1DRule_CCO(qdpts,wghts,ngr,a,b);
  }
  else  if (gridtype=="LU" or gridtype=="LU_N"){
    this->create1DRule_LU(qdpts,wghts,ngr,a,b);
  }
  else  if (gridtype=="HG"){
    this->create1DRule_HG(qdpts,wghts,ngr);
  }
  else  if (gridtype=="JB"){
    this->create1DRule_JB(qdpts,wghts,ngr,a,b);
  }
  else  if (gridtype=="GLG"){
    this->create1DRule_GLG(qdpts,wghts,ngr);
      }
  else  if (gridtype=="SW"){
    this->create1DRule_SW(qdpts,wghts,ngr);
  }
  else  if (gridtype=="pdf"){
    this->create1DRule_pdf(qdpts,wghts,ngr,a,b);
  }
  else  if (gridtype=="GP3"){
    this->create1DRule_pdf(qdpts,wghts,ngr,a,b);
  }
  
  else 
    throw Tantrum("Quad::create1DRule(): Grid type unrecognized! Options are 'CC','CCO','NC','NCO','LU', 'HG', 'JB', 'GLG', 'SW', or 'GP3' ");
  
  return;
}
/**********************************************************************************/
/**********************************************************************************/

void Quad::create1DRule_LU(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b)
{
  qdpts.Resize(ngr,0.e0);
  wghts.Resize(ngr,0.e0);
  

  // Work arrays
  Array1D<double> endpts(2,0.e0); // real array of length two with zeroed values, needed to pass to gaussqC
  Array1D<double> bwork(ngr,0.e0);
  Array1D<double> qdpts_1d(ngr,0.e0);
  int kind=1;
  double alpha=0.0,beta=0.0;

  gq( kind, alpha, beta, qdpts_1d, wghts );

  // Rescale and index
  for(int i=0; i<ngr;i++){ 
    qdpts(i)=a+(b-a)*(qdpts_1d(i)+1.)/2.;
    wghts(i) *= (b-a)/4.; //since the integral is with respect to pdf=1/2
  }    
  
  return;
}
/**********************************************************************************/
void Quad::create1DRule_HG(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr)
{
  qdpts.Resize(ngr,0.e0);
  wghts.Resize(ngr,0.e0);
  

  // Work arrays
  Array1D<double> endpts(2,0.e0); // real array of length two with zeroed values, needed to pass to gaussqC
  Array1D<double> bwork(ngr,0.e0);
  Array1D<double> qdpts_1d(ngr,0.e0);
  int kind=4;
  double alpha=0.0,beta=0.0;

  const double pi = 4.e0*atan(1.e0);
  double spi = sqrt(pi);
  double fac = sqrt(2.0);

  gq( kind, alpha, beta, qdpts_1d, wghts );
  
  // Rescale and index
  for(int i=0; i<ngr;i++){ 
    qdpts(i)=qdpts_1d(i)*fac;
    wghts(i) /= spi;
  }    
  
  return;
}
/**********************************************************************************/
void Quad::create1DRule_NC(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b)
{
  qdpts.Resize(ngr,0.e0);
  wghts.Resize(ngr,0.e0);

  
  if (ngr==1){
    qdpts(0)=0.0;
    wghts(0)=2.0;
  }
  else{
    for (int i=0;i<ngr;i++)
      qdpts(i) =  -1.+2.*double(i)/ double(ngr-1.);
    
    Array1D<double> exact(ngr,0.e0);
    for (int i=0;i<ngr;i=i+2)
      exact(i)=2./(i+1.);
    vandermonde_gq(qdpts,wghts,exact);
    
   }


  // Rescale and index
  for (int i=0;i<ngr;i++){
    qdpts(i) = (b+a)/2.+ qdpts(i)* (b-a)/2.;
    wghts(i) = wghts(i) *(b-a)/4.;  //since the integral is with respect to pdf=1/2
  
    
  }

    return;
}
/**********************************************************************************/
void Quad::create1DRule_NCO(Array1D<double>& qdpts,Array1D<double>& wghts,int ngr, double a, double b)
{
  qdpts.Resize(ngr,0.e0);
  wghts.Resize(ngr,0.e0);
  
  for (int i=0;i<ngr;i++)
    qdpts(i) =  -1.+2.*(i+1.)/ double(ngr+1.);

  
  
  Array1D<double> exact(ngr,0.e0);
  for (int i=0;i<ngr;i=i+2)
    exact(i)=2./(i+1.);
  vandermonde_gq(qdpts,wghts,exact);
  
  // Rescale and index
  for (int i=0;i<ngr;i++){
    qdpts(i) = (b+a)/2.+ qdpts(i)* (b-a)/2.;
    wghts(i) = wghts(i) *(b-a)/4.;  //since the integral is with respect to pdf=1/2

    
  }
  return;
}
/**********************************************************************************/
void Quad::create1DRule_CC(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b)
{
  qdpts.Resize(ngr,0.e0);
  wghts.Resize(ngr,0.e0);
  
  if (ngr==1){
    qdpts(0)=0.0;
    wghts(0)=2.0;
  }
  
  else {
    
    double pi=4.*atan(1.);
    double theta,f;
    
    
    for (int i=0;i<ngr;i++)
      qdpts(i) =  cos( (double) (ngr-1-i) * pi / (double) (ngr-1) );
    
    for (int i=0;i<ngr;i++){
      theta = (double) i * pi / (double) (ngr-1);
      wghts(i) = 1.0;
      
      for (int j=1;j<=( ngr - 1 ) / 2; j++ ){
        if ( 2 * j == ( ngr - 1 ) )
          f = 1.0;
        else
          f = 2.0;
        
        wghts(i)-= f*cos ( 2.0 * ( double ) j * theta ) / ( double ) ( 4 * j * j - 1 );
      }

    }
    
    wghts(0)  /=  (double) ( ngr - 1 );
    for ( int i = 1; i < ngr - 1; i++ )
      wghts(i) *= 2.0 / ( double ) ( ngr - 1 );
    
    wghts(ngr-1) *= 1.0 / ( double ) ( ngr - 1 );
  }
  
  // Rescale and index
  for (int i=0;i<ngr;i++){
    qdpts(i) = (b+a)/2.+ qdpts(i)* (b-a)/2.;
    wghts(i) = wghts(i) *(b-a)/4.;  //since the integral is with respect to pdf=1/2
  
    
  }
  
  return;
}
/**********************************************************************************/
void Quad::create1DRule_CCO(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b)
{
  qdpts.Resize(ngr,0.e0);
  wghts.Resize(ngr,0.e0);

  double pi=4.*atan(1.);
  double theta,f;
  
  for (int i=0;i<ngr;i++)
    qdpts(i) =  cos( (double) (ngr-i) * pi / (double) (ngr+1) );

  
  for (int i=0;i<ngr;i++){
    theta = (double) (ngr-i) * pi / (double) (ngr+1);
    
    wghts(i) = 1.0;
    
    for (int j=1;j<=( ngr - 1 ) / 2; j++ ){
               wghts(i)-= 2.0*cos ( 2.0 * ( double ) j * theta ) / ( double ) ( 4 * j * j - 1 );
    }
    if(ngr%2==1)
      f=ngr;
    else
      f=ngr-1;
    
    wghts(i)-=cos((f+1.)*theta) / f;
    
  }
  
  for ( int i = 0; i < ngr; i++ )
    wghts(i) *= 2.0 / ( double ) ( ngr + 1 );
  
  // Rescale and index
  for (int i=0;i<ngr;i++){
    qdpts(i) = (b+a)/2.+ qdpts(i)* (b-a)/2.;
    wghts(i) = wghts(i) *(b-a)/4.;  //since the integral is with respect to pdf=1/2

 
    
  }

  return;
}
/**********************************************************************************/
void Quad::create1DRule_JB(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b)
{
  qdpts.Resize(ngr,0.e0);
  wghts.Resize(ngr,0.e0);

  // The norm
  double mu0=  pow(2.,alpha_+beta_+1.)*beta(alpha_+1.,beta_+1.);

  gq (5, alpha_, beta_, qdpts, wghts ) ;

  // Rescale and index
  for(int i=0; i<ngr;i++){ 
    qdpts(i)=a+(b-a)*(qdpts(i)+1.)/2.;
    wghts(i) *= (b-a)/(2.*mu0);  
  }    
  
  return;
}
/**********************************************************************************/
void Quad::create1DRule_GLG(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr)
{
  qdpts.Resize(ngr,0.e0);
  wghts.Resize(ngr,0.e0);
  
  gq (6, alpha_, 0.0, qdpts, wghts ) ;

  // Indexing
  for(int i=0; i<ngr;i++){ 
    wghts(i) /= exp(lgamma(alpha_+1));
  }    
  
return;
}
/**********************************************************************************/
void Quad::create1DRule_SW(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr)
{

  qdpts.Resize(ngr,0.e0);
  wghts.Resize(ngr,0.e0);
  
  // Work arrays
  Array1D<double> al(ngr,0.e0);
  Array1D<double> be(ngr, 0.e0);
  
  // \todo check alpha non-negative?    
  double ee=exp(beta_*beta_/2.);
     double eesq = ee*ee ;
  
  for(int i=0; i<ngr;i++){
      al(i)=exp(alpha_)*pow(ee,2.e0*i-1.e0)*((eesq+1.0)*pow(eesq,i)-1.e0);
      be(i)=exp(2.e0*alpha_)*pow(eesq,3.e0*i-2.e0)*(pow(ee,2.e0*(double)(i))-1.e0);
  }
  
  // The norm
  double mu0=1.;
    
  // Computes the quadrature given the recursion coefficients
  gq_gen(al,be,mu0,qdpts, wghts) ;

  // Rescale?
  
return;
}
/**********************************************************************************/
void Quad::create1DRule_pdf(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b)
{

  qdpts.Resize(ngr,0.e0);
  wghts.Resize(ngr,0.e0);

  // Work arrays
  Array1D<double> al(ngr,0.e0);
  Array1D<double> be(ngr, 0.e0);
  
  // Read the recursion coefficients from the data file
  Array2D<double> albe;
  read_datafileVS(albe,"ab.dat");
  if ((int)albe.XSize()<ngr)
  {
    printf("Quad::create1DRule_pdf() : ngr=%d, rows(ab.dat)=%d !\n",ngr,(int) albe.XSize()) ;
    throw Tantrum("The input coefficient file 'ab.dat' has fewer rows than needed");
  }
  for(int i=0; i<ngr;i++){
    al(i)=albe(i,0);
    be(i)=albe(i,1);
  }
  
  // \todo This is problem specific!!!
  double mu0=1.0; //0.5+0.5*erf(3.2/sqrt(2.));
  
  // Computes the quadrature given the recursion coefficients
  gq_gen(al,be,mu0,qdpts, wghts) ;

  
  return;
}
/**********************************************************************************/
void Quad::create1DRule_GP3(Array1D<double>& qdpts,Array1D<double>& wghts, int ngr, double a, double b)
{

  int nqp[7]={1,3,7,15,31,63,127};
 
  if ((ngr<1) || (ngr>7)) {
    printf("Quad::create1DRule_GP3() : ngr=%d !\n",ngr) ;
    throw Tantrum("The above Gauss-Patterson rule is not available!");
  }

  qdpts.Resize(nqp[ngr-1],0.e0);
  wghts.Resize(nqp[ngr-1],0.e0);

  /* Read the quadrature points and weights from data file */
  stringstream fname;
  fname << "gp3rules/gp3_o" << ngr << "_p.dat";

  Array2D<double> din;
  read_datafileVS(din,fname.str().c_str());
  if ((int)din.XSize()!=nqp[ngr])
  {
    printf("Quad::create1DRule_GP3() : ngr=%d, nqp=%d vs %d !\n",ngr,(int) din.XSize(), nqp[ngr]) ;
    throw Tantrum("The the number of quadrature points does not match the expected value");
  }
  for(int i=0; i<(int)din.XSize();i++) qdpts(i) = din(i,0);
  
  fname.str("gp3rules/gp3_o");
  fname << ngr << "_w.dat";

  read_datafileVS(din,fname.str().c_str());
  if ((int)din.XSize()!=nqp[ngr])
  {
    printf("Quad::create1DRule_GP3() : ngr=%d, nqp=%d vs %d !\n",ngr,(int) din.XSize(), nqp[ngr]) ;
    throw Tantrum("The the number of quadrature weights does not match the expected value");
  }
  for(int i=0; i<(int)din.XSize();i++) wghts(i) = din(i,0);
  
  // Rescale
  for(int i=0; i<ngr;i++){ 
    qdpts(i)  = a+(b-a)*(qdpts(i)+1.)/2.;
    wghts(i) *= (b-a)/2.0;  
  }    
  
  return;
}


/**********************************************************************************/
/**********************************************************************************/


void Quad::SetRule()
{

  // Either the full product of the 1D rules
  if (this->fs_type_=="full"){

    QuadRule rule_1d;
    QuadRule rule_cur;
    QuadRule rule_prod;

    for(int id=0;id<this->ndim_;id++){

      Array1D<double> qdpts_1d;
 this->create1DRule(this->grid_types_(id),qdpts_1d,rule_1d.wghts,param_(id),aa_(id),bb_(id));

      array1Dto2D(qdpts_1d,rule_1d.qdpts);

      if(id==0)
        rule_prod=rule_1d;
      else
        this->MultiplyTwoRules(&rule_cur,&rule_1d,&rule_prod);
      
      rule_cur=rule_prod;

    }
     
    rule_=rule_prod;
  }

  
  // or multiply rules in a specific way and combine them to obtain sparse rule
  else if (this->fs_type_=="sparse"){
    
    this->SetLevel(-1);

    // Incremental buildup of levels
    for(int il=0;il<=maxlevel_;il++)
      this->nextLevel();
    
  }
  
  else 
    throw Tantrum("quad.cpp: unknown rule type.");
    
  return;
}

/// \todo Comment this better
void Quad::nextLevel()
{
  
  this->SetLevel(nlevel_+1);
  

   QuadRule rule_level;
   QuadRule rule_cur;
   QuadRule rule_total;
   
   Array2D<int> multiIndexLevel;
   this->getMultiIndexLevel(multiIndexLevel,nlevel_,ndim_);
   
   int nMultiIndicesLevel=multiIndexLevel.XSize();
   
   Array2D<int> multiIndexLevel_npts(nMultiIndicesLevel,ndim_,0);
   
   for(int j=0;j<nMultiIndicesLevel;j++){
     
     QuadRule* rules;
     QuadRule* rules_1;
     QuadRule* srules;
     
     rules=new QuadRule[ndim_];
     rules_1=new QuadRule[ndim_];
     srules=new QuadRule[ndim_];
     
     
     for(int id=0;id<ndim_;id++){
       // multiIndexLevel(j,id)=multiIndex(j+levelCounter_start,id);
       int npts,npts_1;
       
       if (this->growth_rules_(id)==0){
        if ( multiIndexLevel(j,id)==0){
          npts = 1;
          npts_1=0;
        }
        else if ( multiIndexLevel(j,id)==1){
          npts=3;
          npts_1=1;
        }
        else  {
          npts=(int) pow(2,multiIndexLevel(j,id))+1;
          npts_1=(int) pow(2,multiIndexLevel(j,id)-1)+1;
        }
      }
       else if (this->growth_rules_(id)==1){
        npts=(int) pow(2,multiIndexLevel(j,id)+1)-1;
        npts_1=(int) pow(2,multiIndexLevel(j,id))-1;
      }
      Array1D<double> qdpts_1d;
Array1D<int> indices_1d;
 this->create1DRule(this->grid_types_(id),qdpts_1d,rules[id].wghts,npts,aa_(id),bb_(id));
      array1Dto2D(qdpts_1d,rules[id].qdpts);
      
      if(npts_1>0){
        this->create1DRule(this->grid_types_(id),qdpts_1d,rules_1[id].wghts,npts_1,aa_(id),bb_(id));
        array1Dto2D(qdpts_1d,rules_1[id].qdpts);
        this->SubtractTwoRules(&rules[id],&rules_1[id], &srules[id]);
      }
      else
        srules[id]=rules[id];
        
     }//end of id loop

     QuadRule rule_temp;
     
     this->MultiplyManyRules(ndim_,srules,&rule_temp);
       
     if(j==0){
       rule_level=rule_temp;
     }
     else {
       this->AddTwoRules(&rule_cur,&rule_temp,&rule_level);
     }
     rule_cur=rule_level;
     
     delete []rules;
     delete []rules_1;
     delete []srules;
   }
   
   if (nlevel_==0)
     rule_total=rule_level;
   
   else {
     this->AddTwoRules(&rule_,&rule_level,&rule_total);
     // rule_total=*rule_total_p;
   }
   
   rule_=rule_total;
 
   
   this->compressRule();

 
 return;
}




void Quad::getMultiIndexLevel(Array2D<int>& multiIndexLevel, int level,int ndim)
{
  


  int iup=0;
  
  int nup_level=choose(ndim+level-1,level);

 

  multiIndexLevel.Resize(nup_level,ndim,0);
  
  


    if(ndim==1)
      multiIndexLevel(0,0)=level;
    else{

      for (int first=level;first>=0;first--){
        Array2D<int> theRest;
        getMultiIndexLevel(theRest,level-first,ndim-1);
        
        for(int j=0;j<(int)theRest.XSize();j++){
          multiIndexLevel(iup,0)=first;
          for(int id=1;id<ndim;id++){
            multiIndexLevel(iup,id)=theRest(j,id-1);
        }
          iup++;
        }
        
      }

    }


   



  return;
}


// void Quad::ComputeCurrIndicesAtNextLevel(Array2D<int>& indAtNextLevel)
// {
//   int ngrid=grid(nlevel_-1);
//   int nqdpts=rule_.qdpts.XSize();

//  indAtNextLevel.Resize(nqdpts,ndim_,0);
  
//   Array1D<int> projIndex(ndim_,0);


//   for(int i=0;i<nqdpts;i++){
//     for(int id=0;id<this->ndim_;id++)
// projIndex(id)=rule_.indices(i,id);
     

//     for(int id=ndim_-1;id>=0;id--){
//       int pi_next;
//       if(this->growth_rules_(id)==0)
//         {
//           pi_next=2*projIndex(id);
//           if(ngrid==1) pi_next=1;
//         }
//       else if(this->growth_rules_(id)==1)
//         pi_next=2*projIndex(id)+1;
      
//     indAtNextLevel(i,id)=pi_next;
//     }
    
//   }
  
//   return;
// }

void Quad::compressRule()
{
  Array1D<int> ind;
  int nqdpts=this->rule_.qdpts.XSize();
  for(int iq=0;iq<nqdpts;iq++)
    for(int id=0;id<this->ndim_;id++)
      if(fabs(this->rule_.qdpts(iq,id))<1.e-15)
	   this->rule_.qdpts(iq,id)=0.e0;


  Array1D<int> newInd;
 Array1D<int> oldInd;
   for(int i=0;i<nqdpts;i++) oldInd.PushBack(i);

   //  for(int id=this->ndim_-1;id>=0;id--){
   //shell_sort_incr(this->rule_.qdpts,id,newInd,oldInd);
   // }
shell_sort_all(this->rule_.qdpts,newInd,oldInd);
  Array1D<double> tmp;
  tmp=this->rule_.wghts;
  for(int i=0;i<nqdpts;i++){
    this->rule_.wghts(i)=tmp(oldInd(i));
  }
  write_datafile(this->rule_.qdpts,"qdpts.s.dat");

  Array2D<double> qdpts_tmp(1,this->ndim_,0.e0);
  for(int id=0;id<this->ndim_;id++) qdpts_tmp(0,id)=this->rule_.qdpts(0,id);
  Array1D<double> wghts_tmp(1,this->rule_.wghts(0));
  ind.Resize(1,oldInd(0));
  int iqq=0;
  for(int iq=1;iq<nqdpts;iq++){
    
    bool issame=true;
   for(int id=0;id<this->ndim_;id++){
     if(fabs(this->rule_.qdpts(iq,id)-this->rule_.qdpts(iq-1,id))>1.e-15){
       issame=false;
       break;
     }

   }	
   if(!issame){
     Array1D<double> qdpts1d(this->ndim_,0.e0);
     for(int id=0;id<this->ndim_;id++) qdpts1d(id)=this->rule_.qdpts(iq,id);
     paddMatRow(qdpts_tmp,qdpts1d);
     wghts_tmp.PushBack(this->rule_.wghts(iq));
     ind.PushBack(oldInd(iq));
     iqq++;
   }
   else{
     wghts_tmp(iqq) += this->rule_.wghts(iq);
   }
     

  }
  
  this->rule_.qdpts=qdpts_tmp;
  this->rule_.wghts=wghts_tmp;
 	


  return;
}

// int Quad::grid(int il,int growthrule)
// {
//   // Get the number of grid points given the level and using the growth rule
//   int gr;
//   if(growthrule==0){
//     gr=(int) pow(2,il)+1;
//     if(il==0)
//       gr=1;
//   }
//   else if(growthrule==1)
//     gr=(int) pow(2,il+1)-1;
  
//   return gr;
// }

// int Quad::grid(int il)
// {
//   // Get the number of grid points given the level and using the growth rule
//   int gr;
//   if(growth_rule_==0){
//     gr=(int) pow(2,il)+1;
//     if(il==0)
//       gr=1;
//   }
//   else if(growth_rule_==1)
//     gr=(int) pow(2,il+1)-1;
  
//   return gr;
// }

