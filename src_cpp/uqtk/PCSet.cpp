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
/*! \file PCSet.cpp
*/
#include "PCSet.h"
#include "PCBasis.h"
#include <math.h>
#include <sstream>
#include "uqtkconfig.h"
#include "slatec.h"
#include "quad.h"
#include "multiindex.h"
#include "minmax.h"


// Static members need to be pre-declared
int PCSet::next_index_ = 0;
PCSet::OMap_t *PCSet::omap_ = NULL;


PCSet::PCSet(const string sp_type, const int order, const int n_dim, const string pc_type, const double alpha, const double betta): 
  spType_(sp_type), pcType_(pc_type), order_(order), nDim_(n_dim), alpha_(alpha), beta_(betta)
{
  SetVerbosity(0);

  nPCTerms_=computeMultiIndex(nDim_,order_, multiIndex_);

  maxOrdPerDim_.Resize(this->nDim_,this->order_);
  
  Initialize("TotalOrder");

  return;
}

PCSet::PCSet(const string sp_type, const Array1D<int>& maxOrders, const int n_dim, const string pc_type, const double alpha, const double betta): 
  spType_(sp_type), pcType_(pc_type), maxOrders_(maxOrders), nDim_(n_dim), alpha_(alpha), beta_(betta)
{
  SetVerbosity(0);

  nPCTerms_=computeMultiIndexHDMR(nDim_,maxOrders_, multiIndex_);

  order_=maxOrders_(maxIndex(maxOrders_));
  maxOrdPerDim_.Resize(this->nDim_,this->order_);
  
  Initialize("HDMR");

  return;

}

PCSet::PCSet(const string sp_type, const Array2D<int>& customMultiIndex, const string pc_type, const double alpha, const double betta): 
  spType_(sp_type), pcType_(pc_type), nDim_(customMultiIndex.YSize()), alpha_(alpha), beta_(betta)
{

  SetVerbosity(0);

  order_=0;
  this->multiIndex_=customMultiIndex;

  this->nPCTerms_=multiIndex_.XSize();

  for(int i=0;i<this->nPCTerms_;i++){
    int sum=0;
    for(int j=0;j<this->nDim_;j++)
      sum+=customMultiIndex(i,j);
    
    if (order_<=sum)
      order_=sum;
  }

  this->ComputeMaxOrdPerDim();

  Initialize("Custom");

  return;

} 


PCSet::~PCSet()
{
  // Remove this object from the object map
  OMap_t::iterator me = PCSet::omap_->find(this->my_index_);
  if(PCSet::omap_ && (me != PCSet::omap_->end())) PCSet::omap_->erase(me);

  // Destroy the PC basis object
  delete p_basis_;
}


void PCSet::ComputeMaxOrdPerDim()
{
  maxOrdPerDim_.Resize(nDim_,0);
  for(int id=0;id<this->nDim_;id++){
    for(int ipc=0;ipc<this->nPCTerms_;ipc++){
      if(multiIndex_(ipc,id)>maxOrdPerDim_(id))
	maxOrdPerDim_(id)=multiIndex_(ipc,id);
    }
  }

  return;
}


void PCSet::Initialize(const string ordertype)
{
  // Echo the settings if desired
  if(uqtkverbose_>0){
    cout << endl;
    cout << "-----------------------------------------------------" << endl ;
    cout << "Initializing PCSet class with the following settings:" << endl ;
    cout << "# Dim      : " << nDim_   << endl;
    if(!strcmp(ordertype.c_str(),"HDMR")){
      cout << "Customized order for each HDMR-dimensionality: "  << endl;
      for(int i=0;i<(int)maxOrders_.XSize();i++)
      cout << i << "-dim order " << maxOrders_(i) << ", ";
      cout << endl;
    }
    else if(!strcmp(ordertype.c_str(),"Custom")){
      
      cout << "Customized multiIndex with "  << multiIndex_.XSize() << " basis terms."<< endl;
      if(uqtkverbose_>1){
        for(int i=0;i<(int)multiIndex_.XSize();i++){
          for(int j=0;j<(int)multiIndex_.YSize();j++){
           cout << multiIndex_(i,j) << " ";
          }
        cout << endl;
        }
        cout << endl;
      }
    }
    else
      cout << "Order      : " << order_  << endl;
    cout << "Basis type : " << pcType_ << endl;
    cout << "-----------------------------------------------------" << endl ;
  }

  // Add this object to the static object map with
  // an integer handle that will be used to retrieve this object
  // as needed when it is called from Fortran.
  if(! PCSet::omap_) PCSet::omap_ = new PCSet::OMap_t ;
  this->my_index_ = PCSet::next_index_++;
  (*PCSet::omap_)[my_index_] = this;

  // Initialize 1d basis class
  maxorddim_=maxOrdPerDim_(maxIndex(maxOrdPerDim_));
  p_basis_ = new PCBasis(pcType_, alpha_, beta_,maxorddim_);
  pcType_=p_basis_->GetPCType();


  // Get the norms of the multi-D basis functions
  this->EvalNormSq(psiSq_);
 
 
  
  if (!strcmp(spType_.c_str(),"NISP")){
    // Get the basis values at quadrature points (MAYBE have NISP no quad)
    this->InitNISP();
  }
  else if (!strcmp(spType_.c_str(),"NISPnoq")){
  }
  else if (!strcmp(spType_.c_str(),"ISP")){
    if(!strcmp(ordertype.c_str(),"TotalOrder")){
      this->InitISP();
    }
    else{
      throw Tantrum("Intrusive implementation is only avaliable for TotalOrder multiindex construction");
    }
  }
  else {
    string err_message = (string) "PCSet::PCSet(): requested implementation (ISP or NISP) type, " + spType_ +
      (string) ", is not supported";
    throw Tantrum(err_message);
  }


 
  return ;

}






void PCSet::InitISP()
{

  // Initialize the default quadrature points and evaluate basis functions at those
  // locations
  /// \todo Need to find a better way to handle high-dimensional systems
  if (nDim_<=6) {// A temporary solution
    // This may seem like overkill, but in the triple products integrand of order 3*order_ need to be evaluated 
    this->SetQuadRule(pcType_,"full",2*order_+1);
    cout << "Used " << 2*order_+1 << " quadrature points per dimension for initialization." << endl;
  }
  else{
    // The third argument is the level of the sparse quadrature, this should be enough to integrate up to order 3*order_ exactly
    this->SetQuadRule(pcType_,"sparse",order_); 
    //  this->SetQuadRule("CC","sparse",order_); 
    cout << "Used level " << order_ << " sparse quadrature points for initialization." << endl;
  }

 

  // Initialize the factors for products
  this->EvalBasisProd3();

  // Set some defaults
  this->rTolTaylor_    = 1.e-6        ;
  this->maxTermTaylor_ = 500          ; 
  this->SMALL_         = 1.e-15       ;
  this->rTolGMRESDiv_  = 1.e-8        ;
  this->logMethod_     = TaylorSeries ;

  // Set cvode defaults
  this->CVmaxord_      = 8      ; 
  this->CVmaxnumsteps_ = 5000   ;
  this->CVinitstep_    = 2.e-8  ;
  this->CVmaxstep_     = 2.e-1  ;
  this->CVrelt_        = 1.e-10 ;
  this->CVabst_        = 1.e-14 ;

  return;
}

void PCSet::InitNISP()
{

  
 // Initialize the default quadrature points and evaluate basis functions at those
  // locations
  /// \todo Need to find a better way to handle high-dimensional systems
  if (nDim_<=6) {// A temporary solution
    this->SetQuadRule(pcType_,"full",order_+1);
    cout << "Used " << order_+1 << " quadrature points per dimension for initialization." << endl;
  }
  else{
    // The third argument is the level of the sparse quadrature, this should be enough to integrate up to order 3*order_ exactly
    this->SetQuadRule(pcType_,"sparse",order_); 
    //  this->SetQuadRule("CC","sparse",order_); 
    cout << "Used level " << order_ << " sparse quadrature points for initialization." << endl;
  }

 
  
  return;
}


void PCSet::EvalBasisProd3()
{

  // Allocate storage for non-zero <\Psi_i \Psi_j \Psi_k> and their indices
  this->iProd2_.Clear();
  this->iProd2_.Resize(this->nPCTerms_);
  this->jProd2_.Clear();
  this->jProd2_.Resize(this->nPCTerms_);
  this->psiIJKProd2_.Clear();
  this->psiIJKProd2_.Resize(this->nPCTerms_);

  // Evaluate all possible <\Psi_i \Psi_j \Psi_k> that would enter the computation
  // of the k-th PC coefficient in the product of two PCEs. Only retain the ones
  // that are non-zero.
  if(uqtkverbose_>0){
    cout << endl;
    cout << "Computation of <\\Psi_i \\Psi_j \\Psi_k>'s. The number of non-zero entries for" << endl;
  }
  for(int k=0; k < nPCTerms_; k++){
    // Make sure there are no elements in the <\Psi_i \Psi_j \Psi_k> storage vector
    iProd2_(k).Clear();
    jProd2_(k).Clear();
    psiIJKProd2_(k).Clear();
    for(int j=0; j < nPCTerms_; j++){
      for(int i=0; i < nPCTerms_; i++){
        double wrkSum = 0.e0;
        for(int iqp=0; iqp < nQuadPoints_; iqp++){
          wrkSum += psi_(iqp,i)*psi_(iqp,j)*psi_(iqp,k)*quadWeights_(iqp);
        }
        
          

        // All these terms are supposed to be integer-valued actually, but
        // roundoff and inaccuracies from the integration can make them
        // slightly non-zero, so we check to see if they are non-zero
        // by checking whether they are larger than 0.001 in magnitude
        if(fabs(wrkSum) > 0.001){
          // Store the term and its indices
          iProd2_(k).PushBack(i);
          jProd2_(k).PushBack(j);
          psiIJKProd2_(k).PushBack(wrkSum);
        }
      }
    }
    if(uqtkverbose_>0){
      cout << "k = " << k << " : " << psiIJKProd2_(k).Length() << endl;
    }
    if(uqtkverbose_>1){
     	Array1D<double> tmpw;
      tmpw=psiIJKProd2_(k);
      Array1D<int> tmpi, tmpj;
      tmpi=iProd2_(k);
      tmpj=jProd2_(k);
      for (int m=0;m<(int) psiIJKProd2_(k).Length();m++)
        cout << tmpi(m) << " " << tmpj(m) << " " << k << " : " << tmpw(m) << endl;
      cout << endl;
    }   

  }
  cout << endl;
 

  return;
}


void PCSet::SetQuadRule(const string grid_type,const string fs_type,int param)
{
  // Move strings to char*, since Quad needs char* 
  /// \todo Need to improve it
  char grid_type_c[10], fs_type_c[10];
  strcpy(grid_type_c, grid_type.c_str());
  strcpy(fs_type_c, fs_type.c_str());

  Quad spRule(grid_type_c,fs_type_c,nDim_,param,alpha_,beta_); 
  spRule.SetRule();

  // Get the intergation rule information
  spRule.GetRule(this->quadPoints_, this->quadWeights_,this->quadIndices_);
  this->nQuadPoints_ = this->quadPoints_.XSize();
  

  // Get the basis values at the default quadrature points
  this->EvalBasisAtCustPts(this->quadPoints_,this->psi_);
  

  return;
}

void PCSet::SetQuadRule(Quad &quadRule)
{
  
  // Get the intergation rule information
  quadRule.GetRule(this->quadPoints_, this->quadWeights_,this->quadIndices_);
  this->nQuadPoints_ = this->quadPoints_.XSize();

  
  // Get the basis values at the default quadrature points
  this->EvalBasisAtCustPts(this->quadPoints_,this->psi_);
  

 
  return;
}


void PCSet::PrintMultiIndex() const
{
  
  cout << "====================================================" << endl;
  cout << "multi-indices(i,j), with i = 0 ... P, j = 1 ... nDim" << endl;
  cout << "with P = " << nPCTerms_ - 1 << ", and nDim = " << nDim_ << endl;
  cout << "----------------------------------------------------" << endl;
  cout << "i\\j | ";
  for(int idim=0; idim < nDim_; idim++){
    cout << idim + 1  << "  ";
  }
  cout << endl;
  cout << "------";
  for(int idim=0; idim < nDim_; idim++){
    cout << "---";
  }
  cout << endl;
  for(int ip=0; ip < nPCTerms_; ip++){  
    cout << ip << "   | ";
    for(int idim=0; idim < nDim_; idim++){
      cout << multiIndex_(ip,idim) << "  ";
    }
    cout << endl;
  }
  return;
}

void PCSet::PrintMultiIndexNormSquared() const
{

  cout << "===================================================================" << endl;
  cout << "multi-indices(i,j), and \\Psi^2(i) with i = 0 ... P, j = 1 ... nDim" << endl;
  cout << "with P = " << nPCTerms_ - 1 << ", and nDim = " << nDim_ ;
  cout << " (terms up to order " << order_ << ")" << endl;
  cout << "-------------------------------------------------------------------" << endl;
  cout << "i\\j | ";
  for(int idim=0; idim < nDim_; idim++){
    cout << idim + 1  << "  ";
  }
  cout << "\\Psi^2" << endl;
  cout << "------";
  for(int idim=0; idim < nDim_; idim++){
    cout << "---------";
  }
  cout << endl;
  for(int ip=0; ip < nPCTerms_; ip++){  
    cout << ip << "   | ";
    for(int idim=0; idim < nDim_; idim++){
      cout << multiIndex_(ip,idim) << "  ";
    }
    cout << psiSq_(ip) << endl;
  }

}

void PCSet::InitMeanStDv( const double& m, const double& s, double* p ) const
{

  // Check to make sure we have a 1D PCE
  if ( this->nDim_ != 1 ){ 
    string err_message = (string) "PCSet::InitMeanStDv(): Is only implemented for 1D PCEs";
    throw Tantrum(err_message);
  }

  // Check to make sure we have enough terms to set the standard deviation
  if ( this->nPCTerms_ < 2){
    string err_message = (string) "PCSet::InitMeanStDv(): At least 2 terms needed to set a standard deviation";
    throw Tantrum(err_message);
  }

  // Make sure we have a positive standard deviation s 
  if ( s < 0.e0 ){
    string err_message = (string) "PCSet::InitMeanStDv(): input standard deviation is negative";
    throw Tantrum(err_message);
  }

  // Set the mean
  p[0] = m;

  // Set the standard deviation
  p[1] = s/(sqrt(this->psiSq_(1)));

  // Set the remainder of coeffients to 0
  for (int k = 2; k < this->nPCTerms_; k++){
    p[k] = 0.e0;
  }

  return;

}

void PCSet::InitMeanStDv( const double& m, const double& s, Array1D<double>& p ) const
{

  // Check array sizes
  if ( (int) p.Length() != this->nPCTerms_ ){ 
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::InitMeanStDv(): array size does not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }

  // Initialize Mean and Standard Deviation
  this->InitMeanStDv( m, s, p.GetArrayPointer() ) ;

  return;

}

void PCSet::Copy(double* p1, const double* p2) const
{
  // p1 = p2
  for(int ip=0; ip < this->nPCTerms_; ip++){
    p1[ip] = p2[ip];
  }
  return;
}

void PCSet::Copy(Array1D<double>& p1, const Array1D<double>& p2) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) || 
      ( (int) p2.Length() != this->nPCTerms_ ) ){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::Copy(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // p1 = p2
  this->Copy(p1.GetArrayPointer(), p2.GetConstArrayPointer());

  return;
}

void PCSet::Add(const double* p1, const double* p2, double* p3) const
{
  // Add two arrays and return result in p3
  for(int ip=0; ip < this->nPCTerms_; ip++){
    p3[ip] = p1[ip] + p2[ip];
  }
  return;
}

void PCSet::Add(const Array1D<double>& p1, const Array1D<double>& p2, Array1D<double>& p3) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) || 
      ( (int) p2.Length() != this->nPCTerms_ ) || 
      ( (int) p3.Length() != this->nPCTerms_ ) ){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::Add(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Add two arrays and return result in p3
  this->Add(p1.GetConstArrayPointer(), p2.GetConstArrayPointer(), p3.GetArrayPointer());

  return;
}

void PCSet::AddInPlace(double* p1, const double* p2) const
{
  // Add p2 to p1 and return result in p1
  for(int ip=0; ip < this->nPCTerms_; ip++){
    p1[ip] += p2[ip];
  }
  return;
}

void PCSet::AddInPlace(Array1D<double>& p1, const Array1D<double>& p2) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) || 
      ( (int) p2.Length() != this->nPCTerms_ ) ){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::AddInPlace(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Add p2 to p1 and return result in p1
  this->AddInPlace(p1.GetArrayPointer(), p2.GetConstArrayPointer());

  return;
}

void PCSet::Multiply(const double* p1, const double& a, double* p2) const
{
  // Multiply p1 with a and return result in p2
  for(int ip=0; ip < this->nPCTerms_; ip++){
    p2[ip] = p1[ip]*a;
  }
  return;
}

void PCSet::Multiply(const Array1D<double>& p1, const double& a, Array1D<double>& p2) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_) ||
      ( (int) p2.Length() != this->nPCTerms_) ){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::Multiply(): array size does not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Multiply p1 with a and return result in p2
  this->Multiply(p1.GetConstArrayPointer(), a, p2.GetArrayPointer());

  return;
}

void PCSet::MultiplyInPlace(double* p1, const double& a) const
{
  // Multiply p1 with a and return result in p1
  for(int ip=0; ip < this->nPCTerms_; ip++){
    p1[ip] *= a;
  }
  return;
}

void PCSet::MultiplyInPlace(Array1D<double>& p1, const double& a) const
{
  // Check array sizes
  if( (int) p1.Length() != this->nPCTerms_){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::MultiplyInPlace(): array size does not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Multiply p1 with a and return result in p1
  this->MultiplyInPlace(p1.GetArrayPointer(), a);

  return;
}

void PCSet::Subtract(const double* p1, const double* p2, double* p3) const
{
  // p3 = p1 - p2
  for(int ip=0; ip < this->nPCTerms_; ip++){
    p3[ip] = p1[ip] - p2[ip];
  }
  return;
}

void PCSet::Subtract(const Array1D<double>& p1, const Array1D<double>& p2, Array1D<double>& p3) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) || 
      ( (int) p2.Length() != this->nPCTerms_ ) || 
      ( (int) p3.Length() != this->nPCTerms_ ) ){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::Subtract(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // p3 = p1 - p2
  this->Subtract(p1.GetConstArrayPointer(), p2.GetConstArrayPointer(), p3.GetArrayPointer());

  return;
}

void PCSet::SubtractInPlace(double* p1, const double* p2) const
{
  // p1 = p1 - p2
  for(int ip=0; ip < this->nPCTerms_; ip++){
    p1[ip] = p1[ip] - p2[ip];
  }
  return;
}

void PCSet::SubtractInPlace(Array1D<double>& p1, const Array1D<double>& p2) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) || 
      ( (int) p2.Length() != this->nPCTerms_ ) ){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::SubtractInPlace(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // p1 = p1 - p2
  this->SubtractInPlace(p1.GetArrayPointer(), p2.GetConstArrayPointer());

  return;
}

void PCSet::Prod(const double* p1, const double* p2, double* p3) const
{
  // work variables
  int i;
  int j;
  double c;

  // Multiply two arrays and return result in p3
  for(int k=0; k < this->nPCTerms_; k++){
    double tmpSum = 0.e0;
    // Summation over i, j, using only the terms with non-zero Cijk's
    for(int ic=0; ic < (int) (this->psiIJKProd2_(k).Length()); ic++){
      i = this->iProd2_(k)(ic);
      j = this->jProd2_(k)(ic);
      c = this->psiIJKProd2_(k)(ic);
      tmpSum += p1[i] * p2[j] * c;
    }
    p3[k] = tmpSum/this->psiSq_(k);
  }

  return;
}

void PCSet::Prod(const Array1D<double>& p1, const Array1D<double>& p2, Array1D<double>& p3) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) || 
      ( (int) p2.Length() != this->nPCTerms_ ) || 
      ( (int) p3.Length() != this->nPCTerms_ ) ){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::Prod(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Multiply two arrays and return result in p3
  this->Prod(p1.GetConstArrayPointer(), p2.GetConstArrayPointer(), p3.GetArrayPointer());

  return;
}

void PCSet::Polyn(const double* polycf, const int npoly, const double* p1, double* p2) const
{
  // Work variables
  double ptemp[this->nPCTerms_];

  

  double ps[npoly-1];
  for(int i=0;i<npoly-1;i++) ps[i]=polycf[i+1];

  if (npoly>1){
    this->Polyn(ps,npoly-1,p1, ptemp);
    this->Prod(p1,ptemp,p2);   
  }
  

  p2[0]+=polycf[0]; 


  return;
}

void PCSet::Polyn(const Array1D<double>& polycf, const Array1D<double>& p1, Array1D<double>& p2) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) || 
      ( (int) p2.Length() != this->nPCTerms_ ) ){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::Polyn(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Evaluate a polynomial of p1 and return result in p2
  this->Polyn(polycf.GetConstArrayPointer(), polycf.Length(), p1.GetConstArrayPointer(), p2.GetArrayPointer());

  return;
}

void PCSet::PolynMulti(const Array1D<double>& polycf, const Array2D<int>& mindex, const Array2D<double>& p1, Array1D<double>& p2) const
{
  // Work variables
  int nd=mindex.YSize();
  int nterms=mindex.XSize();


  // Check array sizes
  if ( (int) p1.XSize() != this->nPCTerms_ ) {
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::PolynMulti(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  if ( (int) polycf.Length() != nterms ){
    string err_message = (string) "PCSet::PolynMulti(): array sizes do not match" ;
    throw Tantrum(err_message);
  }
  if (  (int) p1.YSize() != nd ){
    string err_message = (string) "PCSet::PolynMulti(): polynomial multiindex array size does not match the number"
      + " of PC coefficient vectors";
    throw Tantrum(err_message);
  }

  p2.Resize(this->nPCTerms_,0.e0);


  // Prepare arrays for recursive calls
  Array1D<int> zeroind, nzeroind;
  Array1D<double> polycf1, polycf2;
  for (int ii=0; ii < nterms; ii++){
    if (mindex(ii,0)==0){
      polycf1.PushBack(polycf(ii));
      zeroind.PushBack(ii);
    }
    else {
      polycf2.PushBack(polycf(ii));     
      nzeroind.PushBack(ii);
    }
  }
  int nz=zeroind.Length();
  int nnz=nzeroind.Length();


  Array2D<int> mindex1(nz,nd-1,0);
  for (int i=0; i<nz; i++)
    for(int j=0; j<nd-1; j++)
      mindex1(i,j)=mindex(zeroind(i),j+1);
  
  Array2D<double> p1_1(this->nPCTerms_,nd-1,0.e0);
  for (int i=0; i<this->nPCTerms_; i++)
    for(int j=0; j<nd-1; j++)
      p1_1(i,j)=p1(i,j+1);


  Array2D<int> mindex2(nnz,nd,0);
  for (int i=0; i<nnz; i++)
    for(int j=0; j<nd; j++)
      mindex2(i,j)=mindex(nzeroind(i),j)-(j==0);

  Array1D<double> pfirst(this->nPCTerms_,0.e0);
  for (int i=0; i<this->nPCTerms_; i++)
    pfirst(i)=p1(i,0);

  // Evaluate recursively

  Array1D<double> ptmp2(this->nPCTerms_,0.e0);
  if (nd==1){
    for (int i=0;i<nz;i++)
      p2(0)+=polycf1(i);
  }
  else {
    if (nz>0)
      this->PolynMulti(polycf1,mindex1,p1_1,p2);
  }


  if (nnz>0){
    Array1D<double> ptmp1;
    this->PolynMulti(polycf2,mindex2,p1,ptmp1);
    this->Prod(pfirst,ptmp1,ptmp2);
  }

  this->AddInPlace(p2,ptmp2);

  return;
}

void PCSet::Exp(const double* p1, double* p2) const
{
  // work variables
  const double maxMult = 1.e3;  // max growth in rel. error w/o stopping the run

  // array to store p1 minus its mean
  double* x = new double[this->nPCTerms_];
  
  x[0] = 0.e0; // random variable with mean zero
  for(int k=1; k < this->nPCTerms_; k++){
    x[k] = p1[k];
  }

  // exp() of the mean of p1
  double expMean = exp(p1[0]);

  //
  // find exp(x) with Taylor series 1 + x + x^2/2! + x^3/3! + ...
  // where each term d_n is computed recursively
  // as d_n = d_{n-1}*x/n
  //

  // Work array for storing the current term in the Taylor series
  double* d = new double[this->nPCTerms_];

  // Set d equal to the first order term, which is x
  for(int k=0; k < this->nPCTerms_; k++){
    d[k] = x[k];
  }

  // Store the first two terms: 1 + x = 1 + d in p2
  p2[0] = 1.e0 + d[0];
  for(int k=1; k < this->nPCTerms_; k++){
    p2[k] = d[k];
  }

  // work variables for progress tracking. Set them to satisfy the while condition initially
  int tOrder = 1;   // We are up to 1st order in the terms
  double rErr = 1.e50; // relative error
  double rErrOld = (maxMult + 1.e0)*rErr; // Old relative error

  // other work vars
  double sc;
  double* dn = new double[this->nPCTerms_]; // term n
  double* fac = new double[this->nPCTerms_]; // factor to multiply d_{n-1} with to get d_n

  while (rErr > this->rTolTaylor_ && 
         tOrder+1 < this->maxTermTaylor_ && 
         rErr <= maxMult*rErrOld){

    rErrOld = rErr;
    tOrder++;

    sc = 1.e0/(double) tOrder;
    for(int k=0; k < this->nPCTerms_; k++){
      fac[k] = x[k]*sc; // fac = x/n
    }

    this->Prod(d,fac,dn); // dn = d*x/n

    rErr = 0.e0;
    for(int k=0; k < this->nPCTerms_; k++){
      p2[k] += dn[k];   // Update p2 with new term
      d[k] = dn[k];     // Save term for next step
      rErr = max(rErr,fabs(dn[k]/p2[0])); // error as mag. of PC coeff. over the mean
    }
  }

  if(uqtkverbose_>0){
    cout << "number of terms in exp(p1) = " << tOrder + 1 << endl;
  }

  // Check on exit criteria to see if something went wrong

  if(tOrder >= this->maxTermTaylor_){
    ostringstream buffer1, buffer2;
    buffer1 << tOrder+1; // The total # of terms is the order + 1
    buffer2 << rErr;     // Relative error at this point
    string err_message = (string) "PCSet::Exp(): Rel. tolerance criterium for Taylor"
      + "series not met after " + buffer1.str() + " terms.\n" 
      + "Relative error at this point is " + buffer2.str() + ".";
    throw Tantrum(err_message);
  }

  if(rErr > maxMult*rErrOld){
    ostringstream buffer1, buffer2;
    buffer1 << tOrder+1; // The total # of terms is the order + 1
    buffer2 << rErr;     // Relative error at this point
    string err_message = (string) "PCSet::Exp(): Taylor series diverging after "
      + buffer1.str() + " terms.\n"
      + "Relative error at this point is " + buffer2.str() + ".";
    throw Tantrum(err_message);
  }

  // multiply the result back with exp(p1(0)) to factor in the exp of the mean
  for(int k=0; k < this->nPCTerms_; k++){
    p2[k] = p2[k]*expMean;
  }

  // clear out work memory
  delete[] x;
  delete[] d;
  delete[] dn;
  delete[] fac;

  return;
}

void PCSet::Exp(const Array1D<double>& p1, Array1D<double>& p2) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) || 
      ( (int) p2.Length() != this->nPCTerms_ ) ){ 
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::Exp(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Take Exp of p1 and return it in p2
  this->Exp(p1.GetConstArrayPointer(), p2.GetArrayPointer());

  return;
}

void PCSet::Log(const double* p1, double* p2) const
{
 
  
  // Computes natural logarithm of PC expansion p1 using Taylor expansion:
  if ( logMethod_ == TaylorSeries )
    LogTaylor(p1,p2) ;
  else if ( logMethod_ == Integration )
    LogInt(p1,p2) ;
  else
    {
      ostringstream buffer ;
      buffer << logMethod_ ; // unknown method for computing natural logarithm
      string err_message = (string) "PCSet::Log(): Unknown method : "
	+ buffer.str() + ", for computing logarithm of a PC expansion ";
      throw Tantrum(err_message);
    }  

  return ;
}

void PCSet::Log(const Array1D<double>& p1, Array1D<double>& p2) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) || 
      ( (int) p2.Length() != this->nPCTerms_ ) ) { 
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::Log(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Take Log of p1 and return it in p2
  this->Log(p1.GetConstArrayPointer(), p2.GetArrayPointer());

  return;
}

void PCSet::Log10(const double* p1, double* p2) const
{

  // Computes logarithm to base 10 of PC expansion p1 and stores the result in p2
  // First use Log() to compute the natural logarithm and then divide it by log(10)

  this->Log(p1,p2);

  double sc = 1.e0 / log( 1.e1 ) ;
  for(int k=0; k < this->nPCTerms_; k++){
    p2[k] *= sc;  
  }

  return;

}

void PCSet::Log10(const Array1D<double>& p1, Array1D<double>& p2) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) || 
      ( (int) p2.Length() != this->nPCTerms_ ) ){ 
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::Log10(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Take Log of p1 and return it in p2
  this->Log10(p1.GetConstArrayPointer(), p2.GetArrayPointer());

  return;
}

void PCSet::RPow(const double* p1, double* p2, const double& a) const
{

  // Compute p2=p1^a as exp(a*log(p1))

  // allocate intermediate work array 
  double* pwrk = new double[this->nPCTerms_];

  // 1) take log of p1
  this->Log(p1,pwrk);
  // 2) multiply result by a
  this->MultiplyInPlace(pwrk,a) ;
  // 3) Compute exp to obtain the final result
  this->Exp(pwrk,p2) ;
  
  delete [] pwrk ;

  return;

}

void PCSet::RPow(const Array1D<double>& p1, Array1D<double>& p2, const double& a) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) || 
      ( (int) p2.Length() != this->nPCTerms_ ) ){ 
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::RPow(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Compute p1^a and return it in p2
  this->RPow(p1.GetConstArrayPointer(), p2.GetArrayPointer(), a);

  return;
}


void PCSet::IPow(const double* p1, double* p2, const int& ia) const
{

  // Set p2 = 1
  p2[0] = 1.0 ;
  for ( int k = 1 ; k < this->nPCTerms_ ; k++ )
    p2[k] = 0.0 ;

  // if ia = 0 -> return
  if ( ia == 0 ) return;

  // Check if power is negative; if yes, make it positive and remember this
  int  ialocal = ia     ;
  bool isneg   = false  ;
  if ( ialocal < 0 ){
    ialocal *= -1    ;
    isneg    =  true ;
  }

  // Check if ia is even/odd and re-initialize p2 appropriately
  if ( ialocal & 1 )
    for ( int k = 0 ; k < this->nPCTerms_ ; k++ )
      p2[k] = p1[k] ;

  // allocate intermediate work arrays and initialize
  double* pwrk1 = new double[this->nPCTerms_];
  double* pwrk2 = new double[this->nPCTerms_];
  for ( int k = 0 ; k < this->nPCTerms_ ; k++ )
    pwrk1[k] = p1[k] ;

  // define additional pointer that will switch between pwrk1 and pwrk2
  double *pw1 = 0, *pw2 = 0, *ptmp = 0 ;
  pw1 = pwrk1 ;
  pw2 = pwrk2 ;
 
  // determine the number of digits in the binary representation of ia
  int ndig ;
  ndig = (int) (log( (double) ia + SMALL ) / log(2.0) ) + 1 ; 
  
  // start loop
  for ( int idig = 1 ; idig < ndig ; idig++ ){
    
    // compute pw2 = pw1^2
    this->Prod(pw1,pw1,pw2) ;
    
    // check if binary representation of ia contains this power  
    if ( (ialocal & (1<<idig) )>>idig ) {
      this->Prod(p2,pw2,pw1) ;
      for ( int k = 0 ; k < this->nPCTerms_ ; k++ )
        p2[k] = pw1[k] ;  
      
    }

    // Switch pointers
    ptmp = pw2  ;
    pw2  = pw1  ;
    pw1  = ptmp ;

  }

  // if power is negative complete the process by taking the inverse:
  // p2 = 1/(p1^|ia|)
  if ( isneg ){
    
    for ( int k = 0 ; k < this->nPCTerms_ ; k++ )
      pwrk1[k] = p2[k] ;  
    
    this->Inv(pwrk1,p2) ;

  }

  delete [] pwrk1 ;
  delete [] pwrk2 ;

  return;

}

void PCSet::IPow(const Array1D<double>& p1, Array1D<double>& p2, const int& ia) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) || 
      ( (int) p2.Length() != this->nPCTerms_ )){ 
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::IPow(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Compute p1^ia and return it in p2
  this->IPow(p1.GetConstArrayPointer(), p2.GetArrayPointer(), ia);

  return;
}


void PCSet::Inv(const double* p1, double* p2) const
{

  // Compute p2=1/p1 

  // allocate and initialize intermediate work array 
  double* pwrk = new double[this->nPCTerms_];
  pwrk[0] = 1.0 ;
  for ( int k = 1 ; k < this->nPCTerms_ ; k++ )
    pwrk[k] = 0.0 ;

  // Compute 1/p1
  this->Div(pwrk,p1,p2) ;
  
  delete [] pwrk ;

  return;

}

void PCSet::Inv(const Array1D<double>& p1, Array1D<double>& p2) const
{
  // Check array sizes
  if( ((int) p1.Length() != this->nPCTerms_ ) || 
      ((int) p2.Length() != this->nPCTerms_ ) ){ 
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::Inv(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Compute p1*a and return it in p2
  this->Inv( p1.GetConstArrayPointer(), p2.GetArrayPointer() );

  return;
}


void PCSet::GMRESMatrixVectorProd(const double* x, const double*a, double* y) const
{
  // Perform matrix vector product as product between two PC variables
  this->Prod(x,a,y);

  return;
}

void PCSet::Div(const double* p1, const double* p2, double* p3) const
{
#define DIV_USE_GMRES
  // Get the division p3 = p1/p2 by solving the system of eqn. p2*p3 = p1

  // Do a test to see whether p1 happens to be zero, in which case p3 = 0
  // (assuming p2 is not zero) (GMRES would choke if p1 is 0)
  double p1_rms = this->GetModesRMS(p1);

  if(p1_rms == 0.e0){
    for(int k=0; k < this->nPCTerms_; k++){ // Set p3 to zero and return
      p3[k] = 0.e0;
    }
    return;
  }

  // Initialize p3 with the right hand side p1
  for(int k=0; k < this->nPCTerms_; k++){
    p3[k] = p1[k];
  }

#ifdef DIV_USE_GMRES

  // Use GMRES to solve the system of equations. GMRES assumes
  // a sparse form is used to store the matrix A, but does not
  // case what format is used, as long as we provide it with
  // a routine to do preconditioning and matrix-vector products
  // with. As the matrix-vector multiplication in this case
  // comes down to a Product between two PC variables,
  // we can just use PCSet::Prod to do the product and
  // access all required info through the data members of PCSet.
  // So we pass p2 as A and say it has nPCTerms_ non-zero elements
  // In this implementation, there is also no preconditioning
  // nor scaling of X and B used.

  // Set up parameters and work space for GMRES
  // The parameters below correspond to the descriptions given in the dgmres.f
  // implementation in the slatec library. For the tolerance parameter, we will
  // use rTolGMRESDiv_. The isym argument will be repurposed to pass an index
  // to a map containing pointers to this PCSet class object (in order to 
  // enable the callbacks from GMRES to the matrix-vector and preconditioning
  // routines in this class).
  int iter;           // number of iterations performed
  int ierr;           // error flag
  double err;         // to hold error estimate of final solution
  int itol = 0;       // convergence test on residual
  int itmax = 0;      // dummy variable
  int iunit = 0;      // no output desired
  int jscal = 0;      // no scaling arrays SB and SX used
  int jpre = 0;       // no preconditioning
  int nrmax = 10;     // max number of restarts in Krylov iteration

  int maxl = 20; // Max dimension of Krylov subspace
  int ligw = 20; // dimension of integer work array
  int lrgw = 1 + this->nPCTerms_ * (maxl + 6) + maxl * (maxl+3); // size of real work array

  double* rgwk = new double[lrgw];  // real work array
  double* rdum = new double[1];     // dummy real array
  int* igwk = new int[ligw];        // integer work array
  int* idum = new int[1];           // dummy integer array
  int* ia = new int[this->nPCTerms_]; // Arrays for matrix structure of A
  int* ja = new int[this->nPCTerms_]; // Arrays for matrix structure of A

  igwk[0] = maxl;
  igwk[1] = maxl;
  igwk[2] = jscal;
  igwk[3] = jpre;
  igwk[4] = nrmax;

  // transfer some arguments to maintain the imposed const constraints on original data
  int nelt = this->nPCTerms_;
  double rtol = this->rTolGMRESDiv_;

  double* p1_w = new double[this->nPCTerms_];
  double* p2_w = new double[this->nPCTerms_];

  for ( int k = 0 ; k < this->nPCTerms_ ; k++ ){
    p1_w[k] = p1[k];
    p2_w[k] = p2[k];
  }

  // Call GMRES to solve the system of equations, and pass the appropriate
  // handles to this class and the static void callback functions for
  // matrix vector multiplications and preconditioning.

  int handle = my_index_;
  f77_matvecprod matvec_callee = GMRESMatrixVectorProdWrapper;
  f77_precond precond_callee = GMRESPreCondWrapper;

  FTN_NAME(dgmres)(&nelt, p1_w, p3, &nelt, ia, ja,
		   p2_w, &handle, matvec_callee, precond_callee,
		   &itol, &rtol, &itmax, &iter, &err, &ierr, &iunit, rdum, rdum,
		   rgwk, &lrgw, igwk, &ligw, rdum, idum);

  if (ierr != 0) {
    ostringstream buffer;
    buffer << ierr;
    string err_message = (string) "PCSet::Div(): error " + buffer.str() + " occurred in GMRES";
    throw Tantrum(err_message);
  }
    
  // Clean up work variables
  delete [] rgwk;
  delete [] rdum;
  delete [] igwk;
  delete [] idum;
  delete [] ia;
  delete [] ja;
  delete [] p1_w;
  delete [] p2_w;

#else /* use Gauss Elimination approach */

  // Allocate work arrays and variables
  Array1D<double> rwrk(this->nPCTerms_,0.e0); // real work array
  Array1D<int> iwrk(this->nPCTerms_,0);       // integer work array
  Array2D<double> a(this->nPCTerms_,this->nPCTerms_,0.e0); // matrix
  int i;
  int j;
  double c;

  // Fill matrix A to represent p2*p3 with p3 unknown
  for(int k=0; k < this->nPCTerms_; k++){
    // Summation over i, j, using only the terms with non-zero Cijk's
    for(int ic=0; ic < this->psiIJKProd2_(k).Length(); ic++){
      i = this->iProd2_(k)(ic);
      j = this->jProd2_(k)(ic);
      c = this->psiIJKProd2_(k)(ic)/this->psiSq_(k);
      a(k,i) +=  p2[j] * c;
    }
  }

  /*
    cout << "Elements a(k,i) are" << endl;
    for(int k=0; k < this->nPCTerms_; k++){
    for(int i=0; i < this->nPCTerms_; i++){
    cout << a(k,i) << " ";
    }
    cout << endl;
    }
  */

  // Solve the system with Gauss Elimination
  int itask = 1;
  int ind = 0;
  int nTerms = this->nPCTerms_;
  FTN_NAME(dgefs)(a.GetArrayPointer(),&nTerms,&nTerms,p3,&itask,&ind,
		  rwrk.GetArrayPointer(),iwrk.GetArrayPointer());
  if(ind < 8){
    cout << "PCSet:Div(), estimated number of accurate digits from dgefs = " << ind << endl;
  }

#endif

  return;
}

void PCSet::Div(const Array1D<double>& p1, const Array1D<double>& p2, Array1D<double>& p3) const
{
  // Check array sizes
  if( ( (int) p1.Length() != this->nPCTerms_ ) ||
      ( (int) p2.Length() != this->nPCTerms_ ) || 
      ( (int) p3.Length() != this->nPCTerms_ ) ){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::Div(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Divide p1 by p2 and return result in p3
  this->Div(p1.GetConstArrayPointer(), p2.GetConstArrayPointer(), p3.GetArrayPointer());

  return;
}

double PCSet::StDv( const double* p ) const
{

  // work variables
  double tsum = 0.e0 ;

  // First calculate the variance
  for ( int k = 1 ; k < this->nPCTerms_ ; k++ ){
    tsum += ( p[k] * p[k] * this->psiSq_(k) ) ;
  }

  // Next take the square root to get standard deviation and return it
  double stdev = sqrt( tsum ) ;

  return stdev ;

}

double PCSet::StDv( const Array1D<double>& p ) const
{

  // Check array sizes
  if ( (int) p.Length() != this->nPCTerms_ ){ 
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::StDv(): array size does not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }

  // Take StDv of p
  double stdv = this -> StDv( p.GetConstArrayPointer() ) ;

  return stdv ;

}

double PCSet::GetModesRMS(const double* p1) const
{
  // Get the rms of PC coefficients.
  double sum = 0;
  for(int k=0; k < this->nPCTerms_; k++){
    sum += p1[k]*p1[k];
  }

  return sqrt(sum/(double) this->nPCTerms_);
}

double PCSet::GetModesRMS(const Array1D<double>& p1) const
{
  // Check array size
  if( (int) p1.Length() != this->nPCTerms_){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::GetModesRMS(): array size does not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }

  // Get the rms of the PC coeffs
  return this->GetModesRMS(p1.GetConstArrayPointer());
}

void PCSet::LogTaylor(const double* p1, double* p2) const
{
  // work variables
  const double maxMult = 0.8;  // max growth in rel. error w/o stopping the run
  double  sc; // used to compute recursive factors

  // mean of p1
  double p1Mean = p1[0];

  if( fabs(p1Mean) < SMALL_ ){
    ostringstream buffer ;
    buffer << p1Mean; // The mean value of PC expansion p1
    string err_message = (string) "PCSet::Log(): Cannot use Taylor series to compute log(p1)\n"
      + "The mean value, " + buffer.str() + ", of PC expansion p1 is too small.";
    throw Tantrum(err_message);
  }

  // array storing (p1-p1Mean)/p1Mean
  double* x = new double[this->nPCTerms_];
  
  sc = 1.e0 / p1Mean ;
  x[0] = 0.e0 ;
  for(int k=1; k < this->nPCTerms_; k++){
    x[k] = p1[k]*sc ;
  }

  // Work array for storing the current term in the Taylor series
  double* d = new double[this->nPCTerms_];

  // work variables for progress tracking. Set them to satisfy the while condition initially
  int    tOrder  = 1                     ; // Taylor expansion contains at least up to 1st order terms
  double rErr    = 1.e50                 ; // relative error
  double rErrOld = (maxMult + 1.e0)*rErr ; // old relative error

  // other work vars
  double* dn  = new double[this->nPCTerms_]; // term n
  double* fac = new double[this->nPCTerms_]; // factor to multiply d_{n-1} with to get d_n

  // Set d equal to the first order term
  for(int k=0; k < this->nPCTerms_; k++){
    d[k] = x[k] ;
  }

  // Store the first two terms: log(p1Mean)+(p1-p1Mean)/p1
  p2[0] = log( p1Mean ) ;
  for(int k=1; k < this->nPCTerms_; k++){
    p2[k] = d[k] ;
  }

  // Start loop to add higher order terms
  while ( rErr     >  this->rTolTaylor_    && 
	  tOrder+1 <  this->maxTermTaylor_ && 
	  rErr     <= maxMult*rErrOld        ){

    rErrOld = rErr;
    tOrder++ ;

    sc = - (double) (tOrder-1) / (double) tOrder ;

    for(int k=0; k < this->nPCTerms_; k++){
      fac[k] = x[k]*sc; // fac = (-(n-1)/n)*(p1-p1Mean)/p1Mean
    }

    this->Prod(d,fac,dn); // dn = d*fac

    rErr = 0.e0;
    for(int k=0; k < this->nPCTerms_; k++){
      p2[k] += dn[k];   // Update p2 with new term
      d[k] = dn[k];     // Save term for next step
      rErr = max(rErr,fabs(dn[k]/p2[0])); // error as mag. of PC coeff. over the mean
    }
  }

  if(uqtkverbose_>0){
    cout << "number of terms in log(p1) = " << tOrder + 1 << endl;
  }

  // Check on exit criteria to see if something went wrong

  if(tOrder >= this->maxTermTaylor_){
    ostringstream buffer1, buffer2;
    buffer1 << tOrder+1; // The total # of terms is the order + 1
    buffer2 << rErr    ; // Relative error at this point
    string err_message = (string) "PCSet::Log(): Rel. tolerance criterium for Taylor"
      + "series not met after " + buffer1.str() + " terms.\n" 
      + "Relative error at this point is " + buffer2.str() + ".";
    throw Tantrum(err_message);
  }

  if(rErr > maxMult*rErrOld){
    ostringstream buffer1, buffer2;
    buffer1 << tOrder+1; // The total # of terms is the order + 1
    buffer2 << rErr;     // Relative error at this point
    string err_message = (string) "PCSet::Log(): Taylor series diverging after "
      + buffer1.str() + " terms.\n"
      + "Relative error at this point is " + buffer2.str() + ".";
    throw Tantrum(err_message);
  }

  // clear out work memory
  delete[] x;
  delete[] d;
  delete[] dn;
  delete[] fac;

  return;
}

void PCSet::LogInt(const double* p1, double* p2) const
{
  // Computes natural logarithm of PC expansion p1 using numerical integration

  // initial condition array
  N_Vector u ;
  u = N_VNew_Serial( this->nPCTerms_ ) ;
  double p1Mean = p1[0];  
  NV_Ith_S(u,0) = log( p1Mean ) ;
  for(int k=1; k < (this->nPCTerms_); k++)
    NV_Ith_S(u,k) = 0.0 ;

  // initialize tolerances
  realtype relT ;
  N_Vector absT ;
  relT = CVrelt_ ;
  absT = N_VNew_Serial( this->nPCTerms_ ) ;   
  for ( int k = 0 ; k < ( this->nPCTerms_ ) ; k++ )
    NV_Ith_S(absT,k) = CVabst_ ;
  
  // initialize integration limits
  realtype tstart,tend,tret ;
  tstart = 0.0 ;
  tend   = 1.0 ;

  // cvode flag
  int cvflag ;

  // initialize f_data (work array to be passed to cvode)
  // position  0                                contains the index of the current PCSet object
  // positions 1...nPCTerms_                    contain x0={p1[0],0,0,...,0}
  // positions (nPCTerms_+1)....(2*nPCTerms_)   contain a ={0,p1[1],p1[2],....,p1[nPCTerms_-1]}
  // positions (2*nPCTerms_+1)....(3*nPCTerms_) contain p[t] = x0+a*t
  // positions (3*nPCTerms_+1)....(4*nPCTerms_) contain a/p[t] (to be computed by LogIntRhs) 
  double *f_data ;
  f_data = new double[4*(this->nPCTerms_)+1];

  f_data[0] = (double) (this->my_index_) ;

  f_data[1] = p1[0] ;
  for (int k = 1 ; k < (this->nPCTerms_) ; k++)
    f_data[k+1] = 0.0 ;

  f_data[(this->nPCTerms_)+1] = 0.0;
  for (int k = 1 ; k < (this->nPCTerms_) ; k++)
    f_data[k+(this->nPCTerms_)+1] = p1[k] ;

  for (int k = 0 ; k < (this->nPCTerms_) ; k++)
    f_data[k+2*(this->nPCTerms_)+1] = 0.0 ;

  for (int k = 0 ; k < (this->nPCTerms_) ; k++)
    f_data[k+3*(this->nPCTerms_)+1] = 0.0 ;

  // Create cvode solver
  void *cvode_mem = NULL;
  cvode_mem = CVodeCreate(CV_ADAMS, CV_NEWTON);
  this->Check_CVflag(cvode_mem, "PCSet::LogInt : CVodeCreate", 0) ;

  /* Allocate memory */
  cvflag = CVodeInit(cvode_mem, &LogIntRhsWrapper, tstart, u);
  this->Check_CVflag(&cvflag, "CVodeInit", 1) ;
  cvflag = CVodeSVtolerances(cvode_mem, relT, absT);
  this->Check_CVflag(&cvflag, "CVodeSVtolerances", 1) ;

  // Set dense solver
  cvflag = CVDense(cvode_mem, (this->nPCTerms_) );
  this->Check_CVflag(&cvflag, "PCSet::LogInt : CVDense", 1) ;

  /* Set work array */
  cvflag = CVodeSetUserData(cvode_mem,  (void *) f_data);
  this->Check_CVflag(&cvflag, "CVodeSetUserData", 1) ;

  // Set maximum order for the integratiom method
  cvflag = CVodeSetMaxOrd(cvode_mem, CVmaxord_);
  this->Check_CVflag(&cvflag, "PCSet::LogInt : CVodeSetMaxOrd", 1) ;

  // Set maximum number of steps
  cvflag = CVodeSetMaxNumSteps(cvode_mem, CVmaxnumsteps_);
  this->Check_CVflag(&cvflag, "PCSet::LogInt : CVodeSetMaxNumSteps", 1) ;

  // Set initial step size
  cvflag = CVodeSetInitStep(cvode_mem, CVinitstep_);
  this->Check_CVflag(&cvflag, "PCSet::LogInt : CVodeSetInitStep", 1) ;

  // Set maximum step size
  cvflag = CVodeSetMaxStep(cvode_mem, CVmaxstep_);
  this->Check_CVflag(&cvflag, "PCSet::LogInt : CVodeSetMaxStep", 1) ;

  // Set stop value
  cvflag = CVodeSetStopTime(cvode_mem, tend);
  this->Check_CVflag(&cvflag, "PCSet::LogInt : CVodeSetStopTime", 1) ;

  cvflag = CVode(cvode_mem, tend, u, &tret, CV_NORMAL);
  this->Check_CVflag(&cvflag, "PCSet::LogInt : CVode", 1) ;

  // Retrieve PC terms
  for (int k = 0 ; k < (this->nPCTerms_) ; k++)
    p2[k] = NV_Ith_S(u,k) ;

  // Destroy temporary memory
  N_VDestroy_Serial( u    ) ;
  N_VDestroy_Serial( absT ) ;

  CVodeFree( &cvode_mem ) ;

  delete [] f_data ;

  return ;

}

int PCSet::LogIntRhs(realtype t, N_Vector y, N_Vector ydot, void *f_data) const
{

  // do nothing
  double *x0,*a,*xat,*abyxat ;

  x0       = & ((double *) f_data) [1] ;
  a        = x0 +(this->nPCTerms_) ;
  xat      = a  +(this->nPCTerms_) ;
  abyxat   = xat+(this->nPCTerms_) ;

  // compute x0+a*t
  for (int i = 0 ; i < (this->nPCTerms_) ; i++)
    xat[i] = x0[i]+a[i]*t ;

  // compute a/(x0+a*t)
  this->Div(a,xat,abyxat) ;

  // send rhs back to cvode
  for (int i = 0 ; i < (this->nPCTerms_) ; i++)
    NV_Ith_S(ydot,i) = abyxat[i] ;

  return ( 0 ) ;

}

void PCSet::Derivative(const double* p1, double* p2) const
{
  // Make sure it is 1-d implementation
  if (this->nDim_ != 1)
    throw Tantrum("PCSet::Derivative(): Derivative computation supports 1d PC only");

  // Compute p2=p1' 

 if(!strcmp(this->pcType_.c_str(),"LU")){

  // allocate and initialize intermediate work array (x)
  double* px = new double[this->nPCTerms_];
  for ( int k = 0 ; k < this->nPCTerms_ ; k++ )
    px[k] = 0.0 ;
  px[1] = 1.0 ;

  // allocate and compute an intermediate work array
  double* pn = new double[this->nPCTerms_];
  for ( int k = 0 ; k < this->nPCTerms_ ; k++ )
    pn[k]= -k*p1[k];

  // allocate and compute an intermediate work array
  double* pnx = new double[this->nPCTerms_];
  this->Prod(pn,px,pnx);

  // allocate and compute an array according to recursive relation 
  double* ptmp = new double[this->nPCTerms_];
  for ( int k = 0 ; k < this->nPCTerms_-1 ; k++ )
    ptmp[k]=(k+1)*p1[k+1];
  ptmp[this->nPCTerms_-1]=0.0;

  this->AddInPlace(pnx,ptmp);

 // allocate and initialize intermediate work array (1-x^2)
  double* p1_x2 = new double[this->nPCTerms_];
  for ( int k = 0 ; k < this->nPCTerms_ ; k++ )
    p1_x2[k] = 0.0 ;
  p1_x2[0] = 2.0/3.0 ;
  p1_x2[2] = -2.0/3.0 ;


  // Compute pnx/p1_x2
  this->Div(pnx,p1_x2,p2) ;
  
  delete [] p1_x2 ;
  delete [] px ;
  delete [] pn ;
  delete [] pnx ;
  delete [] ptmp ;
 }

 else if(!strcmp(this->pcType_.c_str(),"HG")){
   for ( int k = 0 ; k < this->nPCTerms_-1 ; k++ )
     p2[k]=(k+1)*p1[k+1];
   p2[this->nPCTerms_-1]=0.0;

 }

 else{
   string err_message = (string) "PCSet::Derivative(): PC type "+this->pcType_+" is not supported by derivative computation";
   throw Tantrum(err_message);

 }

  return;

}

void PCSet::Derivative(const Array1D<double>& p1, Array1D<double>& p2) const
{
  // Check array sizes
  if( ((int) p1.Length() != this->nPCTerms_ ) || 
      ((int) p2.Length() != this->nPCTerms_ ) ){ 
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::Derivative(): array sizes do not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }
  // Compute p1*a and return it in p2
  this->Derivative( p1.GetConstArrayPointer(), p2.GetArrayPointer() );

  return;
}

int PCSet::Check_CVflag(void *flagvalue, const char *funcname, int opt) const
{

  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if ( ( opt == 0 ) && ( flagvalue == NULL ) ) 
    {
      string err_message = (string) "\nCVODE_ERROR: " 
	+funcname+" failed - returned NULL pointer\n\n" ;
      throw Tantrum(err_message);
    }

  /* Check if flag < 0 */
  else if ( opt == 1 ) 
    {
      errflag = (int *) flagvalue;
      if ( *errflag < 0 ) 
	{
	  ostringstream buffer1 ;
	  buffer1 << (int *) flagvalue;
	  string err_message = (string) "\nCVODE_ERROR: " 
	    +funcname+" failed with flag =" + buffer1.str() + "\n\n" ;
	  throw Tantrum(err_message);
	}
    }

  return ( 0 ) ;

}


void PCSet::SeedBasisRandNumGen(const int& seed) const
{
  // Seed the PC Basis random number generator
  this->p_basis_->SeedRandNumGen(seed);

  return;
}

void PCSet::DrawSampleSet(const Array1D<double>& p, Array1D<double>& samples)
{
  // Check array size
  if( (int) p.Length() != this->nPCTerms_){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::DrawSampleSet: p array size does not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }

  // Draw samples
  this->DrawSampleSet(p.GetConstArrayPointer(), samples.GetArrayPointer(), samples.Length());

  return;
}

void PCSet::DrawSampleSet(const double* p, double* samples, const int& nSamples)
{
  
  // Draw samples
  for (int js = 0 ; js < nSamples ; js++) {
    // Draw a sample of the nDim_ dimensional random number vector xi
    Array1D<double> xia(this->nDim_,0.e0);
    this->p_basis_->GetRandSample(xia);

    // Evaluate the PC expansion with this sample xia
    samples[js] = this->EvalPC(p, xia.GetConstArrayPointer());
  }  

  return;
}

void PCSet::DrawSampleVar(Array2D<double>& samples) const
{
  // Check array size
  if( (int) samples.YSize() != this->nDim_){
    ostringstream buffer;
    buffer << this->nDim_;
    string err_message = (string) "PCSet::DrawSampleVar: samples array size does not match the dimensionality"
      + " of PC (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }

  
  int nSamples=samples.XSize();
  
  for (int js = 0 ; js < nSamples ; js++) {
    // Draw a single nDim_ dimensional random vector xia
    Array1D<double> xia(this->nDim_,0.e0);
    this->p_basis_->GetRandSample(xia);
    for(int idim=0;idim<this->nDim_;idim++)
      samples(js,idim)=xia(idim);
  }
  
  return;
}


double PCSet::EvalPC(const Array1D<double>& p, Array1D<double>& randVarSamples)
{
  // Check array sizes
  if( (int) p.Length() != this->nPCTerms_){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::EvalPC: p array size does not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }

  if( (int) randVarSamples.Length() != this->nDim_){
    ostringstream buffer;
    buffer << this->nDim_;
    string err_message = (string) "PCSet::EvalPC: randVarSamples array size does not match the number"
      + " of PC dimensions (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }

  // Evaluate the PC expansion at the random variables given
  return this->EvalPC(p.GetConstArrayPointer(), randVarSamples.GetConstArrayPointer());
}

double PCSet::EvalPC(const double* p, const double* randVarSamples)
{
  // Evaluate PC expansion at a single point with coefficient vector p

  // Create a single sample in a 2D Array
  Array2D<double> singleSample(1,nDim_,0.e0);
  for(int id=0;id<nDim_;id++)
    singleSample(0,id)=randVarSamples[id];
  // Copy the coefficient array into Array1D
  Array1D<double> pa(nPCTerms_,0.e0);
  for(int ip=0;ip<nPCTerms_;ip++)
    pa(ip)=p[ip];

  // Call the general PC evaluation function for this sample
  Array1D<double> xch(1,0.e0);
  this->EvalPCAtCustPoints(xch,singleSample,pa);

  return xch(0);
}


void PCSet::EvalPCAtCustPoints(Array1D<double>& xch,Array2D<double>& custPoints, Array1D<double>& p)
{
  // Evaluate PC expansion at a set of points with coefficient vector p

  // Check array sizes
  if( (int) custPoints.YSize() != this->nDim_){
    ostringstream buffer;
    buffer << this->nDim_;
    string err_message = (string) "PCSet::EvalPCAtCustPoints: randVarSamples array size does not match the number"
      + " of PC dimensions (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }

  // Check array sizes
  if( (int) p.Length() != this->nPCTerms_){
    ostringstream buffer;
    buffer << this->nPCTerms_;
    string err_message = (string) "PCSet::EvalPCAtCustPoints: p array size does not match the number"
      + " of PC terms (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }


  int nSample=custPoints.XSize();
  xch.Resize(nSample,0.e0);
  Array2D<double> psi;
  this->EvalBasisAtCustPts(custPoints,psi);
  for (int is=0;is<nSample;is++){
    double sum=0.;
    for (int ip=0;ip<this->nPCTerms_;ip++)
      sum += p(ip)*psi(is,ip);
    xch(is)=sum;
  }


  return;
}
 



void PCSet::EvalBasisAtCustPts(const Array2D<double>& custPoints,Array2D<double>& psi)
{
  // Check array sizes
  if((int) custPoints.YSize() != this->nDim_){
    ostringstream buffer;
    buffer << this->nDim_;
    string err_message = (string) "PCSet::EvalBasisAtCustPoints: custPoints array size does not match the number"
      + " of PC dimensions (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }

  int npts=custPoints.XSize();

  // Set array size and fill in with units
  psi.Resize(npts,this->nPCTerms_,1.e0);

  // Evaluate the basis functions at each of the points
  for(int isp=0; isp < npts; isp++){
    for(int id=0;id<nDim_;id++){
      // Define a temporary array to store the basis
      // values in for this point
      Array1D<double> basisVals(maxOrdPerDim_(id)+1);
      p_basis_->EvalBasis(custPoints(isp,id),basisVals);
      for(int ipc=0; ipc < this->nPCTerms_; ipc++){
	psi(isp,ipc)*=basisVals(this->multiIndex_(ipc,id));
      }
    }
   
  }


  return;
}


void PCSet::GalerkProjection(const Array1D<double>& fcn, Array1D<double>& ck)
{

  // Performs Galerkin projection, given function evaluations at quadrature points

  ck.Resize(nPCTerms_,0.e0);

  // Check array sizes
  if( (int) fcn.Length() != this->nQuadPoints_){
    ostringstream buffer;
    buffer << this->nQuadPoints_;
    string err_message = (string) "PCSet::GalerkProjection: fcn array size does not match the number"
      + " of quadrature points (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }

  
  if(uqtkverbose_>1){
    cout << "PCSet::GalerkinProjection() : using " << nQuadPoints_ << " quadrature points" << endl;
  }
  
  // Compute the projection integral
  for(int ipc=0;ipc<nPCTerms_;ipc++){
    
    if(uqtkverbose_>1){
      cout << "PCSet::GalerkinProjection() : Computing " << ipc << " th coefficient out of " << nPCTerms_ << endl;
    }
    
    double num=0.0 ;
    //double den=0.0 ;
    for(int iqp=0;iqp<nQuadPoints_;iqp++){
      num+=fcn(iqp)*this->psi_(iqp,ipc)*this->quadWeights_(iqp);
      // den+=this->psi_(iqp,ipc)*this->psi_(iqp,ipc)*this->quadWeights_(iqp);
    }
    //ck(ipc)=num/den;
     ck(ipc)=num/psiSq_(ipc);
  }

  
  if(uqtkverbose_>1){
    cout << "PCSet::GalerkinProjection() : Done" << endl;
  }
  
  return;

}





void PCSet::GalerkProjectionMC(const Array2D<double>& x, const Array1D<double>& fcn, Array1D<double>& ck)
{
  // Performs Galerkin projection by Monte Carlo integration, given function evaluations at any points

  ck.Resize(nPCTerms_,0.e0);

  int npts = x.XSize();

  // Check array sizes
  if( (int) fcn.Length() != npts){
    string err_message = (string) "PCSet::GalerkProjectionMC: fcn array size does not match the number"
      + " of given points ";
    throw Tantrum(err_message);
  }
  if( (int) x.YSize() != this->nDim_){
    ostringstream buffer;
    buffer << this->nDim_;
    string err_message = (string) "PCSet::GalerkProjectionMC: data array size does not match the number"
      + " of dimensions (" + buffer.str() + ")";
    throw Tantrum(err_message);
  }

  if(uqtkverbose_>1){
    cout << "Galerkin Projection via Monte-Carlo: using " << npts << " points" << endl;
  }
  
  // Compute the projection integral
  for(int ipc=0;ipc<this->nPCTerms_;ipc++){
    
    if(uqtkverbose_>1){
      cout << "Galerkin ProjectionMC: Computing " << ipc << " th coefficient out of " << nPCTerms_ << endl;
    }
    
    Array2D<double> psi;
    this->EvalBasisAtCustPts(x,psi);

    double num=0.0 ;
      
    for(int ip=0;ip<npts;ip++){
      num+=(fcn(ip)*psi(ip,ipc)/npts);
    }
    ck(ipc)=num/psiSq_(ipc);
  }

  if(uqtkverbose_>1){
    cout << "Galerkin ProjectionMC: Done" << endl;
  }

  return;

}



int PCSet::ComputeOrders(Array1D<int>& orders)
{
  int maxOrder=0;
  orders.Resize(this->nPCTerms_,0);

  // Loop over all basis terms
  for(int ipc=0;ipc<this->nPCTerms_;ipc++){
    
    // Add the orders of the current basis term
    for(int id=0;id<this->nDim_;id++)
      orders(ipc) += multiIndex_(ipc,id);
	
    // Check to store the maximal order
    if (orders(ipc)>maxOrder)
      maxOrder=orders(ipc);

  }
 
  return maxOrder;
}


int PCSet::ComputeEffDims(Array1D<int>& effdim)
{
  int maxEffDim=0;
  effdim.Resize(this->nPCTerms_,0);

  // Loop over all basis terms
  for(int ipc=0;ipc<this->nPCTerms_;ipc++){
    
    // Compute the number of dimensions with non-zero order in the current basis term
    for(int id=0;id<this->nDim_;id++)
      effdim(ipc) += (multiIndex_(ipc,id)>0);
	
    // Check and store the maximum dimensionality among all bases
    // Note: this is not the classical definition of effective dimenionality, 
    // since all dimensions can still be involved.
    if (effdim(ipc)>maxEffDim)
      maxEffDim=effdim(ipc);

  }
 
  return maxEffDim;
}

void PCSet::EncodeMindex(Array1D< Array2D<int> >& sp_mindex)
{
 
  // Compute the effective dimensionailities of all terms
  Array1D<int> effdim;
  int maxEffDim=ComputeEffDims(effdim);
        
  // Compute and store the number of terms with a given effective dimensionality
  Array1D<int> num_effdim(maxEffDim+1,0);
  for(int ipc=0;ipc<this->nPCTerms_;ipc++)
    num_effdim(effdim(ipc))+=1;

    
    
  // Resize for storage 
  sp_mindex.Resize(maxEffDim+1);
  for(int i_effdim=0;i_effdim<=maxEffDim;i_effdim++)
    // \todo what if num_effdim==0? (clear instead of resize?) or i_effdim==0? be ready for segfaults!
    sp_mindex(i_effdim).Resize(num_effdim(i_effdim),2*i_effdim);
  Array1D<int> ii(maxEffDim+1,0);
  
    
    
  // Store the multiindices in the sparse format, i.e. 
  // the i-th element of sp_mindex is a 2D int matrix correspoding to bases with i non zero orders
  // each row of that 2D-int matrix corresponds to a basis and has 2*f terms:
  // first f are the indices of the non-zero dimensions, and the second half are the correponding orders.
  // E.g. a term with multiindex (3,0,2,1,2,0) will correspond to a row (0,2,3,4,3,2,1,2) in the 2D array sp_mindex(4)
  for(int ipc=0;ipc<this->nPCTerms_;ipc++){
    int cur_effdim=effdim(ipc);
    
    int jj=0;
    for(int id=0;id<this->nDim_;id++){
        
      if(multiIndex_(ipc,id)){
	sp_mindex(cur_effdim)(ii(cur_effdim),jj)=id;
	sp_mindex(cur_effdim)(ii(cur_effdim),jj+cur_effdim)=multiIndex_(ipc,id);
	jj++;
      }
	    
    }
    ii(cur_effdim)++;

  }
    
    
  return;
}

double PCSet::ComputeMean(Array1D<double>& coef)
{

  // Compute the effective dimensionalities
  Array1D<int> effdim;
  int maxEffDim=ComputeEffDims(effdim);
     
  // Search for the term with effective dimensionality = 0, i.e. the 0th order term
  // Note: most often it is the very first term in the multiindex list, hence this may be an overkill.
  for(int ipc=0;ipc<this->nPCTerms_;ipc++){
    if(effdim(ipc)==0)
      return coef(ipc);
  }

  return 0.e0;
}


double PCSet::ComputeVarFrac(Array1D<double>& coef, Array1D<double>& varfrac)
{
  // Compute the second moment
  varfrac.Resize(this->nPCTerms_,0.e0);
  double var=0.e0;
  for(int ipc=0;ipc<this->nPCTerms_;ipc++){
    varfrac(ipc)=pow(coef(ipc),2.e0)*psiSq_(ipc);
    var+=varfrac(ipc);
  }
  
  // Subtract the mean-squared to compute the variance
  double mean=ComputeMean(coef);
  var-=pow(mean,2.e0);

  // Scale the fractional contributions for all terms
  for(int ipc=0;ipc<this->nPCTerms_;ipc++)
    varfrac(ipc) /= var;


  return var;
}





void PCSet::ComputeMainSens(Array1D<double>& coef, Array1D<double>& mainsens)
{
  // Compute effective dimensionalities for all basis terms
  Array1D<int> effdim;
  int maxEffDim=ComputeEffDims(effdim);
  Array1D<double> varfrac;

  // Compute variance fractions per basis term
  double var=ComputeVarFrac(coef,varfrac);

  mainsens.Resize(this->nDim_,0.e0);
  // Loop over all basis terms
  for(int ipc=0;ipc<this->nPCTerms_;ipc++){
    // Search for univariate terms only
    if (effdim(ipc)==1){
      // Add the variance contributions for the univariate terms
      for(int id=0;id<this->nDim_;id++){
	if ( multiIndex_(ipc,id) != 0 ){
	  mainsens(id)+=varfrac(ipc);
	  break;
	}
      }
    }
  }


  return;
}

void PCSet::ComputeTotSens(Array1D<double>& coef, Array1D<double>& totsens)
{

  // Compute variance fractions
  Array1D<double> varfrac;
  double var=ComputeVarFrac(coef,varfrac);

  totsens.Resize(this->nDim_,0.e0);

  // Loop over all basis terms
  for(int ipc=0;ipc<this->nPCTerms_;ipc++){
    for(int id=0;id<this->nDim_;id++){
      if ( multiIndex_(ipc,id) != 0 ){
	// Add the appropriate variance contribution to the corresponding element in totsens
	totsens(id)+=varfrac(ipc);
      }
    }
    
  }

  return;
}

void PCSet::ComputeJointSens(Array1D<double>& coef, Array2D<double>& jointsens)
{

  // Compute variance fractions
  Array1D<double> varfrac;
  double var=ComputeVarFrac(coef,varfrac);

  jointsens.Resize(this->nDim_,this->nDim_,0.e0);

  
  for(int ipc=0;ipc<this->nPCTerms_;ipc++){
    Array1D<int> nz;
    // Store the order non-zero elements
    for(int id=0;id<this->nDim_;id++)
      if ( multiIndex_(ipc,id) != 0 )
	nz.PushBack(id);

    // Add the appropriate variance contribution to the corresponding element in jointsens
    for(int iz=0;iz<(int) nz.Length();iz++)
      for(int jz=iz+1;jz<(int) nz.Length();jz++)
	jointsens(nz(iz),nz(jz))+=varfrac(ipc);
    
  }

  return;
}


void PCSet::EvalNormSq(Array1D<double>& normsq)
{
  normsq.Resize(nPCTerms_,1.e0);
  
  // Get the 1d norms-squared
  Array1D<double> norms1d;
  p_basis_->Get1dNormsSq(norms1d);

  // FOr each term, multiply appropriate 1d norms-squared
  for(int ipc=0; ipc<nPCTerms_; ipc++)
    for(int id=0; id<nDim_; id++)
      normsq(ipc) *= norms1d(multiIndex_(ipc,id));
  

  return;
}



bool PCSet::IsInDomain(double x)
{
  string type=p_basis_->GetPCType();

  // For LU and JB, the domain is [-1,1]
  if( !strcmp(type.c_str(),"LU") or !strcmp(type.c_str(),"JB") or !strcmp(type.c_str(),"LU_N") )
    if(fabs(x)>1)
      return false;
  
  // For GLG and SW, the domain is [0,\infty)
  if( !strcmp(type.c_str(),"GLG") or !strcmp(type.c_str(),"SW") )
    if(x<0)
      return false;
  
  return true;
}

