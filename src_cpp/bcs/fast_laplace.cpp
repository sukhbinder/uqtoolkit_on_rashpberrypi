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
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "assert.h"
#include <sstream>
#include <fstream>

#include "fast_laplace.h"
#include "uqtktools.h"
#include "ftndefs.h"
#include "lapack.h"


/// \brief The implementation of the Bayesian COmpressive Sensing algorithm using Laplace Priors
/// \note This function has been written relying on the algorithm and MATLAB code presented in
/// http://ivpl.eecs.northwestern.edu/research/projects/bayesian-compressive-sensing-using-laplace-priors
/// and references therein
/// \todo The array manipulations are not optimized - perhaps they need to be reconsidered using, 
/// say, fortran matrix-vector manipulation routines
void FastLaplace(Array2D<double> &PHI, Array1D<double> &y, double &sigma2,
                 double eta, Array1D<double> &lambda_init, 
		 int adaptive, int optimal, double scale, int verbose,
                 Array1D<double> &weights, Array1D<int> &used,
                 Array1D<double> &errbars, Array1D<double> &basis, 
                 Array1D<double> &alpha, double &lambda) 
{
  int n = (int) PHI.XSize() ;
  int m = (int) PHI.YSize() ;

    
  /* find initial alpha */
  Array2D<double> PHIT ;

  transpose(PHI,PHIT);

  Array1D<double> PHIy ;
  prodAlphaMatVec(PHIT, y, 1.0, PHIy);

  Array1D<double> PHI2(m,0.0);
  for (int j = 0; j<m; j++)
    for (int i = 0; i<n; i++)
      PHI2(j) += (PHI(i,j)*PHI(i,j)) ;

  Array1D<double> ratio(m,0.0);
  for (int j = 0; j<m; j++)
    ratio(j) = pow(PHIy(j),2)/PHI2(j);

  int indx ;
  double maxr = maxVal(ratio,&indx) ;
  //  printf("max=%e, indx=%d\n",maxr,indx); 

  Array1D<int>index(1,indx);
  alpha.Resize(1) ;
  alpha(0)= PHI2(index(0))/(maxr-sigma2);
  //printf("alpha=%e\n",alpha(0));

  /* Compute initial mu, Sig, S, Q */
  Array2D<double> phi(n,1,0.0);
  for (int i = 0; i<n; i++ ) phi(i,0) = PHI(i,index(0));
  double Hessian=alpha(0);
  for (int i = 0; i<n; i++ ) Hessian += phi(i,0)*phi(i,0)/sigma2;
  Array2D<double> Sig(1,1,1.0/Hessian) ;
  double dtmp1=Sig(0,0)*PHIy(index(0))/sigma2;
  Array1D<double> mu(1,dtmp1);
  Array1D<double> left(m,0.0),phitmp(n);
  for ( int i=0; i<n; i++ )phitmp(i)=phi(i,0);
  prodAlphaMatVec(PHIT, phitmp, 1.0/sigma2, left);

  Array1D<double> S(m), Q(m);
  for (int j=0; j<m; j++)
  {
    S(j) = PHI2(j)/sigma2-Sig(0,0)*left(j)*left(j);
    Q(j) = PHIy(j)/sigma2-Sig(0,0)*PHIy(index(0))*left(j)/sigma2;
  }
     
  /* Keep track of the positions selected during the iterations */
  Array1D<int> selected;
  selected=index;
  Array1D<int> deleted ;
  Array1D<double> ML(MAX_IT,0) ;

  int count ;
  for ( count=0; count<MAX_IT; count++ )
  {
    if ( (count%10)==0 )
	 cout << "Iteration # " << count << endl;

    Array1D<double> s,q;
    s=S; q=Q;
    for ( int i = 0; i< (int) index.XSize(); i++)
    {
      s(index(i)) = alpha(i)*S(index(i))/(alpha(i)-S(index(i)));
      q(index(i)) = alpha(i)*Q(index(i))/(alpha(i)-S(index(i)));
       
    }
      
     
    if (lambda_init.XSize() == 0)
    {
      double suma=0.0; 
      for ( int i = 0; i< (int) alpha.XSize(); i++) suma += 1.0/alpha(i);
      lambda = 2*( index.XSize() - 1 ) / suma;
    }
    else
      lambda = lambda_init(0);

    Array1D<double> A(m,0.0), B(m,0.0), C(m,0.0);
    Array1D<double> theta(m,0.0), discriminant(m,0.0), nextAlphas(m,0.0) ;
    for ( int i=0; i<m; i++ )
    {
      A(i) = lambda + s(i) - pow(q(i),2);
      B(i) = 2.0*lambda*s(i) + pow(s(i),2);
      C(i) = lambda*pow(s(i),2);
      theta(i) = pow(q(i),2)-s(i);    
      discriminant(i) = pow(B(i),2) - 4.0*A(i)*C(i);
      nextAlphas(i) = (-B(i) - sqrt(discriminant(i)) ) / (2.0*A(i));
    }

    /* choose the next alpha that maximizes marginal likelihood */
    Array1D<double> ml(m,-1.0e20);
    Array1D<int> ig0 ;

    find(theta, lambda, string("gt"), ig0) ;

    /* indices for reestimation */
    Array1D<int> ire,foo,which ;
    intersect(ig0,index,ire,foo,which);

    if (ire.XSize() > 0)
    {
      Array1D<double> Alpha(ire.XSize(),0.0);
      Array1D<double> delta(ire.XSize(),0.0);
      for ( int i=0; i< (int) ire.XSize(); i++ )
      {
        Alpha(i) = nextAlphas(ire(i)); 
        delta(i) = (alpha(which(i))-Alpha(i))/(Alpha(i)*alpha(which(i)));
        double As=Alpha(i) / (Alpha(i) + s(ire(i)));
          double as=alpha(which(i)) / (alpha(which(i)) + s(ire(i)));
          if (As<=0.0) As=1e-20;
          if (as<=0.0) as=1e-20;
        ml(ire(i)) = pow(q(ire(i)),2)/ (Alpha(i) + s(ire(i)))
                     + log(As) - lambda/ Alpha(i) 
                     - pow(q(ire(i)),2)/ (alpha(which(i)) + s(ire(i))) 
                     - log(as) 
                     + lambda/alpha(which(i));
	
          
      }
    }

    /* indices for adding */
    Array1D<int> iad ;

    setdiff(ig0, ire, iad) ;

    //printf("%d %d iad:%d\n",ig0.XSize(), ire.XSize(),(int) iad.XSize());
    if ( iad.XSize() > 0 )
    {   
      Array1D<double> Alpha(iad.XSize(),0.0);
      for ( int i=0; i< (int) iad.XSize(); i++ )
      {
        Alpha(i) = nextAlphas(iad(i));
        ml(iad(i)) = log(Alpha(i) / (Alpha(i) + s(iad(i))) ) 
	  +pow(q(iad(i)),2) / (Alpha(i) + s(iad(i))) - lambda/Alpha(i);
      }
      intersect(deleted, iad, which);

        for ( int i=0; i< (int) which.XSize(); i++ )
        ml(which(i)) = -1.0e20;
        
    }

    Array1D<int> indxM(m,0); for ( int i=0; i<m; i++ ) indxM(i) = i;
    Array1D<int> is0;

      setdiff_s(indxM,ig0,is0);

    // printf("%d %d is0:%d\n",indxM.XSize(), ig0.XSize(), (int) is0.XSize());

    /* indices for deleting */
    Array1D<int> ide;
    intersect(is0,index,ide,foo,which);

    //printf("ide:%d\n",(int) ide.XSize());
    if ( ide.XSize() > 0 )
    {
      if ( index.XSize() == 1 )
	for ( int i=0; i< (int) ide.XSize(); i++ )
	  ml(ide(i)) = -1.0e200;
      else
	for ( int i=0; i< (int) ide.XSize(); i++ )
	  ml(ide(i)) = -pow(q(ide(i)),2) / (alpha(which(i)) + s(ide(i))) 
	    -log( alpha(which(i)) /(alpha(which(i)) + s(ide(i)))) 
	    +lambda/alpha(which(i));
    }

    int idx ;
    ML(count) = maxVal(ml,&idx);
   

      /* check convergence */
    if (count > 1) 
      if ( fabs(ML(count)-ML(count-1)) < fabs(ML(count)-ML(0))*eta )
        break;
    
    /*
      update alphas
      Choose the basis which results in the largest increase in the
      likelihood
    */

    find(index,idx,string("eq"),which);

    //printf("which:%d\n",(int)which.XSize());
    //printf("theta(%d)=%e,lambda=%e\n",idx,theta(idx),lambda);
    
    if ( theta(idx) > lambda )
    {
      if ( which.XSize() > 0 ) 
      {
        /* reestimate a basis */
        double Alpha = nextAlphas(idx);
        double Sigii = Sig(which(0),which(0)); 
        double mui   = mu(which(0)); 
        Array1D<double> Sigi(Sig.XSize(),0.0) ;
	for ( int i=0; i<(int) Sig.XSize(); i++ ) Sigi(i) = Sig(i,which(0));
	double delta = Alpha-alpha(which(0));
	double ki    = delta/(1+Sigii*delta);
        for ( int i = 0; i < (int) mu.XSize() ; i++ )
	  mu(i) = mu(i)-ki*mui*Sigi(i);
        for ( int j = 0; j < (int) Sig.YSize() ; j++ )
	  for ( int i = 0; i < (int) Sig.XSize() ; i++ )
	    Sig(i,j) = Sig(i,j)-ki*Sigi(i)*Sigi(j);
	Array1D<double> tmp1, comm;

        prodAlphaMatVec (phi,Sigi,1.0,       tmp1);
        prodAlphaMatTVec(PHI,tmp1,1.0/sigma2,comm);
          addVecAlphaVecPow(S,ki,    comm,2);
          addVecAlphaVecPow(Q,ki*mui,comm,1);
          alpha(which(0)) = Alpha;
          
        if ( verbose > 0 )
          printf("Reestimate %d..\n",idx);
      }
      else 
      {
        /* add a basis */   
        double Alpha = nextAlphas(idx);
        Array1D<double> phii(PHI.XSize(),0.0);
	for ( int i = 0; i< (int) phii.XSize(); i++ ) phii(i) = PHI(i,idx); 
	double Sigii = 1.0/(Alpha+S(idx)); 
	double mui   = Sigii*Q(idx);
	Array1D<double> comm1,tmp1;

	prodAlphaMatTVec(phi,phii,1.0/sigma2,tmp1);

        prodAlphaMatVec (Sig,tmp1,1.0,       comm1);
	Array1D<double> ei;
	ei=phii;
      //    cout << "hook A 55" << endl;

        
	prodAlphaMatVec(phi,comm1,-1.0,tmp1);// !!! slowdown later
     //     cout << "hook A 57" << " " << phi.XSize() << " X " << phi.YSize() << " " << comm1.Length() << endl;

	addVecAlphaVecPow(ei,1.0,tmp1,1);

	Array1D<double> off(comm1.XSize(),0.0);
	addVecAlphaVecPow(off,-Sigii,comm1,1);

	for ( int j = 0; j< (int) Sig.YSize(); j++ )
	  for ( int i = 0; i< (int) Sig.XSize(); i++ )
            Sig(i,j) += Sigii*comm1(i)*comm1(j) ;
	paddMatColScal(Sig,off,Sigii) ;
//	cout << "hook A 60" << endl;

	for ( int i = 0; i< (int) mu.XSize(); i++ ) mu(i) -= mui*comm1(i);
	mu.PushBack(mui);

	Array1D<double> comm2(m,0.0);
	prodAlphaMatTVec(PHI,ei,1.0/sigma2,comm2);// !!!early slowdown here
 //         cout << "hook A 70" << " " << PHI.XSize() << " X " << PHI.YSize() << " " << ei.Length() << endl;


	addVecAlphaVecPow(S,-Sigii,comm2,2);
  //        cout << "hook A 71" << endl;

	addVecAlphaVecPow(Q,-mui,comm2,1);
 //         cout << "hook A 72" << endl;

	index.PushBack(idx);
	alpha.PushBack(Alpha);

	//paddMatCol(phi,phii); // !!! slowdown later---fixed!!!
          phi.insertCol(phii,phi.YSize());
          
//          cout << "hook A 73" << endl;

          
        if ( verbose > 0 )
          printf("Add %d.. \n",idx) ;
      }
    }
    else
    {
      if ( ( which.XSize() > 0 ) && ( index.XSize()> 1 ) )
      { 
	/* delete a basis */
        deleted.PushBack(idx);
        double Sigii = Sig(which(0),which(0)); 
        double mui   = mu(which(0)); 
        Array1D<double> Sigi ;
	getCol(Sig,which(0),Sigi) ;
	for ( int j = 0; j< (int) Sig.YSize(); j++ )
	  for ( int i = 0; i< (int) Sig.XSize(); i++ )
	    Sig(i,j) -= Sigi(i)*Sigi(j)/Sigii ;
	//	  cout << "hook AAAA" << endl;
	
	delCol(Sig,which(0));
	delRow(Sig,which(0));
	addVecAlphaVecPow(mu,-mui/Sigii,Sigi,1);
	delCol(mu,which(0));
	
	Array1D<double> comm,tmp1 ;
	prodAlphaMatVec(phi,Sigi,1.0/sigma2,tmp1);
        prodAlphaMatTVec(PHI,tmp1,1.0,comm);
	addVecAlphaVecPow(S,1.0/Sigii,comm,2);
	addVecAlphaVecPow(Q,mui/Sigii,comm,1);
        
	delCol(index,which(0));
	delCol(alpha,which(0));
	delCol(phi,  which(0));
	//  cout << "hook BBBB" << endl;

        if ( verbose > 0 )
          printf("Delete %d.. \n",idx);
      }
      else if ((which.XSize() > 0 ) && (index.XSize() == 1))
      {
        printf("Something is wrong, trying to delete the only coefficient\n");
        printf("that has been added.\n");
	break;
      }
    }
    
    selected.PushBack(idx);
    //for ( int j = 0; j< (int) selected.XSize(); j++ )
    //printf("selected(%d):%d\n",j+1,selected(j));    
  }

  //for ( int j = 0; j< (int) mu.XSize(); j++ )
  //printf("mu(%d):%f\n",j+1,mu(j));
  //for ( int j = 0; j< (int) mu.XSize(); j++ )
  //printf("index(%d):%d\n",j+1,index(j));
  weights = mu    ;
  used    = index ;

  /* re-estimated sigma2 */
  double sum1=0.0, sum3 = 0.0;
  for ( int i = 0; i<(int) y.XSize(); i++)
  {
    double sum2 = 0.0 ;
    for ( int j = 0; j<(int) phi.YSize(); j++)
      sum2 += phi(i,j)*mu(j) ;
    sum1 += pow(y(i)-sum2,2) ;
  }
  for ( int i = 0; i<(int) alpha.XSize(); i++)
    sum3 += alpha(i)*Sig(i,i) ;
  sigma2 = sum1/((double) (n-index.XSize())+sum3) ;

  errbars.Resize(Sig.XSize()) ;
  for ( int i = 0; i<(int) errbars.XSize(); i++) errbars(i) = sqrt(Sig(i,i));

  //for ( int i = 0; i<(int) errbars.XSize(); i++) 
  // printf("eb(%2d)=%e\n",i+1,errbars(i));

  /* generate a basis for adaptive CS? */
  //printf("%d %d\n",adaptive,optimal);
  if ( adaptive == 1 )
  {
    if ( optimal == 1 )
    {
      int nSig = Sig.XSize() ;
      double VL,VU ;
      int IL=nSig,IU=nSig;
      double ABSTOL = -1.0 ;
      int    nEigVals, info, ifail ;
      double eigVal ;
      basis.Resize(nSig) ;
      int lWrk = 8*nSig ;
      double *dwrk = new double[lWrk] ;
      int    *iwrk = new int   [lWrk] ;
      FTN_NAME(dsyevx)((char *)"V", (char *)"I", (char *)"L", &nSig, Sig.GetArrayPointer(), &nSig, 
                       &VL, &VU, &IL, &IU, &ABSTOL, 
                       &nEigVals, &eigVal, basis.GetArrayPointer(), 
		       &nSig, dwrk, &lWrk, iwrk, &ifail, &info ) ;
      delete [] dwrk ;
      delete [] iwrk ;
      if ( info != 0 )
	printf("WARNING : DSYEVX failed with info=%d\n",info) ;
    }
    else
    {
      Array2D<double> tmp1;
      prodAlphaMatTMat(phi,phi,1.0/sigma2,tmp1);
      double tmp1m = 0.0 ;
      for ( int i = 0; i < (int) tmp1.XSize() ; i++ )
        tmp1m += tmp1(i,i); 
      tmp1m /= ( (double) tmp1.XSize() ) ;
      Array2D<double> Sig_inv ;
      Sig_inv = tmp1 ;
      for ( int i = 0 ; i < (int) Sig_inv.XSize() ; i++ )
	Sig_inv(i,i) += scale*tmp1m;

      int nSig = Sig_inv.XSize() ;
      double VL,VU ;
      int IL=1,IU=1;
      double ABSTOL = -1.0 ;
      int    nEigVals, info, ifail ;
      double eigVal ;
      basis.Resize(nSig) ;
      int lWrk = 8*nSig ;
      double *dwrk = new double[lWrk] ;
      int    *iwrk = new int   [lWrk] ;
      FTN_NAME(dsyevx)((char *)"V", (char *)"I", (char *)"L", &nSig, Sig_inv.GetArrayPointer(), &nSig, 
               &VL, &VU, &IL, &IU, &ABSTOL, 
               &nEigVals, &eigVal, basis.GetArrayPointer(), 
	       &nSig, dwrk, &lWrk, iwrk, &ifail, &info ) ;
      delete [] dwrk ;
      delete [] iwrk ;
      if ( info != 0 )
	printf("WARNING : DSYEVX failed with info=%d\n",info) ;
    }
    //printf("%d\n",basis.XSize());
    //for ( int i = 0; i<(int) basis.XSize(); i++) 
    //  printf("basis(%2d)=%e\n",i+1,basis(i));

  }
  printf("BCS algorithm converged, # iterations : %d \n",count);
  return ;

}
