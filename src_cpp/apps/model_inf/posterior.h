#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Array1D.h"
#include "Array2D.h"
#include "uqtkmcmc.h"
#include "uqtktools.h"
#include "PCSet.h"

struct modelAux
{
    Array2D<double> domain;
    PCSet* modelpc;//needed only if modeltype="pcx"
    Array2D<double> pccf_all;
    Array1D<double> pccf;
};

struct postAux
{
    Array2D<double> xdata;
    Array2D<double> ydata;
    string modeltype;
    string liktype;
    Array1D<double> likparams;
    int pdim;
    PCSet* pdfpcmodel;//needed only if pdftype="pc"
    modelAux* modelinfo;
};

// Function declarations
void model(Array1D<double>& moutput, Array2D<double>& xdata,Array1D<double>& param, string modeltype, void *modelinfo);


/// \brief Evaluate the log of the posterior for a given set of model
/// and nuisance parameters
double LogPosterior(Array1D<double>& m, void *info);

/// \brief Compute the profile
//void profile_fcn(Array2D<double>& minput,Array1D<double>& modelpar,Array1D<double>& moutput);
//void profile_fcn(Array2D<double>& minput,Array1D<double>& modelpar,Array1D<double>& moutput, modelstruct* pModel);
//void profile_fcn(Array2D<double>& minput,Array1D<double>& modelpar,Array1D<double>& moutput,char* pModel);



double LogPosterior(Array1D<double>& m, void *info)
{
    double pi=4.*atan(1.);
    
    postAux* pinfo = (postAux*) info;
   // Array2D<double> param_samples(pinfo->nsam,pinfo->pdim);
    int ndim = pinfo->xdata.YSize();
    int nx = pinfo->xdata.XSize();
    assert (nx==pinfo->ydata.YSize());
    int neach = pinfo->ydata.XSize();

    
    Array1D<double> moutput;
    model(moutput,pinfo->xdata,m,pinfo->modeltype,(void*)pinfo->modelinfo);

    double var;
    if(pinfo->liktype=="constVar"){
        var=pinfo->likparams(0);
        
    }
  double  loglik=0.0;
 for(int ie=0;ie<neach; ie++){
     for (int ix=0;ix<nx;ix++){
      double err=pinfo->ydata(ix,ie)-moutput(ix);
     loglik  = loglik - 0.5 * log(2*pi) - 0.5*log(var) - pow(err,2)/(2.0*var);
     }

 }
    
    double logprior=0.0;

    double logpost=logprior+loglik;

 return logpost;
 
}



void model(Array1D<double>& moutput, Array2D<double>& xdata,Array1D<double>& param, string modeltype, void *modelinfo)
{
    
    int nx=xdata.XSize();
    int ndim=xdata.YSize();
    int pdim=param.XSize();
    moutput.Resize(nx);
    
    
    
    if (modeltype=="prop"){
        assert(pdim==1);
        assert(ndim==1);
        
            for (int j=0;j<nx;j++)
                moutput(j)=param(0)*xdata(j,0);
    }
    else if (modeltype=="pcpar"){
        modelAux* minfo = (modelAux*) modelinfo;
        assert(pdim==minfo->modelpc->GetNDim());
        
        for (int j=0;j<nx;j++){
            Array1D<double> pccf_this;
            getCol(minfo->pccf_all,j,pccf_this);
            moutput(j)= minfo->modelpc->EvalPC(pccf_this,param);

    
        }
    }
    
 
//    else if (modeltype=="linear"){
//        assert(pdim==2);
//        assert(ndim==1);
//        
//        for (int i=0;i<nsam;i++)
//            for (int j=0;j<nx;j++)
//                model_samples(i,j)=param_samples(i,0)+param_samples(i,1)*xdata(j,0);
//    }
//    
//    else if (modeltype=="pcx"){
//        
//        modelAux* minfo = (modelAux*) modelinfo;
//        
//        assert((pdim+ndim)==minfo->modelpc->GetNDim());
//        Array2D<double> pcinput(nx*nsam,pdim+ndim);
//        for(int i=0;i<nx;i++){
//            for(int k=0;k<nsam;k++){
//                
//                
//                for(int j=0;j<pdim;j++)
//                    pcinput(k+i*nsam,j)=-1.+2.*( param_samples(k,j)-minfo->domain(j,0) ) / (minfo->domain(j,1)-minfo->domain(j,0));
//                
//                for(int j=0;j<ndim;j++)
//                    pcinput(k+i*nsam,pdim+j)=-1.+2.*( xdata(i,j)-minfo->domain(pdim+j,0) ) / (minfo->domain(pdim+j,1)-minfo->domain(pdim+j,0));
//                
//                
//                
//            }
//            
//        }
//        
//        Array1D<double> model_samples_planar;
//        minfo->modelpc->EvalPCAtCustPoints(model_samples_planar,pcinput,minfo->pccf);
//        
//        for (int i=0;i<nsam;i++)
//            for (int j=0;j<nx;j++)
//                model_samples(i,j)=model_samples_planar(i+j*nsam);
//    }
    
    return;
}


//
//void profile_fcn(Array2D<double>& minput,Array1D<double>& modelpar,Array1D<double>& moutput, modelstruct* pModel)
//{
//  int indim=minput.YSize();
//  int nt=minput.XSize();
//
//  int npar=modelpar.XSize();
//  Array2D<double> xdata(1,npar,0.e0);
//
//  for(int i=0;i<npar;i++){
//    double a=(pModel->pdomain)(i,0);
//    double b=(pModel->pdomain)(i,1);
//    xdata(0,i)=(modelpar(i)-(a+b)/2.) / ( (b-a)/2. );
//  }
//
//  moutput.Clear();
//
//  for(int j=0;j<nt;j++){
//    Array1D<double> c_k;
//
//    getRow(pModel->pccf_all,j,c_k);
//
//    Array1D<double> ydata;
//    PCSet pcm("NISPnoq",pModel->mindex,"LEG",0.0,1.0);
//
//    Array2D<int> mi;
//    pcm.GetMultiIndex(mi);
//
//    pcm.EvalPCAtCustPoints(ydata,xdata,c_k);
//
//    moutput.PushBack(ydata(0));
//  }
//
//  return;
//}
//
//void profile_fcn(Array2D<double>& minput,Array1D<double>& modelpar,Array1D<double>& moutput,char* pModel)
//{
//
//int indim=minput.YSize();
//  int nt=minput.XSize();
//
//  moutput.Resize(nt,0.e0);
//
//  write_datafile(minput,"minput.dat");
//  write_datafile_1d(modelpar,"mparam.dat");
//  system(pModel);
//
//
//  read_datafile_1d(moutput,"moutput.dat");
//
//
//  return;
//}
//
//
//
//void profile_fcn(Array2D<double>& minput,Array1D<double>& modelpar,Array1D<double>& moutput)
//{
//  int indim=minput.YSize();
//  assert(indim==1);
//  int nt=minput.XSize();
//
//  moutput.Resize(nt,0.e0);
//
//
//  double lam=modelpar(0);
//
//  for (int it=0;it<nt;it++){
//    moutput(it)=exp(-lam*minput(it,0));
//  }
//
//
//
//  return;
//}
//
//
