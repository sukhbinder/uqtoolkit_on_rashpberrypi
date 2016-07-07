/* infer_model.cpp */
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>

#include "posterior.h"

#include "uqtkmcmc.h"

#include "uqtktools.h"
#include "Array1D.h"
#include "Array2D.h"

using namespace std;

#define XFILE "xdata.dat"
#define YFILE "ydata.dat"
#define MODELSCRIPT "model.x"
#define NMCMC 1000
#define CHST 10


int usage()
{
  printf("usage: infer_model [-h] [-x<xfile>] [-y<yfile>] [-m<model.x>] [-e<exactPostGrid>] [-f<exactPostGrid>] [-n<chainstates>] [-s<chainsteps>]\n");
  printf(" -h                 : print out this help message \n");
  printf(" -x <xfile>         : file path for x values (default=%s) \n", XFILE);
  printf(" -y <yfile>         : file path for y values, i.e. function evaluations at x (default=%s) \n", YFILE);
  printf(" -m <model.x>       : file path to model script (default=%s) \n", MODELSCRIPT);
  printf(" -e <exactPostGrid> : computes the exact posterior predictions at values given by file <exactPostGrid> \n");
  printf(" -f <exactPostGrid> : same as above, but skip the MCMC step\n");
  printf(" -n <nmcmc>   : number of chainstates for MCMC (default=%d) \n", NMCMC);
  printf(" -s <chainsteps>    : number of chainsteps in MCMC (default=%d) \n", CHST);
  printf("================================================================================\n");
  printf("Input:: MCMC chain input file 'mcmc.xml'\n");
  printf("        In case modelx == \"PC\": 'mindex.dat', 'pccf_all.dat', 'pdomain.dat' \n");
  printf("Output:: 'postdens.dat', 'model_MAP.dat', 'model_chain.dat'\n");

  printf("--------------------------------------------------------------------------------\n");
  printf("Comments: None yet.\n");
  printf("Complexity: Not tested yet.\n");
  printf("Todo: -Implement adaptive domain splitting, either here or using a python wrapper.\n");
  printf("Todo: -Set up and point to some test runs.\n");
  printf("================================================================================\n");
  exit(0);
  return 0;
}

/******************************************************************************/
/// \brief  Main program of inferring input parameter given data
int main (int argc, char *argv[]) {

  // Declare main variables
  int c;
  char* xfile;
  char* yfile;
  char* modelscript;
  char* pargridfile;
  bool eflag=false;
  bool fflag=false;
  int nmcmc=NMCMC;
  int chst=CHST;

  // Set the defaults
  xfile=XFILE;
  yfile=YFILE;
  modelscript=MODELSCRIPT;

  // Read the user input
  while ((c=getopt(argc,(char **)argv,"hx:y:m:e:f:n:s:"))!=-1){
    switch (c) {
    case 'h':
      usage();
      break;
    case 'x':
       xfile =  optarg;	
       break; 
    case 'y':
       yfile = optarg;	
       break; 
    case 'm':
       modelscript = optarg;	
       break;
    case 'e':
       eflag=true;
       pargridfile = optarg;	
       break;    
    case 'f':
       fflag=true;
       pargridfile = optarg;	
       break;   
    case 'n':
       nmcmc =  strtol(optarg, (char **)NULL,0);
       break;
    case 's':
       chst =  strtol(optarg, (char **)NULL,0);	
       break;
    default :
      break;
    }
  }
  
  // Printing user input
  printf("\ninfer_model() : parameters ================================= \n");
  printf("X-data file      = %s \n",xfile);
  printf("Y-data file      = %s \n",yfile);
  printf("Model script     = %s \n",modelscript);
  if (nmcmc==0)
    printf("Postprocessing only the MAP value of the chain \n");
  else  
    printf("Postprocessing %d states every %d step from the chain \n",nmcmc,chst);
  

  if (eflag)
    printf("Exact posterior will be computed at points in file %s \n",pargridfile);

  if (fflag)
    printf("No MCMC: only exact posterior will be computed at points in file %s \n",pargridfile);

  
  
 
  
//  // Merge X- and Y- data
//  Array2D<double> xydata;
//  xydata=xdata;
//  xydata.insertCol(ydata,indim);
//
//  // Append mean and variance
//  Array1D<double> ydata_mean;
//  Array1D<double> ydata_std;
//  Array2D<double> ydata_st;
//  
//  for(int k=0;k<nt;k++){
//  Array1D<double> ydata_cur;
//  getRow(ydata, k, ydata_cur);
//  ydata_mean.PushBack(get_mean(ydata_cur));
//  ydata_std.PushBack(get_std(ydata_cur));
// }
//  
// /*
//  read_datafileVS(ydata_st,"output.st.avestd.dat");
//  getCol(ydata_st,0,ydata_mean);
//  getCol(ydata_st,1,ydata_std);
//  */
//
//
//  xydata.insertCol(ydata_mean,indim+neach);
//  xydata.insertCol(ydata_std,indim+neach+1);
//
//  // Append more info (this is an overkill - a whole vector is appended just to pass a single value)
//  Array1D<double> intinfo(nt,neach);
//  xydata.insertCol(intinfo,indim+neach+2);
// 
//  Array2D<double> info(nt,0);
//  // Casting the data and model pointers
//  void* datapointer=(void*) (&xydata);
//  //  void* infopointer=(void*) (&info);

    postAux* postinfo = new postAux;
    
    
    read_datafileVS(postinfo->xdata,xfile);
    read_datafileVS(postinfo->ydata,yfile);

    Array1D<double> likparams(1,0.1);
    postinfo->modeltype="pcpar";
    postinfo->liktype="constVar";
    postinfo->likparams=likparams;

    // Get problem sizes
    int nx=postinfo->xdata.XSize();
    int ndim=postinfo->xdata.YSize();
    int neach=postinfo->ydata.XSize();
    assert(nx==postinfo->ydata.YSize());

    
    // Print info
    printf("Number of points     = %d \n",nx);
    printf("X-dimensionality     = %d \n",ndim);
    printf("Number of replicas   = %d \n",neach);

    
    int pdim=1;
    double gamma=0.1;
    // Set the model
    modelAux* minfo = new modelAux;
    
    Array2D<int> mi_dummy;
    read_datafileVS(mi_dummy,"pmindex.dat");

    PCSet surrmodel_dummy("NISPnoq",mi_dummy,"LU");
    PCSet* pcptr=&surrmodel_dummy;
    pcptr->PrintMultiIndex();
    //minfo->modelpc=pcptr;
    
    if (postinfo->modeltype=="pcpar"){
        read_datafileVS(minfo->domain,"pdomain.dat");
        postinfo->pdim=minfo->domain.XSize();
        Array2D<int> mindex;
        read_datafileVS(mindex,"pmindex.dat");
        PCSet surrmodel("NISPnoq",mindex,"LU");
        //surrmodel_dummy=surrmodel;
            cout << "B" << endl;
        pcptr->PrintMultiIndex();
//        pcptr=&surrmodel;
        minfo->modelpc= &surrmodel_dummy;
        read_datafileVS(minfo->pccf_all,"ppccf.dat");
        assert(minfo->pccf_all.XSize()==surrmodel.GetNumberPCTerms());
        assert(minfo->pccf_all.YSize()==nx);

        
    //DANGER: if surrmodel is defined inside if{}, then posterior model evaluation breaks, the pointer is lost.
   }
//    else if (postinfo->modeltype=="pcxpar"){
//        read_datafileVS(minfo->domain,"xpdomain.dat");
//        postinfo->pdim=minfo->domain.XSize()-ndim;
//        Array2D<int> mindex;
//        read_datafileVS(mindex,"xpmindex.dat");
//        PCSet surrmodel("NISPnoq",mindex,"LU");
//        
//        minfo->modelpc=&surrmodel;
//        minfo->pccf.Resize(surrmodel.GetNumberPCTerms());
//       
//        read_datafile_1d(minfo->pccf,"xppccf.dat");
//    }
  
    cout << "C" << endl;
    minfo->modelpc->PrintMultiIndex();


    int chdim=pdim;
    Array1D<double> chstart(chdim,0.01);
    void* pinfo=(void*) postinfo;
    MCMC mchain(LogPosterior,pinfo);
    
    mchain.setChainDim(chdim);
    mchain.initAMGamma(gamma);
    mchain.initAdaptSteps(nmcmc/10,10,nmcmc);
    
    mchain.printChainSetup();
    cout << "D" << endl;

    mchain.runOptim(chstart);
    cout << "D" << endl;

    if (nmcmc>0)
        mchain.runChain(nmcmc, chstart);
    cout << "===========================================================" << endl;
    Array1D<double> mapparams;
    mchain.getMode(mapparams);
    write_datafile_1d(mapparams,"MAPparams.dat");

    
//  void* modelpointer;
 //   modelstruct modelpc;
//  
//    string modelscript_str(modelscript);
//  if (modelscript_str=="PC"){
//    read_datafileVS(modelpc.mindex,"mindex.dat");
//    read_datafileVS(modelpc.pccf_all,"pccf_all.dat");
//    read_datafileVS(modelpc.pdomain,"pdomain.dat");
//    modelpc.modelname=modelscript;
//   modelpointer=(void*) (&modelpc);
//  }
//  else if  (modelscript_str=="custom"){
// modelpc.modelname=modelscript;
//    modelpointer=(void*) (&modelpc);
//  }
//  else{
//    modelpc.modelname=modelscript;
//    modelpointer=(void*) (&modelpc);
//  }
//  // Initialize the chain
//  MChain aMCMC(datapointer,modelpointer,LogLikelihood);
//
//  // Read the xml input
//  aMCMC.ReadXMLChainInput("mcmc.xml");

 

    
  return 0;
}
