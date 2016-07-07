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
#include <assert.h>

#include "uqtktools.h"
#include "multiindex.h"


int computeNPCTerms(int ndim,int norder)
{
  if (norder==-1)
    return 0;

  int enume=1;
  int denom=1;

  int minNO = min(norder,ndim);

  for(int k=0; k < minNO; k++){
      enume = enume*(norder+ndim-k);
      denom = denom*(k+1);
  }

  int nPCTerms = enume/denom;
  
  return nPCTerms;
}


int computeMultiIndex(int ndim,int norder, Array2D<int>& mi)
{
  if (ndim==0) return 1;

  // Compute the number of PC terms
  int npc=computeNPCTerms(ndim,norder);

  // Initialize multiIndex
  int iup=0;
  int isum=0;

  // Work arrays
  Array1D<int> ic(ndim,0);
  Array1D<int> ict(ndim,0);

  mi.Resize(npc,ndim,0);


  //-----------zero-th order term-------------------------
  iup=0;
  for(int idim=0;idim < ndim;idim++){
    mi(iup,idim)=0;
  }
  if(norder > 0){
  //-----------first order terms---------------------------
    for(int idim=0;idim < ndim;idim++){
      iup++;
      mi(iup,idim)=1; //multiIndex is a kronecker delta
      ic(idim)=1;
    }
  }
  if(norder > 1){
  //-----------higher order terms--------------------------
    for(int iord=2;iord<norder+1;iord++){
      int lessiord=iup;                   //number of terms of order less than iord
      for(int idim=0;idim<ndim;idim++){
        isum=0;
        for(int ii=idim;ii<ndim;ii++){
          isum=isum+ic(ii);
        }
        ict(idim)=isum;
      }

      for(int idim=0;idim < ndim;idim++){
        ic(idim)=ict(idim);
      }

      for(int idimm=0;idimm<ndim;idimm++){
        for(int ii=lessiord-ic(idimm)+1;ii<lessiord+1;ii++){
          iup++;
          for(int idim=0;idim < ndim;idim++){
            mi(iup,idim)=mi(ii,idim);
          }
          mi(iup,idimm)=mi(iup,idimm)+1;
        }
      }
    }
  }

  return npc;
}


int computeNPCTermsHDMR(int ndim,  Array1D<int>& maxorders)
{
  int nhdmr=maxorders.XSize()-1;

  int cnt=0;
  for(int i=0;i<=nhdmr;i++){
    cnt+= ( choose(maxorders(i),i)*choose(ndim,i) );
  }

  return cnt;
}

int computeMultiIndexHDMR(int ndim,Array1D<int>& maxorders, Array2D<int>& mindex)
{
  int nhdmr=maxorders.XSize()-1;
  int npc=computeNPCTermsHDMR(ndim,maxorders);

  mindex.Resize(npc,ndim,0);
  int iup=1;
  for(int i=1; i<=nhdmr;i++){
    Array2D<int> ind;
   

    chooseComb(ndim,i,ind);
   
    if (maxorders(i)<i)
      return iup;
    else{
      Array2D<int> mi;
      computeMultiIndex(i,maxorders(i)-i,mi);
      for(int j=0;j<(int)ind.XSize();j++){  
    for(int k=0; k<(int)mi.XSize();k++){
        for(int ii=0;ii<i;ii++){
          mindex(iup,ind(j,ii))=mi(k,ii)+1;
        }
        iup++;
      }


    }
    }
  }

  assert(iup==npc);
  return npc; //or iup;
}


void decodeMindex(Array1D< Array2D<int> >& sp_mindex, int ndim, Array2D<int>& mindex){
    
    int npc=0;
    for(int i=0;i<(int) sp_mindex.XSize();i++){
        npc+=sp_mindex(i).XSize();
    }
    mindex.Resize(npc,ndim,0);
    int ipc=0;
    for(int i=0;i<(int) sp_mindex.XSize();i++){
        for(int j=0;j<(int) sp_mindex(i).XSize();j++){
        for(int i_effdim=0;i_effdim<i;i_effdim++){
            
        
            mindex(ipc,sp_mindex(i)(j,i_effdim))=sp_mindex(i)(j,i_effdim+i);
        }
            ipc++;
        }
    }

    assert(ipc==npc);
    
    return;
}



bool is_admis(Array1D<int>& mindex_try,Array2D<int>& mindex){
    
    bool admis=true;
    int npc=mindex.XSize();
    int ndim=mindex.YSize();
    
    Array1D<int> tmp;
    tmp=mindex_try;
    
    for(int j=0;j<ndim;j++){
        
        if(mindex_try(j)>0){
            tmp(j)--;
            for(int ipc=0;ipc<npc;ipc++){
                
                Array1D<int> cur_mindex;
                getRow(mindex,ipc,cur_mindex);
                
                if(is_equal(tmp,cur_mindex)){
                    admis=true;
                    break;
                }
                admis=false;
            }
            
            if(admis==false)
                break;
            
            tmp(j)++;
        }
        
    }
    
    return admis;
}


void upOrder(Array2D<int>& mindex,Array2D<int>& new_mindex){
    
    int npc=mindex.XSize();
    int ndim=mindex.YSize();


    
    Array1D<int> orders;
    getOrders(mindex,orders);


    
    new_mindex=mindex;
    
    int imax;
    int maxOrd=maxVal(orders,&imax);
    
    for(int ipc=0;ipc<npc;ipc++){
        if (orders(ipc)==maxOrd){

        int nzind=0;
        while (mindex(ipc,nzind)==0 && nzind<ndim-1) {
            nzind++;
        }
            
            //cout << nzind << endl;
            Array1D<int> new_mindex_try;
            getRow(mindex,ipc,new_mindex_try);
        for(int j=0;j<=nzind;j++){
            new_mindex_try(j)++;


            if(is_admis(new_mindex_try,new_mindex)){
                

                paddMatRow(new_mindex,new_mindex_try);  
                
            }
            
        new_mindex_try(j)--;
        }
        
        }
    }
    
    return;
}

void getOrders(Array2D<int>& mindex,Array1D<int>& orders){
    
int npc=mindex.XSize();
int ndim=mindex.YSize();

    orders.Resize(npc,0);
    
    for (int ipc=0; ipc<npc ; ipc++) {
        int sum=0;
        for(int id=0;id<ndim;id++)
            sum+=mindex(ipc,id);
        orders(ipc)=sum;
    }
return;
}


int get_invmindex(Array1D<int> mi)
{
  int nd=mi.XSize();
  int ss=0;
      for(int id=0;id<nd;id++)
        ss+=mi(id);
      
      int index=computeNPCTerms(nd,ss-1);


      index+=get_invmindex_ord(mi);

      


  return index;
}

int get_invmindex_ord(Array1D<int> mi)
{
  int index=0;
  int nd=mi.XSize();

int ss=0;
      for(int id=0;id<nd;id++)
        ss+=mi(id);

if (nd>1){

      for(int ii=ss;ii>mi(0);ii--)
        index+=computeNPCTerms(nd-2,ss-ii);

             Array1D<int> mic(nd-1,0);
        for(int id=0;id<nd-1;id++)
          mic(id)=mi(id+1);
        
        index+=get_invmindex_ord(mic);
      }


  return index;

}
