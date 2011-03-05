#include <iostream>
#include <cmath>   
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
using namespace std;
using namespace arma;

// using namespace Rcpp;

class CLP {
private: 
  const int ngrid;
  const int nvertex;
  const int nparam;
  const int nmessages;

  mat Param;
  mat Rate;
  //mat Param_Poisson;
  const colvec Data_t;
  const colvec Data_k;
  const mat Data_v;

public:
  CLP(mat Data_in, int nParam_in, int nVertex_in);
  ~CLP();
  void setParam(mat Param_in);
  void setRate(mat Rate_in);
  mat rawSimulator(double myLHS, double myRHS, colvec Init);
  mat rejSimulator();
  colvec getEST(int maxMC);
};

colvec CLP::getEST(int maxMC){
  colvec retSta = zeros(nvertex);
  colvec retSca = zeros(1);
  colvec retOBJ = zeros(nvertex+1);

  colvec Ts;
  colvec dT;
  colvec row_i;
  colvec row_j;

  mat cur_sim;

  for(int mc_itr = 0; mc_itr < maxMC; ++mc_itr){
    cout << "\t MC#= "<< mc_itr << endl;
    cur_sim = rejSimulator();

    Ts = trans(cur_sim.row(nvertex));

    dT = zeros<colvec>(Ts.n_elem);
    for(int itr =1; itr < Ts.n_elem;++itr){
      dT(itr) = Ts(itr) - Ts(itr-1);
      cout << "HERE" << itr << endl;
    }

    
    retSta = retSta + mean(cur_sim.rows(0,nvertex-1),1);

    for(int i =0; i<nvertex-1;i++){	
      for(int j=(i+1);j <nvertex;j++) {
	row_i = trans(cur_sim.row(i));
	row_j = trans(cur_sim.row(j));
	cout <<"HERE"<<endl;        
	retSca = retSca + sum((row_i % row_j + (1-row_i) % (1-row_j)) % dT)/nmessages;
      }
    }
  }

  retSta = retSta/maxMC;
  retSca = retSca/maxMC;
  retOBJ = join_cols(retSta,retSca);

  return(retOBJ);
}

void CLP::setParam(mat Param_in)
{
  Param = Param_in;
}

void CLP::setRate(mat Rate_in)
{
  Rate = Rate_in;
}


CLP::~CLP(){
}


CLP::CLP(mat Data_in, int nParam_in, int nVertex_in): ngrid(10), nvertex(nVertex_in), nparam(nParam_in), nmessages(Data_in.n_rows), Data_t(Data_in.col(0)), Data_k(Data_in.col(nvertex+1)), Data_v(Data_in.cols(1,nvertex))
{
}

mat CLP::rawSimulator(double myLHS, double myRHS, colvec Init) 
{
  const double dT = (myRHS-myLHS)/ngrid;
  double drift = 0;
  double volat = 0; 
  
  rowvec States_T(ngrid+1);
  mat States_VX(nvertex,ngrid+1);
 
  States_T(0) = myLHS;
  
  States_VX(span(0,nvertex-1),0) = Init(span(0,nvertex-1),0);
  for(int itr_t=1;itr_t <= ngrid; ++itr_t) {
    States_T(itr_t) = States_T(itr_t-1) + dT; 
    for(int itr_v=0;itr_v < nvertex; ++itr_v){
      drift = Param(itr_v,1)*(States_VX(itr_v,itr_t-1) - Param(itr_v,0))*dT;
      volat = Param(itr_v,2)*sqrt(abs(States_VX(itr_v,itr_t)*(1-States_VX(itr_v,itr_t-1))*dT))*as_scalar(randn());
      States_VX(itr_v,itr_t) = States_VX(itr_v,itr_t-1)+ drift + volat; 

      if( States_VX(itr_v,itr_t) < 0){
	States_VX(itr_v,itr_t) = 0.001;
      }

      if( States_VX(itr_v,itr_t) > 1){
	States_VX(itr_v,itr_t) = 0.999;
      } 

    }
  }
  return(join_cols(States_VX,States_T));
}

mat CLP::rejSimulator(){
  mat retSTATE = Param.col(0); //columns will grow...
  mat retTIMES = zeros(1); //columns will grow...
  
  double myLHS_cur;
  double myRHS_cur;

  colvec myINIT_cur;
  rowvec myVERTEX_cur;
  int myTOPIC_cur;

  for(int itr_message=0; itr_message < nmessages-1; itr_message++){
    // cout << "Message # " << itr_message+1 << "/" << nmessages << endl;
    myLHS_cur = Data_t(itr_message);
    myRHS_cur = Data_t(itr_message+1);
    myVERTEX_cur = Data_v.row(itr_message+1);
    myTOPIC_cur = Data_k(itr_message+1);
    
    if(itr_message == 0){
      myINIT_cur = Param.col(0);
    } else {
      // otherwise updated at the end of the while-loop below
    }
   
    bool do_more = true;
    
    while(do_more){  
      mat retOBJ_prop = rawSimulator(myLHS_cur,myRHS_cur,myINIT_cur);
      mat State_VX_prop = retOBJ_prop.rows(0,nvertex-1);
      mat State_T_prop = retOBJ_prop.row(nvertex);

      int myRHS_index = retOBJ_prop.n_cols-1;
      double mystop_prob = 0;
      double weight = 0;
      
      for(int i=0;i < (nvertex-1); ++i){ //iterate through vertex
	for(int j = (i+1); j < nvertex; ++j){ //iterate through diff vertex
	  rowvec path_i = State_VX_prop.row(i);
	  rowvec path_j = State_VX_prop.row(j);

	  rowvec path_ij_rate_1 = path_i % path_j;
	  rowvec path_ij_rate_2 = (1 - path_i) % (1 - path_j);

	  mat path_ij_rate_both = join_cols(path_ij_rate_1, path_ij_rate_2);
	  rowvec path_ij_rate_sum = path_ij_rate_1 + path_ij_rate_2;

	  if(myVERTEX_cur(i) == 1 && myVERTEX_cur(j) == 1){
	    weight = path_ij_rate_both(as_scalar(myTOPIC_cur)-1,myRHS_index);
	  } else {
	    weight = 1;
	  }

	  mystop_prob = mystop_prob + log(weight);

	  mystop_prob = mystop_prob - Rate(0,0)*mean(path_ij_rate_sum)*(myRHS_cur-myLHS_cur);
	 
	} // j>i for-ends
      } // i for-ends
      //
      // now the decision time
      if(as_scalar(log(randu(1))) < mystop_prob) {
	do_more = false;
	myINIT_cur = State_VX_prop.col(myRHS_index);	
	retSTATE = join_rows(retSTATE,State_VX_prop.cols(1,myRHS_index));
	retTIMES = join_rows(retTIMES,State_T_prop.cols(1,myRHS_index));
      } else {
	do_more = true;
      }

    } // this finishes one interval ... if not rejected ...
    //    cout << "---- accept " << endl;

    // starting new interval
  }
  return(join_cols(retSTATE,retTIMES));
}

int main(int argc, char **argv) {
  // const int NVERTEX_Normal = 2;
  // const int NVERTEX_Abnorm = 2;
  int NVERTEX = atoi(argv[1]);
  int NPARAM = atoi(argv[2]);
  int maxSearch = atoi(argv[3]);
  int maxMC = atoi(argv[4]);

  mat myParam(NVERTEX,NPARAM);
  mat myRate(1,1);
  
  mat SIMIN;
  SIMIN.load("SIMIN.mat", raw_ascii);
  myParam.col(0) = trans(SIMIN(SIMIN.n_rows - 1, span(0,NVERTEX-1)));
  myParam.col(1) = ones(NVERTEX)*(-1);
  myParam.col(2) = ones(NVERTEX);
  myRate = SIMIN(SIMIN.n_rows-1,NPARAM);

  mat myData;
  myData.load("myData.txt",raw_ascii);

  cout << myData <<endl;
  cout << myParam << endl;
  cout << myRate << endl;

  CLP myCLP(myData, NPARAM, NVERTEX);
  myCLP.setParam(myParam);
  myCLP.setRate(myRate);
 
  mat SaveME(maxSearch,NVERTEX+1);
  mat myITRresult; 
  
  for(int i=0;i < maxSearch; ++i){
    cout << "Search #" << i << endl;

    myITRresult = myCLP.getEST(maxMC);    


    myParam.col(0) = myITRresult(span(0,NVERTEX-1),0);
    myCLP.setParam(myParam);
    myCLP.setRate(ones(1)*myITRresult(NVERTEX,0));
    
    SaveME.row(i) = reshape(myITRresult,1,NVERTEX+1);
    
    cout << SaveME.row(i) << endl;
  }
  
  SaveME.save("SIMOUT.mat",raw_ascii);
  
  return 0;
}
