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
  int nvertex;
  int nparam;
  int nmessages;
  int ngrid;

  mat Param;
  colvec Data_t;
  mat Data_v;
  colvec Data_k;

public:
  CLP(mat Data_in, int nParam_in, int nVertex_in);
  ~CLP();
  void setParam(mat Param_in);
  mat rawSimulator(double myLHS, double myRHS, colvec Init);
  mat rejSimulator();
  colvec getMeans(int maxMC);
};


colvec CLP::getMeans(int maxMC){
  colvec retVec= zeros(nvertex);
  //  colvec retVec= zeros<colvec>(nvertex);
  //  cout << "Test 7" << endl;
 mat cur_sim;

  for(int mc_itr = 0; mc_itr < maxMC; ++mc_itr){
    cur_sim = rejSimulator();
    retVec = retVec + mean(cur_sim(span(0,nvertex-1),span::all),1);
  }
  
  return(retVec);
}

void CLP::setParam(mat Param_in)
{
  Param = Param_in;
}

CLP::~CLP(){
}


CLP::CLP(mat Data_in, int nParam_in, int nVertex_in): ngrid(10)
{
  nvertex = nVertex_in;
  nparam = nParam_in;
  nmessages = Data_in.n_rows;

  Data_t = Data_in(span::all, 0);
  Data_v = Data_in(span::all, span(1,nvertex));
  Data_k = Data_in(span::all, nvertex);

}

mat CLP::rawSimulator(double myLHS, double myRHS, colvec Init) 
{
  const double dT = (myRHS-myLHS)/ngrid;
  double drift = 0;
  double volat = 0; 
  
  rowvec States_T(ngrid+1);
  mat States_VX(nvertex,ngrid+1);
 
  States_T(0) = myLHS;
  
  States_VX(span(0,nvertex-1),0) = Init;
  //cout << "init 1" << endl;
  for(int itr_t=1;itr_t <= ngrid; ++itr_t) {
    States_T(itr_t) = States_T(itr_t-1) + dT; 
    for(int itr_v=0;itr_v < nvertex; ++itr_v){
      drift = Param(itr_v,1)*(States_VX(itr_v,itr_t-1) - Param(itr_v,0))*dT;
      //  cout << "Hohoho" << endl;
      volat = Param(itr_v,2)*sqrt(abs(States_VX(itr_v,itr_t)*(1-States_VX(itr_v,itr_t-1))*dT))*as_scalar(randn());
      //cout << "Hohoho2" << endl;
      States_VX(itr_v,itr_t) = drift + volat; 
      //cout << "Hohoho3" << endl;
    }
  }
  //  cout << "Test 10" << endl;
  return(join_cols(States_VX,States_T));
}

mat CLP::rejSimulator(){
  mat retOBJ = Param(span(),0);

  double myLHS_cur;
  double myRHS_cur;

  colvec myINIT_cur;
  rowvec myVERTEX_cur;
  int myTOPIC_cur;

  for(int itr_message=0; itr_message < nmessages-1; itr_message++){
    //cout << "inside loop1" << endl;
    myLHS_cur = Data_t(itr_message);
    //cout << "inside loop2" << endl;
    myRHS_cur = Data_t(itr_message+1);
    //cout << "inside loop3" << endl;
    // myVERTEX_cur = Data_v.row(itr_message+1,span::all); 
    myVERTEX_cur = Data_v.row(itr_message+1);
    //cout << "inside loop4" << endl;
    myTOPIC_cur = Data_k(itr_message+1);
    //cout << "inside loop5" << endl;

    if(itr_message == 0){
      myINIT_cur = Param(span::all,0);
    } else {
      // otherwise updated at the end of the while-loop below
    }
   
    bool do_more = true;
    
    while(do_more){  
      cout << itr_message << "------ rej" << endl;  

      mat retOBJ_prop = rawSimulator(myLHS_cur,myRHS_cur,myINIT_cur);
      mat State_VX_prop = retOBJ_prop(span(0,nvertex-1),span::all);
      //      mat State_T_prop = retOBJ_prop(nvertex,span::all);

      int myRHS_index = retOBJ_prop.n_cols-1;
      double mystop_prob = 0;
      double weight = 0;
	
      
      for(int i=0;i < (nvertex-1); ++i){ //iterate through vertex
	for(int j = (i+1); j < nvertex; ++j){ //iterate through diff vertex
	  // colvec path_i = State_VX_prop.row(i,span::all);
	  rowvec path_i = State_VX_prop.row(i);
	  cout << path_i << endl;
	  // colvec path_j = State_VX_prop(j,span::all);
	  rowvec path_j = State_VX_prop.row(j);
	  

	  rowvec path_ij_rate_1 = path_i % path_j;
	  rowvec path_ij_rate_2 = (1 - path_i) % (1 - path_j);

	  //	  mat temppath_ij_1 = join_rows(path_i, path_j);
	  // mat temppath_ij_2 = join_rows(1 - path_i, 1 - path_j);
	  // colvec path_ij_rate_1 = prod(temppath_ij_1,1);
	  // colvec path_ij_rate_2 = prod(temppath_ij_2,1);
	  mat path_ij_rate_both = join_cols(path_ij_rate_1, path_ij_rate_2);
	  rowvec path_ij_rate_sum = path_ij_rate_1 + path_ij_rate_2;
	  //	  cout << "Test 11" << endl;
	  // colvec path_ij_rate_sum = sum(path_ij_rate_both,1);

	  //	  cout << Data_k << endl;
	  //	  cout << myVERTEX_cur << endl;
	  if(myVERTEX_cur(i) == 1 && myVERTEX_cur(j) == 1){
	    //cout << myRHS_index << "\t" << myTOPIC_cur;
	    weight = path_ij_rate_both(as_scalar(myTOPIC_cur)-1,myRHS_index);

	    if(weight <= 0){ 
	      weight = 0.01;
	    }
	      
	    // cout << "here ?" << endl;
	  } else {
	    weight = 1;
	  }
	  //  cout << "Before mystop" << endl;
	  cout << weight << endl;
	  cout << myRHS_cur << " " << myLHS_cur << endl;	
	  cout << mean(path_ij_rate_sum) << endl;
	  cout << path_ij_rate_sum << endl;
	  mystop_prob = mystop_prob + log(weight) - mean(path_ij_rate_sum)*(myRHS_cur-myLHS_cur);
	} // j>i for-ends
      } // i for-ends
      cout << mystop_prob << endl;
      // now the decision time
      if(as_scalar(log(randu(1))) < mystop_prob) {
	do_more = false;
	myINIT_cur = State_VX_prop(span::all,myRHS_index);	
	retOBJ_prop = State_VX_prop(span::all,span(1,myRHS_index));
	cout << retOBJ << endl;
	cout << retOBJ_prop << endl;
      } else {
	do_more = true;
      }

    } // this finishes one interval ... if not rejected ...
    cout << "---- accept " << endl;
    
    // starting new interval
  }

  return(retOBJ);
}

int main(int argc, char **argv) {
  // const int NVERTEX_Normal = 2;
  // const int NVERTEX_Abnorm = 2;
  int NVERTEX = atoi(argv[1]);
  int NPARAM = atoi(argv[2]);

  // mat normal_vertex_param(NVERTEX_Normal, NPARAM);
  // mat abnorm_vertex_param(NVERTEX_Abnorm, NPARAM);
  mat myParam(NVERTEX,NPARAM);
  
  mat tmpmat;
  tmpmat.load("SIMIN.mat", raw_ascii);
  myParam.col(0) = trans(tmpmat.row(tmpmat.n_rows - 1));
  myParam.col(1) = ones(NVERTEX)*(-1);
  myParam.col(2) = ones(NVERTEX);
  
  // normal_vertex_param << 0.5 << -1 << 1 << endr << 0.5 << -1 << 1 << endr;
  // abnorm_vertex_param << 0.5 << -1 << 1 << endr << 0.5 << -1 << 1 << endr;
  // cout << "Test 1" << endl;
  
  // myParam = join_cols(normal_vertex_param,abnorm_vertex_param);
  //  cout << "Test 2" << endl;
  mat myData;
  myData.load("myData.txt",raw_ascii);
  //  cout << "Test 4" << endl;
  CLP myCLP(myData, NPARAM, NVERTEX);
  //  cout << "Test 3" << endl;
  
  myCLP.setParam(myParam);
  //  cout << "Test 5" << endl;
  int maxSearch = 100;
  int maxMC = 100;

  mat SaveME(maxSearch,NVERTEX);
  for(int i=0;i < maxSearch; ++i){
    cout << "Search #" << i << endl;
    myParam.col(0) = myCLP.getMeans(maxMC);
    myCLP.setParam(myParam);
    SaveME.row(i) = reshape(myParam.col(0),1,NVERTEX);
    cout << reshape(myParam.col(0),1,NVERTEX) << endl;
  }
  SaveME.save("SIMOUT.mat",raw_ascii);

  return 0;
}
