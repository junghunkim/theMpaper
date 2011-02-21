#include <cmath>   
#include <armadillo>

using namespace arma;

// using namespace Rcpp;
// arma::mat Param = as<arma::mat>( Param_In ) ;
// arma::mat Data = as<arma::mat>( Data_In ) ;
// arma::colvec Init = as<arma::colvec>( Init_In ) ;

class CondtionalLatentPosition {
private: 
  int nvertex;
  int nparam;
  int nmessages;
  int ngrid;

  mat Param;
  colvec Data_t;
  mat Data_s;
  colvec Data_k;

public:
  ConditionalLatentPosition(mat Param_In, colvec Init_In);
  ~ConditionalLatentPosition();
  mat rawSimulator();
  mat rejSimulator();
};

typedef CLP ConditionalLatentPosition;

CLP::ConditionalLatentPosition(mat Data_in, mat Param_in): ngrid(10)
{
  Param = Param_in;
  
  nvertex = Param.n_rows;
  nparam = Param.n_cols;

  nmessages = Data_in.nrows;

  Data_t = Data_in(span::all, 0);
  Data_v = Data_in(span::all, span(1,nvertex));
  Data_k = Data_in(span::all, nvertex+1));
}

mat CLP::rawSimulator(double myLHS, double myRHS, colvec Init) 
{
  const double dT = (myRHS-myLHS)/ngrid;
  double drift = 0;
  double volat = 0; 
  
  rowvec States_T(ngrid+1);
  mat States_VX(nvertex,ngrid+1);
 
  States_T(0) = myLHS;
  States_VX(span::all,0) = Init;

  for(int itr_t=1;itr_t <= ngrid; ++itr_t) {
    States_T(itr_t) = States_T(itr_t-1) + dT; 
    for(int itr_v=0;itr_v < nvertex; ++itr_v){
      drift = Param(itr_v,1)*(States_VX(itr_v,itr_t-1) - Param(itr_v,0))*dT;
      volat = Param(itr_v,2)*sqrt(States_VX(itr_v,itr_t)*(1-States_VX(itr_v,itr_t-1))*dT)*as_scalar(randn());
      States_VX(itr_v,itr_t) = drift + volat; 
    }
  }

  return(join_cols(States_VX,States_T));
}

mat CLP::rejSimulator(){
  colvec retOBJ = Param(span::all,0);

  double myLHS_cur;
  double myRHS_cur;

  colvec myINIT_cur;
  colvec myVERTEX_cur;
  int myTOPIC_cur;

  for(int itr_message=0; itr_message < nmessages-1; itr_message++){
    myLHS_cur = Data_t(itr_message);
    myRHS_cur = Data_t(itr_message+1);
    myVERTEX_cur = Data_v(itr_message+1,span::all); 
    myTOPIC_cur = Data_k(itr_message+1);

    if(itr_message == 0){
      myINIT_cur = Param(span::all,0);
    } else {
      // otherwise updated at the end of the while-loop below
    }
   
    bool do_more = true;
    
    while(do_more){  
      mat retOBJ_prop = rawSimulator(myLHS_cur,myRHS_cur,myINIT_cur);
      mat State_VX_prop = retOBJ_prop(span(0,nvertex-1),span::all);
      mat State_T_prop = retOBJ_prop(nvertex,span::all);

      int myRHS_index = retOBJ_prop.ncol-1;
      double mystop_prob = 0;
      double weight = 0;
	
      
      for(int i=0;i < (nvertex-1); ++i){ //iterate through vertex
	for(int j = (i+1); j < nvertex; ++j){ //iterate through diff vertex
	  colvec path_i = State_VX_prop(i,span::all);
	  colvec path_j = STate_VX_prop(j,span::all);

	  mat path_i_pvec = join_rows(path_i,1-path_i);
	  mat path_j_pvec = join_rows(path_j,1-path_j);
	  
	  mat temppath_ij_1 = join_rows(curpath_i, curpath_j);
	  mat temppath_ij_2 = join_rows(1 - curpath_i, 1 - curpath_j);
	  colvec path_ij_rate_1 = prod(temppath_ij_1,1);
	  colvec path_ij_rate_2 = prod(temppath_ij_2,1);

	  mat path_ij_rate_both = join_rows(path_ij_rate_1, path_ij_rate_2);
	  colvec path_ij_rate_sum = sum(path_ij_rate_both,1);

	  if(myVERTEX_cur(i) == 1 && myVERTEX_cur(j) == 1){
	    weight = path_ij_both(myRHS_index, myTOPIC_cur);
	  } else {
	    weight = 1;
	  }
	  
	  mystop_prob = mystop_prob + log(weight) - mean(path_ij_rate_sum)*(myRHS_cur-myLHS_cur);
	} // j>i for-ends
      } // i for-ends
      
      // now the decision time
      if(log(runif(1)) < mystop_prob) {
	do_more = false;
	myINIT_cur = State_VX_prop(span::all,myRHS_index);	
      } else {
	do_more = true;
      }

    } // this finishes one interval ... if not rejected ...
    
    retOBJ_prop = retOBJ_prop(span::all,span(1,myRHS_index));
    retOBJ = join_rows(retOBJ,retOBJ_prop);
    // starting new interval
  }

  return(retOBJ);
}

int main() {
  const int NVERTEX = 4;
  const int NPARAM = 3;
  const int NGRID = 100;

  arma::mat normal_vertex_param(NVERTEX,NPARAM);
  arma::mat abnorm_vertex_param(NVERTEX,NPARAM);

  normal_vertex_param.row(1);
  normal_vertex_param.row(2);

  abnormal_vertex_param.row(1);
  abnormal_vertex_param.row(2);
  
  arma::mat X = arma::join_rows(normal_vertex_param,abnormal_vertex_param);
  
  CLP myLP(X,);
  arma::rowvec dT(NGRID);

  arma::mat output(NVERTEX,NGRID) = myLP.simulate(dT);

  for(int i=1;i < NVERTEX; i++){
    for(int j=(i+1); j < NVERTEX;j++){
      
    };
  };
  
  return 0;
}
