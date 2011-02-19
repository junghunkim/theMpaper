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
    };
  };

  return(join_rows(States_VX,States_T));
}

mat CLP::rejSimulator(){
  mat retOBJ = zeros(1,ngrid+1);
  mat retOBJ_prop; 
  bool do_more = true;

  double myLHS_cur;
  double myRHS_cur;

  colvec myINIT_cur;
  colvec myVERTEX_cur;
  int myTOPIC_cur;

  for(int itr_message=0; itr_message < nmessages-1; itr_message++){
    myLHS_cur = Data_t(itr_message);
    myRHS_cur = Data_t(itr_message+1);
    myVERTEX_cur = Data_s(itr_message+1,span(1,nvertex)); 
    myTOPIC_cur = Data_s(itr_message+1,nvertex+1);

    if(itr_message == 0){
      myINIT_cur = Param(span::all,0);
    } else {
      // otherwise updated at the end of the while-loop below
    }

    while(do_more){  
      retOBJ_prop = rawSimulator(myLHS_cur,myRHS_cur,myINIT_cur); 

      double mystop_prob = 0;
      
      for(int i=0;i < (nvertex-1); ++i){ //iterate through vertex
	for(int j = (i+1); j < nvertex; ++j){ //iterate through diff vertex
	  colvec curpath_i = retOBJ_prop(i,span::all);
	  colvec curpath_j = retOBJ_prop(j,span::all);

	  mat temppath_i = join_cols(curpath_i,1-curpath_i);
	  mat temppath_j = join_cols(curpath_j,1-curpath_j);
	  
	  mat temppath_ija = join_cols(curpath_i, cur_path_j);
	  mat temppath_ijb = join_cols(1 - curpath_i, 1 - curpath_j);

	  colvec templambda_ij_1 = prod(temppath_ija,1);
	  colvec templambda_ij_2 = prod(temppath_ijb,1);

	  mat templambda_ij_both = join_cols(templambda_ij_1, templambda_ij_2);
	  
	  colvec templambda_ij = sum(templambda_ij_both,1);
	  int nrow = templambda_ij.nelem;

	  double templambda_ij_end; 

	  if(myVERTEX(i) == 1 && myVERTEX(j) == 1){
	    templambda_ij_end = templambda_ij_both(nrow,myTOPIC);
	  } else {
	    templambda_ij_end = 1;
	  }
	  
	  mycumsum = mycumsum +log(templambda_ij_end) + log(exp(-mean(templambda_ij)*(myRHS-myLHS)));
	  
	  if(log(runif(1)) < mycumsum) {
	    do_more = false;
	  } else {
	    do_more = true;
	  }
	  
	  retOBJ = join_rows(retOBJ,sim_out);
	}
      } 
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
