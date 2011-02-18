#include <cmath>   
#include <armadillo>

using namespace arma;
using namespace Rcpp;

arma::mat Param = as<arma::mat>( Param_In ) ;
arma::mat Data = as<arma::mat>( Data_In ) ;
arma::colvec Init = as<arma::colvec>( Init_In ) ;

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

  Data_t = Data_in(span(1,nmessages),0);
  Data_s = Data_in(span(1,nmessages),span(1,nvertex));
  Data_k = Data_in(span(1,nmessages),nvertex+1));
}

mat CLP::rawSimulator(double myLHS, double myRHS, colvec Init) 
{
  double drift_term = 0;
  double volat_term = 0;
  
  mat States_VT(1+nvertex,ngrid+1);
  
  double dT = (myRHS-myLHS)/ngrid;

  States_VT(0,0) = myLHS;
  States_VT(span(1,nvertex),0) = Init;

  for(int itr_t=1;itr_t <= ngrid; itr_t++) {
    States_VT(0,itr_t) = States_VT(0,itr_t-1) + dT; 
    for(int itr_v=1;itr_v <= nvertex; itr_v++){
      drift_term = Param(itr_v,1)*(States_VT(itr_v,itr_t-1) - Param(itr_v,0))*dT;
      volat_term = Param(itr_v,2)*sqrt(States_VT(itr_v,itr_t)*(1-States_VT(itr_v,itr_t-1))*dT)*as_scalar(randn());
      States_VT(itr_v,itr_t) = draft_term + volat_term; 
    };
  };

  return(States_VT);
}

mat CLP::rejSimulator(){
  mat retOBJ = zeros(1,ngrid+1);

  mat sim_out; 
  bool do_more = true;
  colvec myInit;
  double myLHS;
  double myRHS;
  colvec myVERTEX;
  int myTOPIC;

  for(int itr_message=0; itr_message < nmessages; itr_messages++){
    myLHS = Data_t(itr_message);
    myRHS = Data_t(itr_message+1);
    myVERTEX = Data_s(itr_message+1,span(1,nvertex)); 
    myTOPIC = Data_s(itr_message+1,nvertex+1);

    if(itr_message == 0){
      myInit = Param(span(1,nvertex),0);
    } else {
      myInit = sim_out(span(1,nvertex),ngrid);      
    }

    while(do_more){  
      sim_out = rawSimulator(myLHS,myRHS,myInit); 

      double mycumsum = 0;
      
      for(int i=1;i<nvertex; i++){ //iterate through vertex
	for(int j = (i+1); j<= nvertex;j++){ //iterate through diff vertex
	  colvec curpath_i = sim_out(i,span(0,ngrid));
	  colvec curpath_j = sim_out(j,span(0,ngrid));

	  mat temppath_i = join(curpath_i,1-curpath_i);
	  mat temppath_j = join(curpath_j,1-curpath_j);
	  
	  mat temppath_ija = join(curpath_i, cur_path_j);
	  mat temppath_ijb = join(1 - curpath_i, 1 - curpath_j);
	  colvec templambda_ij_1 = prod(temppath_ija, 1);
	  colvec templambda_ij_2 = prod(temppath_ijb,1);

	  mat templambda_ij_both = join(templambda_ij_1, templambda_ij_2);
	  
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
	  
	  retOBJ = join(retOBJ,sim_out);
      }
    } 
  }
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
