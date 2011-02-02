//mysig = signature(param="matrix", init="numeric", LHS="numeric", RHS="numeric")
//myinc = '
#include <cmath> 

arma::colvec x = Rcpp::as<arma::colvec> (x_);
arma::mat Y = Rcpp::as<arma::mat>( Y_ ) ;
arma::colvec z = Rcpp::as<arma::colvec>( z_ ) ;


class LatentPosition {
private: 
  long nvertex;
  int nparam;
  arma::mat Param;
  arma::colvec Init;
  
public:
  LatentPosition(arma::mat Param_In, arma::colvec Init_In);
  arma::mat Simulator(arma::rowvec dT);
};

LatentPosition::Simulator(arma::mat Param_in, arma::colvec Init_in)
{
  Param = Param_in;
  Init = Init_in;
  
  nvertex = Param.n_rows;
  nparam = Param.n_cols;
}

arma::mat Simulator(arma::rowvec dT) 
{
  double drift_term = 0;
  double volat_term = 0;
  
  int ngrid = dT.n_elem;
  
  arma::mat States_VT(nvertex,ngrid);
  arma::colvec States_CU(Init);

  arma::colvec dW(1);
  
  for(int itr_t=0;itr_t < ngrid;itr_t++) {
    for(int itr_v=0;itr_v <nvertex;itr_v++){
      drift_term = Param(itr_v,1)*(States_CU(itr_v) - Param(itr_v,2))*dT(itr_t);
      volat_term = Param(itr_v,3)*sqrt(States_CU(itr_v)*(1-States_CU(itr_v))*dT(itr_t))*as_scalar(arma::randn());
      States_VT(itr_v,itr_t) = draft_term + volat_term; 
    }
    States_CU = States_VT(arma::span::all,itr_t);
  }
}

//'

