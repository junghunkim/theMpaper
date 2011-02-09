#include <cmath>   

arma::mat Param = Rcpp::as<arma::mat>( Param_In ) ;
arma::mat Data = Rcpp::as<arma::mat>( Data_In ) ;
arma::colvec Init = Rcpp::as<arma::colvec>( Init_In ) ;

class CondtionalLatentPosition {
private: 
  long nvertex;
  int nparam;

  arma::mat Param;
  arma::mat Data_t;
  arma::mat Data_s;

public:
  Conditional_LatentPosition(arma::mat Param_In, arma::colvec Init_In);
  ~Conditional_LatentPosition();
  arma::mat rawSimulator();
  arma::mat rejSimulator();
};

LatentPosition::LatentPosition(arma::mat Data_in, arma::mat Param_in)
{
  Param = Param_in;
  
  nvertex = Param.n_rows;
  nparam = Param.n_cols;
  
  nmessages = Data_in.nrows;

  Data_t = Data_in(span(1,nmessages),0);
  Data_s = Data_in(span(1,nmessages),span(1,nvertex));
}

arma::mat LatentPosition::rawSimulator(double myLHS, double myRHS, arma::mat Init, int ngrid=10) 
{
  double drift_term = 0;
  double volat_term = 0;
  
  arma::mat States_VT(1+nvertex,ngrid+1);
  
  const double dT = (myRHS-myLHS)/ngrid;

  States_VT(0,0) = myLHS;
  States_VT(span(1,nvertex),0) = Init;

  for(int itr_t=1;itr_t <= ngrid; itr_t++) {
    States_VT(0,itr_t) = States_VT(0,itr_t-1) + dT; 
    for(int itr_v=1;itr_v <= nvertex; itr_v++){
      drift_term = Param(itr_v,1)*(States_VT(itr_v,itr_t-1) - Param(itr_v,0))*dT;
      volat_term = Param(itr_v,2)*sqrt(States_VT(itr_v,itr_t)*(1-States_VT(itr_v,itr_t-1))*dT)*as_scalar(arma::randn());
      States_VT(itr_v,itr_t) = draft_term + volat_term; 
    };
  };

  return(States_VT);
}

arma::mat Conditional_LatentPosition::rejSimulator(){
  arma::mat sim_out; 
  bool do_more = true;

  for(int itr.message=0;itr.message<nmessages;itr.messages++){
    double myLHS =;
    double myRHS =;
    arma::mat myInit=;
    while(do_more){  
      
      sim_out = Simulate(myLHS,myRHS,myInit); 
      
      for(int i=1;i<nvertex;i++){
	mean(sim_out,dim=1)*arma::sum(dT);
      };
      
      for(int i=1;i<nvertex; i++){ //iterate through vertex
	for(int j = (i+1); j <=nvertex;j++){ //iterate through diff vertex
	  int temp_int = sim_out[i-1,j-1];
	  if(temp_int == 1){ // 1 means topic one 
	    temp_val_2 *= sim_out[i-1,t-1]*sim_out[j-1,t-1];
	  } else if(temp_int == 2) {
	    temp_val_2 *= (1-sim_out[i-1,t-1])*(1-sim_out[j-1,t-1]);
	  } else {
	    temp_val_2 *= 1;
	  };	// then, need to know if two communicated and compute the weight here
	}; 
      };
    };  
    
    if(){
      = T;
    };
  };
  
};

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
  
  LatentPosition myLP(X,);
  arma::rowvec dT(NGRID);

  arma::mat output(NVERTEX,NGRID) = myLP.simulate(dT);

  for(int i=1;i < NVERTEX; i++){
    for(int j=(i+1); j < NVERTEX;j++){
      
    };
  };
  
  return 0;
}
