//mysig = signature(param="matrix", init="numeric", LHS="numeric", RHS="numeric")
//myinc = '
#include <cmath>   

arma::colvec x = Rcpp::as<arma::colvec> (x_);
arma::mat Y = Rcpp::as<arma::mat>( Y_ ) ;
arma::colvec z = Rcpp::as<arma::colvec>( z_ ) ;

class Conditional_LatentPosition:public LatenPosition {
private:
  arma::mat Data;

public:
  void setData(arma::mat Data);
};

Conditional_LatentPosition::Conditional_LatentPosition(arma::mat Data, arma::mat Param_in, arma::colvec Init_in) : LatentPosition(arma::mat Param_in, arma::colvec Init_in)
{
  this.Data = Data;
};

void Conditional_LatentPosition::setData(arma::mat Data) {
  this.Data = Data;
};

arma::mat Conditional_LatentPosition::Simulate(arma::rowvec dT){
  arma::mat sim_out;
  int ngrid = dT.n_elem;
  bool do_more = true;
  double temp_val_1 = 0;
  double temp_val_2 = 0;
  while(do_more){  // need to know boolean type in c++
    sim_out = LatentPosition::Simulate(dT); //need to know scope resolution
    
    for(int i=1;i<nvertex;i++){
      mean(sim_out,dim=1)*arma:sum(dT);
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
    

    if(){
      = T;
    };
  };

};



class LatentPosition {
private: 
  long nvertex;
  int nparam;
  arma::mat Param;
  arma::colvec Init; 
public:
  LatentPosition(arma::mat Param_In, arma::colvec Init_In);
  virtual ~LatentPosition();
  virtual arma::mat Simulator(arma::rowvec dT);
  void setParam(arma::mat Param_in);
  void setInit(arma::mat Init_in);
};

void LatentPosition::setParam(arma::mat Param_in):
{
  this.Param = Param_in;
  this.nparam = Param.n_rows;
};

LatentPosition::LatentPosition(arma::mat Param_in, arma::colvec Init_in)
{
  Param = Param_in;
  Init = Init_in;
  
  nvertex = Param.n_rows;
  nparam = Param.n_cols;
}

arma::mat LatentPosition::Simulator(arma::rowvec dT) 
{
  double drift_term = 0;
  double volat_term = 0;
  
  int ngrid = dT.n_elem;
  
  arma::mat States_VT(nvertex,ngrid);
  arma::colvec States_CU(Init);
  
  for(int itr_t=0;itr_t < ngrid;itr_t++) {
    for(int itr_v=0;itr_v <nvertex;itr_v++){
      drift_term = Param(itr_v,1)*(States_CU(itr_v) - Param(itr_v,2))*dT(itr_t);
      volat_term = Param(itr_v,3)*sqrt(States_CU(itr_v)*(1-States_CU(itr_v))*dT(itr_t))*as_scalar(arma::randn());
      States_VT(itr_v,itr_t) = draft_term + volat_term; 
    }
    States_CU = States_VT(arma::span::all,itr_t);
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
  
  LatentPosition myLP(X,);
  arma::rowvec dT(NGRID);

  arma::mat output(NVERTEX,NGRID) = myLP.simulate(dT);

  for(int i=1;i < NVERTEX; i++){
    for(int j=(i+1); j < NVERTEX;j++){
      
    };
  };
  
  return 0;
}
