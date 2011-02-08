//mysig = signature(param="matrix", init="numeric", LHS="numeric", RHS="numeric")
//myinc = '
#include <cmath>  

arma::mat myParam = Rcpp::as<arma::colvec> (Param_);
arma::mat myData = Rcpp::as<arma::mat>(Data_);
arma::rowvec mydT = Rcpp::as<arma::rowvec>(dT_);
arma::rowvec myInit = myParam.col(1); // 1 is for the mean values
double maxMC = 100;

const int NVERTEX = myParam.n_rows;
const int NPARAM = myParam.n_cols;
const int NGRID = 100;

Conditional_LatentPosition myLP(myData,myParam,myInit);

Rcpp::List myMCresult = Rcpp::List(maxMC);

for(int mymc_itr=0; mymc_itr < maxMC; mymc_itr++){
  myMCresult[mymc_itr] = myLP.simulate();
 };

return myMCresult;

/*       */
class Conditional_LatentPosition:public LatenPosition {
private:
  arma::mat Data; // nrow = ntotal_emails & ncol= (v1,...,vN,time), and for each row, the last entry is 
  
public:
  Conditional_LatentPosition(arma::mat Data, arma::mat Param_in, arma::colvec Init_in);
  ~Conditional_LatentPosition();
};

Conditional_LatentPosition::Conditional_LatentPosition(long mcrep,arma::mat Data, arma::mat Param_in, arma::colvec Init_in) : LatentPosition(arma::mat Param_in, arma::colvec Init_in)
{
  this.Data = Data;//should this be instantiated?
};

arma::mat Conditional_LatentPosition::Simulate(arma::rowvec dT){//for the conditional distribution, lhs,rhs can be varying
  arma::mat sim_out; //nrow == nvertex  && ncol = original ngrid?
  int ngrid = dT.n_elem;
  bool do_more = true;
  double temp_val_1 = 0;
  double temp_val_2 = 0;
  
  while(do_more){  // need to know boolean type in c++
    sim_out = LatentPosition::Simulate(dT); //need to know scope resolution
    
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
    

    if(){
      = T;
    };
  };

};



class LatentPosition {
private: 
  long nvertex;
  int nparam;

  double lhs;
  double rhs;

  arma::mat Param;
  arma::colvec Init; 
public:
  LatentPosition(arma::mat Param_In, arma::colvec Init_In);
  virtual ~LatentPosition();
  virtual arma::mat Simulator(arma::rowvec dT);
  void setParam(arma::mat Param_in);
  void setInit(arma::mat Init_in);
  void setLHSRHS(double lhs, double rhs);
};

void LatentPosition::setLHSRHS(double lhs, double rhs)
{
  this.lhs = lhs;
  this.rhs = rhs;
};

void LatentPosition::setParam(arma::mat Param_in)
{
  this.Param = Param_in;
  this.nparam = Param.n_cols;
};

void LatentPosition::setInit(arma::mat Init_in)
{
  this.Init = Init_in;
  this.nvertex = Init_in.n_rows;
}

LatentPosition::LatentPosition(arma::mat Param_in, arma::colvec Init_in, double lhs_in, double rhs_n)
{
  lhs = lhs_in;
  rhs = rhs_in;

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
