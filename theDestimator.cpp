#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

using namespace std;
using namespace arma;

class SingleEstep {
private:
  const int nvertex;
  
  mat info_param;
  mat rate_param;

public:
  SingleEstep(mat Data_in, int nVertex_n);
  ~SingleEstep();
  mat rawSimulator(double myLHS, double myRHS, colvec Init);
};


class singleMstep {
private:

public:

};


int main(int argc, char **argv){
  int maxEMsteps = atoi(argv[1])

  for(int i =0; i < maxEMsteps; ++i) {
    
  }
}

