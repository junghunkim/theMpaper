#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
  mat B = randn<mat>(4,5);
  
  cout << B << endl;  
  B.rows(1,2) = ones(2,5);
  cout << B << endl;

  B.load("myData.txt");
  
  cout << B(span::all,0)<< endl;
  return 0;
  }
