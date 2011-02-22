#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
  mat A = randn<mat>(3,3);
  mat B = randn<mat>(4,5);

  cout << B << endl;  
  
  A = 1-B;

  cout << A << endl;

  B = join_rows(A,B);

  cout << B << endl;

  B.load("myData.txt");
  
  cout << B(span::all,0)<< endl;
  return 0;
  }
