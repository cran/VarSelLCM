#ifndef ALGO_H
#define ALGO_H

#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Algo{
  public:
  Mat<double> x, priors;
  colvec z, omega;
  double micl;
  int g;
  
  Algo();
  Algo(const S4 & obj);
  ~Algo(){};
  
  double Integre_Complete_Like(vec z);
  void Optimize_model();
  void Optimize_partition();
  void Run();

};
#endif