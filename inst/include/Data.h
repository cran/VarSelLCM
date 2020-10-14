/*
  Cette classe contient les éléments relatifs aux données 

Ces élements sont:
  m_nrows : nombre d'observations
  m_ncols : nombre de variables
*/
#ifndef Data_H
#define Data_H

#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Data{
  public:
  int m_nrows, m_ncols;
  
  Data(){};
  ~Data(){};
  
};
#endif
