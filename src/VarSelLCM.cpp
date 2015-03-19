
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "Algo.h"
using namespace std;
using namespace Rcpp;
using namespace arma;
 
//[[Rcpp::export]]
S4  OptimizeMICL(S4 obj, int optimize){
  
 Algo algorithm(obj); 
 
 if (optimize==1){
   algorithm.Run();
 }


   
  S4 model = obj.slot("model");
  model.slot("omega") = wrap(trans(algorithm.omega));
  S4 partitions = obj.slot("partitions");
  partitions.slot("zOPT") = wrap(trans(algorithm.z));
 S4 criteria = obj.slot("criteria");
 criteria.slot("MICL") = algorithm.micl;

 return obj;
  
}


 
