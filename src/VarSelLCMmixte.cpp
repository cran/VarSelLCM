#include "AlgorithmContinuous.h"
#include "XEMContinuous.h"
#include "XEMContinuous.h"

//[[Rcpp::export]]
S4  OptimizeMICL(S4 reference, StringVector name){
  S4 * reference_p=&reference;
  string namestr = as<std::string>(name);
  DataContinuous * data_p = new DataContinuous(as<S4>(reference.slot("data")));
  AlgorithmContinuous *algo_p = new AlgorithmContinuous(data_p, reference_p);
  algo_p->Run(reference_p);
  XEMContinuous *xem_p  = new XEMContinuous(data_p, reference_p);
  xem_p->Run(); 
  xem_p->Output(reference_p);
  
  return reference;
}
