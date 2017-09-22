#include "AlgorithmMixed.h"
#include "XEMContinuous.h"
#include "XEMCategorical.h"
#include "XEMInteger.h"
#include "XEMMixed.h"
#include "XEMPen.h" 

//[[Rcpp::export]]
S4  OptimizeMICL(S4 reference, StringVector name){
  S4 * reference_p=&reference;
  string namestr = as<std::string>(name); 
  if (namestr == "Continuous"){
    DataContinuous * data_p = new DataContinuous(as<S4>(reference.slot("data")));
    AlgorithmContinuous *algo_p = new AlgorithmContinuous(data_p, reference_p);
    algo_p->Run(reference_p);
    XEMContinuous *xem_p  = new XEMContinuous(data_p, reference_p);
    xem_p->Run(); 
    xem_p->Output(reference_p);
  }else if (namestr == "Integer"){
    DataInteger * data_p = new DataInteger(as<S4>(reference.slot("data")));
    AlgorithmInteger *algo_p = new AlgorithmInteger(data_p, reference_p);
    algo_p->Run(reference_p);
    XEMInteger *xem_p  = new XEMInteger(data_p, reference_p);
    xem_p->Run(); 
    xem_p->Output(reference_p);
 }else if (namestr == "Categorical"){
    DataCategorical * data_p = new DataCategorical(as<S4>(reference.slot("data")));
    AlgorithmCategorical *algo_p = new AlgorithmCategorical(data_p, reference_p);
    algo_p->Run(reference_p);
    XEMCategorical *xem_p  = new XEMCategorical(data_p, reference_p);
    xem_p->Run(); 
    xem_p->Output(reference_p);
  }else if (namestr == "Mixed"){
    DataMixed *data_p = new DataMixed(as<S4>(reference.slot("data")));
    AlgorithmMixed *algo_p = new AlgorithmMixed(data_p, reference_p);
    algo_p->Run(reference_p);
    XEMMixed *xem_p  = new XEMMixed(data_p, reference_p);
    xem_p->Run(); 
    xem_p->Output(reference_p);
  }
  return reference;
}

//[[Rcpp::export]]
S4  ComputeMICL(S4 reference, StringVector name){
  S4 * reference_p=&reference;
  string namestr = as<std::string>(name);
  if (namestr == "Continuous"){
    DataContinuous * data_p = new DataContinuous(as<S4>(reference.slot("data")));
    AlgorithmContinuous *algo_p = new AlgorithmContinuous(data_p, reference_p);
    algo_p->ComputeMICL(reference_p);
  }else if (namestr == "Integer"){
    DataInteger * data_p = new DataInteger(as<S4>(reference.slot("data")));
    AlgorithmInteger *algo_p = new AlgorithmInteger(data_p, reference_p);
    algo_p->ComputeMICL(reference_p);
  }else if (namestr == "Categorical"){
    DataCategorical * data_p = new DataCategorical(as<S4>(reference.slot("data")));
    AlgorithmCategorical *algo_p = new AlgorithmCategorical(data_p, reference_p);
    algo_p->ComputeMICL(reference_p);
  }else if (namestr == "Mixed"){
    DataMixed *data_p = new DataMixed(as<S4>(reference.slot("data")));
    AlgorithmMixed *algo_p = new AlgorithmMixed(data_p, reference_p);
    algo_p->ComputeMICL(reference_p);
  }
  return reference;
}

//[[Rcpp::export]]
S4  OptimizePenLike(S4 reference, double pen){
  S4 * reference_p=&reference;
  XEMPen * xem_p  = new XEMPen(reference_p, pen);
  xem_p->Run(); 
  xem_p->Output(reference_p);
  return reference;
}
