#include "AlgorithmContinuous.h"
#include "XEMContinuous.h"

AlgorithmContinuous::AlgorithmContinuous(const DataContinuous * data, const S4 * reference_p){
  vbleSelec = as<S4>(reference_p->slot("strategy")).slot("vbleSelec");
  if (vbleSelec){
    data_p = data;
    InitCommumParamAlgo(as<S4>(reference_p->slot("model")).slot("g"), as<S4>(reference_p->slot("strategy")).slot("initModel"), data_p->m_nrows, data_p->m_ncols);
    m_integralenondiscrim=ones<vec>(data_p->m_ncols);
    vec tmp;
    for (int j=0; j<data_p->m_ncols; j++){
      tmp = data_p->m_x.col(j);
      m_integralenondiscrim(j)=IntegreOneVariable(tmp(find(data_p->m_notNA.col(j) == 1)), j);
    }
  }
}

double AlgorithmContinuous::IntegreOneVariable(const vec & v, const int & j){
  double output = 0;
  double n = v.n_rows;
  if (n> 0){ 
    double n1 = n + data_p->m_priors(j,3);
    double s1 = sqrt( data_p->m_priors(j,1)*data_p->m_priors(j,1) + var(v) *(n-1)  + pow((data_p->m_priors(j,2) - mean(v)),2) /(1/data_p->m_priors(j,3) + 1/n)  );
    output =  -log(sqrt( M_PI))*n + lgamma((n + data_p->m_priors(j,0))*0.5) - lgamma(data_p->m_priors(j,0)*0.5) +   data_p->m_priors(j,0) * log(data_p->m_priors(j,1)/s1) - n*log(s1) + log(sqrt(data_p->m_priors(j,3) / n1) );
  }
  return output;
}

double AlgorithmContinuous::Integre_Complete_Like_Cand(){
  double outmicl = lgamma(m_g*0.5) - m_g*lgamma(0.5) - lgamma(data_p->m_nrows + m_g*0.5) + sum(m_integralenondiscrim);
  for (int k=0; k<m_g; k++)  outmicl += lgamma(sum(m_zCandCurrent==k) + 0.5);
  vec tmp;
  for (int j=0; j<data_p->m_ncols; j++){
    if (m_omegaCurrent(j)==1){
      tmp = data_p->m_x.col(j);
      for (int k=0; k<m_g; k++)  outmicl +=  IntegreOneVariable(tmp(find(((m_zCandCurrent == k) + data_p->m_notNA.col(j)) == 2)), j);        
      outmicl -=  m_integralenondiscrim(j);
    }
  }
  return outmicl;
}

void AlgorithmContinuous::Optimize_model(){
  m_miclCurrent = lgamma(m_g*0.5) - m_g*lgamma(0.5) - lgamma(data_p->m_nrows + m_g*0.5) + sum(m_integralenondiscrim);
  for (int k=0; k<m_g; k++) m_miclCurrent += lgamma(sum(m_zCandCurrent==k) + 0.5);
  double discrim=0;
  vec tmp;
  for (int j=0; j<data_p->m_ncols; j++){
    tmp = data_p->m_x.col(j);
    discrim =  -  m_integralenondiscrim(j) ;   
    for (int k=0; k<m_g; k++)  discrim +=  IntegreOneVariable(tmp(find(((m_zCandCurrent == k) + data_p->m_notNA.col(j)) == 2)), j);        
    if (discrim > 0){
      m_omegaCurrent(j) = 1;
      m_miclCurrent +=  discrim;   
    }else{
      m_omegaCurrent(j) = 0;
    }
  }
}

void AlgorithmContinuous::zCandInit(){
  XEMContinuous xem(data_p, m_omegaCurrent, m_g);
  xem.Run();
  m_zCandCurrent = xem.FindZMAP();
  m_zStarCurrent = m_zCandCurrent;
}
