#include "AlgorithmInteger.h"
#include "XEMInteger.h"

AlgorithmInteger::AlgorithmInteger(const DataInteger * data, const S4 * reference_p){
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

double AlgorithmInteger::IntegreOneVariable(const vec & v, const int & j){
  double output = 0;
  if (v.n_rows> 0){ 
    double alpha = sum(v) + data_p->m_priors(j,0);
    double beta = v.n_rows + data_p->m_priors(j,1);
    output = data_p->m_priors(j,0) * log(data_p->m_priors(j,1)) - lgamma(data_p->m_priors(j,0)) + lgamma(alpha) - alpha * log(beta);
    for (int i=0; i<v.n_rows; i++)
      output -= lgamma( v(i) + 1 );
  }
  return output;
}

double AlgorithmInteger::Integre_Complete_Like_Cand(){
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

void AlgorithmInteger::Optimize_model(){
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

void AlgorithmInteger::zCandInit(){
  XEMInteger xem(data_p, m_omegaCurrent, m_g);
  xem.Run();
  m_zCandCurrent = xem.FindZMAP();
  m_zStarCurrent = m_zCandCurrent;
}
