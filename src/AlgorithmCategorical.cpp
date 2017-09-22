#include "AlgorithmCategorical.h"
#include "XEMCategorical.h"

AlgorithmCategorical::AlgorithmCategorical(const DataCategorical * data,  const S4 * reference_p){
  vbleSelec = as<S4>(reference_p->slot("strategy")).slot("vbleSelec");
  if (vbleSelec){
    data_p = data;
    InitCommumParamAlgo(as<S4>(reference_p->slot("model")).slot("g"), as<S4>(reference_p->slot("strategy")).slot("initModel"), data_p->m_nprofiles, data_p->m_ncols);
    m_integralenondiscrim=ones<vec>(data_p->m_ncols);
    for (int j=0; j<data_p->m_ncols; j++) m_integralenondiscrim(j)=IntegreOneVariableCategoricalNotDiscrim(j);
  }    
}

double AlgorithmCategorical::IntegreOneVariableCategoricalDiscrim(const int & j){
  int mj = data_p->m_whotakewhat[j].size();
  Mat <double> nkjh = 0.5 * ones<mat>(m_g, mj);
  for (int h=0; h<mj;h++){
    for (int i=0; i<data_p->m_whotakewhat[j][h].n_rows; i++) nkjh(m_zCandCurrent(data_p->m_whotakewhat[j][h](i)), h)+=data_p->m_w(data_p->m_whotakewhat[j][h](i));
  }
  double output = m_g*lgamma(mj * 0.5) - (m_g*mj) * lgamma(0.5);
  double tmpmargin = 0;
  for (int k=0; k<m_g; k++){
    tmpmargin = 0;
    for (int h=0; h<mj; h++){
      output += lgamma(nkjh(k,h));
      tmpmargin += nkjh(k,h);
    }
    output -= lgamma(tmpmargin);
  }
  return output;
}

double AlgorithmCategorical::IntegreOneVariableCategoricalDiscrim(const int & j, const colvec & z){
  int mj = data_p->m_whotakewhat[j].size();
  Mat <double> nkjh = 0.5 * ones<mat>(m_g, mj);
  for (int h=0; h<mj;h++){
    for (int i=0; i<data_p->m_whotakewhat[j][h].n_rows; i++) nkjh(z(data_p->m_whotakewhat[j][h](i)), h)+=data_p->m_w(data_p->m_whotakewhat[j][h](i));
  }
  double output = m_g*lgamma(mj * 0.5) - (m_g*mj) * lgamma(0.5);
  double tmpmargin = 0;
  for (int k=0; k<m_g; k++){
    tmpmargin = 0;
    for (int h=0; h<mj; h++){
      output += lgamma(nkjh(k,h));
      tmpmargin += nkjh(k,h);
    }
    output -= lgamma(tmpmargin);
  }
  return output;
}


double AlgorithmCategorical::IntegreOneVariableCategoricalNotDiscrim(const int & j){
  int mj = data_p->m_whotakewhat[j].size();
  Col <double> nh = 0.5 * ones<vec>(mj);
  for (int h=0; h<mj;h++)  nh(h) += sum(data_p->m_w(data_p->m_whotakewhat[j][h]));
  double output = lgamma(mj * 0.5) - (mj) * lgamma(0.5) -lgamma(sum(nh));
  for (int h=0; h<mj; h++) output += lgamma(nh(h));
  return output;
}

double AlgorithmCategorical::Integre_Complete_Like_Cand(){
  double outmicl = lgamma(m_g*0.5) - m_g*lgamma(0.5) - lgamma(sum(data_p->m_w) + m_g*0.5) + sum(m_integralenondiscrim);
  for (int k=0; k<m_g; k++) outmicl += lgamma(  sum(data_p->m_w(find(m_zCandCurrent==k))) + 0.5);
  for (int j=0; j<data_p->m_ncols; j++){
    if (m_omegaCurrent(j)==1) 
    outmicl +=  IntegreOneVariableCategoricalDiscrim(j) - m_integralenondiscrim(j);        
  }
  return outmicl;
}

void AlgorithmCategorical::Optimize_model(){
  m_miclCurrent = lgamma(m_g*0.5) - m_g*lgamma(0.5) - lgamma(sum(data_p->m_w) + m_g*0.5) + sum(m_integralenondiscrim);
  for (int k=0; k<m_g; k++) m_miclCurrent += lgamma(  sum(data_p->m_w(find(m_zCandCurrent==k))) + 0.5);
  double discrim=0;
  for (int j=0; j<data_p->m_ncols; j++){
    discrim =  IntegreOneVariableCategoricalDiscrim(j) - m_integralenondiscrim(j);   
    if (discrim > 0){
      m_omegaCurrent(j) = 1;
      m_miclCurrent +=  discrim;  
    }else{
      m_omegaCurrent(j) = 0;
    }      
  }
}

void AlgorithmCategorical::zCandInit(){
  XEMCategorical xem(data_p, m_omegaCurrent, m_g);
  xem.Run();
  m_zCandCurrent = xem.FindZMAP();
  m_zStarCurrent = m_zCandCurrent;
}
