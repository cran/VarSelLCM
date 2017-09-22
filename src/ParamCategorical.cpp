#include "ParamCategorical.h"

ParamCategorical::ParamCategorical(){
  this->m_alpha.resize(0);
  this->m_pi = ones<vec>(0);  
}

ParamCategorical::ParamCategorical(const ParamCategorical & param){
  this->m_alpha = param.m_alpha;
  this->m_pi = param.m_pi;  
}

ParamCategorical::ParamCategorical(const DataCategorical * data,  const colvec & omega, const int & g){
  this->m_alpha.resize(sum(omega));
  this->m_pi = ones<vec>(g)/g;  
  if (sum(omega)>0){
    uvec location = find(omega == 1);
    for (int j=0; j<sum(omega); j++){
      m_alpha[j] = randu(g, data->m_nmodalities(location(j)));
      for (int k=0; k<g; k++) m_alpha[j].row(k) = m_alpha[j].row(k) / sum(m_alpha[j].row(k));
    }
  }
}

void ParamCategorical::egalise(const colvec omega){
  for (int j=0; j<omega.n_rows; j++){
    if (omega(j)==0){
      for (int k=1; k<m_alpha[j].n_rows;k++){
        m_alpha[j].row(k)=m_alpha[j].row(0);
      }
    }
  }
}
