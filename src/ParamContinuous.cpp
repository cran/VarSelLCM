#include "ParamContinuous.h"

ParamContinuous::ParamContinuous(){
  this->m_mu = ones<mat>(0,0);
  this->m_sd = m_mu;
  this->m_pi = ones<vec>(0);  
}

ParamContinuous::ParamContinuous(const ParamContinuous & param){
  this->m_mu = param.m_mu;
  this->m_sd = param.m_sd;
  this->m_pi = param.m_pi;  
}

ParamContinuous::ParamContinuous(const DataContinuous * data,  const colvec & omega, const int & g){
  this->m_mu = ones<mat>(g, sum(omega));
  this->m_sd = m_mu;
  this->m_pi = ones<vec>(g)/g;  
  int k=0, li=0;
  if (sum(omega)>0){
    uvec location = find(omega == 1);
    for (int j=0; j<m_mu.n_cols; j++){
      ivec who = randi<ivec>(data->m_nrows, distr_param(0, data->m_nrows -1));
      k=0;
      li=0;
      while (k<g){
        if (data->m_notNA(who(li),location(j)) == 1){
          m_mu(k,j) = data->m_x(who(li),location(j));
          k++;
        }
        li++;
      }
      vec tmp = data->m_x.col(location(j));
      vec keep=tmp(find(data->m_notNA.col(j) == 1));
      m_sd.col(j)=m_sd.col(j)*sqrt(sum(pow(( keep - mean(keep)),2) ) / keep.n_rows);
    }
  }
}
