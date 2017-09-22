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
  if (sum(omega)>0){
  uvec location = find(omega == 1);
  for (int j=0; j<m_mu.n_cols; j++){
    vec tmp = data->m_x.col(location(j));
    vec keep=tmp(find(data->m_notNA.col(j) == 1));
    m_mu.col(j)=m_mu.col(j)*mean(keep);
    m_sd.col(j)=m_sd.col(j)*sqrt(sum(pow(( keep - mean(keep)),2) ) / keep.n_rows);
  }
  this->m_pi = ones<vec>(g)/g;  
  ivec who = randi<ivec>(data->m_nrows, distr_param(0, data->m_nrows -1));
    for (int k=0; k<g; k++){
      for (int j=0; j<m_mu.n_cols; j++){
        if (data->m_notNA(who(k),location(j)) == 1)
          m_mu(k,j) = data->m_x(who(k),location(j));
      }
    }
  }
}


ParamContinuous::ParamContinuous(const DataContinuous * data,  const colvec & omega, const int & g, ivec who){
  this->m_mu = ones<mat>(g, sum(omega));
  this->m_sd = m_mu;
  if (sum(omega)>0){
  uvec location = find(omega == 1);
  for (int j=0; j<m_mu.n_cols; j++){
    vec tmp = data->m_x.col(location(j));
    vec keep=tmp(find(data->m_notNA.col(j) == 1));
    m_mu.col(j)=m_mu.col(j)*mean(keep);
    m_sd.col(j)=m_sd.col(j)*sqrt(sum(pow(( keep - mean(keep)),2) ) / keep.n_rows);
  }
  this->m_pi = ones<vec>(g)/g;  
    for (int k=0; k<g; k++){
      for (int j=0; j<m_mu.n_cols; j++){
        if (data->m_notNA(who(k),location(j)) == 1)
          m_mu(k,j) = data->m_x(who(k),location(j));
      }
    }
  }
}


void ParamContinuous::egalise(const colvec omega){
  for (int j=0; j<m_mu.n_cols; j++){
    if (omega(j)==0){
      for (int k=1; k<m_mu.n_rows;k++){
        m_mu(k,j) = m_mu(0,j);
        m_sd(k,j) = m_sd(0,j);
      }
    }
  }
}
