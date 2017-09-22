#include "ParamInteger.h"

ParamInteger::ParamInteger(){
  this->m_lambda = ones<mat>(0,0);
  this->m_pi = ones<vec>(0);  
}

ParamInteger::ParamInteger(const ParamInteger & param){
  this->m_lambda = param.m_lambda;
  this->m_pi = param.m_pi;  
}

ParamInteger::ParamInteger(const DataInteger * data,  const colvec & omega, const int & g){
  this->m_lambda = ones<mat>(g, sum(omega));
  this->m_pi = ones<vec>(g)/g;  
  if (sum(omega)>0){
    uvec location = find(omega == 1);
    for (int j=0; j<m_lambda.n_cols; j++){
      vec tmp = data->m_x.col(location(j));
      vec keep=tmp(find(data->m_notNA.col(j) == 1));
      m_lambda.col(j)=m_lambda.col(j)*mean(keep);
    }
    ivec who = randi<ivec>(data->m_nrows, distr_param(0, data->m_nrows -1));
    for (int k=0; k<g; k++){
      for (int j=0; j<m_lambda.n_cols; j++){
        if (data->m_notNA(who(k),location(j)) == 1)
        m_lambda(k,j) = data->m_x(who(k),location(j))+0.1;
      }
    }
  }
}


ParamInteger::ParamInteger(const DataInteger * data,  const colvec & omega, const int & g, ivec who){
  this->m_lambda = ones<mat>(g, sum(omega));
  this->m_pi = ones<vec>(g)/g;  
  if (sum(omega)>0){
    uvec location = find(omega == 1);
    for (int j=0; j<m_lambda.n_cols; j++){
      vec tmp = data->m_x.col(location(j));
      vec keep=tmp(find(data->m_notNA.col(j) == 1));
      m_lambda.col(j)=m_lambda.col(j)*(0.01+mean(keep));
    }
    for (int k=0; k<g; k++){
      for (int j=0; j<m_lambda.n_cols; j++){
        if (data->m_notNA(who(k),location(j)) == 1)
        m_lambda(k,j) = data->m_x(who(k),location(j))+0.01;
      }
    }
  }
}

void ParamInteger::egalise(const DataInteger * data,const colvec omega){
  for (int j=0; j<m_lambda.n_cols; j++){
    if (omega(j)==0){
      vec tmp = data->m_x.col(j);
      vec keep=tmp(find(data->m_notNA.col(j) == 1));
      m_lambda.col(j)=ones<vec>(m_lambda.n_rows)*mean(keep);
    }
  }
}
