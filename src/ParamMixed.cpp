#include "ParamMixed.h"

ParamMixed::ParamMixed(){
  this->m_paramContinuous = ParamContinuous();
  this->m_paramInteger = ParamInteger();
  this->m_paramCategorical = ParamCategorical();
  this->m_pi = ones<vec>(0);  
}

ParamMixed::ParamMixed(const ParamMixed & param){
  this->m_paramContinuous = param.m_paramContinuous;
  this->m_paramInteger = param.m_paramInteger;
  this->m_paramCategorical = param.m_paramCategorical;
  this->m_pi = param.m_pi;  
}

ParamMixed::ParamMixed(const DataMixed * data,  const colvec & omega, const int & g){
  this->m_pi = ones<vec>(g)/g;  
  int vu=0;
  ivec who = randi<ivec>(data->m_nrows, distr_param(0, data->m_nrows -1));
  if (data->m_withContinuous){
    m_paramContinuous = ParamContinuous(data->m_continuousData_p, omega.subvec(vu, vu + data->m_continuousData_p->m_ncols - 1), g, who);
    vu += data->m_continuousData_p->m_ncols;
  }
  if (data->m_withInteger){
    m_paramInteger = ParamInteger(data->m_integerData_p, omega.subvec(vu, vu + data->m_integerData_p->m_ncols - 1), g, who);
    vu += data->m_integerData_p->m_ncols;
  }
  if (data->m_withCategorical){
    m_paramCategorical = ParamCategorical(data->m_categoricalData_p, omega.subvec(vu, vu + data->m_categoricalData_p->m_ncols - 1), g);
    vu += data->m_categoricalData_p->m_ncols;
  }
}

void ParamMixed::egalise(const DataMixed * data, const colvec omega){
  int vu=0;
  if (data->m_withContinuous){
    m_paramContinuous.egalise(omega.subvec(vu, vu + data->m_continuousData_p->m_ncols - 1));
    vu += data->m_continuousData_p->m_ncols;
  }
  if (data->m_withInteger){
    m_paramInteger.egalise(data->m_integerData_p, omega.subvec(vu, vu + data->m_integerData_p->m_ncols - 1));
    vu += data->m_integerData_p->m_ncols;
  }
  if (data->m_withCategorical){
    m_paramCategorical.egalise(omega.subvec(vu, vu + data->m_categoricalData_p->m_ncols - 1));
    vu += data->m_categoricalData_p->m_ncols;
  }
}
