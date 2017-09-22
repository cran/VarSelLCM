#include "AlgorithmMixed.h"
#include "XEMMixed.h"

AlgorithmMixed::AlgorithmMixed(const DataMixed * data, const S4 * reference_p){
  vbleSelec = as<S4>(reference_p->slot("strategy")).slot("vbleSelec");
  if (vbleSelec){
    data_p = data;
    InitCommumParamAlgo(as<S4>(reference_p->slot("model")).slot("g"), as<S4>(reference_p->slot("strategy")).slot("initModel"), data_p->m_nrows, data_p->m_ncols);
    if (data_p->m_withContinuous) algoCont_p = new AlgorithmContinuous(data_p->m_continuousData_p, reference_p);
    if (data_p->m_withInteger) algoInte_p = new AlgorithmInteger(data_p->m_integerData_p, reference_p);
    if (data_p->m_withCategorical) algoCate_p = new AlgorithmCategorical(data_p->m_categoricalData_p, reference_p);
  }
}

double AlgorithmMixed::Integre_Complete_Like_Cand(){
  double outmicl = lgamma(m_g*0.5) - m_g*lgamma(0.5) - lgamma(data_p->m_nrows + m_g*0.5);
  for (int k=0; k<m_g; k++)  outmicl += lgamma(sum(m_zCandCurrent==k) + 0.5);
  int loc=0;
  if (data_p->m_withContinuous){
    outmicl += sum(algoCont_p->m_integralenondiscrim);
    vec tmp;
    for (int j=0; j<data_p->m_continuousData_p->m_ncols; j++){
      if (m_omegaCurrent(loc)==1){
        tmp = data_p->m_continuousData_p->m_x.col(j);
        for (int k=0; k<m_g; k++)  outmicl +=  algoCont_p->IntegreOneVariable(tmp(find(((m_zCandCurrent == k) + data_p->m_continuousData_p->m_notNA.col(j)) == 2)), j);        
        outmicl -=  algoCont_p->m_integralenondiscrim(j);
      }
      loc++;
    }    
  }
  if (data_p->m_withInteger){
    outmicl += sum(algoInte_p->m_integralenondiscrim);
    vec tmp;
    for (int j=0; j<data_p->m_integerData_p->m_ncols; j++){
      if (m_omegaCurrent(loc)==1){
        tmp = data_p->m_integerData_p->m_x.col(j);
        for (int k=0; k<m_g; k++)  outmicl +=  algoInte_p->IntegreOneVariable(tmp(find(((m_zCandCurrent == k) + data_p->m_integerData_p->m_notNA.col(j)) == 2)), j);      
        outmicl -=  algoInte_p->m_integralenondiscrim(j);
      }
      loc++;
    }    
  }
  if (data_p->m_withCategorical){
    outmicl += sum(algoCate_p->m_integralenondiscrim);
    for (int j=0; j<data_p->m_categoricalData_p->m_ncols; j++){
      if (m_omegaCurrent(loc)==1) outmicl +=  algoCate_p->IntegreOneVariableCategoricalDiscrim(j, m_zCandCurrent) - algoCate_p->m_integralenondiscrim(j);        
      loc++;
    }  
  }  
  return outmicl;
}

void AlgorithmMixed::Optimize_model(){
  m_miclCurrent = lgamma(m_g*0.5) - m_g*lgamma(0.5) - lgamma(data_p->m_nrows + m_g*0.5);
  for (int k=0; k<m_g; k++) m_miclCurrent += lgamma(sum(m_zCandCurrent==k) + 0.5);
  double discrim=0;
  int loc=0;
  if (data_p->m_withContinuous){
    vec tmp;
    for (int j=0; j<data_p->m_continuousData_p->m_ncols; j++){
      tmp = data_p->m_continuousData_p->m_x.col(j);
      discrim = 0;
      for (int k=0; k<m_g; k++)  discrim +=  algoCont_p->IntegreOneVariable(tmp(find(((m_zCandCurrent == k) + data_p->m_continuousData_p->m_notNA.col(j)) == 2)), j);        
      if (discrim > algoCont_p->m_integralenondiscrim(j)){
        m_omegaCurrent(loc) = 1;
        m_miclCurrent +=  discrim;   
      }else{
        m_omegaCurrent(loc) = 0;
        m_miclCurrent +=  algoCont_p->m_integralenondiscrim(j);
      }
      loc++;
    }
  }
  if (data_p->m_withInteger){
    vec tmp;
    for (int j=0; j<data_p->m_integerData_p->m_ncols; j++){
      tmp = data_p->m_integerData_p->m_x.col(j);
      discrim = 0;
      for (int k=0; k<m_g; k++)  discrim +=  algoInte_p->IntegreOneVariable(tmp(find(((m_zCandCurrent == k) + data_p->m_integerData_p->m_notNA.col(j)) == 2)), j);        
      if (discrim > algoInte_p->m_integralenondiscrim(j)){
        m_omegaCurrent(loc) = 1;
        m_miclCurrent +=  discrim;   
      }else{
        m_omegaCurrent(loc) = 0;
        m_miclCurrent +=  algoInte_p->m_integralenondiscrim(j);
      }
      loc++;
    }
  }  
  if (data_p->m_withCategorical){
    for (int j=0; j<data_p->m_categoricalData_p->m_ncols; j++){
      discrim=  algoCate_p->IntegreOneVariableCategoricalDiscrim(j, m_zCandCurrent);
      if (discrim > algoCate_p->m_integralenondiscrim(j)){
        m_omegaCurrent(loc) = 1;
        m_miclCurrent +=  discrim;   
      }else{
        m_omegaCurrent(loc) = 0;
        m_miclCurrent +=  algoCate_p->m_integralenondiscrim(j);
      }
      loc++;
    }  
  }    
}

void AlgorithmMixed::zCandInit(){
  XEMMixed xem(data_p, m_omegaCurrent, m_g);
  xem.Run();
  m_zCandCurrent = xem.FindZMAP();
  m_zStarCurrent = m_zCandCurrent;
}
