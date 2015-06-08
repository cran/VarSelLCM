/*
Cette classes définie les paramètres pour des données continues

Ces éléments sont:
m_mu : matrice des moyennes
m_sd : matrice des ecarts-types
m_pi : proportions

*/
#ifndef ParamContinuous_H
#define ParamContinuous_H

#include "DataContinuous.h"
#include "Param.h"

class ParamContinuous : public Param{
  public:
  Mat<double> m_mu, m_sd;
    
  ParamContinuous();
  ParamContinuous(const ParamContinuous & param);
  ParamContinuous(const DataContinuous *, const colvec & , const int &);
  ~ParamContinuous(){};
};
#endif