/*
Cette classes définie les paramètres pour des données qualitatives

Ces éléments sont:
m_alpha : liste où m_alpha[j] donne les probabilités associées à chaque modalité pour chaque classe
m_pi : proportions

*/
#ifndef ParamCategorical_H
#define ParamCategorical_H

#include "Param.h"
#include "DataCategorical.h"

class ParamCategorical : public Param{
  public:
  vector< Mat<double> > m_alpha;
  
  ParamCategorical();
  ParamCategorical(const ParamCategorical & param);
  ParamCategorical(const DataCategorical *, const colvec & , const int &);
  ~ParamCategorical(){};
    void egalise(const colvec );

};
#endif
