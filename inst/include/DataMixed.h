/*
  Cette classe contient les éléments relatifs aux données mixtes

Ces élements sont:
*/
#ifndef DataMixed_H
#define DataMixed_H
  
#include "DataContinuous.h"
#include "DataInteger.h"
#include "DataCategorical.h"


class DataMixed : public Data{
  public:
    DataContinuous * m_continuousData_p;
    DataInteger * m_integerData_p;
    DataCategorical * m_categoricalData_p;
    bool m_withContinuous, m_withInteger, m_withCategorical;

  DataMixed();
  DataMixed(const S4 &);
  ~DataMixed(){};
};
#endif