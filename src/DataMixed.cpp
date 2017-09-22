#include "DataMixed.h"

DataMixed::DataMixed(const S4 & obj){
  this->m_nrows = obj.slot("n");
  this->m_ncols = obj.slot("d");
  this->m_withContinuous = obj.slot("withContinuous");
  this->m_withInteger = obj.slot("withInteger");
  this->m_withCategorical = obj.slot("withCategorical");
  if (m_withContinuous) m_continuousData_p = new DataContinuous(as<S4>(obj.slot("dataContinuous")));
  if (m_withInteger) m_integerData_p = new DataInteger(as<S4>(obj.slot("dataInteger")));
  if (m_withCategorical) m_categoricalData_p = new DataCategorical(as<S4>(obj.slot("dataCategorical")));
}
