#include "DataInteger.h"

DataInteger::DataInteger(const S4 & obj){
 this->m_x = as<mat>(obj.slot("data"));
 this->m_priors  = as<mat>(obj.slot("priors")); 
 this->m_nrows = m_x.n_rows;
 this->m_ncols = m_x.n_cols;
 this->m_notNA = as<mat>(obj.slot("notNA"));
}
