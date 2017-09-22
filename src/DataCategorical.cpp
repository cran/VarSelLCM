#include "DataCategorical.h"

DataCategorical::DataCategorical(const S4 & obj){
 this->m_profiles = as<mat>(obj.slot("shortdata"));
 this->m_w  = (as<vec>(obj.slot("weightdata")));
 this->m_nrows = sum(m_w);
 this->m_nmodalities = max(m_profiles);
 // Les données manquantes sont codées par la valeur -1
 m_profiles = m_profiles - 1;
 this->m_ncols = m_profiles.n_cols;
 this->m_nprofiles = m_profiles.n_rows;
 this->m_whotakewhat.resize(m_ncols);
 this->m_dl= trans(m_nmodalities);
 for (int j=0; j<m_ncols; j++){
   m_whotakewhat[j].resize(m_nmodalities(j));
   for (int h=0; h<m_nmodalities(j); h++)  m_whotakewhat[j][h] = find(m_profiles.col(j) == h);
   this->m_dl(j) = m_nmodalities(j) - 1;
 }
}
