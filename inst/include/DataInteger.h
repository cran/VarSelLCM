/*
Cette classe contient les éléments relatifs aux données continues

Ces élements sont:
m_x : matrice des données
m_priors : matrice des hyper-paramètres (1 ligne par variable)
m_notNA : matrice binaire indiquant si l'observation est manquante ou non (attention dans m_x les valeurs
manquantes ont été artificellement remplacées par la valeur 0 pour pouvoir utiliser armadillo)
*/
#ifndef DataInteger_H
#define DataInteger_H

#include "Data.h"

class DataInteger: public Data{
  public:
    Mat<double> m_x, m_priors, m_notNA;
  
  DataInteger(){};
  DataInteger(const S4 &);
  ~DataInteger(){};
  
};
#endif