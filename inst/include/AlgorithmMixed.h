/*
Cette classe est héritière d'Algorithm et permet d'effectuer le choix de modèle au sens du critère MICL pour des données continues

Ces élements sont:
ceux de la classe Algorithm
data_p : pointeur sur le jeux de données mixed

*/
#ifndef AlgorithmMixed_H
#define AlgorithmMixed_H

#include "DataMixed.h"
#include "AlgorithmContinuous.h"
#include "AlgorithmInteger.h"
#include "AlgorithmCategorical.h"


class AlgorithmMixed : public Algorithm{ 
  
  public:
  const DataMixed * data_p;
  AlgorithmContinuous * algoCont_p;
  AlgorithmInteger * algoInte_p;
  AlgorithmCategorical * algoCate_p;
  
  // Constructeur et destructeur par défaut
  AlgorithmMixed(){};
  ~AlgorithmMixed(){};
  
  // Constructeur utilisé qui prend en entrée un pointeur du jeu de données et un pointeur de l'object S4 retourné sous R
  AlgorithmMixed(const DataMixed *, const S4 *);
  // Les trois fonction suivantes sont à redéfinir pour chaque classe héritiaire d'Algorithm
  // Calcul la vraisemblance complétée intégrée pour le modèle m_omegaCurrent et la partition m_zCandCurrent
  virtual double Integre_Complete_Like_Cand();
  // Définit m_omegaCurrent comme le modèle maximisant la vraisemblance complétée intégrée pour la partition m_zCandCurrent
  virtual void Optimize_model();
  // Définit la partition m_zCandCurrent initiale associée au modèle initial
  virtual void zCandInit();
};
#endif