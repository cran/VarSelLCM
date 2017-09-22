/*
Cette classe est héritière d'Algorithm et permet d'effectuer le choix de modèle au sens du critère MICL pour des données qualitatives

Ces élements sont:
ceux de la classe Algorithm
data_p : pointeur sur le jeux de données qualitatif

*/
#ifndef AlgorithmCategorical_H
#define AlgorithmCategorical_H

#include "DataCategorical.h"
#include "Algorithm.h"


class AlgorithmCategorical : public Algorithm{ 
  public:
  const DataCategorical * data_p;
  
  // Constructeur et destructeur par défaut
  AlgorithmCategorical(){};
  ~AlgorithmCategorical(){};
  
  // Constructeur utilisé qui prend en entrée un pointeur du jeu de données et un pointeur de l'object S4 retourné sous R
  AlgorithmCategorical(const DataCategorical *, const S4 *);
  // Permet d'inialiser les éléments de la classe qui sont liées aux données
 // void InitSpecificParamAlgo(const DataCategorical * datapasse);
  // Calcul de la vraisemblance intégrée pour une variable non discriminante (indice en entrée)
  double IntegreOneVariableCategoricalNotDiscrim(const int &);
  // Calcul de la vraisemblance complétée intégrée pour une variable (indice en entrée)
  double IntegreOneVariableCategoricalDiscrim(const int &);
  
  // introduit pour faciliter le cas mixte
  double IntegreOneVariableCategoricalDiscrim(const int &, const colvec &);
  
  // Les trois fonction suivantes sont à redéfinir pour chaque classe héritiaire d'Algorithm
  // Calcul la vraisemblance complétée intégrée pour le modèle m_omegaCurrent et la partition m_zCandCurrent
  virtual double Integre_Complete_Like_Cand() ;
  // Définit m_omegaCurrent comme le modèle maximisant la vraisemblance complétée intégrée pour la partition m_zCandCurrent
  virtual void Optimize_model();
  // Définit la partition m_zCandCurrent initiale associée au modèle initial
  virtual void zCandInit();
};
#endif