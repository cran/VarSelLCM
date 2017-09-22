#include "Algorithm.h"

void Algorithm::InitCommumParamAlgo(const int & g, const int & nbinit, const int & nrows, const int & ncols){
  m_zStarBest = ones<vec>(nrows);
  m_zStarCurrent = ones<vec>(nrows);
  m_zCandCurrent = ones<vec>(nrows);
  m_miclBest = log(0);
  m_miclCurrent = log(0);
  m_g = g;
  m_integralenondiscrim = ones<vec>(0);
  Mat<double> tmp=randu<mat>(ncols, nbinit);
  omegainit = ones<mat>(ncols, nbinit);
  uvec loc;
  Col<double> tmpvec ;
  //omegainit.col(0) = ones<vec>(ncols);
  for (int ini=1; ini<nbinit; ini++){
    loc = find(tmp.col(ini)>0.5);
    tmpvec = zeros<vec>(ncols);
    tmpvec(loc)+=1;
    omegainit.col(ini)=tmpvec;
  }
}

void Algorithm::Run(S4 * output_p){
  if (vbleSelec){
    double prec = log(0);
    m_omegaBest = omegainit.col(0);
    int cvrate=0;
    for (int ini=0; ini<omegainit.n_cols; ini++){
      m_omegaCurrent = omegainit.col(ini);
      if (sum(m_omegaCurrent)==0)
        m_omegaCurrent = ones<vec>(omegainit.n_rows);
      prec = log(0);
      zCandInit();
      m_miclCurrent = Integre_Complete_Like_Cand();
      while (prec < m_miclCurrent){
        // cout << trans(m_omegaCurrent ) << " " << m_miclCurrent << endl;
        prec = m_miclCurrent;
        Optimize_partition();
        Optimize_model();
      }
      if (m_miclCurrent > m_miclBest){
        m_miclBest = m_miclCurrent;
        m_omegaBest = m_omegaCurrent;
        m_zStarBest = m_zStarCurrent;
        cvrate = 1;
      }else if (m_miclCurrent ==  m_miclBest){
        cvrate ++;
      }
    }
    as<S4>(output_p->slot("model")).slot("omega") = wrap(trans(m_omegaBest));
    as<S4>(output_p->slot("partitions")).slot("zOPT") = wrap(trans(m_zStarBest));
    as<S4>(output_p->slot("criteria")).slot("MICL") = m_miclBest;  
    as<S4>(output_p->slot("criteria")).slot("cvrate") = cvrate;  
  }
}

void Algorithm::ComputeMICL(S4 * output_p){
  double prec = log(0);
  Col<double> tmp = as<S4>(output_p->slot("model")).slot("omega");
  //Col<double> tmp2= as<S4>(output_p->slot("partitions")).slot("zMAP");
  m_omegaCurrent = tmp;
  m_omegaBest= tmp;
  for (int ini=0; ini<50; ini++){
    prec = log(0);
    m_omegaCurrent = tmp;
    //m_zCandCurrent = tmp2-1;
    //m_zStarCurrent = m_zCandCurrent;
    zCandInit();
    m_miclCurrent = Integre_Complete_Like_Cand();
    while (prec < m_miclCurrent){
      prec = m_miclCurrent;
      Optimize_partition();
      Optimize_partition();
      Optimize_partition();
    }
    if (m_miclCurrent > m_miclBest){
      m_miclBest = m_miclCurrent;
      m_zStarBest = m_zStarCurrent;
    }
  }
  as<S4>(output_p->slot("partitions")).slot("zOPT") = wrap(trans(m_zStarBest));
  as<S4>(output_p->slot("criteria")).slot("MICL") = m_miclBest;
}

void Algorithm::Optimize_partition(){
  int chgt = m_zCandCurrent.n_rows ;
  double critere_cand=0;
  m_zCandCurrent = m_zStarCurrent;
  while (chgt > 0 ){
    ivec who = randi<ivec>(chgt, distr_param(0,  m_zStarBest.n_rows -1));
    chgt = 0;
    for (int it=0; it < who.n_rows; it++){ 
      for (int k=0; k<m_g; k++){
        m_zCandCurrent( who(it) ) = k;
        critere_cand = Integre_Complete_Like_Cand();
        if (critere_cand > m_miclCurrent){
          m_zStarCurrent(who(it)) = k;
          m_miclCurrent = critere_cand;
          chgt = it ;
        }
      }
      m_zCandCurrent( who(it) ) = m_zStarCurrent( who(it) );
    }
  }
}
