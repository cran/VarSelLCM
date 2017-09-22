#include "XEMInteger.h"


Col<double> dlogPoisson(const Col<double> & x, const Col<double> & o, const double  lambda){
  Col<double>tmpval= -lambda + x*log(lambda);
  for (int i=0; i<tmpval.n_rows;i++) tmpval(i) -= lgamma(x(i)+1);
  if (sum(o)<x.n_rows)  tmpval(find( o == 0)) = zeros<vec>(x.n_rows-sum(o));
  return  tmpval;
};

XEMInteger::XEMInteger(const DataInteger * datapasse, const S4 * reference_p){
  paramEstim = as<S4>(reference_p->slot("strategy")).slot("paramEstim");
  if (paramEstim){
    InitCommumParamXEM(as<S4>(reference_p->slot("model")).slot("omega"), as<S4>(reference_p->slot("model")).slot("g"), as<S4>(reference_p->slot("strategy")));
    InitSpecificParamXEMInteger(datapasse);
  }
}

XEMInteger::XEMInteger(const DataInteger * datapasse, const colvec & omega, const int & g){
  paramEstim = TRUE;
  InitCommumParamXEM(omega, g);  
  InitSpecificParamXEMInteger(datapasse);
}

void XEMInteger::InitSpecificParamXEMInteger(const DataInteger * datapasse){
  data_p = datapasse;
  for (int i=0;i<nbSmall;i++) paramCand.push_back( ParamInteger(data_p, omega, g) );
  tmplogproba = zeros<mat>(data_p->m_nrows, g);
  maxtmplogproba = ones<vec>(data_p->m_nrows);
  rowsums = ones<mat>(data_p->m_nrows);
  m_weightTMP = zeros<vec>(data_p->m_nrows);
}

void XEMInteger::SwitchParamCurrent(int ini){paramCurrent_p = &paramCand[ini];}

void XEMInteger::ComputeTmpLogProba(){
  for (int k=0; k<g; k++){
    tmplogproba.col(k) = zeros<vec>(data_p->m_nrows) + log(paramCurrent_p->m_pi(k));
    for (int j=0; j< sum(omega); j++) tmplogproba.col(k) += dlogPoisson(data_p->m_x.col(location(j)), data_p->m_notNA.col(location(j)), paramCurrent_p->m_lambda(k,j));
  }
}


void XEMInteger::Mstep(){
  paramCurrent_p->m_pi = trans(sum(tmplogproba,0));
  paramCurrent_p->m_pi = paramCurrent_p->m_pi / sum(paramCurrent_p->m_pi);
  for (int k=0; k<g; k++){
    for (int j=0; j< sum(omega); j++){
      m_weightTMP = tmplogproba.col(k) % data_p->m_notNA.col(location(j));
      paramCurrent_p->m_lambda(k,j) = sum(data_p->m_x.col(location(j)) % m_weightTMP ) / sum(m_weightTMP);
    }
  }
}

void XEMInteger::Output(S4 * reference_p){
  if (paramEstim){
    if (m_nbdegenere<nbKeep){
      Mat<double> lambda=ones<mat>(g, data_p->m_ncols);
      int loc=0;
      for (int j=0; j<data_p->m_ncols; j++){
        if (omega(j) == 0){
          vec tmp = data_p->m_x.col(j);
          vec keep = tmp(find(data_p->m_notNA.col(j) == 1));
          lambda.col(j) = lambda.col(j)*mean(keep);
          loglikeoutput += sum(dlogPoisson(data_p->m_x.col(j), data_p->m_notNA.col(j), lambda(0,j)));
        }else{
          lambda.col(j) = paramCurrent_p->m_lambda.col(loc);
          loc ++;
        }
      }
      as<S4>(reference_p->slot("param")).slot("lambda") = wrap(trans(lambda));
      as<S4>(reference_p->slot("criteria")).slot("loglikelihood") = loglikeoutput;
      as<S4>(reference_p->slot("criteria")).slot("degeneracyrate") =  double(m_nbdegenere)/double(nbKeep);
      if (sum(omega)==0){
        as<S4>(reference_p->slot("param")).slot("pi") = wrap(ones<vec>(g) * (1/g));
        tmplogproba = ones<mat>(data_p->m_nrows,g) * (1/g);
      }else{
        as<S4>(reference_p->slot("param")).slot("pi") = wrap(trans(paramCurrent_p->m_pi));
        Estep();
      } 
      as<S4>(reference_p->slot("partitions")).slot("tik") = wrap(tmplogproba);
      as<S4>(reference_p->slot("partitions")).slot("zMAP") = wrap(FindZMAP());
    }else{
      as<S4>(reference_p->slot("criteria")).slot("degeneracyrate") = 1;
    }
    
  }
}
