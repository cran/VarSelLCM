#include "XEMContinuous.h"


Col<double> dlogGaussian(const Col<double> & x, const Col<double> & o, const double  mu, const double  sd){
  Col<double>tmpval= - 0.5*pow((x - mu),2)/pow(sd,2) - log(sd * sqrt( 2*M_PI));
  if (sum(o)<x.n_rows)  tmpval(find( o == 0)) = zeros<vec>(x.n_rows-sum(o));
  return  tmpval;
}

XEMContinuous::XEMContinuous(const DataContinuous * datapasse, const S4 * reference_p){
  paramEstim = as<S4>(reference_p->slot("strategy")).slot("paramEstim");
  if (paramEstim){
    InitCommumParamXEM(as<S4>(reference_p->slot("model")).slot("omega"), as<S4>(reference_p->slot("model")).slot("g"), as<S4>(reference_p->slot("strategy")));
    InitSpecificParamXEMContinuous(datapasse);
  }
}

XEMContinuous::XEMContinuous(const DataContinuous * datapasse, const colvec & omega, const int & g){
  paramEstim = TRUE;
  InitCommumParamXEM(omega, g);  
  InitSpecificParamXEMContinuous(datapasse);
}

void XEMContinuous::InitSpecificParamXEMContinuous(const DataContinuous * datapasse){
  data_p = datapasse;
  for (int i=0;i<nbSmall;i++) paramCand.push_back( ParamContinuous(data_p, omega, g) );
  tmplogproba = zeros<mat>(data_p->m_nrows, g);
  maxtmplogproba = ones<vec>(data_p->m_nrows);
  rowsums = ones<mat>(data_p->m_nrows);
  m_weightTMP = zeros<vec>(data_p->m_nrows);
}

void XEMContinuous::SwitchParamCurrent(int ini){paramCurrent_p = &paramCand[ini];}

void XEMContinuous::ComputeTmpLogProba(){
  for (int k=0; k<g; k++){
    tmplogproba.col(k) = zeros<vec>(data_p->m_nrows) + log(paramCurrent_p->m_pi(k));
    for (int j=0; j< sum(omega); j++) tmplogproba.col(k) += dlogGaussian(data_p->m_x.col(location(j)), data_p->m_notNA.col(location(j)), paramCurrent_p->m_mu(k,j),  paramCurrent_p->m_sd(k,j));
  }
}

void XEMContinuous::Mstep(){
  paramCurrent_p->m_pi = trans(sum(tmplogproba,0));
  paramCurrent_p->m_pi = paramCurrent_p->m_pi / sum(paramCurrent_p->m_pi);
  for (int k=0; k<g; k++){
    for (int j=0; j< sum(omega); j++){
      m_weightTMP = tmplogproba.col(k) % data_p->m_notNA.col(location(j));
      paramCurrent_p->m_mu(k,j) = sum(data_p->m_x.col(location(j)) % m_weightTMP ) / sum(m_weightTMP);
      paramCurrent_p->m_sd(k,j) = sqrt(sum( pow(data_p->m_x.col(location(j)) - paramCurrent_p->m_mu(k,j),2) % m_weightTMP) / sum(m_weightTMP));
    }
    if (any(paramCurrent_p->m_sd.row(k)<0.0001)) {degeneracy=1;}
  }
}

void XEMContinuous::Output(S4 * reference_p){
  if (paramEstim){
    if (m_nbdegenere<nbKeep){
      Mat<double> mu=ones<mat>(g, data_p->m_ncols);
      Mat<double> sd=ones<mat>(g, data_p->m_ncols);
      int loc=0;
      for (int j=0; j<data_p->m_ncols; j++){
        if (omega(j) == 0){
          vec tmp = data_p->m_x.col(j);
          vec keep = tmp(find(data_p->m_notNA.col(j) == 1));
          mu.col(j) = mu.col(j)*mean(keep);
          sd.col(j) = sd.col(j)*sqrt(sum(pow(( keep - mean(keep)),2) ) / keep.n_rows);
          loglikeoutput += sum(dlogGaussian(data_p->m_x.col(j), data_p->m_notNA.col(j), mu(0,j), sd(0,j)));
        }else{
          mu.col(j) = paramCurrent_p->m_mu.col(loc);
          sd.col(j) = paramCurrent_p->m_sd.col(loc);
          loc ++;
        }
      }
      as<S4>(reference_p->slot("param")).slot("mu") = wrap(trans(mu));
      as<S4>(reference_p->slot("param")).slot("sd") = wrap(trans(sd));
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
