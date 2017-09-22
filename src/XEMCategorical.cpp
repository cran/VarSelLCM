#include "XEMCategorical.h"

XEMCategorical::XEMCategorical(const DataCategorical * datapasse, const S4 * reference_p){
  paramEstim = as<S4>(reference_p->slot("strategy")).slot("paramEstim");
  if (paramEstim){
    InitCommumParamXEM(as<S4>(reference_p->slot("model")).slot("omega"), as<S4>(reference_p->slot("model")).slot("g"), as<S4>(reference_p->slot("strategy")));
    InitSpecificParamXEMCategorical(datapasse);
  }
}

XEMCategorical::XEMCategorical(const DataCategorical * datapasse, const colvec & omega, const int & g){
  paramEstim = TRUE;
  InitCommumParamXEM(omega, g);  
  InitSpecificParamXEMCategorical(datapasse);
}

void XEMCategorical::InitSpecificParamXEMCategorical(const DataCategorical * datapasse){
  data_p = datapasse;
  for (int i=0;i<nbSmall;i++) paramCand.push_back( ParamCategorical(data_p, omega, g) );
  tmplogproba = zeros<mat>(data_p->m_nprofiles, g);
  maxtmplogproba = ones<vec>(data_p->m_nprofiles);
  rowsums = ones<mat>(data_p->m_nprofiles);
  tmpval=zeros<vec>(data_p->m_nprofiles);
}

void XEMCategorical::SwitchParamCurrent(int ini){paramCurrent_p = &paramCand[ini];}

void XEMCategorical::ComputeTmpLogProba(){
  for (int k=0; k<g; k++){
    tmpval = zeros<vec>(data_p->m_nprofiles) + log(paramCurrent_p->m_pi(k));
    for (int j=0; j<sum(omega); j++){
      for (int h=0; h<data_p->m_nmodalities(location(j)); h++){
        tmpval(data_p->m_whotakewhat[location(j)][h]) += log(paramCurrent_p->m_alpha[j](k,h));
      }
    }
    tmplogproba.col(k) = tmpval;
  }  
}

double XEMCategorical::ComputeLogLike(){
  ComputeTmpLogProba();
  maxtmplogproba = max(tmplogproba, 1);
  for (int k=0; k<g; k++)tmplogproba.col(k) -= maxtmplogproba;
  tmplogproba = exp(tmplogproba);
  rowsums = sum(tmplogproba,1);
  return sum(maxtmplogproba % data_p->m_w )+ sum((log(rowsums)) %data_p->m_w );
}

void XEMCategorical::Mstep(){
  for (int k=0; k<g; k++) paramCurrent_p->m_pi(k) = sum((tmplogproba.col(k)) % data_p->m_w);
  paramCurrent_p->m_pi = paramCurrent_p->m_pi / sum(paramCurrent_p->m_pi);
  for (int j=0; j< sum(omega); j++){
    for (int h=0; h< data_p->m_nmodalities(location(j)); h++){
      paramCurrent_p->m_alpha[j].col(h) = trans(trans(data_p->m_w(data_p->m_whotakewhat[location(j)][h])) *tmplogproba.rows(data_p->m_whotakewhat[location(j)][h]));
    }
    for (int k=0; k<g; k++){
      paramCurrent_p->m_alpha[j].row(k) = paramCurrent_p->m_alpha[j].row(k)/sum(paramCurrent_p->m_alpha[j].row(k));        
    }
  }
}

void XEMCategorical::Output(S4 * reference_p){
  if (paramEstim){
    vector< Mat<double> >  alpha;
    alpha.resize(data_p->m_ncols);
    int loc=0;
    for (int j=0; j<data_p->m_ncols; j++){
      alpha[j] = zeros<mat>(g, data_p->m_nmodalities(j));
      if (omega(j) == 0){
        vec tmp = zeros<vec>(data_p->m_nmodalities(j));
        for (int h=0; h<data_p->m_nmodalities(j); h++)
        tmp(h) = sum(data_p->m_w(data_p->m_whotakewhat[j][h]));
        tmp = tmp/sum(tmp);
        for (int k=0; k<g; k++) alpha[j].row(k)=trans(tmp);
        for (int h=0; h<data_p->m_nmodalities(j); h++){
          loglikeoutput += sum(data_p->m_w(data_p->m_whotakewhat[j][h])) * log(alpha[j](0,h));
        }
        
      }else{
        alpha[j]=paramCurrent_p->m_alpha[loc];
        loc ++;
      }
    }
    as<S4>(reference_p->slot("param")).slot("alpha") = wrap(alpha);    
    as<S4>(reference_p->slot("criteria")).slot("loglikelihood") = loglikeoutput;
    as<S4>(reference_p->slot("criteria")).slot("degeneracyrate") = m_nbdegenere/nbKeep;
    if (sum(omega)==0){
      as<S4>(reference_p->slot("param")).slot("pi") = wrap(ones<vec>(g) * (1/g));
      tmplogproba = ones<mat>(data_p->m_nrows,g) * (1/g);
    }else{
      as<S4>(reference_p->slot("param")).slot("pi") = wrap(trans(paramCurrent_p->m_pi));
      Estep();
    }    
    as<S4>(reference_p->slot("partitions")).slot("tik") = wrap(tmplogproba);
    as<S4>(reference_p->slot("partitions")).slot("zMAP") = wrap(FindZMAP());
  }
}
