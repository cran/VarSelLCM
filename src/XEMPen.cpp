#include "XEMPen.h"

Col<double> dlogGaussianter(const Col<double> & x, const Col<double> & o, const double  mu, const double  sd){
  Col<double>tmpval= - 0.5*pow((x - mu),2)/pow(sd,2) - log(sd * sqrt( 2*M_PI));
  if (sum(o)<x.n_rows)  tmpval(find( o == 0)) = zeros<vec>(x.n_rows-sum(o));
  return  tmpval;
}

Col<double> dlogPoissonter(const Col<double> & x, const Col<double> & o, const double  lambda){
  Col<double>tmpval= -lambda + x*log(lambda);
  for (int i=0; i<tmpval.n_rows;i++) tmpval(i) -= lgamma(x(i)+1);
  if (sum(o)<x.n_rows)  tmpval(find( o == 0)) = zeros<vec>(x.n_rows-sum(o));
  return  tmpval;
}

XEMPen::XEMPen(const S4 * reference_p, const double pen){
  data_p = new DataMixed(as<S4>(reference_p->slot("data")));
  S4 strat = as<S4>(reference_p->slot("strategy"));
  nbSmall = strat.slot("nbSmall");
  iterSmall = strat.slot("iterSmall");
  nbKeep = strat.slot("nbKeep");
  iterKeep = strat.slot("iterKeep");
  tolKeep = strat.slot("tolKeep");
  loglikepen = ones<vec>(nbSmall) * log(0);
  m_nbdegenere=0;
  degeneracy=0;
  g=as<S4>(reference_p->slot("model")).slot("g");
  iterCurrent = iterSmall;
  m_penalty = pen;
  for (int i=0;i<nbSmall;i++){
    vec V = randu<vec>(data_p->m_ncols);
    omegaCand.push_back( zeros<vec>(data_p->m_ncols) );
    nbparamCand.push_back( zeros<vec>(data_p->m_ncols) );
    omegaCand[i] = omegaCand[i] + (V>0.5);
    int vu=0;
    if (data_p->m_withContinuous){
      nbparamCand[i].subvec(vu, vu + data_p->m_continuousData_p->m_ncols - 1) = 2*(g-1)*omegaCand[i].subvec(vu, vu + data_p->m_continuousData_p->m_ncols - 1) + 2*data_p->m_continuousData_p->m_ncols;
      vu += data_p->m_continuousData_p->m_ncols;
    }
    if (data_p->m_withInteger){
      nbparamCand[i].subvec(vu, vu + data_p->m_integerData_p->m_ncols - 1) = (g-1)*omegaCand[i].subvec(vu, vu + data_p->m_integerData_p->m_ncols - 1) + data_p->m_integerData_p->m_ncols;
      vu += data_p->m_integerData_p->m_ncols;
    }
    if (data_p->m_withCategorical){
      nbparamCand[i].subvec(vu, vu + data_p->m_categoricalData_p->m_ncols - 1) = (g-1)*data_p->m_categoricalData_p->m_dl%omegaCand[i].subvec(vu, vu + data_p->m_categoricalData_p->m_ncols - 1) + data_p->m_categoricalData_p->m_dl;
      vu += data_p->m_categoricalData_p->m_ncols;
    }
    // Genere tout les parametres
    paramCand.push_back( ParamMixed(data_p,  ones<vec>(data_p->m_ncols), g) );
    // Mets a egalite les parametres en fonction de omega
    paramCand[i].egalise(data_p, omegaCand[i]);
  } 
  tmplogproba = zeros<mat>(data_p->m_nrows, g);
  maxtmplogproba = ones<vec>(data_p->m_nrows);
  rowsums = ones<mat>(data_p->m_nrows);
  m_weightTMP=zeros<vec>(data_p->m_nrows);
  m_loglikenondis=zeros<vec>(data_p->m_ncols);
  // partie continue
  int repere=0;
  if (data_p->m_withContinuous){
    munondisc = zeros<vec>(data_p->m_continuousData_p->m_ncols);
    sdnondisc = zeros<vec>(data_p->m_continuousData_p->m_ncols);
    for (int j=0; j< data_p->m_continuousData_p->m_ncols; j++){
      m_weightTMP = data_p->m_continuousData_p->m_notNA.col(j);
      munondisc(j) = sum(data_p->m_continuousData_p->m_x.col(j) % m_weightTMP ) / sum(m_weightTMP);
      sdnondisc(j) = sqrt(sum( pow(data_p->m_continuousData_p->m_x.col(j) - munondisc(j),2) % m_weightTMP) / sum(m_weightTMP));
      m_loglikenondis(repere)=sum(dlogGaussianter(data_p->m_continuousData_p->m_x.col(j), data_p->m_continuousData_p->m_notNA.col(j), munondisc(j),  sdnondisc(j)));
      repere++;
    }
  }
  // partie entiere
  if (data_p->m_withInteger){
    lambdanondisc = zeros<vec>(data_p->m_integerData_p->m_ncols);
    for (int j=0; j< data_p->m_integerData_p->m_ncols; j++){
      m_weightTMP =data_p->m_integerData_p->m_notNA.col(j);
      lambdanondisc(j) = sum(data_p->m_integerData_p->m_x.col(j) % m_weightTMP ) / sum(m_weightTMP);
      m_loglikenondis(repere)=sum(dlogPoissonter(data_p->m_integerData_p->m_x.col(j), data_p->m_integerData_p->m_notNA.col(j), lambdanondisc(j)));
      repere++;
    }
  }  
  // partie qualitative
  if (data_p->m_withCategorical){
    for (int j=0; j< data_p->m_categoricalData_p->m_ncols; j++){
      alphanondisc.push_back( zeros<vec>(data_p->m_categoricalData_p->m_nmodalities(j) ));
      alphanondisc[j]=ones<vec>(data_p->m_categoricalData_p->m_nmodalities(j));
      for (int h=0; h< data_p->m_categoricalData_p->m_nmodalities(j); h++) alphanondisc[j](h) = sum(data_p->m_categoricalData_p->m_w(data_p->m_categoricalData_p->m_whotakewhat[j][h]));
      alphanondisc[j] = alphanondisc[j]/ sum(alphanondisc[j]);
      m_loglikenondis(repere)=0;
      for (int h=0; h<data_p->m_categoricalData_p->m_nmodalities(j); h++) m_loglikenondis(repere) = m_loglikenondis(repere)  + (data_p->m_categoricalData_p->m_whotakewhat[j][h].n_rows)*log(alphanondisc[j](h));
      repere++;
    }
  }
  
}

void XEMPen::Run(){
  // Partie Small EM
  for (int ini=0; ini<nbSmall; ini++){
    SwitchCurrent(ini);
    OneEM();
    loglikepen(ini) = ComputeLoglikepen();
  }
  // On conserve les meilleurs initialisations
  uvec indices = sort_index(loglikepen);
  iterCurrent = iterKeep;
  m_nbdegenere = 0;
  for (int tmp1=0; tmp1<nbKeep; tmp1++){
    SwitchCurrent(indices(nbSmall - tmp1 - 1));
    OneEM();
    loglikepen(indices(nbSmall - tmp1 - 1)) = ComputeLoglikepen();
    m_nbdegenere += degeneracy;
  }
  uword  index;
  double indicebest = (loglikepen).max(index);
  SwitchCurrent(index);
  // important pour obtenir la bonne partition a la sortie
  indicebest = ComputeLoglikepen();
  indices = sort_index(loglikepen);  
}

void XEMPen::SwitchCurrent(int ini){
  omegaCurrent_p = &omegaCand[ini];
  paramCurrent_p = &paramCand[ini];
  nbparamCurrent_p = &nbparamCand[ini];
}

colvec XEMPen::FindZMAP(){
  Col<double> zMAP=ones<vec>(tmplogproba.n_rows);
  uword  index;
  double max_val=0;
  for (int i=0; i<tmplogproba.n_rows; i++){
    max_val = (tmplogproba.row(i)).max(index);
    zMAP(i)=index;
  }
  return zMAP;
}

double XEMPen::ComputeLoglikepen(){
  double output=-99999999999999;
  if (degeneracy==0){
    ComputeTmpLogProba();
    maxtmplogproba = max(tmplogproba, 1);
    if (min(maxtmplogproba) == 0){
      output = log(0);
    }else{
      for (int k=0; k<g; k++) tmplogproba.col(k)-=maxtmplogproba;
      tmplogproba = exp(tmplogproba);
      rowsums = sum(tmplogproba,1);
      output = sum(maxtmplogproba) + sum(log(rowsums));
    }
    // Ajouter la penalite
    double nbparam = g-1 + sum(*nbparamCurrent_p);
    output = output - nbparam*m_penalty;
  }  
  return output;
}

void XEMPen::OneEM(){
  degeneracy=0;
  double loglikepen = ComputeLoglikepen(), prec = -99999999999999;
  int it=0;
  ParamMixed backupparam ;
  Col<double> backupmodel, backupnbparam;
  while ( (it<iterCurrent) && ((loglikepen-prec)>tolKeep) && (degeneracy==0)){
    it ++;
    Estep();
    backupparam = (* paramCurrent_p);
    backupmodel = (* omegaCurrent_p);
    backupnbparam=(*nbparamCurrent_p);
    Mstep();
    prec = loglikepen;
    loglikepen = ComputeLoglikepen();
  }
  // Une verif
  if ((degeneracy==0)&&(prec>(loglikepen+tolKeep))){
     int repere=0;
  } 
}

void XEMPen::ComputeTmpLogProba(){
  for (int k=0; k<g; k++) tmplogproba.col(k) = log(paramCurrent_p->m_pi(k))*ones<vec>(data_p->m_nrows);
  int repere=0;
  if (data_p->m_withContinuous){
    for (int j=0; j< data_p->m_continuousData_p->m_ncols; j++){
      for (int k=0; k<g; k++) tmplogproba.col(k) += dlogGaussianter(data_p->m_continuousData_p->m_x.col(j), data_p->m_continuousData_p->m_notNA.col(j), paramCurrent_p->m_paramContinuous.m_mu(k,j),  paramCurrent_p->m_paramContinuous.m_sd(k,j));
      repere++;
    } 
  }
  if (data_p->m_withInteger){
    for (int j=0; j< data_p->m_integerData_p->m_ncols; j++){
      for (int k=0; k<g; k++) tmplogproba.col(k) += dlogPoissonter(data_p->m_integerData_p->m_x.col(j), data_p->m_integerData_p->m_notNA.col(j), paramCurrent_p->m_paramInteger.m_lambda(k,j));
      repere++;
    } 
  }
  if  (data_p->m_withCategorical){
    for (int j=0; j<data_p->m_categoricalData_p->m_ncols; j++){
      for (int k=0; k<g; k++){
        Col<double> probavari=ones<vec>(data_p->m_nrows);
        for (int h=0; h<data_p->m_categoricalData_p->m_nmodalities(j); h++) probavari(data_p->m_categoricalData_p->m_whotakewhat[j][h]) = probavari(data_p->m_categoricalData_p->m_whotakewhat[j][h])*log(paramCurrent_p->m_paramCategorical.m_alpha[j](k,h));
        tmplogproba.col(k) += probavari;
      }      
      repere++;
    }    
  }
}

void XEMPen::Mstep(){
  paramCurrent_p->m_pi = trans(sum(tmplogproba,0));
  paramCurrent_p->m_pi = paramCurrent_p->m_pi / sum(paramCurrent_p->m_pi);
  // partie continue
  int repere=0;
  if (data_p->m_withContinuous){
    for (int j=0; j< data_p->m_continuousData_p->m_ncols; j++){
      Col<double> tmpmu = paramCurrent_p->m_paramContinuous.m_mu.col(j)*0;
      Col<double> tmpsd = paramCurrent_p->m_paramContinuous.m_sd.col(j)*0;
      double tmploglike=0;
      for (int k=0; k<g; k++){
        m_weightTMP = tmplogproba.col(k) % data_p->m_continuousData_p->m_notNA.col(j);
        tmpmu(k) = sum(data_p->m_continuousData_p->m_x.col(j) % m_weightTMP ) / sum(m_weightTMP);
        tmpsd(k) = sqrt(sum( pow(data_p->m_continuousData_p->m_x.col(j) - tmpmu(k),2) % m_weightTMP) / sum(m_weightTMP));
        tmploglike+= sum(dlogGaussianter(data_p->m_continuousData_p->m_x.col(j), data_p->m_continuousData_p->m_notNA.col(j), tmpmu(k),  tmpsd(k)) % m_weightTMP);
      }
      if (any(tmpsd<0.0001)) {degeneracy=1;}
      if (tmploglike > (m_loglikenondis(repere) + 2*(g-1)*m_penalty)){
        (*omegaCurrent_p)(repere)=1;
        (*nbparamCurrent_p)(repere)=2*g;
        paramCurrent_p->m_paramContinuous.m_mu.col(j) = tmpmu;
        paramCurrent_p->m_paramContinuous.m_sd.col(j) = tmpsd;

      }else{
        (*omegaCurrent_p)(repere)=0;
        (*nbparamCurrent_p)(repere)=2;
        paramCurrent_p->m_paramContinuous.m_mu.col(j) = ones<vec>(g) * munondisc(j);
        paramCurrent_p->m_paramContinuous.m_sd.col(j) = ones<vec>(g) * sdnondisc(j);
      }      
      repere++;
    }
  }
  // partie entiere
  if (data_p->m_withInteger){
    for (int j=0; j< data_p->m_integerData_p->m_ncols; j++){
      Col<double> tmplambda = paramCurrent_p->m_paramInteger.m_lambda.col(j);
      double tmploglike=0;
      for (int k=0; k<g; k++){
        m_weightTMP = tmplogproba.col(k) % data_p->m_integerData_p->m_notNA.col(j);
        tmplambda(k) = sum(data_p->m_integerData_p->m_x.col(j) % m_weightTMP ) / sum(m_weightTMP);
        tmploglike+= sum(dlogPoissonter(data_p->m_integerData_p->m_x.col(j), data_p->m_integerData_p->m_notNA.col(j), tmplambda(k)) % m_weightTMP);
      }
      if  (tmploglike!=tmploglike)  {degeneracy=1;}
      if  ((tmploglike!=tmploglike) || (tmploglike > (m_loglikenondis(repere) + (g-1)*m_penalty))){
        (*omegaCurrent_p)(repere)=1;
        (*nbparamCurrent_p)(repere)=g;
        paramCurrent_p->m_paramInteger.m_lambda.col(j) = tmplambda;
      }else{
        (*omegaCurrent_p)(repere)=0;
        (*nbparamCurrent_p)(repere)=1;   
        paramCurrent_p->m_paramInteger.m_lambda.col(j) =  ones<vec>(g) *lambdanondisc(j);
      } 
      repere++;
    }
  }
  // partie qualitative
  if (data_p->m_withCategorical){
    for (int j=0; j< data_p->m_categoricalData_p->m_ncols; j++){
      Mat<double> tmpalpha = zeros<mat>(g, data_p->m_categoricalData_p->m_nmodalities(j));
      double tmploglike=0;
      for (int h=0; h< data_p->m_categoricalData_p->m_nmodalities(j); h++)       tmpalpha.col(h) = trans(sum(tmplogproba.rows(data_p->m_categoricalData_p->m_whotakewhat[j][h])));
      for (int k=0; k<g; k++){  
        double coeff=sum(tmpalpha.row(k));
        tmpalpha.row(k) = tmpalpha.row(k)/sum(tmpalpha.row(k));
        for (int h=0; h<data_p->m_categoricalData_p->m_nmodalities(j); h++) tmploglike+=log(tmpalpha(k,h)) * tmpalpha(k,h) * coeff;
      }
      //if (tmploglike!=tmploglike) degeneracy=1;
      if ((tmploglike!=tmploglike) || (tmploglike > (m_loglikenondis(repere) + data_p->m_categoricalData_p->m_dl(j)*(g-1)*m_penalty))){
        (*omegaCurrent_p)(repere)=1;
        (*nbparamCurrent_p)(repere)=data_p->m_categoricalData_p->m_dl(j)*g;
        paramCurrent_p->m_paramCategorical.m_alpha[j] = tmpalpha;
        for (int k=0; k<g; k++){if (any(tmpalpha.row(k)==1)) degeneracy=1;}
      }else{
        (*omegaCurrent_p)(repere)=0;
        (*nbparamCurrent_p)(repere)=data_p->m_categoricalData_p->m_dl(j);           
        for (int k=0; k<g; k++) paramCurrent_p->m_paramCategorical.m_alpha[j].row(k) = trans(alphanondisc[j]);
      }
    repere++;
    }
  }  
}


void XEMPen::Output(S4 * reference_p){
  as<S4>(reference_p->slot("model")).slot("omega")= wrap(trans( (*omegaCurrent_p)));
  if (m_nbdegenere<nbKeep){
    // Partie Continue 
    if (data_p->m_withContinuous){
      as<S4>(as<S4>(reference_p->slot("param")).slot("paramContinuous")).slot("mu") = wrap(trans( paramCurrent_p->m_paramContinuous.m_mu));
      as<S4>(as<S4>(reference_p->slot("param")).slot("paramContinuous")).slot("sd") = wrap(trans(paramCurrent_p->m_paramContinuous.m_sd));
    }
    // Partie entiere 
    if (data_p->m_withInteger)  as<S4>(as<S4>(reference_p->slot("param")).slot("paramInteger")).slot("lambda") = wrap(trans(paramCurrent_p->m_paramInteger.m_lambda));
    // partie qualitative
    if (data_p->m_withCategorical) as<S4>(as<S4>(reference_p->slot("param")).slot("paramCategorical")).slot("alpha") = wrap(paramCurrent_p->m_paramCategorical.m_alpha);    
    double nbparam = g-1 + sum(*nbparamCurrent_p);
    as<S4>(reference_p->slot("criteria")).slot("degeneracyrate") = m_nbdegenere/nbKeep;
    as<S4>(reference_p->slot("criteria")).slot("loglikelihood") = max(loglikepen) +  nbparam*m_penalty;;
    as<S4>(reference_p->slot("param")).slot("pi") = wrap(trans(paramCurrent_p->m_pi));
    Estep();
    as<S4>(reference_p->slot("partitions")).slot("tik") = wrap(tmplogproba);
    as<S4>(reference_p->slot("partitions")).slot("zMAP") = wrap(FindZMAP());
  }else{
    as<S4>(reference_p->slot("criteria")).slot("degeneracyrate") = 1;
  }  
}
