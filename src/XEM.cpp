#include "XEM.h"

void XEM::InitCommumParamXEM(const colvec & om, const int & gv){
  nbSmall = 10;
  iterSmall = 20;
  nbKeep = 1;
  iterKeep = 1;
  tolKeep = 0.001;
  m_nbdegenere=0;
  loglikeSmall = ones<vec>(nbSmall) * log(0);
  omega = om;
  g=gv;
  location = find(omega == 1);
  iterCurrent = iterSmall;
  loglikeoutput=log(0);
}

void XEM::InitCommumParamXEM(const colvec & om, const int & gv, const S4 & strategy){
  nbSmall = strategy.slot("nbSmall");
  iterSmall = strategy.slot("iterSmall");
  nbKeep = strategy.slot("nbKeep");
  iterKeep = strategy.slot("iterKeep");
  tolKeep = strategy.slot("tolKeep");
  loglikeSmall = ones<vec>(nbSmall) * log(0);
  m_nbdegenere=0;
  omega = om;
  g=gv;
  location = find(omega == 1);
  iterCurrent = iterSmall;
  loglikeoutput=log(0);
  if (sum(omega) == 0){
    loglikeoutput = 0;
  }
}

void XEM::Run(){
  if (paramEstim){
    // Partie Small EM
    if (sum(omega)>0){
      for (int ini=0; ini<nbSmall; ini++){
        SwitchParamCurrent(ini);
        OneEM();
        loglikeSmall(ini) = ComputeLogLike();
      }
      if (any(loglikeSmall != loglikeSmall))
      loglikeSmall(find((loglikeSmall != loglikeSmall))) = -99999999999999999*ones<vec>(sum(loglikeSmall != loglikeSmall));
      
      if (any(loglikeSmall == log(0)))
      loglikeSmall(find((loglikeSmall == log(0)))) = -99999999999999999*ones<vec>(sum(loglikeSmall == log(0)));
      
      // On conserve les meilleurs initialisations
      uvec indices = sort_index(loglikeSmall);
      iterCurrent = iterKeep;
      m_nbdegenere = 0;
      // if (nbSmall > nbKeep)    loglikeSmall( indices.head(nbSmall - nbKeep) ) = loglikeSmall( indices.head(nbSmall - nbKeep) ) + log(0);
      int degenere = 0;
      for (int tmp1=0; tmp1<nbKeep; tmp1++){
        SwitchParamCurrent(indices(nbSmall - tmp1 - 1));
        OneEM();
        loglikeSmall(indices(nbSmall - tmp1 - 1)) = ComputeLogLike();
        degenere = FiltreDegenere();
        if (degenere==1){
          m_nbdegenere ++;
          loglikeSmall(indices(nbSmall - tmp1 - 1)) = -99999999999999999;
        }
      }
      if (any(loglikeSmall != loglikeSmall))
      loglikeSmall(find((loglikeSmall != loglikeSmall))) = -99999999999999999*ones<vec>(sum(loglikeSmall != loglikeSmall));
               
      uword  index;
      double indicebest = (loglikeSmall).max(index);
      SwitchParamCurrent(index);
      loglikeoutput = ComputeLogLike();
      
      
      degenere = FiltreDegenere();
      if (degenere==1){
        m_nbdegenere ++;
        loglikeoutput = -99999999999999999;
      }
      indices = sort_index(loglikeSmall);
    }
  }
}

colvec XEM::FindZMAP(){
  Col<double> zMAP=ones<vec>(tmplogproba.n_rows);
  uword  index;
  double max_val=0;
  for (int i=0; i<tmplogproba.n_rows; i++){
    max_val = (tmplogproba.row(i)).max(index);
    zMAP(i)=index;
  }
  return zMAP;
}

double XEM::ComputeLogLike(){
  ComputeTmpLogProba();
  maxtmplogproba = max(tmplogproba, 1);
  double output=0;
  if (min(maxtmplogproba) == 0){
    output = log(0);
  }else{
    for (int k=0; k<g; k++) tmplogproba.col(k)-=maxtmplogproba;
    tmplogproba = exp(tmplogproba);
    rowsums = sum(tmplogproba,1);
    output = sum(maxtmplogproba) + sum(log(rowsums));
  }
  return output;
}


void XEM::Estep(){for (int k=0; k<g; k++) tmplogproba.col(k) = tmplogproba.col(k)/rowsums;}

void XEM::OneEM(){
  double loglike = ComputeLogLike(), prec = log(0);
  int it=0;
  while ( (it<iterCurrent) && ((loglike-prec)>tolKeep) ){
    it ++;
    Estep();
    Mstep();
    prec = loglike;
    loglike = ComputeLogLike();
  }
  // Une verif
  //if (prec>(loglike+tolKeep)) cout << "pb EM " << endl;
}
