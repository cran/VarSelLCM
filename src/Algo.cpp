#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h> 

using namespace std;
using namespace arma;
using namespace Rcpp;
#include "Algo.h"

Algo::Algo(const S4 &  obj){
  this->x = as<mat>(obj.slot("data"));
  S4 model = obj.slot("model");
  this->omega = as<vec>(model.slot("omega"));
  this->g = model.slot("g");
  S4 partitions = obj.slot("partitions");
  this->z  = as<vec>(partitions.slot("zOPT")); 
  this->micl = log(0);
  this->priors = as<mat>(obj.slot("priors"));
  micl = Integre_Complete_Like(z);
}




double logdinvgamma(const double & x, const double & alpha, const double  & beta){
  return alpha*log(beta) - lgamma(alpha) - (alpha+1)*log(x) - beta/x ;
}


double IntegreOneVariable(const vec & v, const double & nu, const double & s0, const double & mu,  const double & n0){
  double output = 0;
  double n = v.n_rows;
  if (n> 0){ 
    double n1 = n + n0;
    double theta1 = (n0*mu + n*mean(v))/n1;
    double s1 = sqrt( s0*s0 + var(v) *(n-1)  + pow((mu - mean(v)),2) /(1/n0 + 1/n)  );
    output =  -log(sqrt( M_PI))*n + lgamma((n + nu)*0.5) - lgamma(nu*0.5) +   nu * log(s0/s1) - n*log(s1) + log(sqrt(n0 / n1) );
  }
  return output;
  
}


double Algo::Integre_Complete_Like(vec partition){
  double outmicl = lgamma(g*0.5) - g*lgamma(0.5) - lgamma(x.n_rows + g*0.5);
  for (int k=0; k<g; k++){
    outmicl += lgamma(sum(partition==k) + 0.5);
  }
  for (int j=0; j<x.n_cols; j++){
    if (omega(j)==0){
      outmicl +=  IntegreOneVariable(x.col(j), priors(j,0), priors(j,1), priors(j,2), priors(j,3));
    }else{
      vec tmp = x.col(j);
      for (int k=0; k<g; k++){
        vec sup=(tmp(find(partition == k)));
         // cout << sup.n_rows << endl;
         outmicl +=  IntegreOneVariable(tmp(find(partition == k)), priors(j,0), priors(j,1), priors(j,2), priors(j,3));
        
      }
    }
  }
  return outmicl;
}


void Algo::Optimize_partition(){
  int chgt = x.n_rows ;
  vec z_cand=z;
  double critere_cand;
  while (chgt > 0 ){
    ivec who = randi<ivec>(chgt, distr_param(0, x.n_rows -1));
    chgt = 0;

    for (int it=0; it < who.n_rows; it++){ 
      for (int k=0; k<g; k++){
          z_cand( who(it) ) = k;
          critere_cand = Integre_Complete_Like( z_cand );
          if (critere_cand > micl){
            z(who(it)) = k;
            micl = critere_cand;
            chgt = it ;
        }
      }
      z_cand( who(it) ) = z( who(it) );
    }
  }

}


void Algo::Optimize_model(){
  micl = lgamma(g*0.5) - g*lgamma(0.5) - lgamma(x.n_rows + g*0.5);
  for (int k=0; k<g; k++){
    micl += lgamma(sum(z==k) + 0.5);
  }
  double nondiscrim=0;
  double discrim=0;
  for (int j=0; j<x.n_cols; j++){
    nondiscrim = IntegreOneVariable(x.col(j), priors(j,0), priors(j,1), priors(j,2), priors(j,3));
    discrim = 0;
    vec tmp = x.col(j);
    for (int k=0; k<g; k++){
      discrim += IntegreOneVariable(tmp(find(z == k)), priors(j,0), priors(j,1), priors(j,2), priors(j,3));
    }
    if (discrim > nondiscrim){
      omega(j) = 1;
      micl +=  discrim ;
    }else{
      omega(j) = 0 ;
      micl +=  nondiscrim ;
    }      
  }
}


void Algo::Run(){
  double prec = log(0);
  while (prec < micl){
    prec = micl;
    Optimize_model();
    if (sum(omega)>0){
      Optimize_partition();
    }
    Optimize_model();
  }
}
