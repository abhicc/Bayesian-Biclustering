#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix arr(IntegerMatrix m, int nreps, Function f){
            NumericMatrix mat(nreps+1,nreps+1);
            for(int i=0; i<mat.nrow()-1; i++){
            for(int j=(i+1); j<mat.ncol(); j++){
            IntegerVector v1 = m(_,i), v2 = m(_,j);
            List out = f(v1,v2);
            mat(i,j)=out["Rand"];
            }
            }
            return mat;
            }
// [[Rcpp::export]]
NumericMatrix acr(IntegerMatrix m, int nreps, Function f){
            NumericMatrix mat(nreps+1,nreps+1);
            for(int i=0; i<mat.nrow()-1; i++){
            for(int j=(i+1); j<mat.ncol(); j++){
            IntegerVector v1 = m(i,_), v2 = m(j,_);
            List out = f(v1,v2);
            mat(i,j)=out["Rand"];
            }
            }
            return mat;
            }			
// [[Rcpp::export]]
NumericMatrix aer(IntegerMatrix m, int nreps, Function f){
            NumericMatrix mat(nreps+1,nreps+1);
            for(int i=0; i<mat.nrow()-1; i++){
            for(int j=(i+1); j<mat.ncol(); j++){
            IntegerVector v1 = m(_,i), v2 = m(_,j);
            List out = f(v1,v2);
            mat(i,j)=out["Rand"];
            }
            }
            return mat;
            }
