#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(BH)]]
#define BOOST_DISABLE_ASSERTS
//#include "header.h"
#include "header.h"
Rcpp::List GDHingeCxx(mat theta, mat X,mat Y,int K,double lambda,double v,bool ls, int B, int fam, double eps)
{
	  Hinge *hi=new Hinge;
    int p=X.n_cols;
    int n=X.n_rows;
    mat m;
    double s;
    double bound;
    if(fam==0)
    {
      GradientDescentF0<Hinge> GD(K,hi);
        mat mres=GD(theta(span(0,p-1),0),X,Y,v,lambda,ls,B);
        m=mres(0,span(0,p-1));
        if(1-v*lambda*lambda/(4*n)>0){
        bound=1.0/lambda*(GD.Get_F()+lambda*lambda/(4*n)-0.5*log(1-v*lambda*lambda/(4*n))-log(eps));
        }else{
        bound=-1;
        }
        return Rcpp::List::create(Rcpp::Named("m") = m,Rcpp::Named("bound")=bound);
    }else{
      GradientDescent<Hinge> GD(K,hi);
        mat mres=GD(theta(span(0,p-1),0),theta(p,0),X,Y,v,lambda,ls,B);
        m=mres(0,span(0,p-1));
        s=mres(0,p);
        if(1-v*lambda*lambda/(4*n)>0){
            bound=1.0/lambda*(GD.Get_F()+lambda*lambda/(4*n)-0.5*log(1-v*lambda*lambda/(4*n))-log(eps));
        }else{
            bound=-1;
        }
        return Rcpp::List::create(Rcpp::Named("m") = m,
                              Rcpp::Named("s") = s,
                              Rcpp::Named("bound") = bound
                              );
    }
}
RCPP_MODULE(GDHingeCxx) {
  Rcpp::function( "GDHingeCxx", &GDHingeCxx);
}

Rcpp::List GDHingeAUCCxx(mat theta, mat X,mat Y,int K,double lambda,double v,bool ls, int B, int fam,double eps)
{
	HingeAUC *hi=new HingeAUC(Y);
    int p=X.n_cols;
    int n=X.n_rows;
    mat m;
    double s;
    double bound;
    if(fam==0)
    {
      GradientDescentF0<HingeAUC> GD(K,hi);
        mat mres=GD(theta(span(0,p-1),0),X,Y,v,lambda,ls,B);
        m=mres(0,span(0,p-1));
        if(1-v*lambda*lambda/(2*(n-1))>0){
        bound=1.0/lambda*(GD.Get_F()+lambda*lambda/(2*(n-1))-0.5*log(1-v*lambda*lambda/(2*(n-1)))-log(eps));
        }else{
        bound=-1;
        }
        return Rcpp::List::create(Rcpp::Named("m") = m,Rcpp::Named("bound")=bound);
    }else{
      GradientDescent<HingeAUC> GD(K,hi);
        mat mres=GD(theta(span(0,p-1),0),theta(p,0),X,Y,v,lambda,ls,B);
        m=mres(0,span(0,p-1));
        s=mres(0,p);
        if(1-v*lambda*lambda/(2*(n-1))>0){
        bound=1.0/lambda*(GD.Get_F()+lambda*lambda/(4*n)-0.5*log(1-v*lambda*lambda/(4*n))-log(eps));
        }else{
            bound=-1;
        }
        return Rcpp::List::create(Rcpp::Named("m") = m,
                              Rcpp::Named("s") = s,
                              Rcpp::Named("bound") = bound);
    }
}
  RCPP_MODULE(GDHingeAUCCxx) {
  Rcpp::function( "GDHingeAUCCxx", &GDHingeAUCCxx);
  }
