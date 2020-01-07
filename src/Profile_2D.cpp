#include <cmath>
#include <algorithm>  
#include <vector>
#include <limits>     
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
  double ft_y(double y, double mung, double mutg, double sng, double stg, double pi1, double pi2)
{
    int i;
    double rtval;
    double tmp;
	int integ=200;
    double tmp_pos[integ];

    rtval=0;
    for(i=0;i<integ;i++) tmp_pos[i] = (y)/(double)integ*(i+0.5);

    for(i=0;i<integ;i++)
    {
      tmp= -log(sng)-log(stg) -log(y-tmp_pos[i])-log(tmp_pos[i]);
      tmp= tmp-0.5*pow(log2(tmp_pos[i])-(log2(pi1) + mung), 2.0) /sng/sng;
      tmp= tmp-0.5*pow(log2(y - tmp_pos[i])-(log2(1-pi1-pi2) + mutg), 2.0)/stg/stg;

      tmp= exp(tmp);

      rtval +=tmp;
    }

    if(rtval<=0)
      rtval=1e-313;

    return (log(rtval/(double)integ *y));
  }

// [[Rcpp::export]]
double pf_y(NumericMatrix y, int samp, NumericVector mung, NumericVector mutg, NumericVector sng, NumericVector stg, double pi1)
{
   int nG=y.nrow();
  double new_like;
  double tmp;
  int i;
  new_like=0.0;
    for (i=0; i<nG; i++){
       tmp = ft_y(y(i,samp),mung[i], mutg[i], sng[i], stg[i], pi1, 0.0);
       new_like+= tmp;
    }
  return(new_like*(-1.0));
}

// [[Rcpp::export]]
double pmin_y(double ax, double bx, int samp, NumericMatrix y, NumericVector mung, NumericVector mutg, NumericVector sng, NumericVector stg, double tol)
{
  /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

  /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
  double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

  /*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
  tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

  a = ax;
  b = bx;
  v = a + c * (b - a);
  w = v;
  x = v;

  d = 0.;/* -Wall */
    e = 0.;
  fx = pf_y(y, samp, mung, mutg, sng, stg, x);
  fv = fx;
  fw = fx;
  tol3 = tol / 3.;

  /*  main loop starts here ----------------------------------- */

    for(;;) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;

      /* check stopping criterion */

        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (fabs(e) > tol1) { /* fit parabola */

                              r = (x - w) * (fx - fv);
                            q = (x - v) * (fx - fw);
                            p = (x - v) * q - (x - w) * r;
                            q = (q - r) * 2.;
                            if (q > 0.) p = -p; else q = -q;
                            r = e;
                            e = d;
      }

      if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

                                                      if (x < xm) e = b - x; else e = a - x;
                                                    d = c * e;
      }
      else { /* a parabolic-interpolation step */

               d = p / q;
             u = x + d;

             /* f must not be evaluated too close to ax or bx */

               if (u - a < t2 || b - u < t2) {
                 d = tol1;
                 if (x >= xm) d = -d;
               }
      }

      /* f must not be evaluated too close to x */

        if (fabs(d) >= tol1)
          u = x + d;
      else if (d > 0.)
        u = x + tol1;
      else
        u = x - tol1;

      fu = pf_y(y, samp, mung, mutg, sng, stg, u);

      /*  update  a, b, v, w, and x */

        if (fu <= fx) {
          if (u < x) b = x; else a = x;
          v = w;    w = x;   x = u;
          fv = fw; fw = fx; fx = fu;
        } else {
          if (u < x) a = u; else b = u;
          if (fu <= fw || w == x) {
            v = w; fv = fw;
            w = u; fw = fu;
          } else if (fu <= fv || v == x || v == w) {
            v = u; fv = fu;
          }
        }
    }
  /* end of main loop */

    return x;

}


/* Tumor Component Update */
/*   void gettumor(int genes, int h) // option h = 1 for one component; h = 2 for two
{
    double mu, sigma,upp;
    if(opt_seq == 1){
        upp = 18;
    }else{
        upp = 1500;
    }
    
  if(h == 1)
  {
    mu = tmin_y(1, upp, genes,h, mint, 0.001);
    
    sigma =  tmin_y2(0.0001, 1, genes, mu, tf_y, 0.001);
  }
  p->Tavg[genes] = mu;
  p->Tsigma[genes] = sigma;

} */

// [[Rcpp::export]]
double tf_y(int genes, double mu, double sigma, NumericMatrix y, NumericVector pi1, double mung, double sng)
{
	int intx=y.ncol();
  double new_like;
  double tmp;
  int i;
  new_like=0.0;
     for (i=0; i<intx; i++){
         tmp = ft_y(y(genes,i), mung, mu, sng, sigma, pi1[i], 0);
         new_like += tmp;
     }
  return (new_like*(-1.0));
}

// [[Rcpp::export]]
double tmin_y2(double ax, double bx, int genes, double mu, NumericMatrix y, NumericVector pi1, double mung, double sng, double tol)
{
  /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

  /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
  double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

  /*  eps is approximately the square root of the relative machine precision. */
    eps = std::numeric_limits<double>::epsilon();
  tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

  a = ax;
  b = bx;
  v = a + c * (b - a);
  w = v;
  x = v;

  d = 0.;/* -Wall */
    e = 0.;
  fx = tf_y(genes, mu, x, y, pi1, mung, sng);
  fv = fx;
  fw = fx;
  tol3 = tol / 3.;

  /*  main loop starts here ----------------------------------- */

    for(;;) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;

      /* check stopping criterion */

        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (fabs(e) > tol1) { /* fit parabola */

                              r = (x - w) * (fx - fv);
                            q = (x - v) * (fx - fw);
                            p = (x - v) * q - (x - w) * r;
                            q = (q - r) * 2.;
                            if (q > 0.) p = -p; else q = -q;
                            r = e;
                            e = d;
      }

      if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

                                                      if (x < xm) e = b - x; else e = a - x;
                                                    d = c * e;
      }
      else { /* a parabolic-interpolation step */

               d = p / q;
             u = x + d;

             /* f must not be evaluated too close to ax or bx */

               if (u - a < t2 || b - u < t2) {
                 d = tol1;
                 if (x >= xm) d = -d;
               }
      }

      /* f must not be evaluated too close to x */

        if (fabs(d) >= tol1)
          u = x + d;
      else if (d > 0.)
        u = x + tol1;
      else
        u = x - tol1;

      fu = tf_y(genes, mu, u, y, pi1, mung, sng);
	  //printf("%f\n",fu);

      /*  update  a, b, v, w, and x */

        if (fu <= fx) {
          if (u < x) b = x; else a = x;
          v = w;    w = x;   x = u;
          fv = fw; fw = fx; fx = fu;
        } else {
          if (u < x) a = u; else b = u;
          if (fu <= fw || w == x) {
            v = w; fv = fw;
            w = u; fw = fu;
          } else if (fu <= fv || v == x || v == w) {
            v = u; fv = fu;
          }
        }
    }
  /* end of main loop */

    return x;

}

/* Get tumor */
// [[Rcpp::export]]
double mint(int genes, double mu, NumericMatrix y, NumericVector pi1, double mung, double sng) // opt h =1:one component; h =2:two component
{
  double sigma;
  double tmp;
  
    sigma =  tmin_y2(0.0001, 1, genes, mu, y, pi1, mung, sng, 0.001);
    tmp = tf_y(genes, mu, sigma, y, pi1, mung, sng);
  return (tmp);
}

// [[Rcpp::export]]
double tmin_y(double ax, double bx, int genes, NumericMatrix y, NumericVector pi1, double mung, double sng, double tol)
{
  /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

  /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
  double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

  /*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
  tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

  a = ax;
  b = bx;
  v = a + c * (b - a);
  w = v;
  x = v;

  d = 0.;/* -Wall */
    e = 0.;
  fx = mint(genes, x, y, pi1, mung, sng);
  fv = fx;
  fw = fx;
  tol3 = tol / 3.;

  /*  main loop starts here ----------------------------------- */

    for(;;) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;

      /* check stopping criterion */

        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (fabs(e) > tol1) { /* fit parabola */

                              r = (x - w) * (fx - fv);
                            q = (x - v) * (fx - fw);
                            p = (x - v) * q - (x - w) * r;
                            q = (q - r) * 2.;
                            if (q > 0.) p = -p; else q = -q;
                            r = e;
                            e = d;
      }

      if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

                                                      if (x < xm) e = b - x; else e = a - x;
                                                    d = c * e;
      }
      else { /* a parabolic-interpolation step */

               d = p / q;
             u = x + d;

             /* f must not be evaluated too close to ax or bx */

               if (u - a < t2 || b - u < t2) {
                 d = tol1;
                 if (x >= xm) d = -d;
               }
      }

      /* f must not be evaluated too close to x */

        if (fabs(d) >= tol1)
          u = x + d;
      else if (d > 0.)
        u = x + tol1;
      else
        u = x - tol1;

      fu = mint(genes, u, y, pi1, mung, sng);

      /*  update  a, b, v, w, and x */

        if (fu <= fx) {
          if (u < x) b = x; else a = x;
          v = w;    w = x;   x = u;
          fv = fw; fw = fx; fx = fu;
        } else {
          if (u < x) a = u; else b = u;
          if (fu <= fw || w == x) {
            v = w; fv = fw;
            w = u; fw = fu;
          } else if (fu <= fv || v == x || v == w) {
            v = u; fv = fu;
          }
        }
    }
  /* end of main loop */
    return x;
}
  
// [[Rcpp::export]]
double Loglikelihood_ft_y(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   double L=0.0;
   int S=y.ncol();
   int G=y.nrow();
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   L=L-ft_y(y(i,j), MuN[i], MuT[i], SigmaN[i], SigmaT[i], Pi[j], 0.0);
	   }
   }
   return(L);
}

// [[Rcpp::export]]
double SoftThreshold(double x, double t) {
	double y=0;
	if (x>t){
		y=x-t;
	}
    if (x<(-t)){
		y=x+t;
	}
	return(y);
}

// [[Rcpp::export]]
NumericVector SoftThreshold_vec(NumericVector x, double t) {
	int n=x.size();
	NumericVector y(n);
	for(int i=0; i<n; i++){
		if (x[i]>t){
			y[i]=x[i]-t;
		}
		if (x[i]<(-t)){
			y[i]=x[i]+t;
		}
	}
	return(y);
}

// [[Rcpp::export]]
NumericVector Gt(NumericVector x, NumericVector Dg, double t) {
	int n=x.size();
	NumericVector y(n);
	for(int i=0; i < n; i++){
		y[i] = ( x[i] - SoftThreshold(x[i]-t*Dg[i], t) )/t;
	}
	return(y);
}

// [[Rcpp::export]]
NumericVector Gt_vec(NumericVector x, NumericVector Dg, double t) {
	int n=x.size();
	NumericVector y(n);
	y=1.0/t*(x - SoftThreshold_vec(x-t*Dg, t) );
	return(y);
}

// [[Rcpp::export]]
double Alpha_search_2D(NumericVector x, NumericVector p, int S, int G) {
   double alpha_max=0.0;
   double alpha=1.0;
   for(int j=0; j < S; j++){
	   if(p[j]<0){
		   alpha_max=-x[j]/p[j];
	   }else{
		   alpha_max=(1.0-x[j])/p[j];
	   }
	   if(alpha>alpha_max){
		   alpha=alpha_max;
	   }
   }
   for(int i=0; i < 2*G; i++){
	   if(p[S+i]<0){
		   alpha_max=-x[S+i]/p[S+i];
		   if(alpha>alpha_max){
			   alpha=alpha_max;
		   }
	   }
   }
   return(alpha);
}

// [[Rcpp::export]]
double Alpha_search_Pi_2D(NumericVector x, NumericVector p, int S) {
   double alpha_max=0.0;
   double alpha=1.0;
   for(int j=0; j < S; j++){
	   if(p[j]<0){
		   alpha_max=-x[j]/p[j];
	   }else{
		   alpha_max=(1.0-x[j])/p[j];
	   }
	   if(alpha>alpha_max){
		   alpha=alpha_max;
	   }
   }
   return(alpha);
}

// [[Rcpp::export]]
double Alpha_search_MuT_2D(NumericVector x, NumericVector p, int G) {
   double alpha_max=0.0;
   double alpha=1.0;
   for(int i=0; i < G; i++){
	   if(p[i]<0){
		   alpha_max=-x[i]/p[i];
		   if(alpha>alpha_max){
			   alpha=alpha_max;
		   }
	   }
   }
   return(alpha);
}

// [[Rcpp::export]]
double Alpha_search_SigmaT_2D(NumericVector x, NumericVector p, int G) {
	double alpha_max=0.0;
	double alpha=1.0;
	for(int i=0; i < G; i++){
		if(p[i]<0){
			alpha_max=-x[i]/p[i];
		}else{
			alpha_max=(1.0-x[i])/p[i];
		}
		if(alpha>alpha_max){
			alpha=alpha_max;
		}
	}
	return(alpha);
}

// [[Rcpp::export]]
double Alpha_search_2D_sigma(NumericVector x, NumericVector p, int S, int G) {
   double alpha_max=0.0;
   double alpha=100.0;
   for(int j=0; j < S; j++){
	   if(p[j]<0){
		   alpha_max=-x[j]/p[j];
	   }else{
		   alpha_max=(1.0-x[j])/p[j];
	   }
	   if(alpha>alpha_max){
		   alpha=alpha_max;
	   }
   }
   for(int i=0; i < G; i++){
	   if(p[S+2*i]<0){
		   alpha_max=-x[S+2*i]/p[S+2*i];
		   if(alpha>alpha_max){
			   alpha=alpha_max;
		   }
	   }
	   if(p[S+2*i+1]<0){
		   alpha_max=-x[S+2*i+1]/p[S+2*i+1];
	   }else{
		   alpha_max=(1.0-x[S+2*i+1])/p[S+2*i+1];
	   }
	   if(alpha>alpha_max){
		   alpha=alpha_max;
	   }
   }
   return(alpha);
}


/* 1st order derivative */
// [[Rcpp::export]]
double Unit1(double t, double Mu, double Pit) {
   return( log2(t)-Mu-log2(Pit) ) ;
}

// [[Rcpp::export]]
NumericVector x_update_2D(NumericVector Pi, NumericVector MuT, NumericVector SigmaT, int S, int G){
	NumericVector x(2*G+S);
	for (int j = 0; j < S; j++){
	   x[j]=Pi[j];
   }
   for (int i=0; i < G; i++){
	   x[S+2*i]=MuT[i];
	   x[S+2*i+1]=SigmaT[i];
   }
	return(x);
}

// [[Rcpp::export]]
List x_update_inv_2D(NumericVector x, int S, int G){
	NumericVector Pi(S);
	NumericVector MuT(G);
	NumericVector SigmaT(G);
	for (int j = 0; j < S; j++){
	   Pi[j]=x[j];
	}
	for (int i=0; i < G; i++){   
	   MuT[i]=x[S+2*i];
	   SigmaT[i]=x[S+2*i+1];
    }
	return(List::create(Named("Pi")=Pi,Named("MuT")=MuT,Named("SigmaT")=SigmaT));
}

/* Likelihood */
// [[Rcpp::export]]
double innerCPP_2D(double t, double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT) {
   return(exp( -log(t) - log(y-t) - 0.5*pow( (log2(t)-MuN-log2(Pi))/SigmaN ,2.0) - 0.5*pow( (log2(y-t)-MuT-log2(1.0-Pi))/SigmaT ,2.0) ) ) ;
}

// [[Rcpp::export]]
double inner_trapez_2D(double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT){
	double sum=0.0;
	int integ=200;
	double tmp_pos;
	for(int i=0;i<integ;i++){
		tmp_pos = y/(double)integ*(i+0.5);
		sum = sum + innerCPP_2D(tmp_pos, y, Pi, MuN, MuT, SigmaN, SigmaT);
	}
	sum = sum/SigmaN/SigmaT/(double)integ*y;
	if(sum<=0){
		sum=1e-313;
	}
	return(sum);
}

/* 1st Derivative */
// [[Rcpp::export]]
double DPi_inner_2D(double t, double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT) {
	double Unit1_N1=Unit1(t,MuN,Pi);
	double Unit1_T=Unit1(y-t,MuT,1.0-Pi);
    return( (Unit1_N1/SigmaN/SigmaN/Pi-Unit1_T/SigmaT/SigmaT/(1.0-Pi)) * exp( -log(t) - log(y-t) - 0.5*pow( Unit1_N1/SigmaN ,2.0) - 0.5*pow( Unit1_T/SigmaT ,2.0) ) ) ;
}

// [[Rcpp::export]]
double DPi_outer_2D(double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT){
	double sum=0.0;
	int integ=200;
	double lower=0.0;
	double upper=y;
	double h=(upper-lower)/double(integ);
	for (int i=1;i<=(integ-1);i++){
		sum = sum + DPi_inner_2D(lower+double(i)*h, y, Pi, MuN, MuT, SigmaN, SigmaT);
	}
	sum = h*sum;
	sum = sum/SigmaN/SigmaT/log(2.0);
	return(sum);
}

// [[Rcpp::export]]
NumericMatrix D1f0Pi_func_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes 
   NumericMatrix D1f0Pi1(G,S);
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   D1f0Pi1(i,j)=DPi_outer_2D(y(i,j), Pi[j], MuN[i], MuT[i], SigmaN[i], SigmaT[i]);
	   }
   }
   return(D1f0Pi1);
}

// [[Rcpp::export]]
double DMuT_inner_2D(double t, double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT) {
	double Unit1_N1=Unit1(t,MuN,Pi);
	double Unit1_T=Unit1(y-t,MuT,1.0-Pi);
    return( Unit1_T * exp( -log(t) - log(y-t) - 0.5*pow( Unit1_N1/SigmaN ,2.0) - 0.5*pow( Unit1_T/SigmaT ,2.0) ) ) ;
}

// [[Rcpp::export]]
double DMuT_outer_2D(double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT){
	double sum=0.0;
	int integ=200;
	double lower=0.0;
	double upper=y;
	double h=(upper-lower)/double(integ);
	for (int i=1;i<=(integ-1);i++){
		sum = sum + DMuT_inner_2D(lower+double(i)*h, y, Pi, MuN, MuT, SigmaN, SigmaT);
	}
	sum = h*sum;
	sum = sum/SigmaN/pow(SigmaT,3.0);
	return(sum);
}

// [[Rcpp::export]]
NumericMatrix D1f0MuT_func_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes 
   NumericMatrix D1f0Pi1(G,S);
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   D1f0Pi1(i,j)=DMuT_outer_2D(y(i,j), Pi[j], MuN[i], MuT[i], SigmaN[i], SigmaT[i]);
	   }
   }
   return(D1f0Pi1);
}

// [[Rcpp::export]]
double DSigmaT_inner_2D(double t, double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT) {
	double Unit1_N1=Unit1(t,MuN,Pi);
	double Unit1_T=Unit1(y-t,MuT,1.0-Pi);
    return( (pow(Unit1_T/SigmaT,2.0)-1.0) * exp( -log(t) - log(y-t) - 0.5*pow( Unit1_N1/SigmaN ,2.0) - 0.5*pow( Unit1_T/SigmaT ,2.0) ) ) ;
}

// [[Rcpp::export]]
double DSigmaT_outer_2D(double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT){
	double sum=0.0;
	int integ=200;
	double lower=0.0;
	double upper=y;
	double h=(upper-lower)/double(integ);
	for (int i=1;i<=(integ-1);i++){
		sum = sum + DSigmaT_inner_2D(lower+double(i)*h, y, Pi, MuN, MuT, SigmaN, SigmaT);
	}
	sum = h*sum;
	sum = sum/SigmaN/pow(SigmaT,2.0);
	return(sum);
}

// [[Rcpp::export]]
NumericMatrix D1f0SigmaT_func_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes 
   NumericMatrix D1f0Pi1(G,S);
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   D1f0Pi1(i,j)=DSigmaT_outer_2D(y(i,j), Pi[j], MuN[i], MuT[i], SigmaN[i], SigmaT[i]);
	   }
   }
   return(D1f0Pi1);
}

/* 2nd Order Derivative */

// [[Rcpp::export]]
double D2Pi_inner_2D(double t, double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT, double ln2) {
	double Unit1_N=Unit1(t,MuN,Pi);
	double Unit1_T=Unit1(y-t,MuT,1.0-Pi);
    return( (pow((Unit1_N/SigmaN/SigmaN/Pi-Unit1_T/SigmaT/SigmaT/(1.0-Pi)),2.0) - (ln2+Unit1_N)/pow(Pi*SigmaN,2.0) - (ln2+Unit1_T)/pow((1.0-Pi)*SigmaT,2.0) ) * exp( -log(t) - log(y-t) - 0.5*pow( Unit1_N/SigmaN ,2.0) - 0.5*pow( Unit1_T/SigmaT ,2.0) ) ) ;
}

// [[Rcpp::export]]
double D2Pi_outer_2D(double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT){
	double sum=0.0;
	int integ=200;
	double ln2=1.0/log(2.0);
	double lower=0.0;
	double upper=y;
	double h=(upper-lower)/double(integ);
	for (int i=1;i<=(integ-1);i++){
		sum = sum + D2Pi_inner_2D(lower+double(i)*h, y, Pi, MuN, MuT, SigmaN, SigmaT, ln2);
	}
	sum = h*sum;
	sum = sum/SigmaN/SigmaT*ln2;
	return(sum);
}

// [[Rcpp::export]]
NumericMatrix D2f0Pi_func_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes 
   NumericMatrix D1f0Pi1(G,S);
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   D1f0Pi1(i,j)=D2Pi_outer_2D(y(i,j), Pi[j], MuN[i], MuT[i], SigmaN[i], SigmaT[i]);
	   }
   }
   return(D1f0Pi1);
}

// [[Rcpp::export]]
double D2PiMuT_inner_2D(double t, double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT) {
	double Unit1_N=Unit1(t,MuN,Pi);
	double Unit1_T=Unit1(y-t,MuT,1.0-Pi);
    return( (1.0/(1.0-Pi) + Unit1_T*(Unit1_N/SigmaN/SigmaN/Pi-Unit1_T/SigmaT/SigmaT/(1.0-Pi)) ) * exp( -log(t) - log(y-t) - 0.5*pow( Unit1_N/SigmaN ,2.0) - 0.5*pow( Unit1_T/SigmaT ,2.0) ) ) ;
}

// [[Rcpp::export]]
double D2PiMuT_outer_2D(double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT){
	double sum=0.0;
	int integ=200;
	double lower=0.0;
	double upper=y;
	double h=(upper-lower)/double(integ);
	for (int i=1;i<=(integ-1);i++){
		sum = sum + D2PiMuT_inner_2D(lower+double(i)*h, y, Pi, MuN, MuT, SigmaN, SigmaT);
	}
	sum = h*sum;
	sum = sum/SigmaN/pow(SigmaT,3.0)/log(2.0);
	return(sum);
}

// [[Rcpp::export]]
NumericMatrix D2f0PiMuT_func_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes 
   NumericMatrix D1f0Pi1(G,S);
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   D1f0Pi1(i,j)=D2PiMuT_outer_2D(y(i,j), Pi[j], MuN[i], MuT[i], SigmaN[i], SigmaT[i]);
	   }
   }
   return(D1f0Pi1);
}

// [[Rcpp::export]]
double D2PiSigmaT_inner_2D(double t, double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT) {
	double Unit1_N=Unit1(t,MuN,Pi);
	double Unit1_T=Unit1(y-t,MuT,1.0-Pi);
    return( ((Unit1_N/SigmaN/SigmaN/Pi-Unit1_T/SigmaT/SigmaT/(1.0-Pi))*(pow(Unit1_T/SigmaT,2.0)-1.0) + 2.0*Unit1_T/(1.0-Pi)/SigmaT/SigmaT ) * exp( -log(t) - log(y-t) - 0.5*pow( Unit1_N/SigmaN ,2.0) - 0.5*pow( Unit1_T/SigmaT ,2.0) ) ) ;
}

// [[Rcpp::export]]
double D2PiSigmaT_outer_2D(double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT){
	double sum=0.0;
	int integ=200;
	double lower=0.0;
	double upper=y;
	double h=(upper-lower)/double(integ);
	for (int i=1;i<=(integ-1);i++){
		sum = sum + D2PiSigmaT_inner_2D(lower+double(i)*h, y, Pi, MuN, MuT, SigmaN, SigmaT);
	}
	sum = h*sum;
	sum = sum/SigmaN/pow(SigmaT,2.0)/log(2.0);
	return(sum);
}

// [[Rcpp::export]]
NumericMatrix D2f0PiSigmaT_func_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes 
   NumericMatrix D1f0Pi1(G,S);
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   D1f0Pi1(i,j)=D2PiSigmaT_outer_2D(y(i,j), Pi[j], MuN[i], MuT[i], SigmaN[i], SigmaT[i]);
	   }
   }
   return(D1f0Pi1);
}

// [[Rcpp::export]]
double D2MuT_inner_2D(double t, double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT) {
	double Unit1_N=Unit1(t,MuN,Pi);
	double Unit1_T=Unit1(y-t,MuT,1.0-Pi);
    return( (pow(Unit1_T/SigmaT,2.0) - 1.0) * exp( -log(t) - log(y-t) - 0.5*pow( Unit1_N/SigmaN ,2.0) - 0.5*pow( Unit1_T/SigmaT ,2.0) ) ) ;
}

// [[Rcpp::export]]
double D2MuT_outer_2D(double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT){
	double sum=0.0;
	int integ=200;
	double lower=0.0;
	double upper=y;
	double h=(upper-lower)/double(integ);
	for (int i=1;i<=(integ-1);i++){
		sum = sum + D2MuT_inner_2D(lower+double(i)*h, y, Pi, MuN, MuT, SigmaN, SigmaT);
	}
	sum = h*sum;
	sum = sum/SigmaN/pow(SigmaT,3.0);
	return(sum);
}

// [[Rcpp::export]]
NumericMatrix D2f0MuT_func_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes 
   NumericMatrix D1f0Pi1(G,S);
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   D1f0Pi1(i,j)=D2MuT_outer_2D(y(i,j), Pi[j], MuN[i], MuT[i], SigmaN[i], SigmaT[i]);
	   }
   }
   return(D1f0Pi1);
}

// [[Rcpp::export]]
double D2MuTSigmaT_inner_2D(double t, double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT) {
	double Unit1_N=Unit1(t,MuN,Pi);
	double Unit1_T=Unit1(y-t,MuT,1.0-Pi);
    return( Unit1_T*(pow(Unit1_T/SigmaT,2.0) - 3.0) * exp( -log(t) - log(y-t) - 0.5*pow( Unit1_N/SigmaN ,2.0) - 0.5*pow( Unit1_T/SigmaT ,2.0) ) ) ;
}

// [[Rcpp::export]]
double D2MuTSigmaT_outer_2D(double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT){
	double sum=0.0;
	int integ=200;
	double lower=0.0;
	double upper=y;
	double h=(upper-lower)/double(integ);
	for (int i=1;i<=(integ-1);i++){
		sum = sum + D2MuTSigmaT_inner_2D(lower+double(i)*h, y, Pi, MuN, MuT, SigmaN, SigmaT);
	}
	sum = h*sum;
	sum = sum/SigmaN/pow(SigmaT,4.0);
	return(sum);
}

// [[Rcpp::export]]
NumericMatrix D2f0MuTSigmaT_func_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes 
   NumericMatrix D1f0Pi1(G,S);
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   D1f0Pi1(i,j)=D2MuTSigmaT_outer_2D(y(i,j), Pi[j], MuN[i], MuT[i], SigmaN[i], SigmaT[i]);
	   }
   }
   return(D1f0Pi1);
}

// [[Rcpp::export]]
double D2SigmaT_inner_2D(double t, double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT) {
	double Unit1_N=Unit1(t,MuN,Pi);
	double Unit1_T=Unit1(y-t,MuT,1.0-Pi);
	double powUnit=pow(Unit1_T/SigmaT,2.0);
    return( (powUnit*powUnit - 5.0*powUnit + 2.0) * exp( -log(t) - log(y-t) - 0.5*pow( Unit1_N/SigmaN ,2.0) - 0.5*pow( Unit1_T/SigmaT ,2.0) ) ) ;
}

// [[Rcpp::export]]
double D2SigmaT_outer_2D(double y, double Pi, double MuN, double MuT, double SigmaN, double SigmaT){
	double sum=0.0;
	int integ=200;
	double lower=0.0;
	double upper=y;
	double h=(upper-lower)/double(integ);
	for (int i=1;i<=(integ-1);i++){
		sum = sum + D2SigmaT_inner_2D(lower+double(i)*h, y, Pi, MuN, MuT, SigmaN, SigmaT);
	}
	sum = h*sum;
	sum = sum/SigmaN/pow(SigmaT,3.0);
	return(sum);
}

// [[Rcpp::export]]
NumericMatrix D2f0SigmaT_func_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes 
   NumericMatrix D1f0Pi1(G,S);
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   D1f0Pi1(i,j)=D2SigmaT_outer_2D(y(i,j), Pi[j], MuN[i], MuT[i], SigmaN[i], SigmaT[i]);
	   }
   }
   return(D1f0Pi1);
}

/* Likelihood */
// [[Rcpp::export]]
NumericMatrix f0_func_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes 
   NumericMatrix f0(G,S);
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   f0(i,j)=inner_trapez_2D(y(i,j), Pi[j], MuN[i], MuT[i], SigmaN[i], SigmaT[i]);
	   }
   }
   return(f0);
}

// [[Rcpp::export]]
double Loglikelihood_Pi_2D(NumericMatrix y, double Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT, int j)
{
	double L=0.0;
	int G=y.nrow();
	for (int i = 0; i < G; i++){
		L=L-log(inner_trapez_2D(y(i,j), Pi, MuN[i], MuT[i], SigmaN[i], SigmaT[i]));
	}
	return(L);
}

// [[Rcpp::export]]
double GoldenSection_Loglikelihood_Pi_2D(NumericMatrix y, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT, int j)
{
	double a=0.01;
	double b=0.99;
	double gr=(sqrt(5.0) + 1.0) / 2.0;
	double c = b - (b - a) / gr;
	double d = a + (b - a) / gr;
	int i=1;
	
	while(fabs(c - d) > 1e-5){
		//printf("%f\n",Gfunc_2D_C(c,x,p,S,G,y,MuN,SigmaN));
		if(Loglikelihood_Pi_2D(y,c,MuN,MuT,SigmaN,SigmaT,j) < Loglikelihood_Pi_2D(y,d,MuN,MuT,SigmaN,SigmaT,j)){
			b=d;
		}else{
			a=c;
		}
		c = b - (b - a) / gr;
		d = a + (b - a) / gr;
		i++;
	}
	return((c+d)/2);
}

// [[Rcpp::export]]
double Loglikelihood_MuT_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, double MuT, NumericVector SigmaN, NumericVector SigmaT, int i)
{
	double L=0.0;
	int S=y.ncol();
	for (int j=0; j < S; j++){
		L=L-log(inner_trapez_2D(y(i,j), Pi[j], MuN[i], MuT, SigmaN[i], SigmaT[i]));
	}
	return(L);
}

// [[Rcpp::export]]
double GoldenSection_Loglikelihood_MuT_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector SigmaN, NumericVector SigmaT, int j)
{
	double a=1;
	double b=18;
	double gr=(sqrt(5.0) + 1.0) / 2.0;
	double c = b - (b - a) / gr;
	double d = a + (b - a) / gr;
	int i=1;
	
	while(fabs(c - d) > 1e-5){
		//printf("%f\n",Gfunc_2D_C(c,x,p,S,G,y,MuN,SigmaN));
		if(Loglikelihood_MuT_2D(y,Pi,MuN,c,SigmaN,SigmaT,j) < Loglikelihood_MuT_2D(y,Pi,MuN,d,SigmaN,SigmaT,j)){
			b=d;
		}else{
			a=c;
		}
		c = b - (b - a) / gr;
		d = a + (b - a) / gr;
		i++;
	}
	return((c+d)/2);
}

// [[Rcpp::export]]
double Loglikelihood_SigmaT_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, double SigmaT, int i)
{
	double L=0.0;
	int S=y.ncol();
	for (int j=0; j < S; j++){
		L=L-log(inner_trapez_2D(y(i,j), Pi[j], MuN[i], MuT[i], SigmaN[i], SigmaT));
	}
	return(L);
}

// [[Rcpp::export]]
double GoldenSection_Loglikelihood_SigmaT_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, int j)
{
	double a=0.0001;
	double b=1.0;
	double gr=(sqrt(5.0) + 1.0) / 2.0;
	double c = b - (b - a) / gr;
	double d = a + (b - a) / gr;
	int i=1;
	
	while(fabs(c - d) > 1e-5){
		//printf("%f\n",Gfunc_2D_C(c,x,p,S,G,y,MuN,SigmaN));
		if(Loglikelihood_SigmaT_2D(y,Pi,MuN,MuT,SigmaN,c,j) < Loglikelihood_SigmaT_2D(y,Pi,MuN,MuT,SigmaN,d,j)){
			b=d;
		}else{
			a=c;
		}
		c = b - (b - a) / gr;
		d = a + (b - a) / gr;
		i++;
	}
	return((c+d)/2);
}

/* // [[Rcpp::export]]
double Opt_Loglikelihood_Pi_2D(Function opt, NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT)
{
	int S=y.ncol();
	List out;
	NumericVector bound(2);
	bound[0]=0.01;
	bound[1]=0.99;
	double Pi_out;
	for (int j=0; j < S; j++){
		out=opt(Loglikelihood_Pi_2D, bound, 0.001, y, MuN, MuT, SigmaN, SigmaT, j);
		Pi_out=out[0];
	}
	return(Pi_out);
}*/


// [[Rcpp::export]]
double Loglikelihood_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   double L=0.0;
   int S=y.ncol();
   int G=y.nrow();
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   L=L-log(inner_trapez_2D(y(i,j), Pi[j], MuN[i], MuT[i], SigmaN[i], SigmaT[i]));
	   }
   }
   return(L);
}

// [[Rcpp::export]]
double Loglikelihood_2D_L1(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector Beta, NumericVector SigmaN, NumericVector SigmaT, double lambda) {
   double L=0.0;
   double tmp=0.0;
   int S=y.ncol();
   int G=y.nrow();
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   L=L-log(inner_trapez_2D(y(i,j), Pi[j], MuN[i], MuN[i]+Beta[i], SigmaN[i], SigmaT[i]));
	   }
	   tmp=tmp+fabs(Beta[i]);
   }
   L=L+lambda*tmp;
   return(L);
}

// [[Rcpp::export]]
NumericVector D1Loglikelihood_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes
   NumericMatrix f0=f0_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0Pi=D1f0Pi_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0MuT=D1f0MuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0SigmaT=D1f0SigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   
   NumericVector D1Pi(S);
   NumericVector D1MuT(G);
   NumericVector D1SigmaT(G);
   for (int i = 0; i < G; i++){
	   D1MuT[i]=0.0;
	   D1SigmaT[i]=0.0;
	   for (int j=0; j < S; j++){
		   D1MuT[i]=D1MuT[i]-D1f0MuT(i,j)/f0(i,j);
		   D1SigmaT[i]=D1SigmaT[i]-D1f0SigmaT(i,j)/f0(i,j);
	   }
   }
   for (int j = 0; j < S; j++){
	   D1Pi[j]=0.0;
	   for (int i=0; i < G; i++){
		   D1Pi[j]=D1Pi[j]-D1f0Pi(i,j)/f0(i,j);
	   }
   }
   
   NumericVector D1(2*G+S);
   for (int j = 0; j < S; j++){
	   D1[j]=D1Pi[j];
   }
   for (int i=0; i < G; i++){
	   D1[S+2*i]=D1MuT[i];
	   D1[S+2*i+1]=D1SigmaT[i];
   }
   return(D1);
}

// [[Rcpp::export]]
NumericVector D1Pi_Loglikelihood_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes
   NumericMatrix f0=f0_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0Pi=D1f0Pi_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   
   NumericVector D1Pi(S);
   for (int j = 0; j < S; j++){
	   D1Pi[j]=0.0;
	   for (int i=0; i < G; i++){
		   D1Pi[j]=D1Pi[j]-D1f0Pi(i,j)/f0(i,j);
	   }
   }
   
   return(D1Pi);
}

// [[Rcpp::export]]
NumericVector D1MuTSigmaT_Loglikelihood_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes
   NumericMatrix f0=f0_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0MuT=D1f0MuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0SigmaT=D1f0SigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   
   NumericVector D1MuT(G);
   NumericVector D1SigmaT(G);
   for (int i = 0; i < G; i++){
	   D1MuT[i]=0.0;
	   D1SigmaT[i]=0.0;
	   for (int j=0; j < S; j++){
		   D1MuT[i]=D1MuT[i]-D1f0MuT(i,j)/f0(i,j);
		   D1SigmaT[i]=D1SigmaT[i]-D1f0SigmaT(i,j)/f0(i,j);
	   }
   }
   
   NumericVector D1(2*G);
   for (int i=0; i < G; i++){
	   D1[2*i]=D1MuT[i];
	   D1[2*i+1]=D1SigmaT[i];
   }
   return(D1);
}

// [[Rcpp::export]]
NumericVector D1MuT_Loglikelihood_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes
   NumericMatrix f0=f0_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0MuT=D1f0MuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   
   NumericVector D1MuT(G);
   for (int i = 0; i < G; i++){
	   for (int j=0; j < S; j++){
		   D1MuT[i]=D1MuT[i]-D1f0MuT(i,j)/f0(i,j);
		}
   }
   return(D1MuT);
}

// [[Rcpp::export]]
NumericVector D1SigmaT_Loglikelihood_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes
   NumericMatrix f0=f0_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0SigmaT=D1f0SigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);

   NumericVector D1SigmaT(G);
   for (int i = 0; i < G; i++){
	   D1SigmaT[i]=0.0;
	   for (int j=0; j < S; j++){
		   D1SigmaT[i]=D1SigmaT[i]-D1f0SigmaT(i,j)/f0(i,j);
	   }
   }
   
   return(D1SigmaT);
}

// [[Rcpp::export]]
double log_divide(double x, double f){
	double y=0.0;
	if(x>0.0){
		y=exp(log(x)-log(f));
	}else{
		y=-exp(log(-x)-log(f));
	}
	return(y);
}

// [[Rcpp::export]]
NumericVector D1Loglikelihood_log_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes
   NumericMatrix f0=f0_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0Pi=D1f0Pi_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0MuT=D1f0MuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0SigmaT=D1f0SigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   
   NumericVector D1Pi(S);
   NumericVector D1MuT(G);
   NumericVector D1SigmaT(G);
   for (int i = 0; i < G; i++){
	   D1MuT[i]=0.0;
	   D1SigmaT[i]=0.0;
	   for (int j=0; j < S; j++){
		   D1MuT[i]=D1MuT[i]-log_divide(D1f0MuT(i,j),f0(i,j));
		   D1SigmaT[i]=D1SigmaT[i]-log_divide(D1f0SigmaT(i,j),f0(i,j));
	   }
   }
   for (int j = 0; j < S; j++){
	   D1Pi[j]=0.0;
	   for (int i=0; i < G; i++){
		   D1Pi[j]=D1Pi[j]-log_divide(D1f0Pi(i,j),f0(i,j));
	   }
   }
   
   NumericVector D1(2*G+S);
   for (int j = 0; j < S; j++){
	   D1[j]=D1Pi[j];
   }
   for (int i=0; i < G; i++){
	   D1[S+2*i]=D1MuT[i];
	   D1[S+2*i+1]=D1SigmaT[i];
   }
   return(D1);
}

// [[Rcpp::export]]
NumericMatrix D2Loglikelihood_log_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes
   NumericMatrix Hessian(2*G+S,2*G+S);
   NumericMatrix f0=f0_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0Pi=D1f0Pi_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0MuT=D1f0MuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0SigmaT=D1f0SigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0Pi=D2f0Pi_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0PiMuT=D2f0PiMuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0PiSigmaT=D2f0PiSigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0MuT=D2f0MuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0MuTSigmaT=D2f0MuTSigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0SigmaT=D2f0SigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   
   NumericVector D2Pi(S);
   NumericVector D2MuT(G);
   NumericVector D2SigmaT(G);
   for (int j = 0; j < S; j++){
	   D2Pi[j]=0.0;
	   for (int i=0; i < G; i++){
		   D2Pi[j]=D2Pi[j]+log_divide((log_divide(pow(D1f0Pi(i,j),2.0),f0(i,j))-D2f0Pi(i,j)),f0(i,j));
		   Hessian(j,S+2*i)=log_divide((log_divide(D1f0Pi(i,j)*D1f0MuT(i,j),f0(i,j))-D2f0PiMuT(i,j)),f0(i,j));               //D2PiMuT
		   Hessian(S+2*i,j)=Hessian(j,S+2*i);
		   Hessian(j,S+2*i+1)=log_divide((log_divide(D1f0Pi(i,j)*D1f0SigmaT(i,j),f0(i,j))-D2f0PiSigmaT(i,j)),f0(i,j));       //D2PiSigmaT
		   Hessian(S+2*i+1,j)=Hessian(j,S+2*i+1);
	   }
   }
   for (int i = 0; i < G; i++){
	   D2MuT[i]=0.0;
	   D2SigmaT[i]=0.0;
	   Hessian(S+2*i,S+2*i+1)=0.0;
	   for (int j=0; j < S; j++){
		   D2MuT[i]=D2MuT[i]+log_divide((log_divide(pow(D1f0MuT(i,j),2.0),f0(i,j))-D2f0MuT(i,j)),f0(i,j));
		   D2SigmaT[i]=D2SigmaT[i]+log_divide((log_divide(pow(D1f0SigmaT(i,j),2.0),f0(i,j))-D2f0SigmaT(i,j)),f0(i,j));
		   Hessian(S+2*i,S+2*i+1)=Hessian(S+2*i,S+2*i+1)+log_divide((log_divide(D1f0MuT(i,j)*D1f0SigmaT(i,j),f0(i,j))-D2f0MuTSigmaT(i,j)),f0(i,j));
	   }
	   Hessian(S+2*i+1,S+2*i)=Hessian(S+2*i,S+2*i+1);
   }
   for (int j = 0; j < S; j++){
	   Hessian(j,j)=D2Pi[j];
   }
   for (int i=0; i < G; i++){
	   Hessian(S+2*i,S+2*i)=D2MuT[i];
	   Hessian(S+2*i+1,S+2*i+1)=D2SigmaT[i];
   }
   return(Hessian);
}

// [[Rcpp::export]]
NumericMatrix D2Loglikelihood_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes
   NumericMatrix Hessian(2*G+S,2*G+S);
   NumericMatrix f0=f0_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0Pi=D1f0Pi_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0MuT=D1f0MuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0SigmaT=D1f0SigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0Pi=D2f0Pi_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0PiMuT=D2f0PiMuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0PiSigmaT=D2f0PiSigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0MuT=D2f0MuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0MuTSigmaT=D2f0MuTSigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0SigmaT=D2f0SigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   
   NumericVector D2Pi(S);
   NumericVector D2MuT(G);
   NumericVector D2SigmaT(G);
   for (int j = 0; j < S; j++){
	   D2Pi[j]=0.0;
	   for (int i=0; i < G; i++){
		   D2Pi[j]=D2Pi[j]+(pow(D1f0Pi(i,j),2.0)/f0(i,j)-D2f0Pi(i,j))/f0(i,j);
		   Hessian(j,S+2*i)=(D1f0Pi(i,j)*D1f0MuT(i,j)/f0(i,j)-D2f0PiMuT(i,j))/f0(i,j);               //D2PiMuT
		   Hessian(S+2*i,j)=Hessian(j,S+2*i);
		   Hessian(j,S+2*i+1)=(D1f0Pi(i,j)*D1f0SigmaT(i,j)/f0(i,j)-D2f0PiSigmaT(i,j))/f0(i,j);       //D2PiSigmaT
		   Hessian(S+2*i+1,j)=Hessian(j,S+2*i+1);
	   }
   }
   for (int i = 0; i < G; i++){
	   D2MuT[i]=0.0;
	   D2SigmaT[i]=0.0;
	   Hessian(S+2*i,S+2*i+1)=0.0;
	   for (int j=0; j < S; j++){
		   D2MuT[i]=D2MuT[i]+(pow(D1f0MuT(i,j),2.0)/f0(i,j)-D2f0MuT(i,j))/f0(i,j);
		   D2SigmaT[i]=D2SigmaT[i]+(pow(D1f0SigmaT(i,j),2.0)/f0(i,j)-D2f0SigmaT(i,j))/f0(i,j);
		   Hessian(S+2*i,S+2*i+1)=Hessian(S+2*i,S+2*i+1)+(D1f0MuT(i,j)*D1f0SigmaT(i,j)/f0(i,j)-D2f0MuTSigmaT(i,j))/f0(i,j);
	   }
	   Hessian(S+2*i+1,S+2*i)=Hessian(S+2*i,S+2*i+1);
   }
   for (int j = 0; j < S; j++){
	   Hessian(j,j)=D2Pi[j];
   }
   for (int i=0; i < G; i++){
	   Hessian(S+2*i,S+2*i)=D2MuT[i];
	   Hessian(S+2*i+1,S+2*i+1)=D2SigmaT[i];
   }
   return(Hessian);
}

// [[Rcpp::export]]
List D2Loglikelihood_unit_2D(NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT) {
   int S=y.ncol();  // No. of samples
   int G=y.nrow();  // No. of Genes
   NumericMatrix f0=f0_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0Pi=D1f0Pi_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0MuT=D1f0MuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D1f0SigmaT=D1f0SigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0Pi=D2f0Pi_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0PiMuT=D2f0PiMuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0PiSigmaT=D2f0PiSigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0MuT=D2f0MuT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0MuTSigmaT=D2f0MuTSigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   NumericMatrix D2f0SigmaT=D2f0SigmaT_func_2D(y,Pi,MuN,MuT,SigmaN,SigmaT);
   
   NumericMatrix D2Pi(G,S);
   NumericMatrix D2PiMuT(G,S);
   NumericMatrix D2PiSigmaT(G,S);
   NumericMatrix D2MuT(G,S);
   NumericMatrix D2MuTSigmaT(G,S);
   NumericMatrix D2SigmaT(G,S);
   for (int j = 0; j < S; j++){
	   for (int i=0; i < G; i++){
		   D2Pi(i,j)=(pow(D1f0Pi(i,j),2.0)/f0(i,j)-D2f0Pi(i,j))/f0(i,j);
		   D2PiMuT(i,j)=(D1f0Pi(i,j)*D1f0MuT(i,j)/f0(i,j)-D2f0PiMuT(i,j))/f0(i,j);
		   D2PiSigmaT(i,j)=(D1f0Pi(i,j)*D1f0SigmaT(i,j)/f0(i,j)-D2f0PiSigmaT(i,j))/f0(i,j);
		   D2MuT(i,j)=(pow(D1f0MuT(i,j),2.0)/f0(i,j)-D2f0MuT(i,j))/f0(i,j);
		   D2MuTSigmaT(i,j)=(D1f0MuT(i,j)*D1f0SigmaT(i,j)/f0(i,j)-D2f0MuTSigmaT(i,j))/f0(i,j);
		   D2SigmaT(i,j)=(pow(D1f0SigmaT(i,j),2.0)/f0(i,j)-D2f0SigmaT(i,j))/f0(i,j);
	   }
   }
   return(List::create(Named("D2Pi")=D2Pi,Named("D2PiMuT")=D2PiMuT,Named("D2PiSigmaT")=D2PiSigmaT,Named("D2MuT")=D2MuT,Named("D2MuTSigmaT")=D2MuTSigmaT,Named("D2SigmaT")=D2SigmaT));
}

// [[Rcpp::export]]
double Gfunc_2D_C(double alpha, NumericVector x, NumericVector p, int S, int G, NumericMatrix y,NumericVector MuN, NumericVector SigmaN){
	NumericVector x_mid=x+alpha*p;
	List L=x_update_inv_2D(x_mid,S,G);
	NumericVector Pi_mid=L[0];
	NumericVector MuT_mid=L[1];
	NumericVector SigmaT_mid=L[2];
	double obj_mid=Loglikelihood_2D(y, Pi_mid, MuN, MuT_mid, SigmaN, SigmaT_mid);
	return(obj_mid);
}

// [[Rcpp::export]]
double GoldenLine_search_2D(double alpha, NumericVector x, NumericVector p, NumericMatrix y, NumericVector MuN, NumericVector SigmaN, int S, int G, double steps) {
	double a=0.0;
	double b=alpha;
	double gr=(sqrt(5.0) + 1.0) / 2.0;
	double c = b - (b - a) / gr;
	double d = a + (b - a) / gr;
	int i=1;
	
	while(fabs(c - d) > (1e-5)/steps){
		//printf("%f\n",Gfunc_2D_C(c,x,p,S,G,y,MuN,SigmaN));
		if(Gfunc_2D_C(c,x,p,S,G,y,MuN,SigmaN) < Gfunc_2D_C(d,x,p,S,G,y,MuN,SigmaN)){
			b=d;
		}else{
			a=c;
		}
		c = b - (b - a) / gr;
		d = a + (b - a) / gr;
		i++;
	}
	return((c+d)/2);
}

// [[Rcpp::export]]
double GoldenLine_search_Pi_2D(double alpha, NumericVector x, NumericVector p, NumericMatrix y, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, NumericVector SigmaT, double steps) {
	double a=0.0;
	double b=alpha;
	double gr=(sqrt(5.0) + 1.0) / 2.0;
	double c = b - (b - a) / gr;
	double d = a + (b - a) / gr;
	int i=1;
	
	while(fabs(c - d) > (1e-5)/steps){
		NumericVector Pi_c=x+c*p;
		NumericVector Pi_d=x+d*p;
		if(Loglikelihood_2D(y, Pi_c, MuN, MuT, SigmaN, SigmaT) < Loglikelihood_2D(y, Pi_d, MuN, MuT, SigmaN, SigmaT)){
			b=d;
		}else{
			a=c;
		}
		c = b - (b - a) / gr;
		d = a + (b - a) / gr;
		i++;
	}
	return((c+d)/2);
}

// [[Rcpp::export]]
double GoldenLine_search_MuT_2D(double alpha, NumericVector x, NumericVector p, NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector SigmaN, NumericVector SigmaT, double steps) {
	double a=0.0;
	double b=alpha;
	double gr=(sqrt(5.0) + 1.0) / 2.0;
	double c = b - (b - a) / gr;
	double d = a + (b - a) / gr;
	int i=1;
	
	while(fabs(c - d) > (1e-5)/steps){
		NumericVector MuT_c=x+c*p;
		NumericVector MuT_d=x+d*p;
		if(Loglikelihood_2D(y, Pi, MuN, MuT_c, SigmaN, SigmaT) < Loglikelihood_2D(y, Pi, MuN, MuT_d, SigmaN, SigmaT)){
			b=d;
		}else{
			a=c;
		}
		c = b - (b - a) / gr;
		d = a + (b - a) / gr;
		i++;
	}
	return((c+d)/2);
}

// [[Rcpp::export]]
double GoldenLine_search_SigmaT_2D(double alpha, NumericVector x, NumericVector p, NumericMatrix y, NumericVector Pi, NumericVector MuN, NumericVector MuT, NumericVector SigmaN, double steps) {
	double a=0.0;
	double b=alpha;
	double gr=(sqrt(5.0) + 1.0) / 2.0;
	double c = b - (b - a) / gr;
	double d = a + (b - a) / gr;
	int i=1;
	
	while(fabs(c - d) > (1e-5)/steps){
		NumericVector SigmaT_c=x+c*p;
		NumericVector SigmaT_d=x+d*p;
		if(Loglikelihood_2D(y, Pi, MuN, MuT, SigmaN, SigmaT_c) < Loglikelihood_2D(y, Pi, MuN, MuT, SigmaN, SigmaT_d)){
			b=d;
		}else{
			a=c;
		}
		c = b - (b - a) / gr;
		d = a + (b - a) / gr;
		i++;
	}
	return((c+d)/2);
}
