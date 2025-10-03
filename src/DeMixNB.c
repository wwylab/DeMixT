#include <stdio.h>
#include <R.h>
#include <omp.h>
#include "DeMixNB.h"

double Likelihood_NB(double y_ig, double pi_i, long l_ig, double mu_tg, 
    double phi_tg, double mu_ng, double phi_ng)
{
    double D1 = dnbinom(l_ig, pi_i*mu_tg, phi_tg);
    double D2 = dnbinom((long)floor(y_ig) - l_ig, (1-pi_i)*mu_ng, phi_ng);
    return D1*D2; 
}

void LLK_NB_pi_i(double* pi_i, double* y_i, double* mu_t, double* phi_t,
    double* mu_n, double* phi_n, int* nG, double* out)
{
    double LLk[*nG];
    double LH;
    double lk;
    long l;
    #pragma omp parallel for private(LH, lk, l) schedule(dynamic)
    for(int g = 0; g < *nG; g++)
    {
        LH = 0;
        for(l = 0; l <= floor(y_i[g]); l++)
        {
            lk = Likelihood_NB(y_i[g], *pi_i, l, mu_t[g], phi_t[g], mu_n[g], phi_n[g]);
            LH += lk;
        }
        
        if(LH == 0 || ISNAN(LH) || isinf(LH))
        {
            LLk[g] = -1e-8;
        } 
        else
        {
            LLk[g] = log(LH);
        }
    }
    *out = 0;
    for (int g = 0; g < *nG; g++)
    {
        if (LLk[g] != -1e-8)
        {
           *out += LLk[g] * (-1);
        } 
    }
    return; 
}

void LLK_NB_T_g(double* paras, double* y_g, double* pi,
    double* mu_ng, double* phi_ng, int* nS, double* out)
{
    double mu_tg = paras[0];
    double phi_tg = paras[1];
    double LLk[*nS];
    double LH;
    double lk;
    long l;
    #pragma omp parallel for private(LH, lk, l) schedule(dynamic)
    for(int i = 0; i < *nS; i++)
    {
        LH = 0;
        for(l = 0; l <= floor(y_g[i]); l++)
        {
            lk = Likelihood_NB(y_g[i], pi[i], l, mu_tg, phi_tg, *mu_ng, *phi_ng);
            LH += lk;
        }

        if(LH == 0 || ISNAN(LH) || isinf(LH))
        {
            LLk[i] = -1e-8;
        }
        else
        {
            LLk[i] = log(LH);
        }
    }
    *out = 0;
    for(int i = 0; i < *nS; i++)
    {
        if (LLk[i] != -1e-8)
        {
           *out += LLk[i] * (-1);   
        }
    }
    return;
}

int main()
{
    // Tests to verify that the likelihood calculations are correct
    Rprintf("%f\n", dnbinom(100, 120, 10));
    return 0;
}
