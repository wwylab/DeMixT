#include <math.h>
#include <R.h>

static double stirling_approx(double z)
{
    if(z < 10)
    {
        return lgamma(z+1);
    }
    return z * log(z) - z + 0.5 * (log(M_PI * 2.0) + log(z));
}

static double binomCoef(long k, double r)
{
    //double coef = stirling_approx(r+k-1) - stirling_approx(r-1) - stirling_approx((double) k);
    double coef = lgamma(r + k) - lgamma(r) - lgamma((double) k + 1);
    return coef;
}

static double dnbinom(long x, double mu, double size)
{
    if (x == 0 && size == 0) return 1;
    double binom_coef = binomCoef(x, size);
    double res = binom_coef + ((double) x) * (log(mu) - log(size + mu)) + size * (log(size) - log(size + mu));
    return exp(res);
}
