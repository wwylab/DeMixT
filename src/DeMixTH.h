
#ifndef Zeya_Wang
#define Zeya_Wang
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define min(a,b) \
	({__typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a < _b ? _a : _b;})

// parameter definition
typedef struct param{
	//average parameter
	double *Navg1;
	double *Navg2;
	double *Tavg;
   //variance parameter
	double *Nsigma1;
	double *Nsigma2;
	double *Tsigma;
   //pi
	double *pi1;
	double *pi2;
	double *piT;
	//mle
	double **mle;
    double *obj;
} PARAM;



// Function initialize
void initialSet(PARAM *qq);
void load_data(double *mat1);
void saveFiles(double *put1, double *put3, double *put5, double *put7, double *put9, double *put11, double *obj_out, double *put31, double *put32);


//my new function


void gettumor(int genes, int h); // option h = 1 for one component; h = 2 for two
void getnormal1(int genes, unsigned long *seed, int k, int t, double vap1, double del1);		// opt 0 for Mu_n, 1 for sigma_n
void getnormal2(int genes, unsigned long *seed, int k, int t, double vap1, double del1);

void getpi2(int samp, int h,unsigned long *seed);
double lf1(int opt, int opt2,  double y_sum, double pi_sum, double nval1);
double fmin1(double ax, double bx,int iG, int iS, double y_sum, double pi_sum, double (*f)(int, int, double, double, double), double tol);
double lf2(int opt, int opt2, double nval2);
double fmin2(double ax, double bx,int iG, int iS, double (*f)(int, int, double), double tol);
double pf1(double *y_sum, double pi_sum, double *nval1, double pii1, int *Gid1, int size);
double pmin1(double ax, double bx, double *y_sum, double pi_sum, double *nval1, int *Gid1, int size, double (*f)(double*, double, double*, double, int*, int), double tol);
double pf2(int opt, double *nval1, double *nval2, double pii2, int *Gid2, int size);
double pmin2(double ax, double bx,int iS, double *nval1, double *nval2, int *Gid2, int size, double (*f)(int, double*, double*, double, int*, int), double tol);
void getmle(int opt, int opt2, int h);
///new function to add
double ft_y(double y, double mung, double mutg, double sng, double stg, double pi1, double pi2);
double ft_y_SC(double y, double mung, double mutg, double sng, double stg, double pi1, double pi2);
double ft_y2(double y, double mung1, double mung2, double mutg, double sng1, double sng2, double stg, double pi1, double pi2);
double pf_y(int samp, double pi1);
double pf_y2(int samp, double pi1, double pi2);
double pmin_y(double ax, double bx, int samp, double (*f)(int, double), double tol);
double pmin_y2(double ax, double bx, int samp, double pi2, double (*f)(int, double, double), double tol);
double minpi(int samp, double pi2);

void getpi(int samp, int h);
void getspikeinpi(int samp);
void getpiT(int samp);
double gammaln(double xx);
double gamma(double xx);


double tf_y(int genes, double mu, double sigma);
double tf_y2(int genes, double mu, double sigma);
double mint(int genes, int h, double mu); // opt h =1:one component; h =2:two component
double tmin_y(double ax, double bx, int genes, int h, double (*f)(int, int, double), double tol);
double tmin_y2(double ax, double bx, int genes, double mu, double (*f)(int, double, double), double tol);
double pf_yT(int samp, double pi1, double piT);
double NB_nt(double nt, double y, double mung, double mutg, double neta, double teta, double pi1, double pi2);
double ft_nb(double y, double mung, double mutg, double neta, double teta, double pi1, double pi2);
double ft_nb2(double y, double mung1, double mung2, double mutg, double neta1, double neta2, double teta, double pi1, double pi2);




double nitg_ft_y(double t, double y, double mung, double mutg, double sng, double stg, double pi1);
double min_nitg_ft_y(double ax, double bx, double y, double mung, double mutg, double sng, double stg, double pi1, double (*f)(double, double, double, double, double, double,double), double tol);








///other functions
double pnorm(double tmp1, double *tmp1pos);
double cnorm(double val);
double dnorm(double val, double mu, double sigma);
double sum(double *array, int size);
double cinvgamma(double val);


// Random Number generation function
double kiss(unsigned long *);
double snorm(unsigned long *);
double rnorm(double,double,unsigned long *);
void revsort(double *, int *, int);
int ProbSampleReplace(int, double *, unsigned long *);
double sd(double *, int);
double mean(double *, int);



// Variable allocation, arguements for paramters
PARAM *p;
double **FD, **CD;
double **avgparN, **sigparN, **avgparT, **sigparT;
double  **tmppi1, **tmppi2;
int iteration,iteration1;

int nP,nG, nS, nSp, nHavepi, Cid, nmle;
double vap1, del1, vap1m, del1m;
double pvap1, pdel1, pvap1m, pdel1m;

int fNorm, fNorm1, fNorm2, intx;			// Number of normal samples;
int integ; //number of integration bins
int opt_seq; // indicator of sequencing technique
unsigned long *seed;
double M;


double cnorm(double val){
	double tmpval;
	tmpval = exp(-pow(val, 2)/20000);
	return(tmpval);
}

double cinvgamma(double val){
	double tmpval;
	tmpval = pow(val, -1.001)*exp(-0.001/val);
	return(tmpval);
}


double dnorm(double val, double mu, double sigma){
	double tmpval;
	tmpval = exp(-pow((val - mu), 2)/(2*pow(sigma, 2)))/sigma/pow(2*M_PI, 0.5);
	return(tmpval);
}





double sum(double *array, int size)
{
	int i;
	double sum = 0.0;
	for (i=0;i<size;i++){
		sum += array[i];
	}
	return(sum);
}



double mean(double *data, int n)
{
    return sum(data, n)/((double)n);
}



double sd(double *data, int n)
{
    double mu=0.0, sum_deviation=0.0, tmp;
    int i;
    mu = mean(data, n);
    for(i=0; i<n;i++)
    {
    tmp = pow((data[i] - mu),2);
    sum_deviation += tmp;
    }
    return sqrt(sum_deviation/((double)n-1));
}

double gammaln(double xx)
{
	double x,y,tmp,ser;
	static double cof[14]={57.1562356658629235, -59.5979603554754912, 14.1360979747417471, -0.491913816097620199,
		0.339946499848118887e-4, 0.465236289270485756e-4, -0.983744753048795646e-4, 0.158088703224912494e-3,
		-0.210264441724104883e-3, 0.217439618115212643e-3, -0.164318106536763890e-3, 0.844182239838527433e-4,
		-0.261908384015814087e-4, 0.368991826595316234e-5};
	int j;
	
	y=x=xx;
	tmp=x+5.24218750000000000;
	tmp = (x+0.5)*log(tmp) -tmp;
	ser=0.999999999999997092;
	for (j=0;j<14;j++) ser += cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}

double gamma(double xx)
{
	if(xx==0.0) return 99999.0;
	
	
	return exp(gammaln(xx));
	
}

#endif





