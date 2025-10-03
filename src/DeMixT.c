/***************************************************************************************************************
* DeMixT Latest Version
* By Zeya Wang
* Initial Date 12.13.2014
****************************************************************************************************************/

  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <time.h>
  #include <R.h>
  #include <string.h>
  #include <float.h> /* DBL_EPSILON */
  #include "DeMixTH.h"
  #ifdef _OPENMP
    #include<omp.h> /* parallel computing */
  #endif

int checkopenmp(int* numthread)
{
    #ifdef _OPENMP
    #pragma omp parallel
    #pragma omp master
    *numthread = omp_get_num_threads();
    #else
    #endif
    return 0;
}


void Tdemix(double *data, int *nGroup, int *nsamp, int *ngenes, int *nspikein, int *npi, double *pi01, double *pi02, double *fixpi1, double *fixpi2, double *fixpi3,int *nCid, int *niter, int *ninteg, double *tol, int *thread, double *s0, double *m0, double *output1, double *output3, double *output5, double *output7, double *output9, double *output11, double *obj_out, double *output31, double *output32)
{
    //nCid = 1, we have just one stroma component; =2, we have 2
  int i, j, k, l;

  double *st1_mu, *st2_mu, *st1_sig, *st2_sig;
  double *st1_mu_2, *st2_mu_2, *st1_sig_2, *st2_sig_2;
  double **stroma1, **stroma2;
  double **stroma1_2, **stroma2_2;
  double total_tol;
  double **mixed;
  double *mixedm;
  #ifdef _OPENMP 
    int nthread=*thread;
  #endif
  
    //clock_t start_t, end_t;
    
    //initial value for sigma and mu
    double ss0 = *s0;
    double mm0 = *m0;
    
  //set the number of threads
  #ifdef _OPENMP 
    omp_set_num_threads(nthread);
  #endif
  
	
  //get value from input parameter
  nS=*nsamp;            // Number of Samples
  nG=*ngenes;           // Number of Genes
  nSp=*nspikein;  // Numer of Spike in Samples
  nHavepi =  *npi;      // Have pi or not: 0: Not given pi1 and pi2; 1: given pi1 and pi2; 2: given piT
  Cid = *nCid; //indicator for number of components
  integ = *ninteg;
  total_tol = *tol;
  iteration = *niter; // number of iterations
  //Rprintf("iteration is %d\n", iteration);



  // 1 for component 1, 2 for component 2, 3 for mixed sample
  fNorm1=0;
  fNorm2=0;
  for(j=0;j<nS;j++)
  {
    if(nGroup[j]==1) fNorm1++;
    if(nGroup[j]==2) fNorm2++;
    //Rprintf("%d : %d \n", j, nGroup[j]);
  }

  // Determine the number of tumor sample
  fNorm = fNorm1 + fNorm2;
  intx=nS-fNorm;




  nmle = (Cid+1)*intx;

  p = (PARAM *)calloc(1,sizeof(PARAM));

  // FD, Parameter initialized
  FD = calloc(nS ,sizeof(double *));


  for(j=0;j<nS;j++) FD[j]= calloc(nG, sizeof(double));


  //Rprintf("Setting is over\n");

  // Data transfer
  initialSet(p);
  load_data(data);

  //Rprintf("Loading is over\n");

  tmppi1 =calloc(iteration ,sizeof(double *));  // save pi1 estimates in each iteration
  tmppi2 =calloc(iteration ,sizeof(double *));  // save pi2 estimates in each iteration

  for(j=0;j<iteration;j++)
  {
    tmppi1[j]= calloc(intx, sizeof(double ));
    tmppi2[j]= calloc(intx, sizeof(double ));
  }
  //Rprintf("Loading1 is over\n");

  //calculate mean and sd of normal samples
  stroma1 = calloc(nG, sizeof(double *));
  stroma2 = calloc(nG, sizeof(double *));
  st1_mu = calloc(nG ,sizeof(double));
  st2_mu = calloc(nG ,sizeof(double));
  st1_sig = calloc(nG ,sizeof(double));
  st2_sig = calloc(nG ,sizeof(double));

  stroma1_2 = calloc(nG, sizeof(double *));
  stroma2_2 = calloc(nG, sizeof(double *));
  st1_mu_2 = calloc(nG ,sizeof(double));
  st2_mu_2 = calloc(nG ,sizeof(double));
  st1_sig_2 = calloc(nG ,sizeof(double));
  st2_sig_2 = calloc(nG ,sizeof(double));

  //Rprintf("Loading2 is over\n");

  for(j=0;j<nG;j++) stroma1[j]= calloc(fNorm1, sizeof(double));
  for(j=0;j<nG;j++) stroma2[j]= calloc(fNorm2, sizeof(double));
  for(j=0;j<nG;j++) stroma1_2[j]= calloc(fNorm1, sizeof(double));
  for(j=0;j<nG;j++) stroma2_2[j]= calloc(fNorm2, sizeof(double));

    //add an array to store mixed tumor
    mixed = calloc(nG ,sizeof(double *));
    for(j=0;j<nG;j++) mixed[j]= calloc(intx, sizeof(double));
    mixedm = calloc(nG ,sizeof(double));
  //Rprintf("Loading3 is over\n");

    for(i=0;i<nG;i++)
    {
        for(j=0;j<fNorm1;j++) stroma1[i][j] = log2(FD[j][i]);
        if(Cid == 2){
            for(j=0;j<fNorm2;j++) stroma2[i][j] = log2(FD[j+fNorm1][i]);
        }
        for(j=0;j<intx;j++) mixed[i][j] = FD[j+fNorm][i];
    }

    for(i=0;i<nG;i++)
    {
        st1_mu[i] = mean(stroma1[i], fNorm1);
        st1_sig[i] = sd(stroma1[i], fNorm1);
        mixedm[i] = mean(mixed[i], intx);
        if(Cid == 2){
            st2_mu[i] = mean(stroma2[i], fNorm2);
            st2_sig[i] = sd(stroma2[i], fNorm2);
        }
    }


    for(i=0;i<nG;i++) free(stroma1[i]);
  for(i=0;i<nG;i++) free(stroma2[i]);
  for(i=0;i<nG;i++) free(stroma1_2[i]);
  for(i=0;i<nG;i++) free(stroma2_2[i]);

  free(stroma1);
  free(stroma2);
  free(stroma1_2);
  free(stroma2_2);

    double Tavgtmp;
    for(j=0;j<nG;j++)
    {
        
        
        p->Navg1[j]=st1_mu[j];
        if(Cid == 2) p->Navg2[j]=st2_mu[j];
        p->Nsigma1[j]= pow(st1_sig[j], 2.0);
        if(Cid == 2) p->Nsigma2[j]= pow(st2_sig[j], 2.0);
        //p->Tsigma[j]=0.25;
        p->Tsigma[j] = ss0;
        
        if(Cid == 1){
            Tavgtmp = 2*mixedm[j] - exp(p->Navg1[j]*log(2.0) + p->Nsigma1[j]/2.0*pow(log(2.0), 2.0));
        }else if(Cid == 2){
            Tavgtmp = 3*mixedm[j] - exp(p->Navg1[j]*log(2.0) + p->Nsigma1[j]/2.0*pow(log(2.0), 2.0)) - exp(p->Navg2[j]*log(2.0) + p->Nsigma2[j]/2.0*pow(log(2.0), 2.0));
        }
        
        if(Tavgtmp <= 0) Tavgtmp = exp(10*log(2.0) + p->Tsigma[j]/2.0*pow(log(2.0), 2.0));
        p->Tavg[j] = (log(Tavgtmp) - p->Tsigma[j]/2.0*pow(log(2.0), 2.0))/log(2.0);
        if(mm0 > 0.0) p->Tavg[j] = mm0;
        //p->Tavg[j] = 10.0;
        //Rprintf("gene avg for %d is %lf\n", j, p->Tavg[j]);
    }
	
	//Rprintf("Initial pi0 \n");
	for(k=0;k<intx;k++)
    {
      p->pi1[k]= pi01[k];
	  //Rprintf("%15.3f \t",  p->pi1[k], '\n');
      if(Cid == 2) p->pi2[k]= pi02[k];
    }
	

    if(nHavepi==1)
  {
    for(k=0;k<intx;k++)
    {
      p->pi1[k]= fixpi1[k];
      if(Cid == 2) p->pi2[k]= fixpi2[k];
    }
  }else if(nHavepi == 2){
    for(k=0;k<intx;k++) p->piT[k] = fixpi3[k];
  }

  //Rprintf("Loading is over\n");

  free(st1_mu);
  free(st2_mu);
  free(st1_sig);
  free(st2_sig);

  //begin our iteration

//  for(j=0;j<nG;j++)
  //  {
  //    Rprintf("nsigma1 for %d is %lf\n",j, p->Nsigma1[j]);
  //    Rprintf("nsigma2 for %d is %lf\n",j, p->Nsigma2[j]);
  //    Rprintf("navg1 for %d is %lf\n",j, p->Navg1[j]);
  //    Rprintf("navg2 for %d is %lf\n",j, p->Navg2[j]);
//    }


  CD = calloc(intx ,sizeof(double *));
  for(j=0;j<intx;j++) CD[j]= calloc(nG, sizeof(double));
  //Rprintf("Iteration is starting\n");

  avgparT = calloc(iteration, sizeof(double *));
  sigparT = calloc(iteration, sizeof(double *));
  for(j = 0; j < iteration; j++) avgparT[j] = calloc(nG, sizeof(double));
  for(j = 0; j < iteration; j++) sigparT[j] = calloc(nG, sizeof(double));

    for(k = 0; k < iteration; k++){
        for(j = 0; j < nG; j++){
            avgparT[k][j] = 0.0;
            sigparT[k][j] = 0.0;
        }
        for(l = 0; l < intx; l++){
            tmppi1[k][l] = 0.0;
            tmppi2[k][l] = 0.0;
        }
        p->obj[k] = 0;
    }


  iteration1=0;
  double obj_old_all=0.0;

    for(i=0;i<iteration;i++)
    {
		//if(i==0) start_t=clock();
        
        if (nHavepi != 1)  Rprintf("Iteration %d: updating parameters\n", i+1);
        #ifdef _OPENMP
         #pragma omp parallel for //openmp
        for(j=0;j<nG;j++)
        {
            gettumor(j, Cid);
            avgparT[i][j] = p->Tavg[j];
            sigparT[i][j] = p->Tsigma[j];
        }
         #else
        for(j=0;j<nG;j++)
        {
            gettumor(j, Cid);
            avgparT[i][j] = p->Tavg[j];
            sigparT[i][j] = p->Tsigma[j];
        }
         #endif
        
      if (nHavepi != 1)  Rprintf("Iteration %d: updating proportions\n", i+1);
        //updating pi value
        if(nHavepi==0)
        {
            #ifdef _OPENMP
            // multithreaded OpenMP version of code
            #pragma omp parallel for //openmp
            for(l=0;l<intx;l++) 
			{
				if(l < intx - nSp)
				{
					//Rprintf("Spike in %d \n", nSp);
					getpi(l, Cid);
					//Rprintf("Number of Spike in %d \n", l+1);
				}else{
					getspikeinpi(l);
				}			
			}
           #else
            //single-threaded version of code
            for(l=0;l<intx;l++) 
			{
				if(l < intx - nSp)
				{
					//Rprintf("Spike in %d \n", nSp);
					getpi(l, Cid);
					//Rprintf("Number of Spike in %d \n", l+1);
				}else{
					getspikeinpi(l);
				}	
			}
          #endif
        }else if(nHavepi == 2){
          #ifdef _OPENMP
          #pragma omp parallel for //openmp
            for(l=0;l<intx;l++) getpiT(l); // this is to estimate pi1 and pi2 given piT in the three component case
          #else
            for(l=0;l<intx;l++) getpiT(l);
           #endif
        }
        
        for(j=0;j<intx;j++)
        {
            tmppi1[i][j] = p->pi1[j];
	    if (nHavepi != 1) {};  //Rprintf("%15.3f \t",  p->pi1[j]);
            if(Cid == 2)
            {
                tmppi2[i][j] = p->pi2[j];
		if (nHavepi != 1){};  //Rprintf(" %15.3f \t",  p->pi2[j]);
            }
            if (nHavepi != 1) {};  //Rprintf("\n");
        }
        
        
        // Total time expectation depending on the first five iterations
        if(i==1)
        {
            //end_t=clock();
            //Rprintf("================================================================\n");
            //Rprintf("Deconvolution is estimated to be finished in at most %lf hours\n", (double)(end_t-start_t)/CLOCKS_PER_SEC*((double)iteration/3600.0));
            //Rprintf("================================================================\n");
        }

        // objective function 2D
        double ttmp=0.0, tttmp;
        for(j=0;j<nG;j++)
        {
            for(l=0;l<intx;l++)
            {
				if(Cid == 2){
				 tttmp = ft_y2(FD[l+fNorm][j], p->Navg1[j], p->Navg2[j], p->Tavg[j], p->Nsigma1[j], p->Nsigma2[j], p->Tsigma[j], p->pi1[l], p->pi2[l]);
				  ttmp = ttmp - tttmp;
				}else{
				 tttmp = ft_y(FD[l+fNorm][j], p->Navg1[j], p->Tavg[j], p->Nsigma1[j], p->Tsigma[j], p->pi1[l], p->pi2[l]);
			    ttmp = ttmp - tttmp;
			    //Rprintf("obj teset is %f, %d\n", tttmp, j);
			    //Rprintf("obj val is %f\t%f\t%f\t %f\t%f\t%f\t%f\n",FD[l+fNorm][j], p->Navg1[j], p->Tavg[j], p->Nsigma1[j], p->Tsigma[j], p->pi1[l], p->pi2[l]);
			    
				}

            }
        }

        if(i>0){
            obj_old_all=p->obj[i-1];
        }


        p->obj[i] = ttmp;
//Rprintf("obj teset is %f\n", ttmp);
        iteration1++;
        //Rprintf("Obj error is %f\n", fabs(ttmp - obj_old_all));

        if(fabs(ttmp - obj_old_all)<total_tol*fabs(obj_old_all)){
            Rprintf("Break at %d\n", iteration1);
            break;
        }

    }

    //final pi value

    
    if(nHavepi==0)
    {
      #ifdef _OPENMP
       #pragma omp parallel for //openmp
        for(l=0;l<intx;l++)
		{
			if(l < intx - nSp)
				{
					//Rprintf("Spike in %d \n", nSp);
					getpi(l, Cid);
					//Rprintf("Number of Spike in %d \n", l+1);
				}else{
					getspikeinpi(l);
				}	
		}			
       #else
        for(l=0;l<intx;l++) 
		{
			if(l < intx - nSp)
				{
					//Rprintf("Spike in %d \n", nSp);
					getpi(l, Cid);
					//Rprintf("Number of Spike in %d \n", l+1);
				}else{
					getspikeinpi(l);
				}	
		}
	
       #endif
    }else if(nHavepi == 2){
       #ifdef _OPENMP
    #pragma omp parallel for //openmp
        for(l=0;l<intx;l++) getpiT(l); // this is to estimate pi1 and pi2 given piT in the three component case
   #else
        for(l=0;l<intx;l++) getpiT(l);
   #endif
        
    }
    
    
    // deconvolve mle

    for(j=0;j<nG;j++)
    {
        for(k=0;k<intx;k++)
        {
            getmle(k, j, Cid);
            CD[k][j] = p->mle[j][k]; //deconvolved value
        }

    }

    fflush(NULL);
    saveFiles(output1, output3, output5, output7, output9, output11, obj_out, output31, output32);

    //Rprintf("Step 4: Save files done : \n");
    //free variable
    free(p->Navg1);
    free(p->Navg2);
    free(p->Tavg);
    free(p->Nsigma1);
    free(p->Nsigma2);
    free(p->Tsigma);
    free(p->pi1);
    free(p->pi2);
    free(p->obj);

    for(i=0;i<nG;i++)
    {
        free(p->mle[i]);
    }

    free(p->mle);
    free(p);
}





/************************************************************************************
 * Function to load data and save data
************************************************************************************/
  // Data transfer
void load_data(double *mat1)
{
  int j,k;

  for(k=0;k<nS;k++)
  {
    for(j=0;j<nG;j++)
    {
      FD[k][j]=  mat1[nG*k+j];
    }
  }
  if (nHavepi != 1) Rprintf("There are  %d normals and %d tumors\n", fNorm,intx);

}




void initialSet(PARAM *qq)
{
    int i;
    qq->Navg1 =calloc(nG ,sizeof(double ));
    qq->Navg2 =calloc(nG ,sizeof(double ));
    qq->Tavg =calloc(nG ,sizeof(double ));
    qq->Nsigma1 =calloc(nG ,sizeof(double ));
    qq->Nsigma2 =calloc(nG ,sizeof(double ));
    qq->Tsigma = calloc(nG ,sizeof(double ));
    qq->pi1 =calloc(intx ,sizeof(double ));
    qq->pi2 =calloc(intx ,sizeof(double ));
    qq->piT = calloc(intx ,sizeof(double ));
    qq->mle =calloc(nG ,sizeof(double *));
    qq->obj =calloc(iteration ,sizeof(double ));
    for(i=0;i<nG;i++) qq->mle[i]=calloc(nmle, sizeof(double ));
}


// Deconvolved expressions
void saveFiles(double *put1, double *put3, double *put5, double *put7, double *put9, double *put11, double *obj_out, double *put31, double *put32)
{
    int i, j;
    //pi1 and pi2 save
    if(Cid == 2) //tow components
    {
        for(i=0;i<intx;i++)
        {
            put1[i+intx*0] = p->pi1[i];
            put1[i+intx*1] = p->pi2[i];
        }
    }else{
        for(i=0;i<intx;i++) put1[i+intx*0] = p->pi1[i];
    }
    // pure tumor expression save.
    for(j=0;j<nG;j++)
    {
        for(i=0;i<intx;i++)
        {
            put3[j*intx+i] =  CD[i][j];
			put31[j*intx+i] = p->mle[j][i+intx];
			if(Cid == 2) put32[j*intx+i] = p->mle[j][i+2*intx];
        }
    }
    //mu and sigma save for T
    for(j=0;j<nG;j++)
    {
        for(i=0;i<iteration;i++)
        {
            put5[j*iteration+i] = avgparT[i][j];
        }
    }

    for(j=0;j<nG;j++)
    {
        for(i=0;i<iteration;i++)
        {
            put7[j*iteration+i] = sigparT[i][j];
        }
    }

    for(j=0;j<intx;j++)
    {
        for(i=0;i<iteration;i++)
        {
            put9[j*iteration+i] = tmppi1[i][j];
        }
    }


    for(j=0;j<intx;j++)
    {
        for(i=0;i<iteration;i++)
        {
            put11[j*iteration+i] = tmppi2[i][j];
        }
    }
    for(i=0;i<iteration;i++)
    {
        obj_out[i] = p->obj[i];
    }
}



/************************************************************************************
  * get tumor function
*************************************************************************************/
/*Tumor Component Update*/
void gettumor(int genes, int h) // option h = 1 for one component; h = 2 for two
{
    double mu, sigma,upp,low,upps;
    
    upp = 33;
    if(nHavepi == 1) {
        upps = 100;
    }else{
        upps = 25;
    }
    low = 0;
    double obj_old, obj_new;
    if(h == 1)
    {
        obj_old = mint(genes, h, p->Tavg[genes]);
        mu = tmin_y(low, upp, genes, h, mint, 0.001);
        obj_new = mint(genes, h, mu);
        //if((obj_new < obj_old) || iteration1 == 0)
        p->Tavg[genes] = mu;
        
        obj_old = tf_y(genes, p->Tavg[genes], p->Tsigma[genes]);
        sigma =  tmin_y2(0.0001, upps, genes, p->Tavg[genes], tf_y, 0.0001);
        obj_new = tf_y(genes, p->Tavg[genes], sigma);
        //if((obj_new < obj_old) || iteration1 == 0)
        p->Tsigma[genes] = sigma;
        
        //Rprintf("obj_old is %lf %d\t", obj_old, genes);
        //Rprintf("obj_new is %lf %d\t", obj_new, genes);
        
    }else{
        obj_old = mint(genes, h, p->Tavg[genes]);
        mu = tmin_y(low, upp, genes,h, mint, 0.001);
        obj_new = mint(genes, h, mu);
        //if((obj_new < obj_old) || iteration1 == 0)
        p->Tavg[genes] = mu;
        
        obj_old = tf_y2(genes, p->Tavg[genes], p->Tsigma[genes]);
        sigma =  tmin_y2(0.0001, upps, genes, p->Tavg[genes], tf_y2, 0.0001);
        obj_new = tf_y2(genes, p->Tavg[genes], sigma);
        //if((obj_new < obj_old) || iteration1 == 0)
        p->Tsigma[genes] = sigma;
    }
    
}

double tf_y(int genes, double mu, double sigma)
{
    double new_like;
    double tmp;
    int i;
    new_like=0.0;

    for (i=0; i<intx; i++){
        tmp = ft_y(FD[i+fNorm][genes],p->Navg1[genes], mu, p->Nsigma1[genes], sigma, p->pi1[i], 0.0);
        new_like+= tmp;
    }

    return (new_like*(-1.0));
}

double tf_y2(int genes, double mu, double sigma)
{
    double new_like;
    double tmp;
    int i;
    new_like=0.0;
    for (i=0; i<intx; i++){
        tmp = ft_y2(FD[i+fNorm][genes], p->Navg1[genes], p->Navg2[genes], mu, p->Nsigma1[genes], p->Nsigma2[genes], sigma, p->pi1[i], p->pi2[i]);
        new_like+= tmp;
    }

    return (new_like*(-1.0));
}


double mint(int genes, int h, double mu) // opt h =1:one component; h =2:two component
{
  double sigma;
  double tmp,tmp0;
  double upps;
  if(nHavepi == 1) {
      upps = 100;
    }else{
        upps = 25;
    }
  if(h == 1)
  {
    sigma =  tmin_y2(0.0001, upps, genes, mu, tf_y, 0.0001);
    tmp = tf_y(genes, mu, sigma);
    tmp0 = tf_y(genes, mu, p->Tsigma[genes]);
    if(tmp0 < tmp) {
      sigma = p->Tsigma[genes];
      tmp = tmp0;
    }
    
  }else{
    sigma =  tmin_y2(0.0001, upps, genes, mu, tf_y2, 0.0001);
    tmp = tf_y2(genes, mu, sigma);
    tmp0 = tf_y2(genes, mu, p->Tsigma[genes]);
    if(tmp0 < tmp) {
      sigma = p->Tsigma[genes];
      tmp = tmp0;
    }
    
  }

  return (tmp);
}




double tmin_y(double ax, double bx, int genes, int h, double (*f)(int, int, double), double tol)
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
    fx = (*f)(genes, h, x);
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

        fu = (*f)(genes, h, u);

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




double tmin_y2(double ax, double bx, int genes, double mu, double (*f)(int, double, double), double tol)
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
    fx = (*f)(genes, mu, x);
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

        fu = (*f)(genes, mu, u);

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


/*************************************************************************************
 *a function to search the mode for better numerical integration
 *************************************************************************************/
double nitg_ft_y(double t, double y, double mung, double mutg, double sng, double stg, double pi1){

    double tmp;
    tmp = -0.5*log(sng)-0.5*log(stg) -log(y-t)-log(t);
    tmp = tmp-0.5*pow(log2(t)-(log2(pi1) + mung), 2.0) /sng;
    tmp = tmp-0.5*pow(log2(y-t)-(log2(1-pi1) + mutg), 2.0)/stg;
    tmp = exp(tmp);
    return -tmp;
}

double min_nitg_ft_y(double ax, double bx, double y, double mung, double mutg, double sng, double stg, double pi1, double (*f)(double, double, double, double, double, double,double), double tol)
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
    fx = (*f)(x, y, mung, mutg, sng, stg, pi1);
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

        fu = (*f)(u, y, mung, mutg, sng, stg, pi1);
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







/************************************************************************************
 * getpi function
 *************************************************************************************/



/*optimize pi for updating*/
//for calculating pi1


/* For update ratio calculation, we need to fix position at the same number*/
double ft_y(double y, double mung, double mutg, double sng, double stg, double pi1, double pi2)
{
  int i;
  double tmode;
  double ltmode, rtmode;
  double rtval=0.0;
  double tmp;
  double tmp_pos[integ];
  double cutd;
  int integ1, integ2;
  double rtval1=0.0, rtval2=0.0;
  
  integ1 = integ/6; //small bin size
  integ2 = integ - integ1; //large bin size
  double tmp_pos1[integ1], tmp_pos2[integ2];
  
  //this step is for finding the mode to strengthen numerical integration
  if(fabs(mutg - mung) > 4.0){
    tmode =  min_nitg_ft_y(0.1, y-0.1, y, mung, mutg, sng, stg, pi1, nitg_ft_y, 0.01);
    ltmode = tmode;
    rtmode = y - tmode;
    //we cut the interval to two segmentations
    if(ltmode - rtmode > 0){
      cutd = max(y/6, tmode - 1.5*(y - tmode));
      for(i=0;i<integ1;i++) tmp_pos1[i] = (cutd)/(double)integ1*(i+0.5);
      for(i=0;i<integ1;i++)
      {
        tmp = -0.5*log(sng)-0.5*log(stg) -log(y-tmp_pos1[i])-log(tmp_pos1[i]);
        tmp = tmp-0.5*pow(log2(tmp_pos1[i])-(log2(pi1) + mung), 2.0) /sng;
        tmp = tmp-0.5*pow(log2(y - tmp_pos1[i])-(log2(1-pi1-pi2) + mutg), 2.0)/stg;
        tmp = exp(tmp);
        rtval1 +=tmp;
      }
      if(rtval1<=0) rtval1 = 1e-313;
      rtval1 = rtval1/(double)integ1*cutd;
      for(i=0;i<integ2;i++) tmp_pos2[i] = (y - cutd)/(double)integ2*(i+0.5) + cutd;

      for(i=0;i<integ2;i++)
      {
        tmp = -0.5*log(sng)-0.5*log(stg) -log(y-tmp_pos2[i])-log(tmp_pos2[i]);
        tmp = tmp-0.5*pow(log2(tmp_pos2[i])-(log2(pi1) + mung), 2.0) /sng;
        tmp = tmp-0.5*pow(log2(y - tmp_pos2[i])-(log2(1-pi1-pi2) + mutg), 2.0)/stg;
        tmp = exp(tmp);
        rtval2 +=tmp;
      }
      if(rtval2<=0) rtval2 = 1e-313;
      rtval2 = rtval2/(double)integ2*(y - cutd);
      return(log(rtval1 + rtval2));
    }else{
      cutd = min(tmode + 1.5*tmode, y/6.0*5.0);
      //
      for(i=0;i<integ2;i++) tmp_pos2[i] = (cutd)/(double)integ2*(i+0.5);
      for(i=0;i<integ2;i++)
      {
        tmp = -0.5*log(sng)-0.5*log(stg) -log(y-tmp_pos2[i])-log(tmp_pos2[i]);
        tmp = tmp-0.5*pow(log2(tmp_pos2[i])-(log2(pi1) + mung), 2.0) /sng;
        tmp = tmp-0.5*pow(log2(y - tmp_pos2[i])-(log2(1-pi1-pi2) + mutg), 2.0)/stg;
        tmp = exp(tmp);
        rtval1 +=tmp;
      }
      if(rtval1<=0) rtval1 = 1e-313;
      rtval1 = rtval1/(double)integ2*cutd;
      for(i=0;i<integ1;i++) tmp_pos1[i] = (y - cutd)/(double)integ1*(i+0.5) + cutd;
      for(i=0;i<integ1;i++)
      {
        tmp = -0.5*log(sng)-0.5*log(stg) -log(y-tmp_pos1[i])-log(tmp_pos1[i]);
        tmp = tmp-0.5*pow(log2(tmp_pos1[i])-(log2(pi1) + mung), 2.0) /sng;
        tmp = tmp-0.5*pow(log2(y - tmp_pos1[i])-(log2(1-pi1-pi2) + mutg), 2.0)/stg;
        tmp = exp(tmp);
        rtval2 +=tmp;
      }
      if(rtval2<=0) rtval2 = 1e-313;
      rtval2 = rtval2/(double)integ1*(y - cutd);
      return(log(rtval1 + rtval2));
      
    }
    
  }else{
    for(i=0;i<integ;i++) tmp_pos[i] = (y)/(double)integ*(i+0.5);
    for(i=0;i<integ;i++)
    {
      tmp= -0.5*log(sng)-0.5*log(stg) -log(y-tmp_pos[i])-log(tmp_pos[i]);
      tmp= tmp-0.5*pow(log2(tmp_pos[i])-(log2(pi1) + mung), 2.0) /sng;
      tmp= tmp-0.5*pow(log2(y - tmp_pos[i])-(log2(1-pi1-pi2) + mutg), 2.0)/stg;
      tmp= exp(tmp);
      rtval +=tmp;
    }
    if(rtval<=0) rtval=1e-313;
    return(log(rtval/(double)integ *y));
  }

}


double ft_y_SC(double y, double mung, double mutg, double sng, double stg, double pi1, double pi2)
{
    int i;
    double rtval;
    double tmp;
    double tmp_pos[integ];

    rtval=0;
    for(i=0;i<integ;i++) tmp_pos[i] = (y)/(double)integ*(i+0.5);

    for(i=0;i<integ;i++)
    {
        tmp= -log(y-tmp_pos[i])-log(tmp_pos[i]);
        tmp= tmp-0.5*pow(log2(tmp_pos[i])-(log2(pi1) + mung), 2.0) /sng;
        tmp= tmp-0.5*pow(log2(y - tmp_pos[i])-(log2(1-pi1-pi2) + mutg), 2.0)/stg;

        tmp= exp(tmp);

        rtval +=tmp;
    }

    if(rtval<=0)
        rtval=1e-313;

    return (log(rtval/(double)integ *y));
}

double pf_y(int samp, double pi1)
{
    double new_like;
    double tmp;
    int i;
    new_like=0.0;
    for (i=0; i<nG; i++){
        tmp = ft_y(FD[samp+fNorm][i],p->Navg1[i], p->Tavg[i], p->Nsigma1[i], p->Tsigma[i], pi1, 0);
        new_like+= tmp;
    }

    return (new_like*(-1.0));
}

double pmin_y(double ax, double bx, int samp, double (*f)(int, double), double tol)
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
    fx = (*f)(samp, x);
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

        fu = (*f)(samp, u);

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


//another function to get pi
void getpi(int samp, int h)  	// option h = 1 for 1 component, 2 for two component
{
    double pii1, pii2;
    double upp;
    pii1 = 0.0;
    pii2= 0.0;
    
    double obj_old, obj_new;
    if(h == 1)
    {
        obj_old = pf_y(samp, p->pi1[samp]);
        pii1 =  pmin_y(0.01, 0.99, samp, pf_y, 0.0001);
        obj_new = pf_y(samp, pii1);
        //if((obj_new < obj_old) || iteration1 == 0)
        p->pi1[samp] = pii1;
        
    }else{ //two component
        obj_old = pf_y2(samp, p->pi1[samp], p->pi2[samp]);
        pii2 = pmin_y(0.01, 0.99, samp, minpi, 0.0001);
        upp = 1 - pii2;
        pii1 = pmin_y2(0.01, upp, samp, pii2, pf_y2, 0.0001);
        obj_new = pf_y2(samp, pii1, pii2);
        //if((obj_new < obj_old) || iteration1 == 0){
        p->pi1[samp] = pii1;
        p->pi2[samp] = pii2;
        //}
    }
    
}

// another function to get spike in pi
void getspikeinpi(int samp)
{
	double pii1;
	double obj_old, obj_new;
	
	obj_old = pf_y(samp, p->pi1[samp]);
	pii1 = 0.99;
	obj_new = pf_y(samp, pii1);
	p->pi1[samp] = pii1;
}

//another function to get pi given piT
void getpiT(int samp)  	// for two component given piT
{
    double pii1, pii2, piiT;
    double upp;
    pii1 = 0.0;
    pii2= 0.0;
    double obj_old, obj_new;
    piiT = p->piT[samp];
    obj_old = pf_yT(samp, p->pi1[samp], piiT);
    //two component
    upp = 1 - piiT;
    pii1 = pmin_y2(0.01, upp, samp, piiT, pf_yT, 0.0001);
    pii2 = 1 - piiT - pii1;
    obj_new = pf_yT(samp, pii1, piiT);
    //if((obj_new < obj_old) || iteration1 == 0){
    p->pi1[samp] = pii1;
    p->pi2[samp] = pii2;
    //}
    
}



// outer integration
double ft_y2(double y, double mung1, double mung2, double mutg, double sng1, double sng2, double stg, double pi1, double pi2)
{
    int i;
    double rtval;
    double tmp;
    double res;
    double tmp_pos[integ];

    rtval=0;
    for(i=0;i<integ;i++) tmp_pos[i] = (y)/(double)integ*(i+0.5);

    for(i=0;i<integ;i++)
    {
        res =y-tmp_pos[i];
        tmp = -log(tmp_pos[i])-0.5*pow(log2(tmp_pos[i])-(log2(pi2) + mung2), 2.0) /sng2;
        tmp += ft_y_SC(res, mung1, mutg, sng1, stg, pi1, pi2);
        tmp = exp(tmp);
        rtval +=tmp;
    }
    rtval=rtval/sqrt(sng1*sng2*stg);

    if(rtval<=0)
        rtval=1e-313;

    return (log(rtval/(double)integ*y));
}


double pf_y2(int samp, double pi1, double pi2)
{
    double new_like;
    double tmp;
    int i;
    new_like=0.0;

    for (i=0; i<nG; i++){
        tmp = ft_y2(FD[samp+fNorm][i], p->Navg1[i], p->Navg2[i], p->Tavg[i], p->Nsigma1[i], p->Nsigma2[i], p->Tsigma[i], pi1, pi2);
        new_like+= tmp;
    }


    return (new_like*(-1.0));
}


//the joint likelihood function given piiT
double pf_yT(int samp, double pi1, double piT)
{
    double new_like;
    double tmp;
    double pi2;
    int i;
    new_like=0.0;
    pi2 = 1 - piT -pi1;

    for (i=0; i<nG; i++){
        tmp = ft_y2(FD[samp+fNorm][i], p->Navg1[i], p->Navg2[i], p->Tavg[i], p->Nsigma1[i], p->Nsigma2[i], p->Tsigma[i], pi1, pi2);
        new_like+= tmp;
    }
    return (new_like*(-1.0));
}

double pmin_y2(double ax, double bx, int samp, double pi2, double (*f)(int, double, double), double tol)
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
    fx = (*f)(samp, x, pi2);
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

        fu = (*f)(samp, u, pi2);

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

//get pi1 function value given pi2

double minpi(int samp, double pi2)
{
    double pi1;
    double tmp;
    double upp = 1 - pi2;
    pi1 =  pmin_y2(0.01, upp, samp, pi2, pf_y2, 0.0001);
    tmp = pf_y2(samp, pi1, pi2);
    return (tmp);
}

/************************************************************************************
 *mle optimization function
 ************************************************************************************/
//optimize function
void getmle(int opt, int opt2, int h)
{
    double tmp, tmp1, tmp2;
    double upp;
    double y_sum, pi_sum;
    if(h == 1) //one component
    {
        upp = (FD[opt+fNorm][opt2] - 1 + p->pi1[opt])/p->pi1[opt];
        y_sum = FD[opt+fNorm][opt2];
        tmp = fmin1(1, upp, opt2, opt, y_sum, 1, lf1, 0.0001);
        p->mle[opt2][opt+intx] = tmp;
        p->mle[opt2][opt] = (FD[opt+fNorm][opt2] - p->pi1[opt]*tmp )/(1 - p->pi1[opt]);
    }else{
        upp = (FD[opt+fNorm][opt2] - 1 + p->pi1[opt] + p->pi2[opt])/p->pi2[opt];
        tmp2 = fmin2(1, upp, opt2, opt, lf2, 0.0001);
        p->mle[opt2][opt+2*intx] = tmp2;
        upp = (FD[opt+fNorm][opt2] - 1 + p->pi1[opt] + p->pi2[opt] - p->pi2[opt]*tmp2)/p->pi1[opt];
        y_sum = FD[opt+fNorm][opt2] - tmp2*p->pi2[opt];
        pi_sum = 1 - p->pi2[opt];
        tmp1 = fmin1(1, upp, opt2, opt, y_sum, pi_sum, lf1, 0.0001);
        p->mle[opt2][opt+intx] = tmp1;
        p->mle[opt2][opt] = (FD[opt+fNorm][opt2] - p->pi1[opt]*tmp1 - p->pi2[opt]*tmp2)/(1 - p->pi1[opt] - p->pi2[opt]);
    }
}
// kernel optimize function for just component 1 normal and tumor

double lf1(int opt, int opt2,  double y_sum, double pi_sum, double nval1)
{
    double new_like;
    double tmp;
    new_like=0.0;

    new_like = -log(nval1) - 0.5*log(p->Nsigma1[opt2]) - 0.5*(log2(nval1) - p->Navg1[opt2])*(log2(nval1) - p->Navg1[opt2])/p->Nsigma1[opt2];
    tmp= (y_sum - p->pi1[opt]*nval1)/(pi_sum - p->pi1[opt]);
    new_like += -log(tmp) -0.5*log(p->Tsigma[opt2]) - 0.5*(log2(tmp) - p->Tavg[opt2])*(log2(tmp) - p->Tavg[opt2])/p->Tsigma[opt2];
    return (new_like*(-1.0));
}

double fmin1(double ax, double bx,int iG, int iS, double y_sum, double pi_sum, double (*f)(int, int, double, double, double), double tol)
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
    fx = (*f)(iS, iG, y_sum, pi_sum, x);
    //Rprintf("optimized parameters are iS: %d, iG: %d, y_sum: %lf, pi_sum: %lf, x: %lf, fx:%lf \n",iS, iG, y_sum, pi_sum, x, fx);

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

        fu = (*f)(iS, iG, y_sum, pi_sum, u);
        //Rprintf("optimized parameters are iS: %d, iG: %d, y_sum: %lf, pi_sum: %lf, x: %lf, fx:%lf \n",iS, iG, y_sum, pi_sum, x, fu);

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


// wrapped optimize function for component 2

double lf2(int opt, int opt2, double nval2)
{
    double tmp;
    double upp;
    double y_sum, pi_sum;
    double nval1;
    double new_like;
    new_like=0.0;
    upp = (FD[opt+fNorm][opt2] - 1 + p->pi1[opt] + p->pi2[opt] - p->pi2[opt]*nval2)/p->pi1[opt];
    y_sum = FD[opt+fNorm][opt2] - nval2*p->pi2[opt];
    pi_sum = 1 - p->pi2[opt];
    nval1 = fmin1(1, upp, opt2, opt, y_sum, pi_sum, lf1, 0.0001);

    new_like = -log(nval1) - 0.5*log(p->Nsigma1[opt2]) - 0.5*(log2(nval1) - p->Navg1[opt2])*(log2(nval1) - p->Navg1[opt2])/p->Nsigma1[opt2];
    new_like += -log(nval2) - 0.5*log(p->Nsigma2[opt2]) - 0.5*(log2(nval2) - p->Navg2[opt2])*(log2(nval2) - p->Navg2[opt2])/p->Nsigma2[opt2];
    tmp= (FD[opt+fNorm][opt2] - p->pi1[opt]*nval1 - p->pi2[opt]*nval2)/(1 - p->pi1[opt] - p->pi2[opt]);
    new_like += -log(tmp) -0.5*log(p->Tsigma[opt2]) - 0.5*(log2(tmp) - p->Tavg[opt2])*(log2(tmp) - p->Tavg[opt2])/p->Tsigma[opt2];
    return (new_like*(-1.0));
}


double fmin2(double ax, double bx,int iG, int iS, double (*f)(int, int, double), double tol)
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
    fx = (*f)(iS, iG, x);
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

        fu = (*f)(iS, iG, u);

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




