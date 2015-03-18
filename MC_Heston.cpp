//***** ORIGINAL CODE OBTAINED FROM: http://www.math.nyu.edu/fellows_fin_math/gatheral/h2.c    *****//
// CODE HAS BEEN ALTERED TO THE MILSTEIN DISCRETISATION OF VOL SDE, TO ALLOW FOR JUMPS AND TO PROVIDE CONFIDENCE INTERVALS //

/* Official Heston Monte Carlo C - Code,  Fall 2002
   with impl. vol. calculator (due to Chamadia) */
/* compile:    cc -lm h2.c */

#include <stdio.h>
#include <math.h>
#include <ctime>
#include "stdlib.h"

#define Pi 3.14159265358979323846
#define max(A,B) ((A) > (B) ? (A) : (B) )

/* To generate two independent normal variables */

void boxmuller(double * B1, double * B2);

/* To generate the Poisson random variable */

void poisson(int*,double,double);

/* For the impl. vol. calculator */

double N(double z);
double BSPrice(double S, double X, double r, double sigma, double time);
double imp_vol(double spot, double strike, double dt, double target_price);


FILE *fp;


int main()
{
	double v,mirror_v, v0, vbar, lambda, eta, rho, Kmin, Kmax, K, dK, T, dt, lambdaJ, alpha, delta;
        	double S, mirror_S, S0, drift, sqrt_one_minus_rho_square, payoff, standard_deviation,
        	impl_vol;

        double B1, B2, mirror_B1, mirror_B2;

	double *expectation, *second_moment;
	
	int N, n, M, jump_flag, flag, MC_flag, std_dev_flag, i, j, k;



	
	
        /* ====== Hard-wire parameter-values ... */
        
        T=1; S0=1; v0=0.04; vbar=0.04; lambda=1.15; eta=0.39; rho = -0.64;
        
		lambdaJ = 0.1308; delta = 0.0967; 
		alpha =-0.1151;

		if (jump_flag==1) 
			drift = 0;// -lambdaJ*(exp(alpha+delta*delta*0.5)-1);
		else
			drift=0;

        Kmin = 0.8; Kmax = 1.2;  M=5;
        
        n =150;   N =30000;     /* # time steps, # number paths   */
        
		jump_flag = 0;   /*	0 = SV model, 1 = SVJ model */

        flag = 0;        /*  0 = absorbing, 1 = reflecting */
        
        MC_flag = 0;      /* use antithetic MC   */
        
        std_dev_flag = 0;    /* do not compute std.dev  */

       

        
        /*  ====== ... or enter them by hand:  
        
        printf("Enter the maturity: \n");
	scanf("%lf",&T);
	
	printf("Enter the following parameters: S0, v, vbar, lambda, eta, rho:\n");
        scanf("%lf %lf %lf %lf %lf %lf", &S0, &v0, &vbar, &lambda, &eta, &rho);

     
	printf("Enter the minimal strike: \n");
        scanf("%lf",&Kmin); 

	printf("Enter the maximal strike: \n");
        scanf("%lf",&Kmax); 

	printf("Enter the number of strikes you want: \n");
        scanf("%d", &M);   					*/ 
              
	printf("How many time steps  ?\n");
	scanf("%d",&n);
	
	printf("How many simulations ?\n");  	
        scanf("%d",&N);

	printf("Do you want no jumps (0) or to include jumps (1) ?\n");
        scanf("%d",&jump_flag); 
        
	printf("Do you want absorbing (0) or reflecting (1) ?\n");
        scanf("%d",&flag); 
        
	printf("Do you want antithetic (0) or standard MC (1) ?\n");
        scanf("%d",&MC_flag); 
        
        printf("Do you want to compute standard deviation (0=no, 1=yes) ? \n");
        scanf("%d",&std_dev_flag);
        

	/* ======= Compute dK, dt (strike resp. time-steps  ...) and other aux. stuff */
 
   
	if (M>1) 
	    dK=(Kmax-Kmin)/(M-1);
	else 
	    dK=0;

	dt =T/n;
	
	sqrt_one_minus_rho_square = sqrt(1-rho*rho);

       
	
    

	/* ======= Memory allocation and other initializations  */

     

        expectation = (double *) malloc( M*sizeof(double));
        
        if (std_dev_flag != 0)      /* allocate only if you want to compute std.deviation */
           second_moment = (double * ) malloc( M*sizeof(double));
        
 	if (expectation == NULL || ((std_dev_flag != 0) && second_moment == NULL)) 
	{ printf("memory allocation problem ");
	  exit(1);
        }
        
        for (j=0;j<M;j++)
	{	expectation[j]=0; 
		
		if (std_dev_flag != 0)
		   second_moment[j]=0; 
        }		
 	//srand(time(0));
	srand(10);


	/* ======= Beginning of the main loop over paths */

	for (k=0;k<N;k++)
	{

	  /* === local initialization */
	
          v=v0;
	  S=S0;
          
          mirror_v=v0;
          mirror_S=S0;

          K=Kmin;

   	  /* === discretization of the stochastics processes */

          for (i=0;i<n;i++)
	  {
	    
	    
	    boxmuller(&B1,&B2);
   /* get two indep. normals */ 

	    B2=rho*B1 + sqrt_one_minus_rho_square*B2;     /* correlation */

	 
	                 /* compute increments of stoch. processes */
                  
	    S+=drift*S*dt+sqrt(v*dt)*S*B1;
	    
//		Euler:  v+=lambda*(vbar-v)*dt+eta*sqrt(v*dt)*B2;
	    v+=lambda*(vbar-v)*dt+eta*sqrt(v*dt)*B2+eta*eta*dt*(B2*B2-1)/4;  //Milstein
                       /* absorbing or reflecting depending on the flag */
		 
            if (v<0) 
               v=-flag*v;
            
               
            if (MC_flag == 0)     /* create anti-paths ... */
            { mirror_B1 = - B1; mirror_B2 = - B2;
              mirror_S+=drift*S*dt+sqrt(mirror_v*dt)*mirror_S*mirror_B1;               

//		Euler:  mirror_v+=lambda*(vbar-mirror_v)*dt+eta*sqrt(mirror_v*dt)*mirror_B2;
		mirror_v+=lambda*(vbar-mirror_v)*dt+eta*sqrt(mirror_v*dt)*mirror_B2
								+eta*eta*dt*(mirror_B2*mirror_B2-1)/4; //Milstein

             if (mirror_v<0) 
                 mirror_v=-flag*mirror_v;
            }
          }

          /* === At maturity, we compute the pay-offs for the different strikes 
             and add to the expectations.
 */ 
/* Include jumps */
		if (jump_flag !=0)
		{
			int P;
			poisson(&P,lambdaJ,T);
			for(i=0;i<P;i=i+2)
			{
				boxmuller(&B1,&B2);
				S*=exp(alpha+delta*B1);
				if (i+1<P) S*=exp(alpha+delta*B2);
			}
			if(MC_flag==0)
			{
				poisson(&P,lambdaJ,T);
				for(i=0;i<P;i=i+2)
				{
					boxmuller(&B1,&B2);
					mirror_S*=exp(alpha+delta*B1);
					if (i+1<P)	mirror_S*=exp(alpha+delta*B2);
				}
			}			
		}
      
          for (j=0;j<M;j++)
	  {
	     if (MC_flag != 0)
	        payoff = max(S-K,0);  
	     else 
	        payoff = (max(S-K,0) + max(mirror_S-K,0))/2;   

	     expectation[j]+=payoff;
	     if (std_dev_flag != 0) 
	       second_moment[j] += payoff*payoff;        /* so that we can estimate variance of MC  */
	     
             K+=dK;
	  }

	}

    
	
    /* ===== We write the results to a data file and print them on the screen */

    if ((fp = fopen("data.dat","w"))==NULL) 
    {
         puts("cannot open file\n");
         exit(1);
    } 
    
    K=Kmin; 
    
    printf("Heston 2002: #timesteps = %ld,  #paths = %ld \n",n,N);

   

   if (MC_flag != 0)  printf ("Standard MC, "); 

   else  printf("Antithetic MC"); 

   if (flag == 0) printf (" absorbing \n");

   else printf("reflecting \n");

   

   
    for (j=0;j<M;j++)

   {   

       expectation[j]/=N;
   /* divide total sum by number of paths and get expectation */
        impl_vol = imp_vol(S0, K, T, expectation[j]);

        if (std_dev_flag != 0) 
        {
           standard_deviation = sqrt((second_moment[j])/N - (expectation[j])*(expectation[j]));
           fprintf(fp, "strike: %lf  value: %lf  std.dev: %lf CI: ( %lf, %lf ) \n", K, expectation[j],standard_deviation,
									expectation[j]-2*standard_deviation/sqrt((double)N), expectation[j]+2*standard_deviation/sqrt((double)N)); 
           printf( "strike: %lf  value: %lf  std.dev: %lf CI: ( %lf, %lf ) \n", K, expectation[j],standard_deviation,
									expectation[j]-2*standard_deviation/sqrt((double)N), expectation[j]+2*standard_deviation/sqrt((double)N)); 
        }
        else
        {  fprintf(fp, "strike: %lf  value: %lf <=> impl. vol %lf \n", K, expectation[j], impl_vol);  
           printf("strike: %lf  value: %lf <=> impl.vol %lf \n", K, expectation[j], impl_vol);  
        }
     	
     	K+=dK;	
    }



    fclose (fp);

    return(0);
}

/********************************************************************/
void boxmuller(double * B1, double * B2)

/********************************************************************/

{ /* thanks to G.Ciresi for the reference: Num. Recipes in C, 1992  */
     
  double v1,v2,rsq,fac;

 
    
  do 
  { 
    v1 = (rand() / (double) 32767)*2-1; 
    v2 = (rand() / (double) 32767)*2-1; 
    
    rsq = v1*v1 + v2*v2;
  }
    while (rsq >= 1.0 || rsq == 0.0);
    
  fac = sqrt ( -2.0 * log(rsq)/rsq);
   
  *B1 = v1*fac;
  *B2 = v2*fac;


  
} 

void poisson(int *P, double l, double T)
{
	double t = 0, U = (rand() / (double) RAND_MAX);

	*P=0;
	t-=log(U)/l;
	while (t<T)
	{
		*P+=1;
		U = (rand() / (double) RAND_MAX);
		t-=log(U)/l;
	}
}
  
double BSPrice(double S, double X, double r, double sigma, double time)
{ 
    double time_sqrt = sqrt(time);
    double d1 = (log(S/X)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt; 
    double d2 = d1-(sigma*time_sqrt);
    double c = S * N(d1) - X * exp(-r*time) * N(d2);
    return c;
}

double N(double z) {
    double b1 =  0.31938153; 
    double b2 = -0.356563782; 
    double b3 =  1.781477937;
    double b4 = -1.821255978;
    double b5 =  1.330274429; 
    double p  =  0.2316419; 
    double c2 =  0.3989423; 

    double a, t, b, n;	    
	
    if (z >  6.0) { return 1.0; }; // this guards against overflow 
    if (z < -6.0) { return 0.0; };
    
    a=fabs(z); 
    t = 1.0/(1.0+a*p); 
    b = c2*exp((-z)*(z/2.0)); 
    n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t; 
    
    n = 1.0-b*n; 
    if ( z < 0.0 ) n = 1.0 - n; 
    return n;
}


double imp_vol(double spot, double strike, double dt, double target_price)
{
  int i = 1;
  double vol_up = 1.50;
  double vol_down = 0;
  double vol_guess;
  double tolerance = target_price * 0.001;
  double price;

  int max_iteration = 100;
  	
  while((i++) <= max_iteration)
    {
      vol_guess = (vol_up +vol_down)/2;
      price = BSPrice(spot, strike, 0.0, vol_guess, dt);   
 
      if(fabs(price-target_price) < tolerance)
        break;

      if (target_price > price)
        vol_down = vol_guess;
      else
        vol_up = vol_guess;
    }
  
  return vol_guess;
}
        
         






   