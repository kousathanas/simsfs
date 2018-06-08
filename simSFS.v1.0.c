#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>

#define maxnd 10000
#define H 0.5

/*
  PROGRAM: simSFS
  AUTHOR: Athanasios Kousathanas
  VERSION: 1.0

  Generate allele frequency vectors (AFVs) using Wright-Fisher transition matrix methods and sample
  to obtain site frequency spectra (SFS). Code implements methods described by Keightley and Eyre-Walker(2007).

  Features:
  Demography: 1-epoch or 2-epoch population changes
  Selection: only single s

  Compilation:
  use makefile or compile as
  gcc -O2 -o simSFS simSFS.v1.0.c -lm -lgsl -lgslcblas -w
*/


////////////////////////////////////////
//calculate allele frequency change, delta(q)
////////////////////////////////////////
double selfn(double q, double s, double h)
{
  double num, denom;
  num = s*q*(1-q)*(q + h*(1-2*q));
  denom = 1 + 2.0*h*s*q*(1-q) + s*q*q;
  if (denom== 0) denom = .00000000000001;
  return( q + (num/denom));
}
////////////////////////////////////////
//initialise rows of matrix (used by setupmatrix)
////////////////////////////////////////
setuprow(double **a, double p, int row, int n){
  double z;         /*probability of zero failures*/
  double x, y, temp1, temp2;
  int j;

  temp1 = 1.0 - p;
  if (temp1 == 0.0)  temp1 = .000000000001;            /*prevent zero divide*/
  temp1 = p/temp1;       /*First prevent negative log*/
  if (temp1 <=0)   temp1 = .000000000001;
  x = log(temp1);    /*used for later multiplication of series*/
  temp2 = 1.0-p;
  if (temp2<=0)   temp2 = .000000000001;
  z = (double)n*log(temp2);
  if (z > 150)   z = 150;
  a[row][0] = exp(z);
  y = a[row][0];        /*add up for later*/
  for (j = 1;j<=n; j++)
    {

      z = z + log(((double)n + 1.0 - (double)j)/(double)j) + x;
      /*             cumulative binomial coeff.*/
      a[row][j] = exp(z);
      y = y + a[row][j];          /*add up*/
    }
  if (y == 0.0)   y = .0000000000001;
  for (j=0; j<=n; j++) a[row][j] = a[row][j]/y;
}
////////////////////////////////////////
//initialise transition matrix
////////////////////////////////////////
setupmatrix(double **a, int n1, int n2, double s, double h){
  int i;
  double p;


  for (i = 1; i<=n1 - 1; i++)             /*do variable rows*/
    {
      p = (double)i/(double)n1;
      p = selfn(p, s,h);
      setuprow(a, p, i, n2);
    }
  for (i = 0; i<=n2; i++)
    {              /*do constant edges*/
      a[0][i] = 0;
      a[n1][i] = 0;
    }
  a[0][0] = 1.0;                          /*do absorbing corners*/
  a[n1][n2] = 1.0;

}
////////////////////////////////////////
//Perform matrix inversion
////////////////////////////////////////
matrixinvert(int N, double **a,double *inverter){

  int i, j, lotkin_signum;

  gsl_matrix  *lotkin_a,*lotkin_inv;
  gsl_vector *x, *lotkin_b, *lotkin_x;

  gsl_permutation *lotkin_perm;

  /* allocate a, x, b */
  lotkin_a = gsl_matrix_alloc(N-1, N-1);
  lotkin_inv = gsl_matrix_alloc(N-1, N-1);
  x = gsl_vector_alloc(N-1);
  lotkin_b = gsl_vector_alloc(N-1);

  lotkin_x = gsl_vector_alloc(N-1);

  gsl_matrix *identity=gsl_matrix_alloc(N-1, N-1);


  for (i = 0; i<=N-2; i++)
    {

      for (j = 0; j<=N-2; j++)
	{
          gsl_matrix_set( identity, i, i,1);
	}

    }

  for (i = 0; i<=N-2; i++)
    {
      for (j = 0; j<=N-2; j++)
	{
	  gsl_matrix_set(lotkin_a,i,j,gsl_matrix_get(identity,i,j)-a[i+1][j+1]);
	}

    }


  /* LU decomposition and forward&backward substition */

  //gsl_blas_dgemv(CblasNoTrans, 1.0, lotkin_a, x, 0.0, lotkin_b);

  lotkin_perm = gsl_permutation_alloc(N-1);
  gsl_linalg_LU_decomp(lotkin_a, lotkin_perm, &lotkin_signum);
  gsl_linalg_LU_solve(lotkin_a, lotkin_perm, lotkin_b, lotkin_x);
  gsl_linalg_LU_invert (lotkin_a,lotkin_perm, lotkin_inv);

  /*apparently the inversion adds +1 to the first element of the vector
    so I substract it*/

  gsl_matrix_set(lotkin_inv,0,0,gsl_matrix_get(lotkin_inv,0,0)-1);



  for (i=0;i<=N-2;i++){

    inverter[i+1]=gsl_matrix_get(lotkin_inv,0,i);

  }


  gsl_matrix_free(lotkin_a);
  gsl_matrix_free(lotkin_inv);
  gsl_vector_free(x);
  gsl_vector_free(lotkin_b);
  gsl_vector_free(lotkin_x);
  gsl_matrix_free(identity);
  gsl_permutation_free(lotkin_perm);

}
////////////////////////////////////////
//perform TM iteration for t generations
////////////////////////////////////////
tmiterate(double **a,double *mut1,
	  int t, int n1, int n2,int decay,double *sumf){

  int k, i, j = 0;
  double z;

  double * steadystatefreq = (double*) calloc (maxnd+1, sizeof(double));
  double * mut2 = (double*) calloc (maxnd+1, sizeof(double));

  for (i=0; i<=n1; i++) steadystatefreq[i] = 0;
  for (k = 1; k<=t; k++)
    {                 /*t iterations*/

      /* Increment steadystatefreq. distribution*/
      for (i = 0; i<=n1; i++)
	{
	  steadystatefreq[i] = steadystatefreq[i] + mut1[i];
	}

      /* Perform matrix multiplication*/
      for (i = 0; i<=n2; i++)
	{
	  z = 0;
	  for (j = 0; j<=n1; j++)
	    {
	      z = z + a[j][i]*mut1[j];
	    }
	  mut2[i] = z;
	}

      /* Copy result of multiplication*/
      for (i=0; i<=n2; i++) {mut1[i] = mut2[i];

	if (decay==0){
	  sumf[i]+=mut1[i];}
	if(decay==1){sumf[i]=mut1[i];}
      }


    }

  free(steadystatefreq);
  free(mut2);
}

eqf_using_matrix_inversion(int n1,double s,double **a,double *egf_out){


  setupmatrix(a, n1, n1, s, H);
  matrixinvert(n1,a,egf_out);


}
////////////////////////////////////////
//scale AFV for mutations eliminated by selection/not mutated
////////////////////////////////////////

egf_scaling_s(int n1,double *v_s_in, double *v_s_out,double f0){
  int n1d=2*n1;
  double sumx,sumy;
  int i;
  sumx=0;
  sumy=0;

  for (i=1;i<n1d;i++){
    sumx+=v_s_in[i];
  }

  for (i=1;i<n1d;i++){
    v_s_out[i]/=sumx;
    sumy+=v_s_out[i];
  }

  v_s_out[0]=(1-sumy)*(1-f0);
}
////////////////////////////////////////
//scale AFV for sites not mutated
////////////////////////////////////////
egf_scaling_f0(int n1,double *fv,double f0){

  int i;
  int n1d=2*n1;
  for (i=1;i<n1d;i++){
    fv[i]*=(1-f0);
  }
  fv[0]+=f0;
}
////////////////////////////////////////
//binomial sampling from AFV to generate SFS
////////////////////////////////////////
binomial_sampling(int n1,int n1b,int sample, double *invy,int *discrete,gsl_rng *rgen){

  int s1,success;
  int i;

  double prob;

  /*
    based on equilibrium frequency vector (use it as a probability vector),
    I randomly generate numbers from 0 to n1d-1 for sample sites.
  */

  gsl_ran_discrete_t *r= gsl_ran_discrete_preproc (n1,invy);
  for (i=0;i<sample;i++){
    s1=gsl_ran_discrete (rgen, r);
    prob=(double)(s1)/n1;

    success=gsl_ran_binomial(rgen,prob,n1b+1);
    discrete[success]++;
  }
  gsl_ran_discrete_free(r);
}

output_sfs_to_file(int n1,int *sfs1,int *sfs2,char *filename,long seed){
  int i;
  FILE *file;

  file = fopen(filename,"w"); 
 fprintf(file,"#seed used: %d\n",seed);
 fprintf(file,"#neutral\n");
  for (i=0;i<=n1;i++){
    fprintf(file,"%d ",sfs1[i]);
  }
 fprintf(file,"\n#selected\n");


  for (i=0;i<=n1;i++){
    fprintf(file,"%d ",sfs2[i]);
  }

  fclose(file);
}
////////////////////////////////////////
//Calculate weighted average of two vectors
// Normally vectors w(s) x(s) by N1 and N2 respectively
////////////////////////////////////////

vector_average(int n1,int n2,double *fv1,double *fv2, double *fv_out){
  int i;
  int n1d=n1*2;
  int n2d=n2*2;
  for (i=0;i<=n2d;i++){

    fv_out[i]=(n1*fv1[i]+n2*fv2[i])/(n1+n2);

  }
}
//////////////////////////////////////
//Obtain AFV for 1-epoch model
////////////////////////////////////////
calculate_FV_one_epoch(int n1,double s,double *FV){
  int i=0;
  int n1d=2*n1;
  double **FVW = calloc(maxnd+1, sizeof(double *));
 
  for(i = 0; i < maxnd+1; i++){
    FVW[i] = calloc(maxnd+1,  sizeof(double));}
  eqf_using_matrix_inversion(n1d,s,FVW,FV);
  FV[1]+=1;
  for(i = 0; i < maxnd+1; i++)
    free(FVW[i]);
  free(FVW);
}
////////////////////////////////////////
//Obtain AFV for 2-epoch model
////////////////////////////////////////
calculate_FV_two_epoch(int n1,int n2,int t,double s,double *FV){
  int i=0;
  double **egf_0 = calloc(maxnd+1, sizeof(double *));
  for(i = 0; i < maxnd+1; i++){
    egf_0[i] = calloc(maxnd+1,  sizeof(double));
  }  
         
  double * FVW = (double*) calloc (maxnd+1, sizeof(double));
  double * FVX = (double*) calloc (maxnd+1, sizeof(double));
  double * startingfreq = (double*) calloc (maxnd+1, sizeof(double));

  int n1d=n1*2;
  int n2d=n2*2;

  /*
    Contributions of mutations to AFV before the change in population size.
  */

  //obtain AFV at equilibrium

  eqf_using_matrix_inversion(n1d,s,egf_0,FVW);
  FVW[1]+=1.0;//cumulative vector includes 1st mutant
  //perform change in population size N1->N2
  setupmatrix(egf_0,n1d,n2d,s,H);
  tmiterate(egf_0,FVW,1, n1d, n2d,1,FVW);
  setupmatrix(egf_0,n2d,n2d,s,H);
  tmiterate(egf_0,FVW,t-1, n2d, n2d,1,FVW);

  //dumpvector(FVW,0,n2,"FVW");
  /*
    Contributions of mutations to AFV after the change in population size.
  */
  setupmatrix(egf_0,n2d,n2d,s,H);
  startingfreq[1]=1.0;
  tmiterate(egf_0,startingfreq, t-1, n2d, n2d,0,FVX);
  FVX[1]+=1.0;//cumulative vector includes 1st mutant
  //dumpvector(FVX,0,n2,"FVX");

  /*
    averaging of vectors FVW (N1) and FVX(N2)
  */
  vector_average(n1,n2,FVW,FVX,FV);


  for(i = 0; i < maxnd+1; i++){
    free(egf_0[i]);
  }

  free(egf_0);
  free(FVX);
  free(FVW);
  free(startingfreq);
}
////////////////////////////////////////
//print given matrix
////////////////////////////////////////

dumpmatrix(double **m, int s1,
	   int rows, int s2, int cols, char *s){
  int i, j;
  printf("\n%s\n", s);
  for (i = s1; i<=rows; i++)
    {
      for (j = s2; j<=cols; j++)
	{
	  printf(" %2.3f", m[i][j]);
	}
      printf("\n");
    }
  printf("\n");
}
////////////////////////////////////////
//print given vector
////////////////////////////////////////
dumpvector(double *v, int min, int max, char *s){
  int i, control = 0;
  //   printf("\n");
  printf("%s", s); printf("\n");
  for (i = min; i<=max; i++)
    {
      printf("%3d %lf ", i, v[i]);
      control++;
      if (control==5)
	{
	  control = 0;
	  printf("\n");
	}
    }
  printf("\n");
}




/*
////////////////////////////////////////
MAIN START
////////////////////////////////////////
*/

int main(int argc, char **argv)
{

  //initialise input variables and use default values
  int N1=100,n=10,sampleS=1000,sampleN=1000;
  int N2=N1;
  int t=N1;
  double f0=0.9,S=0;
  char outfile[maxnd+1]="sfs.out";
  long seed=time (NULL) * getpid();

  //other variables
  int option_index = 0;
  int c=0;
  static int verbose_flag;

  const gsl_rng_type * T;
  T = gsl_rng_taus;
  gsl_rng  *rgen =gsl_rng_alloc(T);
  gsl_rng_set(rgen,seed);

  int * discrete0 = (int*) calloc (maxnd+1, sizeof(int));
  int * discrete1 = (int*) calloc (maxnd+1, sizeof(int));

  double * FV0 = (double*) calloc (maxnd+1, sizeof(double));
  double * FVS = (double*) calloc (maxnd+1, sizeof(double));


  //parse input
  while (1)
    {

      static struct option long_options[] =
	{
	  {"h",no_argument,  &verbose_flag, 1},
	  {"N1",     required_argument,0, 'a'},
	  {"N2",  required_argument,0, 'b'},
	  {"n",  required_argument, 0, 'c'},
	  {"f0",  required_argument, 0, 'd'},
	  {"s",    required_argument, 0, 'e'},
	  {"LS",    required_argument, 0, 'f'},
	  {"LN",    required_argument, 0, 'g'},
	  {"t",    required_argument, 0, 'h'},
	  {"o",    required_argument, 0, 'i'},
	  {"seed",    required_argument, 0, 'j'},
	  {NULL, 0, 0, 0}
	};
      /* getopt_long stores the option index here. */
      c = getopt_long_only (argc, argv, "",
			    long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
	break;

      switch (c)
	{
	case 0:
	  /* If this option set a flag, do nothing else now. */
	  if (long_options[option_index].flag != 0)
	    {
	      printf("\
option\tdefault_value\tdescription\n\
N1\t100\tsize N1\n\
N2\t100\tsize N2\n\
t\t100\ttime since size change\n\
n\t10\tallele sample size\n\
f0\t0.9\t1-f0 is proportional to mutation rate\n\
s\t0\tselection coeffient\n\
LS\t1000\tno. neutral sites\n\
LN\t1000\tno. selected sites\n\
o\tsfs.out\toutput file\n\
seed\ttime*pid\tset seed for random generator\n\
");
	    }
	  break;
	  printf ("option %s", long_options[option_index].name);
	  if (optarg)
	    printf (" with arg %s", optarg);
	  printf ("\n");
	  break;

	case 'a':
	  N1=atoi(optarg);
	  break;
	case 'b':
	  N2=atoi(optarg);
	  break;
	case 'c':
	  n=atoi(optarg);
	  break;
	case 'd':
	  f0=atof(optarg);
	  break;
	case 'e':
	  S=atof(optarg);
	  break;
	case 'f':
	  sampleS=atoi(optarg);
	  break;
	case 'g':
	  sampleN=atoi(optarg);
	  break;
	case 'h':
	  t=atoi(optarg);
	  break;
	case 'i':
	  strncpy(outfile,optarg,maxnd+1);
	  break;
	case 'j':
	  seed=atoi(optarg);
	  gsl_rng_set(rgen,seed);
	  break;
	case '?':
	  abort();
	  break;
	default:
	  abort ();
	}
    }

  if (optind < argc)
    {
      printf ("non-option ARGV-elements: ");
      while (optind < argc)
	printf ("%s ", argv[optind++]);
      putchar ('\n');
    }

  if (N1>1000||N2>1000){printf ("Use N1,N2 <1000\n"); abort();}
  /*
    calculate allele frequency vectors for neutral (FV0)
    and selected sites (FVS)
  */
  if (N1==N2){
    calculate_FV_one_epoch(N1,0.0,FV0);
    calculate_FV_one_epoch(N1,S,FVS);
  }else{
    calculate_FV_two_epoch(N1,N2,t,0.0,FV0);
    calculate_FV_two_epoch(N1,N2,t,S,FVS);
  }

  /*
    scale FVS to include the frequency of sites at which mutations
    have been eliminated by selection or have not experienced a mutation
  */
  egf_scaling_s(N2,FV0,FVS, f0);
  egf_scaling_f0(N2,FVS,f0);

  /*
    scale FV0 to account for sites that have never experienced a mutation.
  */
  egf_scaling_s(N2,FV0, FV0,f0);
  egf_scaling_f0(N2,FV0,f0);


  /*sampling
    sample sampleN neutral sites from FV0
    sample sampleS selected sites from FVS
  */

  binomial_sampling(N2,n,sampleN,FV0,discrete0,rgen);
  binomial_sampling(N2,n,sampleS, FVS,discrete1,rgen);

  output_sfs_to_file(n,discrete0,discrete1,outfile,seed);

  //free vectors
  free(FV0);
  free(FVS);

  free(discrete0);
  free(discrete1);

  gsl_rng_free (rgen);

  return 0;
}

/*
////////////////////////////////////////
END OF MAIN
////////////////////////////////////////
*/
