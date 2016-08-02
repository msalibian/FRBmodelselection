

// file FRB-model-selection.c
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* #define Psi_constant 4.685061 */
/* #define T_Chi_constant 1.987967 */ /* BP: 40% */
/* #define T_Chi_constant 1.54764 */ /* BP: 50% */
#define NA 99999.99
#define EPSILON 1e-20
#define EPS 1e-10
#define ZERO 1e-10
#define INFI 1e+20
#define MAX_ITER_FAST_S 500
#define MAX_NO_RESAMPLES 500
#define MAX_ITER_FIND_SCALE 200
#define TOL_INVERSE ZERO


/* FILE:  regSrand.c
 * Robust bootstrap for regression estimates
 */

double R_rlm_rand(double *X, double *y, int *N, int *P,
		int *Boot_Samp, int *Nres,
		int *M, int *size_boot, double *ours, double *full,
		double *Beta_m, double *Beta_s, double *Scale,
		int *Seed, int *calc_full,
		double *C, double *Psi_c, int *max_it,
		int *converged_mm,
		int *groups, int *n_group, int *k_fast_s)
{
void initialize_mat(double **a, int n, int m);
void initialize_vec(double *a, int n);
void R_S_rlm(double *X, double *y, int *n, int *P,
		int *nres, int *max_it,
		double *SCale, double *beta_s, double *beta_m,
		int *converged_mm,
		int *seed_rand, double *C, double *Psi_c,
		int *Groups, int *N_group, int *K_fast_s);
double Psi_reg(double,double);
double Psi_reg_prime(double,double);
double Chi_prime(double,double);
double Chi(double,double);
void sampler_i(int, int, int *);
int inverse(double **,double **, int);
void matias_vec_vec(double **, double *, double *, int);
void scalar_mat(double **, double, double **, int, int);
void scalar_vec(double *, double, double *, int);
void sum_mat(double **,double **, double **, int, int);
void sum_vec(double *, double *, double *, int);
void dif_mat(double **, double **, double **, int , int );
void dif_vec(double *, double *, double *, int);
void mat_vec(double **, double *, double *, int, int);
void mat_mat(double **, double **, double **, int, int, int);
// void disp_vec(double *, int);
// void disp_mat(double **, int, int);
// void disp_mat_i(int **, int, int);
// void disp_vec(double *, int);
/* double **xb; */
double *Xb, **xb;
int **boot_samp;
double **x, **x2, **x3, **x4, *beta_m, *beta_s,*beta_aux;
double *Fi, *res, *res_s, *w, *ww, dummyscale, scale;
double *v, *v2, *v_aux, *yb; // , timefinish, timestart;
double u,u2,s,c,Psi_constant;
// double test_chi=0, test_psi=0;
int n,p,m,seed; // ,*indices;
int nboot=*size_boot;
// int fake_p = 0;
register int i,j,k;
setbuf(stdout,NULL);
c = *C; Psi_constant = *Psi_c;
n = *N; p = *P; m = *M; seed = *Seed;
boot_samp = (int **) malloc(m * sizeof(int*) );
for(i=0;i<m;i++)
	boot_samp[i] = (int*) malloc(nboot *sizeof(int));
// indices = (int *) malloc( n * sizeof(int) );
v = (double *) malloc( p * sizeof(double) );
v2 = (double *) malloc( p * sizeof(double) );
v_aux = (double *) malloc( p * sizeof(double) );
yb = (double *) malloc( n * sizeof(double) );
Xb = (double*) malloc( n * p * sizeof(double) );
x =  (double **) malloc ( n * sizeof(double *) );
xb =  (double **) malloc ( n * sizeof(double *) );
Fi  = (double *) malloc ( n * sizeof(double) );
res = (double *) malloc ( n * sizeof(double) );
res_s = (double *) malloc ( n * sizeof(double) );
ww  = (double *) malloc ( n * sizeof(double) );
w   = (double *) malloc ( n * sizeof(double) );
x2 = (double **) malloc ( p * sizeof(double *) );
x3 = (double **) malloc ( p * sizeof(double *) );
x4 = (double **) malloc ( p * sizeof(double *) );
beta_aux = (double *) malloc( p * sizeof(double) );
beta_m = (double *) malloc( p * sizeof(double) );
beta_s = (double *) malloc( p * sizeof(double) );
for(i=0;i<n;i++) {
	x[i] =  (double*) malloc (p * sizeof(double) );
	xb[i] =  (double*) malloc ((p+1) * sizeof(double) );
	};
for(i=0;i<p;i++) {
	x2[i] = (double*) malloc (p * sizeof(double) );
	x3[i] = (double*) malloc (p * sizeof(double) );
	x4[i] = (double*) malloc (p * sizeof(double) );
};
/* copy X into x for easier handling */
for(i=0;i<n;i++)
        for(j=0;j<p;j++)
                x[i][j]=X[j*n+i];
/* calculate robust regression estimates */

for(i=0;i<m;i++)
	for(j=0;j<nboot;j++)
		boot_samp[i][j]=Boot_Samp[j*m+i]-1;

R_S_rlm(X, y, N, P, Nres, max_it, &scale, Beta_s, Beta_m,
		converged_mm, &seed, &c,
		Psi_c, groups, n_group, k_fast_s);

*Scale = scale;
/* get M-fitted values in Fi */
mat_vec(x,Beta_m,Fi,n,p);
/* get residuals of M-est in res */
dif_vec(y,Fi,res,n);
/* get S-fitted values in res_s */
mat_vec(x,Beta_s,res_s,n,p);
/* get residuals of S-est in res_s */
dif_vec(y,res_s,res_s,n);
/* set auxiliary matrices to zero */

initialize_mat(x3, p, p);
initialize_mat(x4, p, p);
initialize_vec(v, p);
u2 = 0.0;
/* calculate correction matrix */

for(i=0;i<n;i++) {
	u = res[i]/scale ;
	w[i]  = Psi_reg(u,Psi_constant)/res[i];
        matias_vec_vec(x2,x[i],x[i],p);
	scalar_mat(x2,Psi_reg_prime(u,Psi_constant),
                x2,p,p);
        sum_mat(x3,x2,x3,p,p);
        matias_vec_vec(x2,x[i],x[i],p);
        scalar_mat(x2,w[i],x2,p,p);
        sum_mat(x4,x2,x4,p,p);
	scalar_vec(x[i],Psi_reg_prime(u,Psi_constant)*u,v_aux,p);
	sum_vec(v,v_aux,v,p);
	u2 += Chi_prime(u, c) * u;
};

/* scalar_vec(v, .5 * (double) (n-p) * scale / u2 , v, p);  */
scalar_vec(v, .5 * (double) n * scale / u2 , v, p);
inverse(x3,x2,p);
mat_mat(x2,x4,x3,p,p,p);
mat_vec(x2,v,v2,p,p);
scalar_mat(x3,scale,x3,p,p);
/* the correction matrix is now in x3 */
/* the correction vector is now in v2 */

/* start the bootstrap replications */
for(i=0;i<m;i++) {
	/* change the seed! */
	++seed;
	// sampler_i(n,nboot,indices);
	// for(j=0;j<nboot; j++)
	// 	indices[j]=boot_samp[i][j];
	/* get pseudo observed y's */
	for(j=0;j<nboot;j++) /* xb[j][p] = */
			yb[j] = y[boot_samp[i][j]];
	for(j=0;j<nboot;j++)
		for(k=0;k<p;k++) {
			// xb[j][k] = x[boot_samp[i][j]][k];
			// Xb[k*nboot+j] = X[k*n + indices[j]];
			Xb[k*nboot+j] = x[boot_samp[i][j]][k];
			xb[j][k] = Xb[k*nboot+j];
		};

	/* calculate full bootstrap estimate */

	if( *calc_full == 1 )
	R_S_rlm(Xb,yb,&nboot,P,Nres,max_it,&dummyscale,
			beta_s,beta_m,converged_mm,&seed,&c,
			Psi_c, groups, n_group, k_fast_s);

/* void R_S_rlm(double *X, double *y, int *n, int *P,
		int *nres, int *max_it,
		double *SCale, double *beta_s, double *beta_m,
		int *converged_mm,
		int *seed_rand, double *C, double *Psi_c,
		int *Groups, int *N_group, int *K_fast_s) */

	/*	double *C, double *Psi_c, int *max_it,
		int *groups, int *n_group, int *k_fast_s); */

	// HERE

	/* disp_mat(xb, nboot,p); */
	// disp_vec(yb,nboot);
	// Rprintf("\nfull scale: %f", dummyscale);

	/* calculate robust bootsrap */

	scalar_vec(v,0.0,v,p);	 	/* v <- 0 */
	scalar_mat(x2,0.0,x2,p,p);	/* x2 <- 0 */
	s = 0.0;
	for(j=0;j<nboot;j++) {
		scalar_vec(xb[j],yb[j]*w[boot_samp[i][j]],v_aux,p);
		sum_vec(v,v_aux,v,p);
		matias_vec_vec(x4,xb[j],xb[j],p);
		scalar_mat(x4,w[boot_samp[i][j]],x4,p,p);
		sum_mat(x2,x4,x2,p,p);
		s += Chi(res_s[boot_samp[i][j]] / scale , c);
	};
	/* s = s * scale / .5 / (double) (nboot - p)  ;  */
	s = s * scale / .5 / (double) n;
	inverse(x2,x4,p);		/* x4 <- x2^-1 */
	mat_vec(x4,v,v_aux,p,p);	/* v_aux <- x4 * v */
	dif_vec(v_aux,Beta_m,v_aux,p); 	/* v_aux <- v_aux - beta_m */
	/* v has the robust bootstrapped vector, correct it */
	mat_vec(x3,v_aux,v,p,p);	/* v <- x3 * v_aux */
	scalar_vec(v2,s-scale,v_aux,p);
	sum_vec(v_aux,v,v,p);

	/* store the betas (splus-wise!) */
	for(j=0;j<p;j++) {
		ours[j*m+i]=v[j];
		if( *calc_full == 1 )
			// full[j*m+i]=beta_m[j]-Beta_m[j];
			full[j*m+i]=beta_m[j];
	};
};
for(i=0;i<m;i++)
	free(boot_samp[i]);
free(boot_samp);
for(i=0;i<n;i++) {
	free(x[i]);
	free(xb[i]);
	};
for(i=0;i<p;i++) {
	free(x2[i]);
	free(x3[i]);
	free(x4[i]);
	};
free(x) ;free(x2);free(xb);
free(x3);free(x4);
free(beta_aux);free(beta_m);free(beta_s);
free(w);free(ww);free(Fi);free(res);
free(v);free(v2);free(v_aux);free(yb);
free(res_s);
free(Xb);
return(0);
}

/* FILE:  regSfixed.c
 * Robust bootstrap for regression estimates
 */

double R_rlm_fixed(double *X, double *y, int *N, int *P, int *Nres,
		int *M, int *size_boot, double *ours, double *full,
		double *Beta_m, double *Beta_s, double *Scale, int *Seed,
		int *calc_full,
		double *C, double *Psi_c, int *max_it,
		int *converged_mm,
		int *groups, int *n_group, int *k_fast_s)
{
void initialize_mat(double **a, int n, int m);
void initialize_vec(double *a, int n);
void R_S_rlm(double *X, double *y, int *n, int *P,
		int *nres, int *max_it,
		double *SCale, double *beta_s, double *beta_m,
		int *converged_mm,
		int *seed_rand, double *C, double *Psi_c,
		int *Groups, int *N_group, int *K_fast_s);
double Psi_reg(double,double);
double Psi_reg_prime(double,double);
double Chi_prime(double,double);
double Chi(double,double);
void sampler_i(int, int, int *);
int inverse(double **,double **, int);
void matias_vec_vec(double **, double *, double *, int);
void scalar_mat(double **, double, double **, int, int);
void scalar_vec(double *, double, double *, int);
void sum_mat(double **,double **, double **, int, int);
void sum_vec(double *, double *, double *, int);
void dif_mat(double **, double **, double **, int , int );
void dif_vec(double *, double *, double *, int);
void mat_vec(double **, double *, double *, int, int);
void mat_mat(double **, double **, double **, int, int, int);
// void disp_vec(double *, int);
double **x, **x2, **x3, **x4, *beta_m, *beta_s,*beta_aux;
double *Fi, *res, *res_s, *w, *ww, dummyscale, scale;
double *v, *v2, *v_aux, *yb; // , timefinish, timestart;
double u,u2,s,c,Psi_constant;
int n,p,m,seed,*indices;
int nboot = *size_boot;
register int i,j;

setbuf(stdout,NULL);
c = *C; Psi_constant = *Psi_c;
n = *N; p = *P; m = *M; seed = *Seed;
indices = (int *) malloc( n * sizeof(int) );
v = (double *) malloc( p * sizeof(double) );
v2 = (double *) malloc( p * sizeof(double) );
v_aux = (double *) malloc( p * sizeof(double) );
yb = (double *) malloc( n * sizeof(double) );
x =  (double **) malloc ( n * sizeof(double *) );
Fi  = (double *) malloc ( n * sizeof(double) );
res = (double *) malloc ( n * sizeof(double) );
res_s = (double *) malloc ( n * sizeof(double) );
ww  = (double *) malloc ( n * sizeof(double) );
w   = (double *) malloc ( n * sizeof(double) );
x2 = (double **) malloc ( p * sizeof(double *) );
x3 = (double **) malloc ( p * sizeof(double *) );
x4 = (double **) malloc ( p * sizeof(double *) );
beta_aux = (double *) malloc( p * sizeof(double) );
beta_m = (double *) malloc( p * sizeof(double) );
beta_s = (double *) malloc( p * sizeof(double) );
for(i=0;i<n;i++)
	x[i] =  (double*) malloc ((p+1) * sizeof(double) );
for(i=0;i<p;i++) {
	x2[i] = (double*) malloc (p * sizeof(double) );
	x3[i] = (double*) malloc (p * sizeof(double) );
	x4[i] = (double*) malloc (p * sizeof(double) );
};

/* copy X into x for easier handling */
for(i=0;i<n;i++)
        for(j=0;j<p;j++)
                x[i][j]=X[j*n+i];

/* calculate robust regression estimates */
R_S_rlm(X,y,N,P,Nres,max_it,&scale,Beta_s,Beta_m,
		converged_mm,&seed,C,
		Psi_c, groups, n_group, k_fast_s);
*Scale = scale;

/* get M-fitted values in Fi */
mat_vec(x,Beta_m,Fi,n,p);

/* get residuals of M-est in res */
dif_vec(y,Fi,res,n);

/* get S-fitted values in res_s */
mat_vec(x,Beta_s,res_s,n,p);

/* get residuals of S-est in res_s */
dif_vec(y,res_s,res_s,n);

/* set auxiliary matrices to zero */

//
// scalar_mat(x3,0.0,x3,p,p);
// scalar_mat(x4,0.0,x4,p,p);
// scalar_vec(v,0.0,v,p);
//

initialize_mat(x3,p,p);
initialize_mat(x4,p,p);
initialize_vec(v,p);

u2 = 0.0;

/* calculate correction matrix */
for(i=0;i<n;i++) {
	u = res[i]/scale ;
	w[i]  = Psi_reg(u,Psi_constant)/res[i];
        matias_vec_vec(x2,x[i],x[i],p);
	scalar_mat(x2,Psi_reg_prime(u,Psi_constant),
                x2,p,p);
        sum_mat(x3,x2,x3,p,p);
        matias_vec_vec(x2,x[i],x[i],p);
        scalar_mat(x2,w[i],x2,p,p);
        sum_mat(x4,x2,x4,p,p);
	scalar_vec(x[i],Psi_reg_prime(u,Psi_constant)*u,v_aux,p);
	sum_vec(v,v_aux,v,p);
	u2 += Chi_prime(u,c) * u;
};
scalar_vec(v, .5 * (double) n * scale / u2 , v, p);
inverse(x3,x2,p);
mat_mat(x2,x4,x3,p,p,p);
mat_vec(x2,v,v2,p,p);
scalar_mat(x3,scale,x3,p,p);
/* the correction matrix is now in x3 */
/* the correction vector is now in v2 */

/* disp_mat(x3, p, p);
disp_vec(v2, p);
*/

/* start the bootstrap replications */
for(i=0;i<m;i++) {
	/* change the seed! */
	++seed;
	sampler_i(n,nboot,indices);

	/* get pseudo observed y's */
	for(j=0;j<nboot;j++) yb[j] = x[j][p] = Fi[j] + res[indices[j]];

	/* calculate full bootstrap estimate */

	if( *calc_full == 1 ) {
	R_S_rlm(X,yb,&nboot,P,Nres,max_it,&dummyscale,beta_s,beta_m,
			converged_mm,&seed,C,
			Psi_c, groups, n_group,k_fast_s);
	if( dummyscale == 0.0 )
		for(j=0;j<p;j++) beta_m[j] = NA;
	};

	/* calculate the robust bootstrap estimate */
	scalar_vec(v,0.0,v,p);	 	/* v <- 0 */
	scalar_mat(x2,0.0,x2,p,p);	/* x2 <- 0 */
	s = 0.0;
	for(j=0;j<nboot;j++) {
		scalar_vec(x[j],yb[j]*w[indices[j]],v_aux,p);
		sum_vec(v,v_aux,v,p);
		matias_vec_vec(x4,x[j],x[j],p);
		scalar_mat(x4,w[indices[j]],x4,p,p);
		sum_mat(x2,x4,x2,p,p);
		s += Chi(res_s[indices[j]] / scale , c);
	};
	/* s = s * scale / .5 / (double) (nboot-p); */
	s = s * scale / .5 / (double) n;
	inverse(x2,x4,p);		/* x4 <- x2^-1 */
	mat_vec(x4,v,v_aux,p,p);	/* v_aux <- x4 * v */
	dif_vec(v_aux,Beta_m,v_aux,p); 	/* v_aux <- v_aux - beta_m */

	/* v has the robust bootstrapped vector, correct it */
	mat_vec(x3,v_aux,v,p,p);	/* v <- x3 * v_aux */
	scalar_vec(v2,s-scale,v_aux,p);
	sum_vec(v_aux,v,v,p);
	/* store the betas */
	for(j=0;j<p;j++) {
		ours[j*m+i]=v[j];
		if( *calc_full == 1 )
			full[j*m+i]=beta_m[j] - Beta_m[j];
	};
};

for(i=0;i<n;i++)
	free(x[i]);
for(i=0;i<p;i++) {
	free(x2[i]);
	free(x3[i]);
	free(x4[i]);
	};
free(x) ;free(x2);
free(x3);free(x4);
free(beta_aux);free(beta_m);free(beta_s);
free(w);free(ww);free(Fi);free(res);
free(v);free(v2);free(v_aux);free(yb);
free(indices);
return(0);
}


int lu(double **a,int *P, double *x)
{
int *pp,p;
register int i,j,k;
double *kk,s;
p = *P;
if ((pp = (int *) malloc(p*sizeof(int)))==NULL)
	{ // Rprintf("\nNot enough memory in LU\n");
	  exit(1); }
/* pp vector storing the permutations */
for(j=0;j<p;j++)   /* cols */
{ pp[j]=j;
  for(i=j;i<p;i++)   /* filas */
	if ( fabs( a[i][j] ) > fabs( a[pp[j]][j] ) )
		pp[j]=i;
  if ( pp[j] != j )       /* permuto las filas cambiando los punt */
	{ kk=a[j];
	  a[j]=a[pp[j]];
	  a[pp[j]]=kk;
	};
  /* salida si el sistema resulta singular (det=0)
   * se detecta si el pivote (j,j) es cero  */
/*  if ( a[j][j] == 0 ) {   free(pp);
				return(1);
				}; */
    if ( fabs(a[j][j]) < TOL_INVERSE ) {   free(pp);
				return(1);
				};
  for(k=(j+1);k<p;k++)
	a[k][j] = a[k][j] / a[j][j];
  for(k=(j+1);k<p;k++)
	for(i=(j+1);i<p;i++)
		a[k][i] = a[k][i] - a[k][j] * a[j][i];

};    /* cierra el for de j */
for(i=0;i<p;i++)
	{ s=0.0;
	  for(j=0;j<i;j++)
	    s += a[i][j] * x[j];
	    x[i] = a[i][p] - s;          /* y[i]=a[i][p] */
	};
for(i=(p-1);i>=0;i--)
	{ s=0;
	  for(j=(i+1);j<p;j++)
	    s += a[i][j] * x[j];
	  x[i] = (x[i] - s) / a[i][i];
	  };
free(pp);
return(0);
}


double Chi_prime(double x, double c)
{
/* //
// // Tukey's bisquare loss function
// */
double t;
if( fabs(x) > c ) return(0.0);
else { t = x / c ;
	return( 6.0*t*(1 - t*t) * (1-t*t) / c );
	}
}


double Chi(double x, double c)
{
/* //
// // Tukey's bisquare loss function
// */
double t;
if( fabs(x) > c ) return(1.0);
else { t = x / c;
	return( 3.0*t*t - 3.0*t*t*t*t + t*t*t*t*t*t );
	}
}


double loss_Tukey(double x, double c)
{
if( fabs(x/c) < 1 )
	return( (x/c)*(x/c)/2. *
			( 1 - (x/c)*(x/c) +
		(x/c)*(x/c)*(x/c)*(x/c)/3. ) );
else
	return( 1. / 6. );
}

double Loss_Tukey(double *x, int n, double c)
{
double loss_Tukey(double,double);
double s=0;
register int i;

for(i=0;i<n;i++) s += loss_Tukey(x[i],c);

return(s);
}

double Psi_reg(double x, double c)
{
if (fabs(x)>c) return(0.0);
else	return( x / c * (1.0-(x/c)*(x/c))*
		(1.0-(x/c)*(x/c))  );
}


double Psi_reg_prime(double x, double c)
{
if (fabs(x)>c) return(0.0);
else	return( ( 1.0 - (x/c)*(x/c) ) *
		( 1.0 - 5.0 * x * x / c / c ) / c );
}


/* ----------------------------------------------------------- */


double kthplace(double *a, int n, int k)
{
int jnc,j;
int l,lr;
double ax,w;
k--;
l=0;
lr=n-1;
while (l<lr)
	{ ax=a[k];
	  jnc=l;
	  j=lr;
	  while (jnc<=j)
		{ while (a[jnc] < ax) jnc++;
		  while (a[j] > ax) j--;
		  if (jnc <= j)
			{ w=a[jnc];
			  a[jnc]=a[j];
			  a[j]=w;
			  jnc++;
			  j--;
			};
		};
	  if (j<k) l=jnc;
	if (k<jnc) lr=j;
	};
return(a[k]);
}

/* ----------------------------------------------------------- */

void sampler_i(int n, int m, int *x)
{
/* function to get a random sample of
 * m indices from (0 to n-1)
 * *x receives the output
 * rand() returns an integer between 0 and RAND_MAX
 */
int i;
for(i=0;i<m;i++)
	x[i] = (int) ( (double) rand() / RAND_MAX * (double) (n-1) );
}

void sample_n_outof_N(int n, int N, int *x)
{
/* function to get a random sample of size n
 * of the indices (0 to N) WITHOUT replication
 * *x receives the output
 * rand() returns an integer between 0 and RAND_MAX
 */
register int i,j,cand,flag;
if( N < n ) {
	// Rprintf("\nCant get %d out of %d without replication\n", n, N);
	for(i=0;i<n;i++) x[i] = i;
} else {
for(i=0;i<n;i++) {
	flag=1;
	while (flag==1) {
		flag=0;
		cand = (int) ( (double) rand() / RAND_MAX *
					(double) N );
		for(j=0;j<i;j++)
			if( cand==x[j] ) flag=1;
		};
	x[i]=cand;
	};
}
}

/* ----------------------------------------------------------- */

double variance(double *x,int n)
{
double mean(double *, int);
double s=0, x_bar;
register int i;
x_bar = mean(x,n);
for(i=0;i<n;i++) s += (x[i]-x_bar) * (x[i]-x_bar);
return( s / (double) (n-1) );
}

double mean(double *x, int n)
{
double s=0;
register int i;
for(i=0;i<n; s += x[i++]);
return( s / (double) n);
}



double norm_diff(double *x, double *y, int n)
{
double s=0;
register int i;
for(i=0;i<n;i++)
	s += (x[i]-y[i])*(x[i]-y[i]);
return( sqrt(s) );
}

/*
 *        MATRIX / VECTOR OPERATIONS
 */

void sum_mat(double **a, double **b, double **c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<m;j++)
		c[i][j] = a[i][j] + b[i][j];
}

void matias_vec_vec(double **a, double *v1, double *v2, int n)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<n;j++)
		a[i][j] = v1[i] * v2[j];
/* could modify it to take advantage of symmetry */
}

void scalar_mat(double **a, double b, double **c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
        for(j=0;j<m;j++)
		c[i][j]  = b * a[i][j];
}

void scalar_vec(double *a, double b, double *c, int n)
{
register int i;
for(i=0;i<n;i++)
	c[i]  = b * a[i];
}

double vecprime_vec(double *a, double *b, int n)
{
register int i;
double s = 0.0;
for(i=0;i<n;i++) s += a[i] * b[i];
return(s);
}

void sum_vec(double *a, double *b, double *c, int n)
{
register int i;
for(i=0;i<n;i++) c[i] = a[i] + b[i];
}

void dif_vec(double *a, double *b, double *c, int n)
{
register int i;
for(i=0;i<n;i++) c[i] = a[i] - b[i];
}

void dif_mat(double **a, double **b, double **c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<m;j++) c[i][j] = a[i][j] - b[i][j];
}

void mat_vec(double **a, double *b, double *c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(c[i]=0,j=0;j<m;j++) c[i] += a[i][j] * b[j];
}

void mat_mat(double **a, double **b, double **c, int n,
		int m, int l)
{
register int i,j,k;
for(i=0;i<n;i++)
	for(j=0;j<l;j++) {
	c[i][j] = 0;
	for(k=0;k<m;k++) c[i][j] += a[i][k] * b[k][j];
	};
}

// void disp_vec(double *a, int n)
// {
// register int i;
// printf("\n");
// for(i=0;i<n; i++) printf("%lf ",a[i]);
// printf("\n");
// }


int inverse(double **a, double **b, int n)
{
int lu(double **, int *, double *);
void mat_vec(double **, double *, double *, int, int);
// void disp_vec(double *, int);
register int i,j,k;
double **c, *e;
c = (double **) malloc( n * sizeof(double *));
e = (double *) malloc( n * sizeof(double));
for(i=0;i<n;i++) c[i] = (double *) malloc ( (n+1) * sizeof(double) );
for(i=0;i<n;i++) {   /* i-th column */

for(j=0;j<n;j++)
	for(k=0;k<n;k++) c[j][k] = a[j][k];
for(j=0;j<i;j++) c[j][n] = 0.0;
c[i][n] = 1.0;
for(j=i+1;j<n;j++) c[j][n] = 0.0;
if( lu(c,&n,e) == 1) {
	for(i=0;i<n;i++) free(c[i]);
	free(c);free(e);
	return(1);
	};
for(j=0;j<n;j++) b[j][i] = e[j] ;
};
for(i=0;i<n;i++) free(c[i]);
free(c);free(e);
return(0);
}

/* ********************************* */


void initialize_mat(double **a, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<m;j++)
		a[i][j] = 0.0;
}

void initialize_vec(double *a, int n)
{
register int i;
for(i=0;i<n;i++) a[i] = 0.0;
}
/*
void disp_vec(double *a, int n)
{
register int i;
Rprintf("\n");
for(i=0;i<n; i++) Rprintf("%lf ",a[i]);
Rprintf("\n");
}

void disp_vec_i(int *a, int n)
{
register int i;
Rprintf("\n");
for(i=0;i<n; i++) Rprintf("%d ",a[i]);
Rprintf("\n");
}

void disp_mat(double **a, int n, int m)
{
register int i,j;
for(i=0;i<n;i++) {
Rprintf("\n");
for(j=0;j<m;j++) Rprintf("%10.8f ",a[i][j]);
};
Rprintf("\n");
}

void disp_mat_i(int **a, int n, int m)
{
register int i,j;
for(i=0;i<n;i++) {
Rprintf("\n");
for(j=0;j<m;j++) Rprintf("%d ",a[i][j]);
};
Rprintf("\n");
}
*/

void R_S_rlm(double *X, double *y, int *n, int *P,
		int *nres, int *max_it,
		double *SCale, double *beta_s, double *beta_m,
		int *converged_mm,
		int *seed_rand, double *C, double *Psi_c,
		int *Groups, int *N_group, int *K_fast_s)
/* x <- design, y <- response, n <- #obs
   p <- dimension, nres <- #resampling */
{
void fast_s_large_n(double *X, double *y,
 		int *nn, int *pp, int *NN, int *K,
		int *ggroups, int *nn_group,
		int *bbest_r, double *bb, double *rrhoc,
		double *bbeta, double *sscale);
void fast_s(double *X, double *y,
		int *nn, int *pp, int *NN, int *K,
		int *bbest_r, double *bb,
		double *rrhoc, double *bbeta, double *sscale);
int rwls(double **, int, int, double *, double *, double,
		double, int, double);
int rwls_chi(double **, int, int, double *, double *, double *,
		double,double);
// void disp_vec(double *, int);
void sample_n_outof_N(int, int, int *);
int lu(double **,int *, double *);
register int i,j; /* ,k; */
double best_s=NA, **x_samp, *cand_beta, **x; /* scale */
double *resid; /* ,s;  */
double *beta, *temp1, *temp2;
int N,p,*b_i; /* , zeroes;  */
int bbest_r = 2, ggroups = *Groups,
	nn_group = *N_group, k_fast_s = *K_fast_s;
double b = 0.5, rrhoc = *C;
// c = *C;
N = *n; p = *P;
srand((long)*seed_rand);
x = (double **) malloc( N * sizeof(double*) );
for(i=0;i<N;i++)
        x[i]= (double *) malloc( (p+1) * sizeof(double) );
b_i = (int *) malloc( p * sizeof(int) );
cand_beta = (double *) malloc( p * sizeof(double) );
beta = (double *) malloc( p * sizeof(double) );
resid = (double *) malloc( N * sizeof(double) );
x_samp = (double **) malloc( p * sizeof(double*) );
temp1 = (double*) malloc( N * sizeof(double) );
temp2 = (double*) malloc( N * sizeof(double) );
for(i=0;i<p;i++)
	x_samp[i] = (double *) malloc( (p+1) * sizeof(double) );
/* rearranges X into a matrix of n x p */
for(i=0;i<N;i++) {
        for(j=0;j<p;j++)
                x[i][j]=X[j*N+i];
        x[i][p]=y[i]; /* don't really use x[][p] but it's
                       * convenient when calling rwls */
        };

*nres = 500; /* need less re-sampling candidates! */
if( *n > 2000 )
	fast_s_large_n(X, y, n, P, nres, &k_fast_s,
			&ggroups, &nn_group, &bbest_r,
			&b, &rrhoc, beta_s, &best_s);
else
	fast_s(X, y, n, P, nres, &k_fast_s,
			&bbest_r, &b, &rrhoc, beta_s, &best_s);
if(fabs(best_s) < ZERO) {
			 *SCale = 0.0;
			for(i=0;i<p;i++) 
			        beta_m[i]=beta_s[i];
			free(b_i);free(temp1); free(temp2);
			for(i=0;i<p;i++)
        			free(x_samp[i]);
			free( x_samp ); free(resid);
			free( cand_beta );
			free(beta);
			for(i=0;i<N;i++)
        			free(x[i]);
			free(x);
			return;
};

*SCale = best_s;
/* starting from the S-estimate (beta_s), use
 * irwls to compute the M-estimate (beta_m)  */
*converged_mm = 1; /* converged by default */
if ( rwls(x,N,p,beta_m,beta_s,best_s,EPS,*max_it,*Psi_c) == 1 )  {
	/* Rprintf("\nRWLS did not converge!\n");  */
	for(i=0;i<p;i++) beta_m[i]=beta_s[i];
	*converged_mm = 0;  /* rwls failed to converge */
};
for(i=0;i<p;i++)
	free(x_samp[i]);
free( x_samp ); free(resid); free( cand_beta ); free(beta);
for(i=0;i<N;i++)
	free(x[i]);
free(x);free(temp1); free(temp2);free(b_i);
}


/*
//
// 2004 / 5 -- Matias Salibian-Barrera & Victor Yohai
// Department of Statistics, University of British Columbia
// matias@stat.ubc.ca
// Department of Mathematics, University of Buenos Aires
// vyohai@uolsinectis.com.ar
//
//
// Reference: A fast algorithm for S-regression estimates,
// 2005, Salibian-Barrera and Yohai.


// This function implements the "large n" strategy
*/

void fast_s_large_n(double *X, double *y,
		int *nn, int *pp, int *NN, int *K,
		int *ggroups, int *nn_group,
		int *bbest_r, double *bb, double *rrhoc,
		double *bbeta, double *sscale)
{
/* *X = the design matrix as a vector (as sent by R)
// *y = the response vector
// *nn = the length of y
// *pp = the number of columns in X
// *NN = number of re-sampling candidates to be
//       used in each partition
// *bbest_r = no. of best candidates to be iterated
//            further
// *bb = right-hand side of S-equation (typically 1/2)
// *rrhoc = tuning constant for Tukey's bi-square loss
//          (this should be associated with *bb)
// *bbeta = final estimator
// *sscale = associated scale estimator
// *ggroups = number of groups in which to split the
//            random subsample
// *nn_group = size of each of the (*ggroups) groups
//             to use in the random subsample
*/
void reset_mat(double **a, int n, int m);
double loss_rho(double *r, double scale, int n, int p, double rhoc);
double find_scale(double *r, double b, double rhoc,
			double initial_scale, int n, int p);
// void disp_mat(double **a, int n, int m);
// void disp_vec(double *a, int n);
void fast_s_with_memory(double **x, double *y,
		int *nn, int *pp, int *NN, int *K,
		int *bbest_r, double *bb,
		double *rrhoc,
		double **best_betas, double *best_scales);
void sample_n_outof_N(int n, int N, int *x);
void refine_fast_s(double **x, double *y, double *weights,
			int n, int p, double *res,
			double *tmp, double *tmp2,
			double **tmp_mat, double **tmp_mat2,
			double *beta_cand, int kk,
			int conv, double b, double rhoc, double *is, double *beta_ref,
			double *scale);
int find_max(double *a, int n);
register int i,j,k,k2;
int n = *nn, p = *pp, kk = *K, *indices;
int groups = *ggroups, n_group = *nn_group, best_r = *bbest_r;
double **best_betas, *best_scales;
double **final_best_betas, *final_best_scales;
double **x, **xsamp, *ysamp, *res, sc, *beta_ref;
double *tmp, *tmp2, **tmp_mat, **tmp_mat2, *weights;
double best_sc, worst_sc, b = *bb, rhoc = *rrhoc, aux;
int pos_worst_scale, conv;
res = (double *) malloc( n * sizeof(double) );
weights = (double *) malloc( n * sizeof(double) );
tmp  = (double *) malloc( n * sizeof(double) );
tmp2 = (double *) malloc( n * sizeof(double) );
tmp_mat  = (double **) malloc( p * sizeof(double *) );
tmp_mat2 = (double **) malloc( p * sizeof(double *) );
for(i=0;i<p;i++) {
	tmp_mat[i] = (double *) malloc( p * sizeof(double) );
	tmp_mat2[i] = (double *) malloc( (p+1) * sizeof(double) );
};
beta_ref = (double *) malloc( p * sizeof(double) );
final_best_betas = (double **) malloc( best_r * sizeof( double * ) );
for(i=0; i < best_r; i++)
	final_best_betas[i] = (double *) malloc(
				p * sizeof(double) );
final_best_scales = (double *) malloc( best_r * sizeof(double) );
k = best_r * groups;
best_betas = (double **) malloc( k * sizeof( double * ) );
best_scales = (double *) malloc( k * sizeof( double ) );
for(i=0; i < k; i++)
		best_betas[i] = (double*) malloc( p * sizeof(double) );
x = (double**) malloc( n * sizeof(double *) );
for(i=0; i<n; i++)
	x[i] = (double*) malloc( p * sizeof(double) );
k = n_group * groups;
indices = (int *) malloc( k * sizeof(int) );
xsamp = (double**) malloc( k * sizeof(double *) );
ysamp = (double*) malloc( k * sizeof(double) );
for(i=0;i<k;i++)
	xsamp[i] = (double*) malloc( p * sizeof(double) );
for(i=0;i<n;i++)
	for(j=0;j<p;j++)
		x[i][j]=X[j*n+i];
/* assume that n > 2000
// k = n_group * groups
// set the seed
*/
srand((long)37);
/* get a sample of k indices */
sample_n_outof_N(k, n-1, indices);
/* get the sampled design matrix and response */
for(i=0;i<k;i++) {
	for(j=0;j<p;j++)
		xsamp[i][j] = x[indices[i]][j];
		ysamp[i] = y[indices[i]];
};
/* now we go through the groups and get the
// *bbest_r best betas for each group
*/
for(i=0; i<groups; i++) {
	fast_s_with_memory(xsamp+i*n_group, ysamp+i*n_group,
				&n_group, pp, NN, K,
				bbest_r, bb, rrhoc,
				best_betas+i*best_r,
				best_scales+i*best_r);
};
/* now  iterate (refine) these "best_r * groups"
// best betas in the (xsamp,ysamp) sample
// with kk C-steps
// and keep only the "best_r" best ones
*/
best_sc = INFI;
pos_worst_scale = conv = 0;
for(i=0; i < best_r; i++)
	final_best_scales[i] = INFI;
worst_sc = INFI;
/* set the matrix to zero */
reset_mat(final_best_betas, best_r, p);
k = n_group * groups;
for(i=0; i< (best_r * groups) ; i++) {
	refine_fast_s(xsamp, ysamp, weights, k, p, res,
			tmp, tmp2, tmp_mat, tmp_mat2,
			best_betas[i], kk, conv, b, rhoc,
			best_scales+i, beta_ref,
			&sc);
	if ( loss_rho(res, worst_sc, k, p, rhoc) < b )  {
		/* scale will be better */
		sc = find_scale(res, b, rhoc, sc, k, p);
		k2 = pos_worst_scale;
		final_best_scales[ k2 ] = sc;
		for(j=0;j<p;j++)
			final_best_betas[k2][j] = beta_ref[j];
		pos_worst_scale = find_max(final_best_scales, best_r);
		worst_sc = final_best_scales[pos_worst_scale];
	};
};
/* now iterate the best "best_r"
// betas in the whole sample (until convergence if possible)
*/
best_sc = INFI;
conv = 1;
for(i=0; i<best_r; i++) {
	refine_fast_s(x, y, weights, n, p, res,
			tmp, tmp2, tmp_mat, tmp_mat2,
			final_best_betas[i], kk, conv, b, rhoc,
			final_best_scales+i, beta_ref,
			&aux);
	if(aux < best_sc) {
			*sscale = best_sc = aux;
			for(j=0;j<p;j++)
				bbeta[j] = beta_ref[j];
	};
};
/* Done. Now clean-up. */
for(i=0;i<n;i++) free(x[i]); free(x);
free(best_scales);
k = best_r * groups;
for(i=0;i<k;i++) free( best_betas[i] );
free(best_betas); free(indices); free(ysamp);
k = n_group * groups;
for(i=0;i<k;i++) free(xsamp[i]);
free(xsamp); free(tmp); free(tmp2);
for(i=0;i<p;i++) {
	free(tmp_mat[i]);
	free(tmp_mat2[i]);
};
free(tmp_mat); free(tmp_mat2); free(weights);
for(i=0;i<best_r;i++)
	free(final_best_betas[i]);
free(final_best_betas);
free(final_best_scales);
free(res);
free(beta_ref);

}

void fast_s_with_memory(double **x, double *y,
		int *nn, int *pp, int *NN, int *K,
		int *bbest_r, double *bb,
		double *rrhoc,
		double **best_betas, double *best_scales)
{
/*
// same as fast_s, but it returns the best.r best
// betas, and their associated scales
// useful for the adjustment for large "n"
//
// x an n x p design matrix (including intercept if appropriate)
// y and n vector
// *nn = n, *pp = p
// *NN = number of re-sampling candidates to be taken
// *K = number of refining steps for each candidate
// *bbest_r = number of (refined) to be retained for
// 					full iteration
// 	*bb = right-hand side of the S-equation
// 	*rrhoc  = tuning constant of the \rho function
// 	*bbeta  = returning fast-S estimator
// 	*sscale = returning associated residual scale
*/
void refine_fast_s(double **x, double *y, double *weights,
			int n, int p, double *res,
			double *tmp, double *tmp2,
			double **tmp_mat, double **tmp_mat2,
			double *beta_cand, int kk,
			int conv, double b, double rhoc,
			double *is,
			double *beta_ref, double *scale);
// void disp_vec(double *a, int n);
// void disp_mat(double **a, int n, int m);
int find_max(double *a, int n);
double find_scale(double *r, double b, double rhoc,
			double initial_scale, int n, int p);
double loss_rho(double *r, double s, int n, int p, double rhoc);
double vecprime_vec(double *a, double *b, int n);
void sample_n_outof_N(int n, int N, int *x);
int lu(double **a,int *P, double *x);
register int i,j,k,no_resamples;
int n = *nn, p = *pp, Nres = *NN, kk = *K;
int *b_i, flag, conv;
double **x_samp, *beta_cand, *beta_ref, *res, aux;
double b = *bb, rhoc = *rrhoc, sc, worst_sc = INFI;
double *weights;
int best_r = *bbest_r, pos_worst_scale;
double *tmp, *tmp2, **tmp_mat2, **tmp_mat;

for(i=0;i<best_r;i++) {
		best_scales[i] = INFI;
};
pos_worst_scale = 0;
res       = (double *) malloc( n * sizeof(double) );
tmp       = (double *) malloc( n * sizeof(double) );
tmp2      = (double *) malloc( n * sizeof(double) );
weights   = (double *) malloc( n * sizeof(double) );
beta_cand = (double *) malloc( p * sizeof(double) );
beta_ref  = (double *) malloc( p * sizeof(double) );
b_i       = (int *) malloc( n * sizeof(int) );
x_samp    = (double **) malloc( n * sizeof(double*) );
tmp_mat   = (double **) malloc( p * sizeof(double*) );
tmp_mat2  = (double **) malloc( p * sizeof(double*) );
for(i=0;i<n;i++) {
	x_samp[i]  = (double *) malloc( (p+1) * sizeof(double) );
}
for(i=0;i<p;i++) {
	tmp_mat[i] = (double *) malloc( p * sizeof(double) );
	tmp_mat2[i] = (double *) malloc( (p+1) * sizeof(double) );
};
/* flag for refine(), conv == 0 means do k refining steps
// conv == 1 means refine until convergence
*/
conv = 0;
aux = -1.0;
/* resampling approximation  */
for(i=0;i<Nres;i++) {
	flag = 1;
	/* find a candidate */
	no_resamples = 0;
	while( flag == 1) {
		if( (++no_resamples) > MAX_NO_RESAMPLES ) {
			// Rprintf("\nToo many singular resamples\nAborting\n\n");
			return;
		};
		/* take a sample of the indices  */
		sample_n_outof_N(p,n-1,b_i);
		/* build the submatrix */
		for(j=0;j<p;j++) {
			for(k=0;k<p;k++)
				x_samp[j][k]=x[b_i[j]][k];
			x_samp[j][p]=y[b_i[j]];
			};
		/* solve the system, lu = 1 means
		// matrix is singular
		*/
		flag = lu(x_samp,pp,beta_cand);
	};
	refine_fast_s(x, y, weights, n, p, res,
			tmp, tmp2, tmp_mat, tmp_mat2,
			beta_cand, kk, conv, b, rhoc,
			&aux, beta_ref, &sc);
	if ( loss_rho(res, worst_sc, n, p, rhoc) < b )  {
		/* scale will be better */
		sc = find_scale(res, b, rhoc, sc, n, p);
		k = pos_worst_scale;
		best_scales[ k ] = sc;
		for(j=0;j<p;j++)
			best_betas[k][j] = beta_ref[j];
		pos_worst_scale = find_max(best_scales, best_r);
		worst_sc = best_scales[pos_worst_scale];
	};

};
/* this function returns all the best_betas
// and best_scales
*/
free(tmp); free(tmp2);
free(res); free(weights); free(beta_cand);
free(beta_ref); free(b_i);
for(i=0; i<n; i++) {
	free(x_samp[i]);
};
for(i=0;i<p;i++) {
	free(tmp_mat[i]);
	free(tmp_mat2[i]);
};
free(x_samp); free(tmp_mat); free(tmp_mat2);
}

void fast_s(double *X, double *y,
		int *nn, int *pp, int *NN, int *K,
		int *bbest_r, double *bb,
		double *rrhoc, double *bbeta, double *sscale)
{
/*
// X an n x p design matrix (including intercept if appropriate)
// y and n vector
// *nn = n, *pp = p
// *NN = number of re-sampling candidates to be taken
// *K = number of refining steps for each candidate
// *bbest_r = number of (refined) to be retained for
// 					full iteration
// 	*bb = right-hand side of the S-equation
// 	*rrhoc  = tuning constant of the \rho function
// 	*bbeta  = returning fast-S estimator
// 	*sscale = returning associated residual scale
*/
void refine_fast_s(double **x, double *y, double *weights,
			int n, int p, double *res,
			double *tmp, double *tmp2,
			double **tmp_mat, double **tmp_mat2,
			double *beta_cand, int kk,
			int conv, double b, double rhoc,
			double *is,
			double *beta_ref, double *scale);
// void disp_vec(double *a, int n);
// void disp_mat(double **a, int n, int m);
int find_max(double *a, int n);
double find_scale(double *r, double b, double rhoc,
			double initial_scale, int n, int p);
double loss_rho(double *r, double s, int n, int p, double rhoc);
double vecprime_vec(double *a, double *b, int n);
void sample_n_outof_N(int n, int N, int *x);
int lu(double **a,int *P, double *x);
register int i,j,k;
int n = *nn, p = *pp, Nres = *NN, kk = *K, no_resamples;
int *b_i, flag, conv;
double **x, **x_samp, *beta_cand, *beta_ref, *res, aux;
double b = *bb, rhoc = *rrhoc, sc, worst_sc = INFI;
double **best_betas, *best_scales, *weights;
int best_r = *bbest_r, pos_worst_scale;
double *tmp, *tmp2, **tmp_mat2, **tmp_mat, best_sc;
best_betas = (double **) malloc( best_r * sizeof(double*) );
best_scales = (double *) malloc( best_r * sizeof(double) );
for(i=0;i<best_r;i++) {
		best_betas[i] = (double*) malloc( p * sizeof(double) );
		best_scales[i] = INFI; };
pos_worst_scale = 0;
res       = (double *) malloc( n * sizeof(double) );
tmp       = (double *) malloc( n * sizeof(double) );
tmp2      = (double *) malloc( n * sizeof(double) );
weights   = (double *) malloc( n * sizeof(double) );
beta_cand = (double *) malloc( p * sizeof(double) );
beta_ref  = (double *) malloc( p * sizeof(double) );
b_i       = (int *) malloc( n * sizeof(int) );
x         = (double **) malloc( n * sizeof(double*) );
x_samp    = (double **) malloc( n * sizeof(double*) );
tmp_mat   = (double **) malloc( p * sizeof(double*) );
tmp_mat2  = (double **) malloc( p * sizeof(double*) );
for(i=0;i<n;i++) {
	x[i]       = (double *) malloc( p * sizeof(double) );
	x_samp[i]  = (double *) malloc( (p+1) * sizeof(double) );
};
for(i=0;i<p;i++) {
	tmp_mat[i] = (double *) malloc( p * sizeof(double) );
	tmp_mat2[i] = (double *) malloc( (p+1) * sizeof(double) );
};
for(i=0;i<n;i++)
        for(j=0;j<p;j++)
                x[i][j]=X[j*n+i];
/* set the seed  */
srand((long)37);
/* flag for refine(), conv == 0 means do k refining steps
// conv == 1 means refine until convergence
*/
conv = 0;
aux = -1.0;
/* resampling approximation  */

for(i=0;i<Nres;i++) {
	flag = 1;
	/* find a candidate */
	no_resamples=0;
	while( flag == 1) {
		if( (++no_resamples) > MAX_NO_RESAMPLES ) {
			// Rprintf("\nToo many singular resamples\nAborting\n\n");
			return;
		};
		/* take a sample of the indices  */
		sample_n_outof_N(p,n-1,b_i);
		/* build the submatrix */
		for(j=0;j<p;j++) {
			for(k=0;k<p;k++)
				x_samp[j][k]=x[b_i[j]][k];
			x_samp[j][p]=y[b_i[j]];
			};
		/* solve the system, lu = 1 means
		// matrix is singular
		*/
		flag = lu(x_samp,pp,beta_cand);
	};
	// Rprintf("\n");
	// for(j=0;j<p;j++) Rprintf("%d ", b_i[j]);
	// Rprintf("\n");
	/* improve the re-sampling candidate */
	refine_fast_s(x, y, weights, n, p, res,
			tmp, tmp2, tmp_mat, tmp_mat2,
			beta_cand, kk, conv, b, rhoc,
			&aux, beta_ref, &sc);
	if( fabs(sc) < ZERO) {
		*sscale = sc;
		for(j=0;j<p;j++) bbeta[j] = beta_cand[j];
		free(best_scales); free(tmp); free(tmp2);
		free(res); free(weights); free(beta_cand);
		free(beta_ref); free(b_i);
		for(i=0;i<best_r;i++)
			free(best_betas[i]);
		free(best_betas);
		for(i=0; i<n; i++) {
			free(x[i]);
			free(x_samp[i]);
		};
		for(i=0; i<p; i++) {
			free(tmp_mat[i]);
			free(tmp_mat2[i]);
		};
		free(x); free(x_samp); free(tmp_mat); free(tmp_mat2);
		return;
	};
	if ( loss_rho(res, worst_sc, n, p, rhoc) < b )  {
		/* scale will be better */
		sc = find_scale(res, b, rhoc, sc, n, p);
		k = pos_worst_scale;
		best_scales[ k ] = sc;
		for(j=0;j<p;j++)
			best_betas[k][j] = beta_ref[j];
		pos_worst_scale = find_max(best_scales, best_r);
		worst_sc = best_scales[pos_worst_scale];
	};
};
/* now look for the very best */
best_sc = INFI;
conv = 1;
for(i=0; i<best_r; i++) {
	refine_fast_s(x, y, weights, n, p, res,
			tmp, tmp2, tmp_mat, tmp_mat2,
			best_betas[i], kk, conv, b, rhoc,
			best_scales+i, beta_ref,
			&aux);
	if(aux < best_sc) {
			*sscale = best_sc = aux;
			for(j=0;j<p;j++)
				bbeta[j] = beta_ref[j];
	};
};
free(best_scales); free(tmp); free(tmp2);
free(res); free(weights); free(beta_cand);
free(beta_ref); free(b_i);
for(i=0;i<best_r;i++)
	free(best_betas[i]);
free(best_betas);
for(i=0; i<n; i++) {
	free(x[i]);
	free(x_samp[i]);
};
for(i=0; i<p; i++) {
	free(tmp_mat[i]);
	free(tmp_mat2[i]);
};
free(x); free(x_samp); free(tmp_mat); free(tmp_mat2);
}

void refine_fast_s(double **x, double *y, double *weights,
			int n, int p, double *res,
			double *tmp, double *tmp2,
			double **tmp_mat, double **tmp_mat2,
			double *beta_cand, int kk,
			int conv, double b, double rhoc,
			double *is,
			double *beta_ref, double *scale)
{
/*
// weights = vector of length n
// res = vector of length n
// x = matrix with the data
// y = vector with responses
// tmp = aux vector of length n
// tmp2 = aux vector of length n
// tmp_mat = aux matrix p x p
// tmp_mat2 = aux matrix p x (p+1)
*/
void fast_s_irwls(double **x, double *y,
		double *weights, int n, int p, double *beta_ref,
		double **tmp_mat, double *tmp, double *tmp2);
double norm_diff(double *x, double *y, int n);
double norm(double *x, int n);
int lu(double **a,int *P, double *x);
void get_weights_rhop(double *r, double s, int n,
		double rhoc, double *w);
void r_sum_w_x(double **x, double *w, int n, int p,
			double *tmp,
			double *sum);
void r_sum_w_x_xprime(double **x, double *w, int n, int p,
			double **tmp, double **ans);
double loss_rho(double *r, double scale, int n, int p, double rhoc);
double MAD(double *a, int n, int center, double *tmp,
			double *tmp2);
double vecprime_vec(double *a, double *b, int n);
register int i,j;
int zeroes=0;
double initial_scale = *is, s0;

for(j=0;j<n;j++)
		if( fabs(res[j] = y[j] - vecprime_vec(x[j], beta_cand, p))
				< ZERO ) zeroes++;
/* if "perfect fit", return it with a 0 assoc. scale */
/* if( zeroes > (((double)n + (double)p)/2.) ) */
if( zeroes > ((double)n /2.) )
{
	// Rprintf("\nToo many zeroes, refine_fast_s\n");
	for(i=0;i<p;i++) beta_ref[i] = beta_cand[i];
	*scale = 0.0;
	return;
};

if( initial_scale < 0.0 )
	initial_scale = MAD(res, n, 0, tmp, tmp2);

s0 = initial_scale;
if( conv > 0 )
		kk = MAX_ITER_FAST_S;
if(kk > 0) {
for(i=0; i < kk; i++) {

	/* one step for the scale */
	s0 = s0 * sqrt( loss_rho(res,
					s0, n, p, rhoc) / b );
	/* compute weights for IRWLS */
	get_weights_rhop(res, s0, n, rhoc, weights);
	/* compute the matrix for IRWLS */
	r_sum_w_x_xprime(x, weights, n, p, tmp_mat,
				tmp_mat2);
	/* compute the vector for IRWLS */
	for(j=0; j<n; j++)
			weights[j] = weights[j] * y[j];
	r_sum_w_x(x, weights, n, p, tmp, tmp2);
	for(j=0; j<p; j++)
		tmp_mat2[j][p] = tmp2[j];
	/* solve the system for IRWLS */
	lu(tmp_mat2, &p, beta_ref);
	/* check for convergence? */
	if(conv > 0) {
		if(norm_diff(beta_cand, beta_ref, p) /
					norm(beta_cand, p) < EPSILON ) {
			// Rprintf("\nRelative norm less than EPSILON\n");
			break;
		};
	};
	for(j=0;j<n;j++)
		res[j] = y[j] - vecprime_vec(x[j], beta_ref , p);
	for(j=0; j<p; j++)
		beta_cand[j] = beta_ref[j];
};
};
*scale = s0;
}


void fast_s_irwls(double **x, double *y,
		double *weights, int n, int p, double *beta_ref,
		double **tmp_mat, double *tmp, double *tmp2)
{
void mat_prime_vec(double **a, double *b, double *c, int n, int m);
void mat_prime_mat_w(double **a, double *w, double **c,
		int n, int m);
register int i;
for(i=0;i<n;i++) tmp[i] = weights[i] * y[i];
mat_prime_vec(x, tmp, tmp2, n, p);
mat_prime_mat_w(x, weights, tmp_mat, n, p);
for(i=0;i<p;i++) tmp_mat[i][p] = tmp2[i];
lu(tmp_mat, &p, beta_ref);
}

void get_weights_rhop(double *r, double s,
		int n,
		double rhoc, double *w)
{
register int i;
double a;
for(i=0;i<n;i++) {
	a = r[i] / s / rhoc;
	if( fabs(a) > 1 )
			w[i] = 0;
	else
			w[i] = (1. - a*a) * (1. - a*a);
};
}

double find_scale(double *r, double b, double rhoc,
			double initial_scale, int n, int p)
{
double loss_rho(double *r, double scale, int n, int p, double rhoc);
int max_it = MAX_ITER_FIND_SCALE, it = 0;
double e = 1, scale;

while( (++it < max_it) && (fabs(e) > EPSILON) )
{
	scale = initial_scale * sqrt(
			loss_rho(r, initial_scale, n, p, rhoc) /
					b ) ;
	e = fabs( scale / initial_scale - 1);
	initial_scale = scale;
};
return(scale);
}

int find_max(double *a, int n)
{
register int i;
int k=0;
double tmp = a[0];
if(n==1) return(0);
else {
	for(i=1;i<n;i++)
		if(a[i] > tmp) {
				tmp = a[i];
				k = i;
		};
};
return(k);
}

void r_sum_w_x(double **x, double *w, int n, int p,
			double *tmp,
			double *sum)
{
/*
// given a matrix x (n x p) and a vector w of n
// weights, it computes the vector
// \sumin w_i x_i
// need space for p doubles in *tmp
*/
void scalar_vec(double *a, double b, double *c, int n);
void sum_vec(double *a, double *b, double *c, int n);
void reset_vec(double *a, int n);
register int i;


reset_vec(sum, p);

for(i=0; i<n; i++) {
	scalar_vec(x[i], w[i], tmp, p);
	sum_vec(sum , tmp , sum, p);
};


}


void r_sum_w_x_xprime(double **x, double *w, int n, int p,
			double **tmp, double **ans)
{
/*
// given a matrix x (n x p) and a vector w of n
// weights, it computes the matrix
// \sumin w_i x_i x_i'
// need space for p x p "doubles" in tmp
*/
void sum_mat(double **a, double **b, double **c, int n, int m);
void matias_vec_vec(double **a, double *v1, double *v2, int n);
void scalar_mat(double **a, double b, double **c, int n, int m);
void reset_mat(double **a, int n, int m);
register int i;

reset_mat(ans, p, p);

for(i=0; i<n; i++) {
		matias_vec_vec(tmp, x[i], x[i], p);
		scalar_mat(tmp, w[i], tmp, p, p);
		sum_mat(ans, tmp, ans, p, p);
};

}


double loss_rho(double *r, double scale, int n, int p, double rhoc)
{
double Chi(double x, double c);
register int i;
double s = 0;
for(i=0;i<n;i++)
		s += Chi(r[i]/scale, rhoc);
/* return(s / ( (double) n - (double) p ) ); */
return(s / (double) n );
}

double MAD(double *a, int n, int center, double *b,
			double *tmp)
{
/* if center == 0 then do not center */
double median_abs(double *, int , double *);
double median(double *,int, double *);
int i;
double med,q;
if(center > 0) {
	med = median(a,n,b);
} else {
	med = 0.0;
};
for(i=0;i<n;i++)
	b[i] = a[i] - med;
q = median_abs(b,n,tmp) * 1.4826;
return(q);
}

double median(double *x, int n, double *aux)
{
double kthplace(double *,int,int);
double t;
register int i;
for(i=0;i<n;i++) aux[i]=x[i];
if ( (n/2) == (double) n / 2 )
	t = ( kthplace(aux,n,n/2) + kthplace(aux,n,n/2+1) ) / 2 ;
else	t = kthplace(aux,n, n/2+1 ) ;
return(t);
}

double median_abs(double *x, int n, double *aux)
{
double kthplace(double *,int,int);
double t;
register int i;
for(i=0;i<n;i++) aux[i]=fabs(x[i]);
if ( (n/2) == (double) n / 2 )
	t = ( kthplace(aux,n,n/2) + kthplace(aux,n,n/2+1) ) / 2 ;
else 	t = kthplace(aux,n, n/2+1 ) ;
return(t);
}

int rwls(double **a, int n, int p,
			double *estimate,
			double *i_estimate,
			double scale, double epsilon,
			int max_it, double Psi_constant
			)
{
/* a <- matrix n x (p+1) (n rows and p+1 columns) where the data's stored
 * res <- vector n of residuals
 * b <- auxiliar matrix to store A'A | A'c
 */
int lu(double **, int *, double *);
// void disp_vec(double*,int);
double norm_diff(double *, double *, int);
double Psi_reg(double, double);
double Loss_Tukey(double*, int, double);

double **b,s,*beta1, *beta2, *beta0, *weights, *resid;
double r; // ,loss1,loss2,lambda;
int iterations=0; //, iter_lambda;
register int i,j,k;
if ( (b = (double **) malloc ( p * sizeof(double *) ) )==NULL )
	{// Rprintf("\nRun out of memory in rwls\n");
		exit(1); };
for (i=0;i<p;i++)
	if ( (b[i] = (double *) malloc ( (p+1) * sizeof(double) ) )==NULL )
		{// Rprintf("\nRun out of memory in rwls\n");
			exit(1); };
beta1 = (double *) malloc( p * sizeof(double) );
beta2 = (double *) malloc( p * sizeof(double) );
beta0 = (double *) malloc( p * sizeof(double) );
weights = (double *) malloc( n * sizeof(double) );
resid = (double *) malloc( n * sizeof(double) );
for(i=0;i<p;i++)
	beta2[i] = (beta1[i]=i_estimate[i]) + 1;
/* main loop */
while( (norm_diff(beta1,beta2,p) > epsilon) &&
	( ++iterations < max_it ) ) {
for(i=0;i<n;i++) {
	s=0;
	for(j=0;j<p;j++)
		s += a[i][j] * beta1[j];
	r = a[i][p]- s;
	if(fabs(r/scale)<1e-7)
		weights[i] = 1.0 / scale / Psi_constant;
                else
        	weights[i] = Psi_reg(r/scale, Psi_constant) / (r/scale);
};
for(j=0;j<p;j++) beta2[j]=beta1[j];
// /* get the residuals and loss for beta2 */
// for(i=0;i<n;i++)
// 	{ s = 0;
// 	for(j=0;j<p;j++)
// 		s += a[i][j] * beta2[j];
// 	resid[i] = (a[i][p] - s)/scale;
// 	};
// loss2 = Loss_Tukey(resid,n,Psi_constant);
/* S+ version of the following code
 * A <- matrix(0, p, p)
 * Y <- rep(0, p)
 * for(i in 1:n) {
 * A <- A + w[i] * a[i,0:(p-1)] %*% t(a[i,0:(p-1)])
 * Y <- Y + w[i] * a[i,0:(p-1)] * a[i,p]
 * }
 * beta1 <- solve(A, Y)
 */
for(j=0;j<p;j++)
	for(k=0;k<=p;k++)  {
		b[j][k]=0.0;
		 for(i=0;i<n;i++)
			b[j][k] += a[i][j] * a[i][k] * weights[i];
		};
lu(b,&p,beta1);
/* is beta1 good enough? */
/* get the residuals and loss for beta1 */
// for(i=0;i<n;i++)
// 	{ s = 0;
// 	for(j=0;j<p;j++)
// 		s += a[i][j] * beta1[j];
// 	resid[i] = (a[i][p] - s)/scale;
// 	};
// loss1 = Loss_Tukey(resid,n,Psi_constant);
// for(j=0;j<p;j++) beta0[j] = beta1[j];
// lambda = 1.;
// iter_lambda=0;
// while( ( loss1 > loss2 ) ) {
// 	Rprintf("%f - %f\n", loss1, loss2);
// 	lambda /= 2.;
// 	for(j=0;j<p;j++)
// 		beta0[j] = (1 - lambda) * beta2[j] + lambda * beta1[j];
// 	/* get the residuals and loss for beta0 */
// 	for(i=0;i<n;i++)
// 		{ s = 0;
// 		for(j=0;j<p;j++)
// 			s += a[i][j] * beta0[j];
// 		resid[i] = a[i][p] - s;
// 		};
// 	loss1 = Loss_Tukey(resid,n,Psi_constant);
// 	if( ++iter_lambda > 10) {
// 		/* Rprintf("\nStuck in local search. Rwls. ");
// 		Rprintf("%f - %f\n",loss1,loss2);  */
// 		loss1 = loss2; /* force the exit */
// 		for(j=0;j<p;j++) beta0[j] = beta2[j];
// 		/* return(1); */
// 	};
// }; /* end while(loss2 <= loss1 ) */
// for(j=0;j<p;j++) beta1[j] = beta0[j];
}; /* end while(norm_diff(...)   */
// for(j=0;j<p;j++) estimate[j]=beta0[j];
for(j=0;j<p;j++) estimate[j]=beta1[j];
free(weights);free(beta1);free(beta2);
free(beta0);free(resid);
for(i=0;i<p;i++) free(b[i]);
free(b);
if( iterations == max_it )
	return 1;
	else
	return 0;
}

void reset_mat(double **a, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<m;j++)
		a[i][j] = 0.0;
}

double norm(double *x, int n)
{
double s = 0;
register int i;
for(i=0; i<n; i++) s += x[i] * x[i];
return(sqrt(s));
}

void mat_prime_vec(double **a, double *b, double *c, int n, int m)
{
register int i,j;
for(i=0;i<m;i++)
	for(c[i]=0,j=0;j<n;j++) c[i] += a[j][i] * b[j];
}

void mat_prime_mat_w(double **a, double *w, double **c,
		int n, int m)
{
register int i,j,k;
for(i=0;i<m;i++) {
	for(j=0;j<m;j++) {
		c[i][j] = 0;
		for(k=0;k<n;k++)
			c[i][j] += a[k][i] * w[k] * a[k][j];
	};
};
}

void reset_vec(double *a, int n)
{
register int i;
for(i=0;i<n;i++) a[i] = 0.0;
}


