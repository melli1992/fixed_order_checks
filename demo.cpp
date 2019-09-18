/*
	demo-c.c
		test program for the Cuba library
		last modified 13 Mar 15 th
*/

#include  <iostream>
#include <stdlib.h>
#include <math.h>
#include "cuba.h"

using namespace std;

static inline cubareal Sq(cubareal x) {
  return x*x;
}



static int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

#define x xx[0]
#define y xx[1]
#define z xx[2]
#define zz xx[3]
#define f ff[0]
#define f2 ff[1]

#ifndef FUN
#define FUN 3
#endif
#define rsq (Sq(x) + Sq(y) + Sq(z))

#if FUN == 1
  f = sin(x)*cos(y)*exp(z);
#elif FUN == 2
  f = 1/(Sq(x + y) + .003)*cos(y)*exp(z);
#elif FUN == 3
  f = 1/(3.75 - cos(M_PI*x) - cos(M_PI*y) - cos(M_PI*z));
  f2 = sin(x)*cos(y)*exp(z);
#elif FUN == 4
  f = fabs(rsq - .125);
#elif FUN == 5
  f = exp(-rsq);
#elif FUN == 6
  f = 1/(1 - x*y*z + 1e-10);
#elif FUN == 7
  f = sqrt(fabs(x - y - z));
#elif FUN == 8
  f = exp(-x*y*z);
#elif FUN == 9
  f = Sq(x)/(cos(x + y + z + 1) + 5);
#elif FUN == 10
  f = (x > .5) ? 1/sqrt(x*y*z + 1e-5) : sqrt(x*y*z);
#else
  f = (rsq < 1) ? 1 : 0;
#endif

  return 0;
}

/*********************************************************************/

int NDIM; //number of dimensions in integral
int NCOMP = 2; //number of components in integral
// integrand should be a function int integrand(ndim, x, ncomp, f, userdata, nvec, core)
// x(ndim, nvec), f(ncomp, nvec)
// nvec = number of samples x that is received, fill array f with integrand values
// returned value of integrand is irrelevant, unless it is -999, then should the integration be aborted
// nvec, userdata and core are optional in integrand function, not in the cuhre and vegas functions
void* USERDATA = NULL; 
int NVEC = 1;
double EPSREL = 1.e-3;
double EPSABS = 1.e-12;
int VERBOSE = 0; //bit 0-1 is verbosity (so 0,1,2,3) , bit 2 = samples, bit 3 can improve smoothing, bit 4 = state file, bit 5 = vegas specific, bit 8 and higher is random number
int LAST = 4;
int SEED = 0; //other choices of seed
int MINEVAL = 0; //min number of evaluations integrand
int MAXEVAL = 50000; //max number of evaluations integrand
char* STATEFILE = NULL; //filename for storing internal state
int* SPIN = NULL;

//Vegas specific 
int NSTART = 1000; //number evaluations to start with
int NINCREASE = 500; //increase in number of integrands
int NBATCH = 1000; //batch size for sampling
int GRIDNO = 0; //keep grids during one integratoin for the next one (if integrands are similar) if gridno > 0


//Cuhre
int KEY = 0; //cubature rule of degree?  

int main() {
  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
  NDIM = 3;
  printf("-------------------- Vegas test --------------------\n");

  Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral, error, prob);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);

  printf("\n-------------------- Cuhre test --------------------\n");

  Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);

  return 0;
}
