#include <cmath>
#include <string>  
#include "DY-num-v01.h"
#include "LHAPDF/LHAPDF.h"
#include "cuba.h"

using namespace std;
using namespace LHAPDF;

#define pbunits 0.389379323*pow(10.,9)
#define fbunits 0.389379323*pow(10.,12)

#define PI 3.141592653589793
#define PI2 9.869604401089358

#define nf 5.
#define CF 4./3.
#define CA 3.
#define Nc 3.

#define zeta2 1.6449340668482264
#define zeta3 1.2020569031595942

#define E 2.718281828459045

#define EulerGamma 0.5772156649015329
#define EEulerGamma 1.781072417990198

/////////////////////////////////////////////////////////////////////////////////////

// parameters for cuba

#define NDIM1 1
#define NDIM2 2
#define NDIM3 3
#define NDIM4 4
#define NCOMP 1
#define NVEC 1
#define NMIN 2
#define USERDATA NULL
#define EPSREL 1.0e-3
#define EPSABS 1.0e-6
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 2000000

#define NSTART 50000
#define NINCREASE 10000
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL

#define NNEW 100000
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0 

/////////////////////////////////////////////////////////////////////////////////////
///
/// Drell-Yan invariant mass distribution
///
/// Amsterdam, 12/04/2019
///
/////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////
///
/// tau * sigma_v, eq. (2.2) and (A.1) of Hamberg at al, Nucl.Phys. B359 (1991) 343-405
///
/////////////////////////////////////////////////////////////////////////////////////

double Prefactor(double &S, double &Q){

// "pbunits" or "fbunits": fixes the normalisation to picobarn or femtobarn

   double aEM = 1./(132.507);

   double r = pbunits*4.*PI*aEM*aEM/3./Nc/Q/Q/S;

   return (r);
}

/////////////////////////////////////////////////////////////////////////////////////
///
/// Construct luminosity
///
/////////////////////////////////////////////////////////////////////////////////////

// luminosity

double lum(double t, double y, double &mu, LHAPDF::PDF* pdf){
   

   //define electric charge vector 
   double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.};
  
   double r = 0.;

   // (q qbar + qbar q) luminosity for DY (lepton pair produced through a photon)

   for(int i = 1; i <= 5; i++){
      r += -eq[i-1]*eq[i-1]*(pdf->xfxQ(i+6,exp(t*log(y)), mu)*pdf->xfxQ(i-6,y*exp(-t*log(y)), mu)+pdf->xfxQ(i-6,exp(t*log(y)), mu)*pdf->xfxQ(i+6,y*exp(-t*log(y)), mu))*log(y)/y;
   }

   return r;
}

// luminosity and derivatives (not needed at the moment)

double* lumi(double t, double y, double &mu, LHAPDF::PDF* pdf){

   // Not used at the moment 
   
   double delta = y*0.00000001;

   double lumi0 = lum(t,y,mu,pdf);
   double lumidp = lum(t,y+delta,mu,pdf); 
   double lumidm = lum(t,y-delta,mu,pdf);   

   double lumid2p = lum(t,y+2.0*delta,mu,pdf); 
   double lumid2m = lum(t,y-2.0*delta,mu,pdf);  

   double* r = new double[3];

   r[0] = lumi0;                                                                // function
   r[1] = (-lumid2p+8.*lumidp-8*lumidm+lumid2m)/(12.*delta);                    // first derivative
   r[2] = (-lumid2p+16.*lumidp-30.*lumi0+16*lumidm-lumid2m)/(12.*delta*delta);  // second derivative

  return r;
}

/////////////////////////////////////////////////////////////////////////////////////
///
/// LO and NLO integrand - terms proportional to delta(1-z)
///
/////////////////////////////////////////////////////////////////////////////////////

// z independent terms

double fsigmaD(double t, double aux, double &Q, double &S, double &muf, double &mur, LHAPDF::PDF* pdf, int ord){  

   // aux is an auxiliariy parameter to make the integration two-dimensional
   // the cuba libraries can integrate only multi-dimensional functions
   
   double tau = Q*Q/S;
   double Qsq = Q*Q;
   double mursq = mur*mur;

   double as = alphasPDF(mur);

   double pre = Prefactor(S,Q);

   double lumit = lum(t,tau,muf,pdf);

   // This is the part of Delta_{V+S} proportional to delta(1-z) 
   // in eq. (B.3) of Hamberg at al, Nucl.Phys. B359 (1991) 343-405
   double DeltaDnlo = (as/4./PI)*CF*(6.*log(Qsq/mursq)+8.*zeta2-16.);
   
   double r;
 
    if      (ord==0){
	 
       r = pre*lumit;
    }
    else if (ord==1){
    
       r = pre*(1.+DeltaDnlo)*lumit;
    }
    else if (ord==2){
    
       r = 0.;
    }
	return r;
}


/////////////////////////////////////////////////////////////////////////////////////
///
/// NLO integrand - z dependent terms 
///
/////////////////////////////////////////////////////////////////////////////////////

double fsigmaZ(double t, double z, double &Q, double &S, double &muf, double &mur, LHAPDF::PDF* pdf, int ord){  

   double tau = Q*Q/S;
   double Qsq = Q*Q;
   double mursq = mur*mur;

   double as = alphasPDF(mur);

   double pre = Prefactor(S,Q);

   double lumiz = lum(t,tau/z,muf,pdf);
   double lumit = lum(t,tau,muf,pdf);

   // This is the part of Delta_{V+S} proportional to plus distributions 
   // in eq. (B.3) + Delta_{H} in eq. (B.4) of Hamberg at al, Nucl.Phys. B359 (1991) 343-405
   double DeltaZnlo = (as/4./PI)*CF*((16.*log(1.-z)/(1.-z)+8.*log(Qsq/mursq)/(1.-z))*(lumiz/z-lumit) 
				 +(-4.*(1.+z)*log(Qsq/mursq)-8.*(1.+z)*log(1.-z)-4.*(1.+z*z)/(1.-z)*log(z))*lumiz/z);

   double r;
	
    if      (ord==0){
	 
       r = 0.;
    }
    else if (ord==1){
    
       r = pre*DeltaZnlo;
    }
    else if (ord==2){
    
       r = 0.;
    }
    return r;
}


/////////////////////////////////////////////////////////////////////////////////////
///
/// LO and NLO invariant mass distribution - terms proportional to delta(1-z)
///
/////////////////////////////////////////////////////////////////////////////////////

double fmassIntD(double t, double aux, void *params){

   double Q = static_cast<double *>(params)[0];
   double S = static_cast<double *>(params)[1];
   double muf = static_cast<double *>(params)[2];  
   double mur = static_cast<double *>(params)[3]; 
   LHAPDF::PDF pdf = static_cast<LHAPDF::PDF *>(params)[4];
   int ord = static_cast<double *>(params)[5];  

   double r = fsigmaD(t,aux,Q,S,muf,mur,pdf,ord); 

   return (r);
}

static int massIntD(const int *ndim, const double xx[], const int *ncomp, double ff[], void *params){

   #define x1 xx[0]
   #define x2 xx[1]
   #define f ff[0]

   f = fmassIntD(x1, x2, params);
  
   return 0;
}

double* vegasMassD(double &Q, double &S, double &muf, double &mur, LHAPDF::PDF* pdf, int ord, int verbose){

   double params[6];
  
   params[0] = Q;
   params[1] = S;
   params[2] = muf;  
   params[3] = mur; 
   params[4] = pdf;
   params[5] = ord;
  
   int comp, nregions, neval, fail;
   int spin = -1; 
   double integral[NCOMP], error[NCOMP], prob[NCOMP];

   const char *env = getenv("CUBAVERBOSE");
   if( env ) verbose = atoi(env);

   Vegas(NDIM2, NCOMP, massIntD, params, 
      NVEC, EPSREL, EPSABS, verbose, SEED,
      MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
      GRIDNO, STATEFILE, &spin, &neval, &fail, integral, error, prob);

   double* r = new double[2];

   r[0] = integral[0];
   r[1] = error[0];

   return(r);
}

double* suaveMassD(double &Q, double &S, double &muf, double &mur, LHAPDF::PDF* pdf, int ord, int verbose){

   double params[6];
  
   params[0] = Q;
   params[1] = S;
   params[2] = muf;  
   params[3] = mur; 
   params[4] = pdf;
   params[5] = ord;
  
   int comp, nregions, neval, fail;
   int spin = -1; 
   double integral[NCOMP], error[NCOMP], prob[NCOMP];

   const char *env = getenv("CUBAVERBOSE");
   if( env ) verbose = atoi(env);

   Suave(NDIM2, NCOMP, massIntD, params,
      NVEC, EPSREL, EPSABS, verbose, SEED,
      MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS, STATEFILE, 
      &spin, &nregions, &neval, &fail, integral, error, prob);

   double* r = new double[2];

   r[0] = integral[0];
   r[1] = error[0];

   return(r);
}

double* cuhreMassD(double &Q, double &S, double &muf, double &mur, LHAPDF::PDF* pdf, int ord, int verbose){

   double params[6];
  
   params[0] = Q;
   params[1] = S;
   params[2] = muf;  
   params[3] = mur; 
   params[4] = pdf;
   params[5] = ord;

   int comp, nregions, neval, fail;
   int spin = -1; 
   double integral[NCOMP], error[NCOMP], prob[NCOMP];

   const char *env = getenv("CUBAVERBOSE");
   if( env ) verbose = atoi(env);

   Cuhre(NDIM2, NCOMP, massIntD, params,
      NVEC, EPSREL, EPSABS, verbose,
      MINEVAL, MAXEVAL, KEY, STATEFILE, &spin,
      &nregions, &neval, &fail, integral, error, prob);

   double* r = new double[2];

   r[0] = integral[0];
   r[1] = error[0];

   return(r);
}

/////////////////////////////////////////////////////////////////////////////////////
///
/// NLO invariant mass distribution - z dependent terms 
///
/////////////////////////////////////////////////////////////////////////////////////

double fmassIntZ(double t, double z, void *params){

   double Q = static_cast<double *>(params)[0];
   double S = static_cast<double *>(params)[1];
   double muf = static_cast<double *>(params)[2];  
   double mur = static_cast<double *>(params)[3]; 
   LHAPDF::PDF* pdf = static_cast<double *>(params)[4];
   int ord = static_cast<double *>(params)[5];  
  
   double tau = Q*Q/S;
                 
   // exponential mapping of z: j = Jacobian

   double r;  
   double j = -log(tau)*pow(E,log(tau)*z);

   r = j*fsigmaZ(t,pow(E,log(tau)*z),Q,S,muf,mur,pdf,ord); 
   
   return (r);
}

static int massIntZ(const int *ndim, const double xx[], const int *ncomp, double ff[], void *params){

   #define x1 xx[0]
   #define x2 xx[1]
   #define f ff[0]

   f = fmassIntZ(x1, x2, params);
  
   return 0;
}
/*
double* vegasMassZ(double &Q, double &S, double &muf, double &mur, int &member, int ord, int verbose){

   double params[6];
  
   params[0] = Q;
   params[1] = S;
   params[2] = muf;  
   params[3] = mur; 
   params[4] = member;
   params[5] = ord;
  
   int comp, nregions, neval, fail;
   int spin = -1; 
   double integral[NCOMP], error[NCOMP], prob[NCOMP];

   const char *env = getenv("CUBAVERBOSE");
   if( env ) verbose = atoi(env);

   Vegas(NDIM2, NCOMP, massIntZ, params, 
      NVEC, EPSREL, EPSABS, verbose, SEED,
      MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
      GRIDNO, STATEFILE, &spin, &neval, &fail, integral, error, prob);

   double* r = new double[2];

   r[0] = integral[0];
   r[1] = error[0];

   return(r);
}

double* suaveMassZ(double &Q, double &S, double &muf, double &mur, int &member, int ord, int verbose){

   double params[6];
  
   params[0] = Q;
   params[1] = S;
   params[2] = muf;  
   params[3] = mur; 
   params[4] = member;
   params[5] = ord;

   int comp, nregions, neval, fail;
   int spin = -1; 
   double integral[NCOMP], error[NCOMP], prob[NCOMP];

   const char *env = getenv("CUBAVERBOSE");
   if( env ) verbose = atoi(env);

   Suave(NDIM2, NCOMP, massIntZ, params,
      NVEC, EPSREL, EPSABS, verbose, SEED,
      MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS, STATEFILE, 
      &spin, &nregions, &neval, &fail, integral, error, prob);

   double* r = new double[2];

   r[0] = integral[0];
   r[1] = error[0];

   return(r);
}

double* cuhreMassZ(double &Q, double &S, double &muf, double &mur, int &member, int ord, int verbose){

   double params[6];
  
   params[0] = Q;
   params[1] = S;
   params[2] = muf;  
   params[3] = mur; 
   params[4] = member;
   params[5] = ord;

   int comp, nregions, neval, fail;
   int spin = -1; 
   double integral[NCOMP], error[NCOMP], prob[NCOMP];

   const char *env = getenv("CUBAVERBOSE");
   if( env ) verbose = atoi(env);

   Cuhre(NDIM2, NCOMP, massIntZ, params,
      NVEC, EPSREL, EPSABS, verbose,
      MINEVAL, MAXEVAL, KEY, STATEFILE, &spin,
      &nregions, &neval, &fail, integral, error, prob);

   double* r = new double[2];

   r[0] = integral[0];
   r[1] = error[0];

   return(r);
}

/////////////////////////////////////////////////////////////////////////////////////
///
/// invariant mass - sum D and Z terms
///
/////////////////////////////////////////////////////////////////////////////////////

double* fMass(double &Q, double &S, double &muf, double &mur, int &member, int ord, string cuba_method){

   int verbose = 1;

   std::string comp_vegas = "vegas"; 
   std::string comp_suave = "suave"; 
   std::string comp_cuhre = "cuhre"; 

   double* r = new double[2];

     if (cuba_method.compare(comp_vegas) == 0){
     
       double* D = vegasMassD(Q,S,muf,mur,member,ord,verbose);
       double* Z = vegasMassZ(Q,S,muf,mur,member,ord,verbose);

       r[0] = D[0]+Z[0];
       r[1] = sqrt(D[1]*D[1]+Z[1]*Z[1]);

       return (r);      
       }
     else if (cuba_method.compare(comp_suave) == 0){
     
       double* D = suaveMassD(Q,S,muf,mur,member,ord,verbose);
       double* Z = suaveMassZ(Q,S,muf,mur,member,ord,verbose);                      

       r[0] = D[0]+Z[0];
       r[1] = sqrt(D[1]*D[1]+Z[1]*Z[1]);

       return (r);
       }
     else if (cuba_method.compare(comp_cuhre) == 0){
     
       double* D = cuhreMassD(Q,S,muf,mur,member,ord,verbose);
       double* Z = cuhreMassZ(Q,S,muf,mur,member,ord,verbose);

       r[0] = D[0]+Z[0];
       r[1] = sqrt(D[1]*D[1]+Z[1]*Z[1]);

       return (r);
     }
} 


*/
