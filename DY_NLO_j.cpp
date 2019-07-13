#include <cmath>
#include <iostream>
#include <fstream>
#include "cuba.h"
#include "LHAPDF/LHAPDF.h"
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_dilog.h>

using namespace std;
using namespace LHAPDF;

// Integration parameters for CUBA

#define NDIM1 1
#define NDIM2 2
#define NDIM3 3
#define NDIM4 4
#define NCOMP 3
#define NVEC 1
#define NMIN 2
#define USERDATA NULL
#define EPSREL 1.0e-3
#define EPSABS 1.0e-8
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 2000000

#define NSTART 10000
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

// Further definitions
#define PI 3.141592653589793
#define pi 3.141592653589793
#define zeta2 pow(pi,2)/6.
#define CF 4./3.
#define Nc 3.
#define pbunits 0.38937966*pow(10.,9.)


double Q, S, tau, muf, mur, as;

double prefactor(double totMass, double phMass){
    // Given in units of pb/GeV^2, appropriate for dsigma/dQ^2

    double aem = 1./132.507;
    double r =  pbunits*4.*pi*pow(aem,2)/(3.*Nc*pow(phMass,2)*totMass);
    return r;
}

double PDFprod(double a, double b, double scale) {
    double prod=0.;
    double Qf[6] = {-1./3.,2./3.,-1./3.,2./3.,-1./3.,2./3.};
    for(int i=1; i<=5; i++){
        prod += pow(Qf[i-1],2)*(xfx(a,scale)[6+i]/a*xfx(b,scale)[6-i]/b + xfx(a,scale)[6-i]/a*xfx(b,scale)[6+i]/b);
    }
    return prod;
}

double partonLum(double x, double xmin, double xmax, double scale) {
    double jacx = xmax - xmin;
    return 1/(xmin/jacx+x)*PDFprod(xmin+jacx*x,xmin/(xmin+jacx*x),scale);
}


// Functional form for the integrand; note that both the arrays xx and ff may be multidimensional.
// Since the integral over z runs from tau to 1, rather than from 0 to 1, we need a "fixing term" (Dpart2) on top of the typical g(z)[f(z)-f(1)] interal (Dpart1).
// Remaining Log terms are kept in Lterms

static int Dpart1(const int *ndim, const cubareal xx[],
                     const int *ncomp, cubareal ff[], void *userdata) {
    #define z_var xx[1]
    #define x_var xx[0]
    #define dpart1 ff[0]

    double zmax = 1.;
    double zmin = tau;
    double jacz = zmax-zmin;
    double z = zmin + jacz* z_var;
    // NOTE: is the scale M in Hamberg et al. the factorziation scale or renormalization scale?
    // Add for full result: 8*log(pow(Q,2)/pow(muf,2))/(1-z)+
    dpart1 = prefactor(S,Q)*as/(4*pi)*CF*(16*log(1-z)/(1-z))*jacz*(partonLum(x_var, tau/z, 1, Q)/z - partonLum(x_var, tau, 1, Q));
    return 0;
}

static int Dpart2(const int *ndim, const cubareal xx[],
                  const int *ncomp, cubareal ff[], void *userdata) {
#define z_var xx[1]
#define x_var xx[0]
#define dpart2 ff[0]

    double zmax = tau;
    double zmin = 0.;
    double jacz = zmax-zmin;
    double z = zmin + jacz* z_var;
    // Add for full result: 8*log(pow(Q,2)/pow(muf,2))/(1-z)+
    dpart2 = prefactor(S,Q)*as/(4*pi)*CF*(16*log(1-z)/(1-z))*jacz*(-partonLum(x_var, tau, 1, Q));
    return 0;
}

static int Lterms(const int *ndim, const cubareal xx[],
                  const int *ncomp, cubareal ff[], void *userdata) {
#define z_var xx[1]
#define x_var xx[0]
#define lterms ff[0]

    double zmax = 1.;
    double zmin = tau;
    double jacz = zmax-zmin;
    double z = zmin + jacz* z_var;

    lterms = prefactor(S,Q)*as/(4.*pi)*CF*(-8*(1+z)*log(1-z)-4*(1+pow(z,2))/(1-z)*log(z))*jacz*partonLum(x_var, tau/z, 1., Q)/z;
    return 0;
}


static int deltaZ(const int *ndim, const cubareal xx[],
                 const int *ncomp, cubareal ff[], void *userdata) {
#define x_var xx[0]
#define delta ff[0]

    // Add for full result: 6*log(pow(Q,2)/pow(muf,2))
    delta = prefactor(S,Q)*(1+ as/(4*pi)*CF*(8*zeta2 - 16))*partonLum(x_var, tau, 1, Q);
    return 0;
}

static int subLP(const int *ndim, const cubareal xx[],
                 const int *ncomp, cubareal ff[], void *userdata) {
#define z_var xx[1]
#define x_var xx[0]
#define nlp ff[0]
#define nnlp ff[1]
#define n3lp ff[2]

    double zmax = 1.;
    double zmin = tau;
    double jacz = zmax-zmin;
    double z = zmin + jacz* z_var;
    // NOTE: is the scale M in Hamberg et al. the factorziation scale or renormalization scale?
    nlp = prefactor(S,Q)*as/(4*pi)*CF*(8.-16.*log(1-z)-8.*log(pow(Q,2)/pow(muf,2)))*jacz*partonLum(x_var, tau/z, 1, Q)/z;
    nnlp = prefactor(S,Q)*as/(4*pi)*CF*(-4.+8.*log(1-z)+4.*log(pow(Q,2)/pow(muf,2)))*(1-z)*jacz*partonLum(x_var, tau/z, 1, Q)/z;
    n3lp = prefactor(S,Q)*as/(4*pi)*CF*8./3.*pow((1-z),2)*jacz*partonLum(x_var, tau/z, 1, Q)/z;
    return 0;
}

// Derivative method: derivative d(Phi(z)/z)/dz

double partonLumz(double x, double scale, double var){
    double plz = 1./(var)*partonLum(x, tau/(var), 1., scale);
    return plz;
}

double partonLumDer(double x, double z, double scale, double eps) {
    double f  = (partonLumz(x,scale,z+eps)-partonLumz(x,scale,z))/eps;
    return f;
}

static int derNLO(const int *ndim, const cubareal xx[],
                 const int *ncomp, cubareal ff[], void *userdata) {
#define z_var xx[1]
#define x_var xx[0]
#define dernlo ff[0]

    double zmax = 1.;
    double zmin = tau;
    double jacz = zmax-zmin;
    double z = zmin + jacz* z_var;

    dernlo = prefactor(S,Q)*(1 + as/(4*pi)*CF*(8*zeta2 - 16 + 8*pow(log(1-z),2) - 2*(1-z)*(-7 - z + 2*(3+z)*log(1-z))
 - 5 + 4*z + pow(z,2) - 2*z*(2+z)*log(z) + 8*gsl_sf_dilog(1-z)))*jacz*partonLumDer(x_var, z, Q, 0.00001);
    return 0;
}


int main(){

// Load PDFset
    initPDFSet("NNPDF31_nnlo_as_0118", LHGRID);
//    initPDFSet("MSTW2008nnlo90cl", LHGRID);

    usePDFMember(0);

// Set parameters
    Q = 500.;
    S = pow(13000.,2);
    tau = pow(Q,2)/S;
    muf = Q;
    mur = Q;
    as = alphasPDF(Q);
    cout << prefactor(S,Q) << endl;
// Additional parameters for CUBA
    int comp, nregions, neval, fail, verbose;
    verbose = 0;
    int spin = -1;
    double integral[NCOMP], error[NCOMP], prob[NCOMP];
    double params[1];

// Integration

    Vegas(NDIM2, NCOMP, Dpart1  , params,
          NVEC, EPSREL, EPSABS, verbose, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, STATEFILE, &spin, &neval, &fail, integral, error, prob);
    double Dpart1res = integral[0];
    double Dpart1er = error[0];

    Vegas(NDIM2, NCOMP, Dpart2  , params,
          NVEC, EPSREL, EPSABS, verbose, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, STATEFILE, &spin, &neval, &fail, integral, error, prob);
    double Dpart2res = integral[0];
    double Dpart2er = error[0];

    double Dtermsres = Dpart1res + Dpart2res;
    double Dtermser = Dpart1er + Dpart2er;

    Vegas(NDIM2, NCOMP, Lterms  , params,
          NVEC, EPSREL, EPSABS, verbose, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, STATEFILE, &spin, &neval, &fail, integral, error, prob);
    double Ltermsres = integral[0];
    double Ltermser = error[0];

    Vegas(NDIM1, NCOMP, deltaZ  , params,
          NVEC, EPSREL, EPSABS, verbose, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, STATEFILE, &spin, &neval, &fail, integral, error, prob);
    double zdelta = integral[0];
    double deltaer = error[0];

    Vegas(NDIM2, NCOMP, subLP  , params,
          NVEC, EPSREL, EPSABS, verbose, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, STATEFILE, &spin, &neval, &fail, integral, error, prob);
    double NLP = integral[0];
    double NNLP = integral[1];
    double N3LP = integral[2];
    double NLPer = error[0];
    double NNLPer = error[1];
    double N3LPer = integral[2];

    // Derivative approach; full NLO result

    Vegas(NDIM2, NCOMP, derNLO, params,
          NVEC, EPSREL, EPSABS, verbose, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, STATEFILE, &spin, &neval, &fail, integral, error, prob);
    double derApp = integral[0];
    double derApper = error[0];





// Output

    cout << "===========================================================================" << endl;
    cout << endl << "dsigma/dQ^2 for DY @ NLO" << endl <<endl;
    cout << "===========================================================================" << endl;

    cout << endl << "PARAMETERS" << endl ;
    cout << "..........................................................................." << endl;
    cout << "S = " << S << " GeV^2" << endl<< "Q = " << Q << " GeV" << endl << "tau = " << tau << endl << "tau*sigma_0 = "<< prefactor(S,Q) << endl;

    cout << "---------------------------------------------------------------------------" << endl;

    cout << endl << "RESULTS" << endl;
    cout << "..........................................................................." << endl;
    cout << "Logs and distributions: " << Dtermsres+Ltermsres << " pb/GeV^2 (" << Dtermser+Ltermser <<")" << endl <<endl;
    cout << "Delta function terms  : " << zdelta << " pb/GeV^2 (" << deltaer << ")" << endl;
    cout << "..........................................................................." << endl;
    cout << endl << "Total direct integration : " << Dtermsres+Ltermsres+zdelta << " pb/GeV^2 (" << Dtermser+Ltermser+deltaer << ")" << endl;
    cout << endl << "Total derivative approach: " << derApp << " pb/GeV^2 (" << derApper << ")" << endl;

    cout << "---------------------------------------------------------------------------" << endl;

    cout << endl << "POWER EXPANSION" << endl;
    cout << "..........................................................................." << endl;
    cout << "LP: " << Dtermsres  << " pb/GeV^2 (" << Dtermser  << ")"  << endl;
    cout << "LP (+delta): " << Dtermsres + zdelta << " pb/GeV^2 (" << Dtermser + deltaer << ")"  << "  ; fractional: "<<(Dtermsres + zdelta)/(Dtermsres+Ltermsres+zdelta) << endl;
    cout << "NLP        : " << NLP << " pb/GeV^2 (" << NLPer << ")"  << "  ; fractional: "<<(NLP)/(Dtermsres+Ltermsres+zdelta) << endl;
    cout << "NNLP       : " << NNLP << " pb/GeV^2 (" << NNLPer << ")"  << " ; fractional: "<<(NNLP)/(Dtermsres+Ltermsres+zdelta) << endl;
    cout << "N3LP       : " << NNLP << " pb/GeV^2 (" << N3LPer << ")"  << " ; fractional: "<<(N3LP)/(Dtermsres+Ltermsres+zdelta) << endl;
    cout << "---------------------------------------------------------------------------" << endl;

    cout << endl << "CUMULATIVE" << endl;
    cout << "..........................................................................." << endl;
    cout << "LP  : " << Dtermsres + zdelta  << " pb/GeV^2 " << " ; fractional: "<<(Dtermsres + zdelta)/(Dtermsres+Ltermsres+zdelta) << endl;
    cout << "NLP : " << NLP + Dtermsres + zdelta << " pb/GeV^2 " << " ; fractional: "<<(NLP + Dtermsres + zdelta)/(Dtermsres+Ltermsres+zdelta) << endl;
    cout << "NNLP: " << NNLP + NLP + Dtermsres + zdelta << " pb/GeV^2 " << " ; fractional: "<<(NNLP + NLP + Dtermsres + zdelta)/(Dtermsres+Ltermsres+zdelta) << endl;
    cout << "N3LP: " << N3LP + NNLP + NLP + Dtermsres + zdelta << " pb/GeV^2 " << " ; fractional: "<<(N3LP + NNLP + NLP + Dtermsres + zdelta)/(Dtermsres+Ltermsres+zdelta) << endl;

    cout <<endl << "---------------------------------------------------------------------------" << endl;

    return 0;
}
