#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "clooptools.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "k_factors_higgs.h"
using namespace std;

////////////////////////////////////////////////////////////
///
/// contains all K factors for higgs
/// split up in LO, NLO and power corrections
/// also have two routines: either fitted pdfs or real ones
///
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
/// couplings
////////////////////////////////////////////////////////////
//double lambdaSM = 3.*mH2/mZ2;
//double lambdahhh = 3.*cos(2.*alpha)*sin(beta+alpha)+3.*ep;



/////////////////////////////////////////////////////////////////////////////////////////
/// LO
/////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
/// prefactor as given in 1.309.6594, eqn 6 and 7
////////////////////////////////////////////////////////////////////////////////////////////////////

//Fdelta, has right limit for mt > 1.000 GeV (-> 2./3)

complex<double> Cab(double s, double mQ2){
	return C0(0,0,s,mQ2,mQ2,mQ2);
}

complex<double> Cac(double t, double mH2, double mQ2){
	return C0(0,mH2,t,mQ2,mQ2,mQ2);
}

complex<double> Cad(double u, double mH2, double mQ2){
	return C0(0,mH2,u,mQ2,mQ2,mQ2);
}

complex<double> Cbc(double u, double mH2, double mQ2){
	return C0(0,mH2,u,mQ2,mQ2,mQ2);
}

complex<double> Cbd(double t, double mH2, double mQ2){
	return C0(0,mH2,t,mQ2,mQ2,mQ2);
}

complex<double> Ccd(double s, double mH2, double mQ2){
	return C0(mH2,mH2,s,mQ2,mQ2,mQ2);
}

complex<double> Dabc(double s, double u, double mH2, double mQ2){
	return D0(0,0,mH2,mH2,s,u,mQ2,mQ2,mQ2,mQ2);
}

complex<double> Dacb(double t, double u, double mH2, double mQ2){
	return D0(0,mH2,0,mH2,t,u,mQ2,mQ2,mQ2,mQ2);
}

complex<double> Dbac(double s, double t, double mH2, double mQ2){
	return D0(0,0,mH2,mH2,s,t,mQ2,mQ2,mQ2,mQ2);
}

complex<double> Fdelta(double s, complex<double> cab){
	return (2.*mt2*(2. + 4.*cab*mt2 - cab*s))/s;
}
// https://arxiv.org/pdf/1603.00385.pdf formula 19 with TF = 1./2. includedx
double Ftriangle_approx(double s){
	return 1./2.*(4./3.+7./90.*s/mt2+1./126.*pow(s,2)/pow(mt2,2)+13./12600.*pow(s,3)/pow(mt2,3)+8./51975.*pow(s,4)/pow(mt2,4));
}


complex<double> Fbox(double s, double t, double u, complex<double> cab, complex<double> cac, complex<double> cad, complex<double> cbc, complex<double> cbd,complex<double> ccd, complex<double> dabc,complex<double> dbac, complex<double> dacb){
	return (2.*mt2*(2. + 4.*cab*mt2 + (dabc + dacb + dbac)*mt2*(-2.*mH2 + 8.*mt2 - s) - (2.*(mH2 - 4.*mt2)*(cac*(mH2 - 1.*t) + cbc*(mH2 - 1.*u)))/s + (dacb*(mH2 - 4.*mt2)*(pow(mH2,2) - t*u))/s))/s;
}
// https://arxiv.org/pdf/1603.00385.pdf formula 20 with TF = 1./2. included, but see fig3: poor approximation
// pT2 = (t*u-pow(mH2,2))/s
// also you see an extra factor of 1./2. there in the full cross section (11), which is definately due to identical particles
double Fbox_approx(double s,double pT2){
	return 1./2.*(-4./3.-7./14.*mH2/mt2-(45.*pow(mH2,2)-14.*mH2*s+6.*pow(s,2))/(315.*pow(mt2,2))+(13./630*pT2*s/pow(mt2,2))
		-(780.*pow(mH2,3)-620.*pow(mH2,2)*s+355.*mH2*pow(s,2)-16.*pow(s,3))/(18900.*pow(mt2,3))-pT2*(11.*pow(s,2)-36.*mH2*s)/(1890.*pow(mt2,3))
		-(2400.*pow(mH2,4)-3480.*pow(mH2,3)*s+2955.*pow(mH2,2)*pow(s,2)-704.*mH2*pow(s,3)+120.*pow(s,4))/(207000.*pow(mt2,4))
		+ pT2*s*(114.*pow(mH2,2)-85.*mH2*s+16.*pow(s,2)-8.*pT2*s)/(10395.*pow(mt2,4)));

}


complex<double> Gbox(double s, double t, double u, complex<double> cab, complex<double> cac, complex<double> cad, complex<double> cbc, complex<double> cbd,complex<double> ccd, complex<double> dabc,complex<double> dbac, complex<double> dacb){
	return (1.*mt2*(2.*(-ccd + (dabc + dacb + dbac)*mt2)*(-2.*mH2 + 8.*mt2 + s) - 2.*(cab*s + cac*(-mH2 + t) + cbc*(-mH2 + u)) + (1.*(-1.*dbac*s*t*(pow(mH2,2) + t*(-8.*mt2 + t)) + (-2.*mH2 + 8.*mt2 + s)*(ccd*s*(-4.*mH2 + s) + cab*s*(-2.*mH2 + s) + 2.*cac*(mH2 - 1.*t)*t + 2.*cbc*(mH2 - 1.*u)*u) - 1.*dabc*s*u*(pow(mH2,2) + u*(-8.*mt2 + u))))/(-pow(mH2,2) + t*u)))/s;
}
// https://arxiv.org/pdf/1603.00385.pdf formula 21 with TF = 1./2. included
double Gbox_approx(double s,double pT2){
	return 1./2.*pT2/mt2*(-11./45.-(62.*mH2-5.*s)/(630.*mt2)-(400.*pow(mH2,2)-156.*mH2*s+49.*pow(s,2))/(12600.*pow(mt2,2))+103./18900.*pT2*s/pow(mt2,2)-(980.*pow(mH2,3)-867.*pow(mH2,2)*s+469.*mH2*pow(s,2)-34.*pow(s,3))/(103950.*pow(mt2,3))+pT2*s*(24.*mH2-7.*s)/(4950.*pow(mt2,3)));

}

double dihiggs_LO_factor(double scale2, double ctheta){
  //ctheta from -1 to 1, scale2 from (2mH)^2 to S
	double Cbox = 1;
	double Cdelta = 3.*mH2/(scale2-mH2);

	double s = scale2;
	double t = -1./2.*(scale2-2.*mH2-sqrt(pow((scale2-2.*mH2),2)-4.*mH2*mH2)*ctheta);
	double u = 2.*mH2-t-s;
	complex<double> cab = Cab(s,mt2);
	complex<double> cac = Cac(t,mH2,mt2);
	complex<double> cad = Cad(u,mH2,mt2);
	complex<double> cbc = cad; //Cbc(u,mH2,mt2);
	complex<double> cbd = cac; //Cbd(t,mH2,mt2);
	complex<double> ccd = Ccd(s,mH2,mt2);
	complex<double> dabc = Dabc(s,u,mH2,mt2);
	complex<double> dacb = Dacb(t,u,mH2,mt2);
	complex<double> dbac = Dbac(s,t,mH2,mt2);
	complex<double> fdel = Fdelta(s,cab);
	complex<double> fbox = Fbox(s, t, u, cab, cac, cad, cbc, cbd, ccd, dabc, dbac, dacb);
	complex<double> gbox = Gbox(s, t, u, cab, cac, cad, cbc, cbd, ccd, dabc, dbac, dacb);
	// extra 0.5 for idential final states!
	double result = 0.5*fbunits*0.5*sqrt(((pow(scale2-mH2-mH2,2)-4.*mH2*mH2)))*pow(GF,2)*pow(alphas_muR,2)/(256.*pow(2.*M_PI,3))*(norm(Cdelta*fdel+Cbox*fbox)+norm(Cbox*gbox));
	clearcache();
	return result;
}


double dihiggs_LO_factor_approx(double scale2, double ctheta){
  //ctheta from -1 to 1, scale2 from (2mH)^2 to S
	double Cbox = 1;
	double Cdelta = 3.*mH2/(scale2-mH2);
	complex<double> fdel = 2./3.;
	complex<double> fbox = -2./3.;
	complex<double> gbox = 0.;
	// approximation of Fbox and Gbox, implementation does not work yet
	double result = fbunits*0.5*0.5*scale2*sqrt(1.-4.*mH2/scale2)*pow(GF,2)*pow(alphas_muR,2)/(256.*pow(2.*M_PI,3))*(norm(Cdelta*fdel+Cbox*fbox)+norm(Cbox*gbox));
	clearcache();
	return result;
}
