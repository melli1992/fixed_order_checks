#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include "parameters.h"
#include "resum_functions.h"

//////////////////////////////////////////////////////////
///
/// contains all the resummation functions up to NNLL
/// note that no final state resummation effects are here
/// so code is only for processes with IS partons
/// some constants that are known to resum (e.g. in DY)
/// are also included
///
//////////////////////////////////////////////////////////
using namespace std;

// LP LL function h0 (or g1) hep-ph/0306211 eqn 39 (note that gammaE is not part there of lambda and the factor 2 (as they have 2*h0 = g1!) or 1905.11771 eqn 6
// (checked with mathematica), seems that the I*0.0 is another side of the cut taken by mathematica and by c++ (maybe this is bad, have to find out)
complex<double> h0(double A1,complex<double>lambda){
	return A1/(2.*M_PI*pow(b0,2))*(2.*lambda+(1.-2.*lambda)*log(1.-2.*lambda));
}
// LP LL (quark)
complex<double> h0q(complex<double>lambda){
	return h0(A1q,lambda);
}
// LP LL (gluon)
complex<double> h0g(complex<double>lambda){
	return h0(A1g,lambda);;
}

/////////////////////////////////////////////////////////////
// NLP LL function
complex<double> h0NLP(double A1, complex<double> N, complex<double>lambda){
	return -A1/(2.*M_PI*b0)*(log(1.-2.*lambda))/N;
}
// NLP LL (quark)
complex<double> h0qNLP(complex<double> N, complex<double>lambda){
	return h0NLP(A1q, N, lambda);
}
// NLP LL (gluon)
complex<double> h0gNLP(complex<double> N, complex<double>lambda){
	return h0NLP(A1g, N, lambda);
}


/////////////////////////////////////////////////////////////
// LP NLL function h1 (or g2) hep-ph/0306211 (note factor 2 and gammaE) eqn 40 or 1905.11771 eqn 61
// (checked with mathematica)
complex<double> h1(double A1,double A2,complex<double>lambda){
	return 1./(2.*M_PI*b0)*(-A2/(M_PI*b0)+A1*log(Q2/muR2))*(2.*lambda+log(1.-2.*lambda))
	+ A1*b1/(2.*M_PI*pow(b0,3))*(2.*lambda+log(1.-2.*lambda)+1./2.*pow(log(1.-2.*lambda),2))
	- A1/(M_PI*b0)*lambda*log(Q2/muF2);
}
// LP NLL (quark)
complex<double> h1q(complex<double>lambda){
	return h1(A1q,A2q,lambda);
}
// LP NLL (gluon)
complex<double> h1g(complex<double>lambda){
	return h1(A1g,A2g,lambda);
}

////////////////////////////////////////////////////////////////////////
// LP NNLL function g3 hep-ph/0306211 eqn 41 (checked with mathematica)
complex<double> h2(double A1,double A2,double A3,complex<double>lambda){
	return 2.*A1/M_PI*zeta2*lambda/(1.-2.*lambda)+A1*pow(b1,2)/(2.*M_PI*pow(b0,4)*(1.-2.*lambda))*(2.*pow(lambda,2)+2.*lambda*log(1.-2.*lambda)+1./2.*pow(log(1.-2.*lambda),2))
	+A1*b2/(2.*M_PI*pow(b0,3))*(2.*lambda+log(1.-2.*lambda)+2.*pow(lambda,2)/(1.-2.*lambda))+A3/(pow(M_PI,3)*pow(b0,2))*pow(lambda,2)/(1.-2.*lambda)
	-A2*b1/(2.*pow(M_PI,2)*pow(b0,3))*(1./(1.-2.*lambda))*(2.*lambda+log(1.-2.*lambda)+2.*pow(lambda,2))-A2/(pow(M_PI,2)*b0)*lambda*log(Q2/muF2)-A1/(2.*M_PI)*lambda*pow(log(Q2/muF2),2)
	+A1/M_PI*lambda*log(Q2/muR2)*log(Q2/muF2)+1./(1.-2.*lambda)*(A1*b1/(2.*M_PI*pow(b0,2))*(2.*lambda+log(1.-2.*lambda))-2*A2/(pow(M_PI,2)*b0)*pow(lambda,2))*log(Q2/muR2)+A1/M_PI*pow(lambda,2)/(1.-2.*lambda)*pow(log(Q2/muR2),2);
}
// LP NNLL (quark)
complex<double> h2q(complex<double>lambda){
	return h2(A1q,A2q,A3q,lambda);
}
// LP NNLL (gluon)
complex<double> h2g(complex<double>lambda){
	return h2(A1g,A2g,A3g,lambda);
}

// wide angle contribution
complex<double> wideangle(double D2,complex<double>lambda){
	return -D2/(pow(M_PI,2)*b0)*lambda/(1.-2.*lambda);
}


//constants (DY case in MSbar scheme - hep-ph/0508284 Eq. 3.13 and 4.6)
double FDY(){
	return -alphas_muR/(M_PI)*CF*(3./2.*zeta2)+pow(alphas_muR/M_PI,2)*(CA*CF*(607./324.-469./144.*zeta2+1./4.*pow(zeta2,2)-187./72*zeta3)+(-41./162.+35./72.*zeta2+17./36.*zeta3)*nF*CF);
}

///////////////////////////////////////////////////////////////////
// matching functions to NLO and NNLO (checked with mathematica)
// need the extra M_gammaE in front

// (di-)higgs specific
complex<double> NLOmatch_higgs(complex<double> N){
	return 1. + 2.*alphas_muR*(A1g*log(N*exp(M_gammaE))*(ISNLP/N+ISLL*log(N*exp(M_gammaE))-ISNLL*log(Q2/muF2)))/(M_PI);
}
complex<double> NNLOmatch_higgs(complex<double> N){
	return 1. + 2.*alphas_muR*(A1g*log(N*exp(M_gammaE))*(ISNLP/N+ISLL*log(N*exp(M_gammaE))-ISNLL*log(Q2/muF2)))/(M_PI)+
	pow(alphas_muR,2)*log(N*exp(M_gammaE))*((6.*pow(A1g,2)*log(N*exp(M_gammaE))*pow(ISNLP/N + ISLL*log(N*exp(M_gammaE)) - ISNLL*log(Q2/muF2),2) + (-3.*D2higgs*ISNNLL+6.*A2g*ISNLL*log(N*exp(M_gammaE)) + 6.*ISNLP*A1g*b0*log(N*exp(M_gammaE))*M_PI/N +
	4.*A1g*b0*ISLL*pow(log(N*exp(M_gammaE)),2)*M_PI + 2.*A1g*b0*ISNNLL*pow(M_PI,3) - 3.*ISNNLL*A1g*b0*M_PI*pow(log(Q2/muF2),2) -
       6.*ISNLL*A1g*b0*log(N*exp(M_gammaE))*M_PI*log(Q2/muR2) - 6.*ISNNLL*log(Q2/muF2)*(A2g - A1g*b0*M_PI*log(Q2/muR2))))/(3*pow(M_PI,2)));
}

// THESE FUNCTIONS NEED TO BE CHECKED!
//https://arxiv.org/pdf/hep-ph/0306211.pdf eqn 44 with gammaE->0 (we don't need this one)
double Cgg1_higgs()
{
	return pow(M_PI,2)+11./2.+6.*zeta2+(33.-2.*nF)/6.*log(muR2/muF2);
}
double Cgg2_higgs()
{
	return 133./12.*pow(M_PI,2)-5.*nF*pow(M_PI,2)/18.+29./20.*pow(M_PI,4)+22.*zeta3-4.*nF*zeta3/3.+3.*pow(M_PI,2)*pow(log(Q2/muF2),2)
		-1./4.*(33.-2.*nF)*pow(M_PI,2)*log(Q2/muR2)+log(Q2/muF2)*(11./2.*pow(M_PI,2)-nF*pow(M_PI,2)/3.-72.*zeta3)
		+ 11399./144.+133./2.*zeta2-9./20.*pow(zeta2,2)-165./4.*zeta3+(19./8.+2./3.*nF)*log(Q2/mt2)+nF*(-1189./144.-5./3.*zeta2+5./6.*zeta3)
		+pow(33.-2.*nF,2)/48.*pow(log(muF2/muR2),2)-18.*zeta2*pow(log(Q2/muF2),2)+(169./4.+171./2.*zeta3-19./6.*nF+(33.-2.*nF)*zeta2)*log(Q2/muF2)
		+(-465./8.+13./3.*nF-3./2.*(33.-2.*nF)*zeta2)*log(Q2/muR2);
}
//https://arxiv.org/pdf/1807.03704.pdf eqn 11
double Cgg1_dihiggs(double Q2)
{
	complex<double> CLO = 3.*mH2/(Q2-mH2+I*mH*GammaH)-1.;
	double sigma1fin_sigma0 = 1./norm(CLO)*(11.*norm(CLO)+4./3.*real(CLO));
	return CA*pow(M_PI,2)*4./3.+sigma1fin_sigma0;
}
//https://arxiv.org/pdf/1807.03704.pdf eqn 11
// see also  https://arxiv.org/pdf/1505.07122.pdf eqn 14-16
// and https://arxiv.org/pdf/1305.5206.pdf eqn 12, 13 for I2, V2
// and https://arxiv.org/pdf/1408.2422.pdf for R2, F2 (full scale dependence, checked that F2 is the same, R2 cannot be checked)
// check whether this is now fully correct
/*double Cgg2_dihiggs(double Q2, double ctheta)
{
	double t = -1./2.*(Q2-2.*mH2-sqrt(Q2*(Q2-4.*mH2))*ctheta);
	double u = -1./2.*(Q2-2.*mH2+sqrt(Q2*(Q2-4.*mH2))*ctheta);
	complex<double> CLO = 3.*mH2/(Q2-mH2+I*mH*GammaH)-1.;
	double sigma1fin_sigma0 = 1./norm(CLO)*(11.*norm(CLO)+4./3.*real(CLO));
	double tplus = -1./2.*(Q2-2.*mH2-Q*sqrt(Q2-4.*mH2));
	double tmin = -1./2.*(Q2-2.*mH2+Q*sqrt(Q2-4.*mH2));
	double Lm = log(muR2/mt2);
	double Ls = log(muR2/Q2);
	double Lu = log(muR2/(-u));
	double Lt = log(muR2/(-t));
	double V2 = 1./pow(3.*Q2*t*u,2)*(pow(mH2,4)*pow(t+u,2)-2.*mH2*mH2*t*u*pow(t+u,2)+pow(t*u,2)*(4.*Q2+pow(t+u,2)));
	double I2 = 4.*M_PI*(1.+2.*pow(mH2,2)/pow(Q2,2)*log((mH2-t)*(mH2-u)/(t*u));
	double F2 = pow(CA,2)*(23827./648.-83./6.*zeta2-253./36.*zeta3+5./8.*zeta4+7./2.*Lm+89./3.*Ls+121./12.*pow(Ls,2))
				+9.*pow(CF,2)+CA*CF*(-145./6.-11./2.*Lm-11.*Ls)+pow(nF,2)*pow(TF,2)*(4./3.*pow(Ls,2)-22./9.*zeta2)
				-5./24.*CA-1./3.*CF-1./3.*nF*TF*CA*(2255./54.+40.*Ls+22.*pow(Ls,2)-217./6.*zeta2+49./3.*zeta3)
				-1./3.*nF*TF*CF*(41.-12.*Lm-24.*zeta3);
	double R2 = -7.*pow(CA,2)+11.*CA*CF-8.*nF*CF*TF+1./3.*CA*(476./9.+11./3.*(4.*Ls+Lt+Lu+4.*mH2/Q2)
			-8.*CF-4./9.*TF*nF*(10./3.+4.*Ls+Lt+Lu))
			-CA/3.*(1.+2.*pow(mH2,2)/pow(Q2,2))*(2.*li2(1.-pow(mH2,2)/(t*u))+4.*li2(mH2/t)+4.*li2(mH2/u)+4.*log(1.-mH2/t)*log(-mH2/t)
					+4.*log(1.-mH2/u)*log(-mH2/u)-8.*zeta2-pow(log(t/u),2));
	double sigma2fin_sigma0 = 1./norm(CLO)*1./(tplus-tmin)*(norm(CLO)*F2+real(CLO)*R2+im(CLO)*I2+V2);
	return pow(CA,2)*(-55./36.*zeta3+607./81.+67./16.*pow(M_PI,2)+91./144.*pow(M_PI,4))+CA*nF*(5.*zeta3/18.-82./81-5./8.*pow(M_PI,2))+pow(b0,2)*11./3.*pow(M_PI,4)+CA*sigmafin1/sigma0*(4./3.*pow(M_PI,2))+sigmafin2/sigma0;
}*/

//DY specific
complex<double> NLOmatch_DY(complex<double> N){
	return 1. + 2.*alphas_muR*(A1q*log(N*exp(M_gammaE))*(ISNLP/N+ISLL*log(N*exp(M_gammaE))-ISNLL*log(Q2/muF2)))/(M_PI);
}
complex<double> NNLOmatch_DY(complex<double> N){
	return 1. + 2.*alphas_muR*(A1q*log(N*exp(M_gammaE))*(ISNLP/N+ISLL*log(N*exp(M_gammaE))-ISNLL*log(Q2/muF2)))/(M_PI)+
	pow(alphas_muR,2)*log(N*exp(M_gammaE))*((6.*pow(A1q,2)*log(N*exp(M_gammaE))*pow(ISNLP/N + ISLL*log(N*exp(M_gammaE)) - ISNLL*log(Q2/muF2),2) + (-3.*D2DY*ISNNLL+6.*A2q*ISNLL*log(N*exp(M_gammaE)) + 6.*ISNLP*A1q*b0*log(N*exp(M_gammaE))*M_PI/N +
	4.*A1q*b0*ISLL*pow(log(N*exp(M_gammaE)),2)*M_PI + 2.*A1q*b0*ISNNLL*pow(M_PI,3) - 3.*ISNNLL*A1q*b0*M_PI*pow(log(Q2/muF2),2) -
       6.*ISNLL*A1q*b0*log(N*exp(M_gammaE))*M_PI*log(Q2/muR2) - 6.*ISNNLL*log(Q2/muF2)*(A2q - A1q*b0*M_PI*log(Q2/muR2))))/(3*pow(M_PI,2)));
}
