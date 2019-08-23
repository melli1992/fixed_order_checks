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
complex<double> NLOmatch_higgs(complex<double> N){
	return 1. + 2.*alphas_muR*(A1g*log(N*exp(M_gammaE))*(ISNLP/N+ISLL*log(N*exp(M_gammaE))-ISNLL*log(Q2/muF2)))/(M_PI);
}
complex<double> NNLOmatch_higgs(complex<double> N){
	return 1. + 2.*alphas_muR*(A1g*log(N*exp(M_gammaE))*(ISNLP/N+ISLL*log(N*exp(M_gammaE))-ISNLL*log(Q2/muF2)))/(M_PI)+
	pow(alphas_muR,2)*log(N*exp(M_gammaE))*((6.*pow(A1g,2)*log(N*exp(M_gammaE))*pow(ISNLP/N + ISLL*log(N*exp(M_gammaE)) - ISNLL*log(Q2/muF2),2) + (-3.*D2higgs*ISNNLL+6.*A2g*ISNLL*log(N*exp(M_gammaE)) + 6.*ISNLP*A1g*b0*log(N*exp(M_gammaE))*M_PI/N + 
	4.*A1g*b0*ISLL*pow(log(N*exp(M_gammaE)),2)*M_PI + 2.*A1g*b0*ISNNLL*pow(M_PI,3) - 3.*ISNNLL*A1g*b0*M_PI*pow(log(Q2/muF2),2) - 
       6.*ISNLL*A1g*b0*log(N*exp(M_gammaE))*M_PI*log(Q2/muR2) - 6.*ISNNLL*log(Q2/muF2)*(A2g - A1g*b0*M_PI*log(Q2/muR2))))/(3*pow(M_PI,2)));
}

complex<double> NLOmatch_DY(complex<double> N){
	return 1. + 2.*alphas_muR*(A1q*log(N*exp(M_gammaE))*(ISNLP/N+ISLL*log(N*exp(M_gammaE))-ISNLL*log(Q2/muF2)))/(M_PI);
}
complex<double> NNLOmatch_DY(complex<double> N){
	return 1. + 2.*alphas_muR*(A1q*log(N*exp(M_gammaE))*(ISNLP/N+ISLL*log(N*exp(M_gammaE))-ISNLL*log(Q2/muF2)))/(M_PI)+
	pow(alphas_muR,2)*log(N*exp(M_gammaE))*((6.*pow(A1q,2)*log(N*exp(M_gammaE))*pow(ISNLP/N + ISLL*log(N*exp(M_gammaE)) - ISNLL*log(Q2/muF2),2) + (-3.*D2DY*ISNNLL+6.*A2q*ISNLL*log(N*exp(M_gammaE)) + 6.*ISNLP*A1q*b0*log(N*exp(M_gammaE))*M_PI/N + 
	4.*A1q*b0*ISLL*pow(log(N*exp(M_gammaE)),2)*M_PI + 2.*A1q*b0*ISNNLL*pow(M_PI,3) - 3.*ISNNLL*A1q*b0*M_PI*pow(log(Q2/muF2),2) - 
       6.*ISNLL*A1q*b0*log(N*exp(M_gammaE))*M_PI*log(Q2/muR2) - 6.*ISNNLL*log(Q2/muF2)*(A2q - A1q*b0*M_PI*log(Q2/muR2))))/(3*pow(M_PI,2)));
}
