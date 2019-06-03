#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "deriv_pdf.h"
using namespace std;

//////////////////////////////////////////////////////////
///
/// contains all K factors for drell yan 
/// split up in LP, NLP, LO, NLO
///
//////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
/// prefactor as given in Nucl.Phys. B359 (1991) 343-405, (A.1), alpha should be alphaEM
/// and Q^4 is actually S^2 * Q^2 (see e.g. universality paper https://arxiv.org/pdf/1706.04018.pdf
////////////////////////////////////////////////////////////////////////////////////////////////////
double LO_factor(){
	return pbunits*4.*M_PI*alphaEM*alphaEM/(3.*Q2*S2)*1./CA;
}

///////////////////////////////////////////////
///
/// next is the gg NLO channel (int and nonint)
///
///////////////////////////////////////////////

////////////////////////////////////////////////
/// LP NLO contribution of DY integrated over z'
////////////////////////////////////////////////
double k_zint_NLO_gg_dy_LP(double z){
	return CF*alphas_Q/(4.*M_PI)*(8*pow(log(1.-z),2));
}


/////////////////////////////////////////////////
/// NLP NLO contribution of DY integrated over z'
/////////////////////////////////////////////////
double k_zint_NLO_gg_dy_NLP(double z){
	return CF*alphas_Q/(4.*M_PI)*(-8.*(1.-z)*(2.*log(1.-z)-3.));
}


//////////////////////////////////////////////////
/// NNLP NLO contribution of DY integrated over z'
//////////////////////////////////////////////////
double k_zint_NLO_gg_dy_NNLP(double z){
	return CF*alphas_Q/(M_PI)*(pow(1.-z,2.)*(log(1.-z)-1.));
}


//////////////////////////////////////////////////
/// const NLO contribution of DY integrated over z'
//////////////////////////////////////////////////
double k_zint_NLO_gg_dy_delta(){
	return CF*alphas_Q/(4.*M_PI)*2.*(2.*M_PI*M_PI/3. - 8.);
}     


///////////////////////////////////////////////////
/// total NLO contribution of DY integrated over z'
///////////////////////////////////////////////////
double k_zint_NLO_gg_dy_exact(double z){
	return CF*alphas_Q/(4.*M_PI)*2.*(-4.*gsl_sf_dilog(z)+2.*(pow(z,2)+2.*z-2.*log(z)-3.)*log(1.-z)+1./2.*(1.-z)*(z+9.)+4.*pow(log(1.-z),2.)-z*(z+2.)*log(z)+2.*pow(M_PI,2.)/3.);
}


///////////////////////////////////////////////////////////
/// Nucl.Phys. B359 (1991) 343-405 - B.3 - delta(1-z) piece
///////////////////////////////////////////////////////////
double k_NLO_dy_gg_delta(){
	return 6.*log(Q2/muF2)+8*zeta2-16.;
}

/////////////////////////////////////////////////////////////////
/// Nucl.Phys. B359 (1991) 343-405 - B.3+B.4 - non constant piece
/////////////////////////////////////////////////////////////////
double k_NLO_dy_gg_nonconst(double z){
	double minz = 1.-z; //define 1-z because we use it so much
	return 16.*log(minz)/minz + 8.*log(Q2/muF2)/minz-4.*(1.+z)*log(Q2/muF2)-8.*(1+z)*log(minz)-4.*(1+pow(z,2))/minz*log(z);
}

////////////////////////////////////////////////////////
/// the vegas routine function
/// still think whether this is the good thing to do!
/// as you integrate over plus distribution
/// this is probably wrong this way, need to split it up
////////////////////////////////////////////////////////
double vegas_k_NLO_dy_gg(double *k, size_t dim, void *params){
	(void)(dim);
	lumni_params * lp = (lumni_params *)params;
	return k_NLO_dy_gg_nonconst(lp->z)+k_NLO_dy_gg_delta();
}

///////////////////////////////////////////////
///
/// next is the qg NLO channel (int and nonint)
///
///////////////////////////////////////////////

/////////////////////////////////////////////////
/// NLP NLO contribution of DY integrated over z'
/////////////////////////////////////////////////
double k_zint_NLO_qg_dy_NLP(double z){
	return TF*alphas_Q/(4.*M_PI)*(-2*(-1 + z)*(-2 + log((Q*pow(-1. + z,2.))/muF)));
}

//////////////////////////////////////////////////
/// NNLP NLO contribution of DY integrated over z'
//////////////////////////////////////////////////
double k_zint_NLO_qg_dy_NNLP(double z){
	return TF*alphas_Q/(4.*M_PI)*((-pow(-1. + z,2.))*(-7. + 2.*log((Q*pow(-1. + z,2))/muF)));
}

//////////////////////////////////////////////////////
/// exact NLO qg contribution of DY integrated over z'
//////////////////////////////////////////////////////
double k_zint_NLO_qg_dy_exact(double z){
	return alphas_Q*TF/(36.*M_PI)*((z-1.)*(6.*(-2.*pow(z,2)+z-2.)*log(Q2/muF2)+25.*pow(z,2)+z-2.)+12.*(z*((3.-2.*z)*z-3.)+2.)*log(1.-z)+6.*z*(z*(2.*z-3.)+3.)*log(z));
}


//////////////////////////////////////////////
/// LP contribution for NLO DY
//////////////////////////////////////////////
double vegas_sig_LP_1(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*alphas_Q/(4*M_PI)*CF*16.*log(1.-z)/(1.-z)*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z - (1.-tau)*pdf_sum_qq_charge_weighted(tau+k[1]*(1.-tau),tau));
}
//////////////////////////////////////////////////////////////
/// correction term stemming from not integrating from 0 to 1
//////////////////////////////////////////////////////////////
double vegas_sig_LP_correction(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0]; // from 0 to tau
	double jacobian = 1.-tau; //needed to transform the boundary dependent terms
	double x = tau+k[1]*jacobian;
	return LO_factor()*alphas_Q/(4*M_PI)*CF*16.*log(1.-z)/(1.-z)*jacobian*(-pdf_sum_qq_charge_weighted(x,tau));
	
	//return log(1.-z)/(1.-z)*jacobian*(x*(z-tau) - x*(1.-tau)); // test function
}

////////////////////////////////////////
/// NLP contribution DY nlo
////////////////////////////////////////
double vegas_sig_NLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*alphas_Q/(4*M_PI)*CF*(8.-16.*log(1-z)-8.*log(Q2/muF2))*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}


////////////////////////////////////////
/// NNLP contribution DY nlo
////////////////////////////////////////
double vegas_sig_NNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*alphas_Q/(4*M_PI)*CF*((-4.+8.*log(1-z)+4.*log(Q2/muF2))*(1-z))*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}


////////////////////////////////////////////////////////////////
/// constant delta contribution (contains LO)
/// already integrated over z
/// so make sure here no z integral
////////////////////////////////////////////////////////////////
double vegas_LO(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double jacobian = 1.;//-tau; //needed to transform the boundary dependent terms
	double x = /*tau+*/k[0]*jacobian;
	return LO_factor()*(1.)*(jacobian*pdf_sum_qq_charge_weighted(/*tau+*/k[0]*jacobian,tau));
}



////////////////////////////////////////////////////////////////
/// constant delta contribution (contains alphas*delta(1-z))
/// already integrated over z
/// so make sure here no z integral
////////////////////////////////////////////////////////////////
double vegas_sig_delta(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double jacobian = 1.;//-tau; //needed to transform the boundary dependent terms
	double x = /*tau+*/k[0]*jacobian;
	return LO_factor()*((alphas_Q/(4.*M_PI)*CF*(8.*zeta2 - 16.)))*(jacobian*pdf_sum_qq_charge_weighted(/*tau+*/k[0]*jacobian,tau));
}



double vegas_sum_pdf(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double jacobian = 1.;//-tau; //needed to transform the boundary dependent terms
	double x = /*tau+*/k[0]*jacobian;
	return (jacobian*pdf_sum_qq_charge_weighted(/*tau+*/k[0]*jacobian,tau));
}

///////////////////////////////////////////////////////////////////
/// full contribution DY nlo (without the plus distribution term)
///////////////////////////////////////////////////////////////////
double vegas_sig_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*alphas_Q/(4*M_PI)*CF*(-8.*(1.+z)*log(1.-z)-4.*(1.+pow(z,2))/(1.-z)*log(z))*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}

double vegas_sig_LP_int(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0]; // from tau to 1
	double x = k[1]; // from 0 to 1
	return LO_factor()*alphas_Q/(4*M_PI)*CF*16.*1./2.*pow(log(1.-z),2)*derivative_qq_pdf_jac(x, z, tau); // jacobian needs to be in the derivative
}

double vegas_sig_NLP_int(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0]; // from tau to 1
	double x = k[1]; // from 0 to 1
	return LO_factor()*alphas_Q/(4*M_PI)*CF*(-16.*(1.-z)*(log(1.-z)-3./2))*derivative_qq_pdf_jac(x, z, tau); // jacobian needs to be in the derivative
}

double vegas_sig_NNLP_int(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0]; // from tau to 1
	double x = k[1]; // from 0 to 1
	return LO_factor()*alphas_Q/(4*M_PI)*CF*(4.*pow(1.-z,2)*(log(1.-z)-1.))*derivative_qq_pdf_jac(x, z, tau); // jacobian needs to be in the derivative
}

double vegas_sig_full_int(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0]; // from tau to 1
	double x = k[1]; // from 0 to 1
	return LO_factor()*alphas_Q/(4*M_PI)*CF*2.*(-4.*gsl_sf_dilog(z)+2.*(pow(z,2)+2.*z-2.*log(z)-3.)*log(1.-z)+1./2.*(1.-z)*(z+9.)+4.*pow(log(1.-z),2.)-z*(z+2.)*log(z)+2.*pow(M_PI,2.)/3.)*derivative_qq_pdf_jac(x, z, tau); // jacobian needs to be in the derivative
}

double vegas_qg_full_int(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0]; // from tau to 1
	double x = k[1]; // from 0 to 1
	return LO_factor()*k_zint_NLO_qg_dy_exact(z)*derivative_qg_pdf_jac(x, z, tau); // jacobian needs to be in the derivative
}

double vegas_qg_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*alphas_Q/(4*M_PI)*TF*(2.*(1.+2.*pow(z,2)-2.*z)*(log(pow(1.-z,2))-log(z)+log(Q2/muF2))+1.-7.*pow(z,2)+6.*z)*(jacobian*pdf_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}

double vegas_sig_testfunction1(double *k, size_t dim, void *params){
	(void)(dim);
	lumni_params * lp = (lumni_params *)params;
	double z = lp->z;
	double eps = 0.0001;
	double zpeps = z + eps;
	double zmeps = z - eps;
	double jacp = 1.-tau/zpeps, jacm = 1.-tau/zmeps;
	double xpeps = tau/zpeps+k[0]*jacp, xmeps = tau/zmeps+k[0]*jacm;
	return derivative_qq_pdf_jac(k[0], z, tau);
	//return (jacp*xpeps*(zpeps-tau)-jacm*xmeps*(zmeps-tau))/(2*eps); // jacobian needs to be in the derivative
}
double vegas_sig_testfunction2(double *k, size_t dim, void *params){
	(void)(dim);
	lumni_params * lp = (lumni_params *)params;
	double z = lp->z;
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms [be sure to integrate from 0 to 1]
	double x = tau/z+k[0]*jacobian;
	//return jacobian*x*(z-tau); 
	return jacobian*pdf_sum_qq_charge_weighted(x, tau/z)/z;
}
