#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "k_factors_dy.h"
using namespace std;

//////////////////////////////////////////////////////////
///
/// contains all K factors for drell yan 
/// split up in LO, NLO and the power expansion
///
//////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////
/// LO
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
/// prefactor as given in Nucl.Phys. B359 (1991) 343-405, (A.1), alpha should be alphaEM
/// (tau*sigma0!)
////////////////////////////////////////////////////////////////////////////////////////////////////
double DY_LO_factor(){
	return pbunits*4.*M_PI*alphaEM*alphaEM/(3.*Q2*S2)*1./CA;
}

////////////////////////////////////////////////////////////////
/// constant delta contribution
////////////////////////////////////////////////////////////////
double vegas_DY_LO(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double result = DY_LO_factor()*(pdf_sum_qqbar_charge_weighted(k[0],tau));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_LO_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double result = DY_LO_factor()*real(fit_sum_qqbar_charge_weighted(k[0],tau));
	if (isnan(result)){return 0;}
	else{return result;}
}


/////////////////////////////////////////////////////////////////////////////////////////
/// NLO
/////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////
/// qqbar channel
///////////////////////////

/////////// integration routines

double vegas_DY_NLO_qqbar_zdep(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; 
	double x = tau/z+k[1]*jacobian;
	double result = DY_LO_factor()*(DY_NLO_qqbar_reg(z)*jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z+DY_NLO_qqbar_plus(z)*(jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z-(1.-tau)*real(pdf_sum_qqbar_charge_weighted(tau+k[1]*(1.-tau),tau))));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NLO_qqbar_zdep_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double z = k[0];
	double jacobian = 1.-tau/z; 
	double x = tau/z+k[1]*jacobian;
	double result = DY_LO_factor()*(DY_NLO_qqbar_reg(z)*jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z+DY_NLO_qqbar_plus(z)*(jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z-(1.-tau)*real(fit_sum_qqbar_charge_weighted(tau+k[1]*(1.-tau),tau))));
	
	
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NLO_qqbar_zdepcorr(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double z = k[0]; // from 0 to tau
	double jacobian = 1.-tau; 
	double x = tau+k[1]*jacobian;
	double result =  DY_LO_factor()*(DY_NLO_qqbar_plus(z))*jacobian*(-pdf_sum_qqbar_charge_weighted(x,tau));
	
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NLO_qqbar_zdepcorr_fit(double *k, size_t dim, void *params){ /*correction from not integrating from 0 to 1 but tau to 1*/
	(void)(dim);
	(void)(params);

	double z = k[0]; // from 0 to tau
	double jacobian = 1.-tau; 
	double x = tau+k[1]*jacobian;
	
	double result =  DY_LO_factor()*(DY_NLO_qqbar_plus(z))*jacobian*(-real(fit_sum_qqbar_charge_weighted(x,tau)));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NLO_qqbar_delta(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double result = DY_LO_factor()*DY_NLO_qqbar_delta()*real(pdf_sum_qqbar_charge_weighted(k[0],tau));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NLO_qqbar_delta_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	
	double result = DY_LO_factor()*DY_NLO_qqbar_delta()*real(fit_sum_qqbar_charge_weighted(k[0],tau));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NLO_qqbar_LP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NLO_qqbar_plus(z))*(jacobian*pdf_sum_qqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z - (1.-tau)*pdf_sum_qqbar_charge_weighted(tau+k[1]*(1.-tau),tau));
}
double vegas_DY_NLO_qqbar_LP_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NLO_qqbar_plus(z))*(jacobian*real(fit_sum_qqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z))/z - (1.-tau)*real(fit_sum_qqbar_charge_weighted(tau+k[1]*(1.-tau),tau)));
}
double vegas_DY_NLO_qqbar_power(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NLO_qqbar_expansion(z, fp->power))*(jacobian*pdf_sum_qqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_DY_NLO_qqbar_power_fit(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NLO_qqbar_expansion(z, fp->power))*(jacobian*real(fit_sum_qqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z))/z);
	
}


/////////////////////////////////////////////////////////////////
/// Nucl.Phys. B359 (1991) 343-405 - B.4 - non constant piece
/////////////////////////////////////////////////////////////////
double DY_NLO_qqbar_reg(double x){
	return (alphas_muR*CF*(-4*(1 + x)*log(Q2/muF2) - 8*(1 + x)*log(1 - x) + (4*(1 + pow(x,2))*log(x))/(-1 + x)))/(4.*M_PI);
}
///////////////////////////////////////////////////////////
/// Nucl.Phys. B359 (1991) 343-405 - B.3 - plus dist
///////////////////////////////////////////////////////////
double DY_NLO_qqbar_plus(double x){
	return (2*alphas_muR*CF*(2*log(1.-x)/(1.-x) + log(Q2/muF2)/(1.-x)))/M_PI;
}
///////////////////////////////////////////////////////////
/// Nucl.Phys. B359 (1991) 343-405 - B.3 - delta(1-z) piece
///////////////////////////////////////////////////////////
double DY_NLO_qqbar_delta(){
	return (alphas_muR*CF*(-24 + 2*pow(M_PI,2) + 9*log(Q2/muF2)))/(6.*M_PI);
}
// power expansions
double DY_NLO_qqbar_expansion(double x, int power){
	if(power==1){
		return (-2*alphas_muR*CF*(-1 + log(Q2/muF2) + 2*log(1 - x)))/M_PI;
	}
	if(power==2){
		return -((alphas_muR*CF*(-1 + x)*(-1 + log(Q2/muF2) + 2*log(1 - x)))/M_PI);
	}
	if(power==3){
		return (2*alphas_muR*CF*pow(-1 + x,2))/(3.*M_PI);
	}
	if(power==4){
		return -(alphas_muR*CF*pow(-1 + x,3))/(3.*M_PI);
	}
	if(power==5){
		return (7*alphas_muR*CF*pow(-1 + x,4))/(30.*M_PI);
	}
	if(power==6){
		return (-11*alphas_muR*CF*pow(-1 + x,5))/(60.*M_PI);
	}
	if(power==7){
		return (16*alphas_muR*CF*pow(-1 + x,6))/(105.*M_PI);
	}
	if(power==8){
		return (-11*alphas_muR*CF*pow(-1 + x,7))/(84.*M_PI);
	}
	if(power==9){
		return (29*alphas_muR*CF*pow(-1 + x,8))/(252.*M_PI);
	}
	if(power==10){
		return (-37*alphas_muR*CF*pow(-1 + x,9))/(360.*M_PI);
	}
	if(power==-1){
		return DY_NLO_qqbar_reg(x)-(-(alphas_muR*CF*(-7353 + 17287*x - 42482*pow(x,2) + 65038*pow(x,3) - 73142*pow(x,4) + 58570*pow(x,5) - 32570*pow(x,6) + 11974*pow(x,7) - 2621*pow(x,8) + 259*pow(x,9) + 2520*(1 + x)*log(Q2/muF2) + 5040*(1 + x)*log(1 - x)))/(2520.*M_PI));
	}
	if(power==-2){
		return -(alphas_muR*CF*(-7353 + 17287*x - 42482*pow(x,2) + 65038*pow(x,3) - 73142*pow(x,4) + 58570*pow(x,5) - 32570*pow(x,6) + 11974*pow(x,7) - 2621*pow(x,8) + 259*pow(x,9) + 2520*(1 + x)*log(Q2/muF2) + 5040*(1 + x)*log(1 - x)))/(2520.*M_PI);
	}
}


///////////////////////////
/// qg channel
///////////////////////////


/////////// integration routines

double vegas_DY_NLO_qg_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NLO_qg_full(z))*(jacobian*pdf_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_DY_NLO_qg_full_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NLO_qg_full(z))*(jacobian*real(fit_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z))/z);
}
double vegas_DY_NLO_qg_power(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NLO_qg_expansion(z, fp->power))*(jacobian*pdf_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_DY_NLO_qg_power_fit(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NLO_qg_expansion(z, fp->power))*(jacobian*real(fit_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z))/z);
}

///// NLO functions
double DY_NLO_qg_full(double x){
	return (alphas_muR*TF*(1 + 6*x - 7*pow(x,2) + 2*(1 - 2*x + 2*pow(x,2))*(log(Q2/muF2) + 2*log(1 - x) - log(x))))/(4.*M_PI);
}
double DY_NLO_qg_expansion(double x, int power){		
	if(power==1){
		return (alphas_muR*TF*(log(Q2/muF2) + 2*log(1 - x)))/(2.*M_PI);
	}
	if(power==2){
		return (alphas_muR*TF*(-1 + x)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)))/(2.*M_PI);
	}
	if(power==3){
		return (alphas_muR*TF*pow(-1 + x,2)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)))/(2.*M_PI);
	}
	if(power==4){
		return (-2*alphas_muR*TF*pow(-1 + x,3))/(3.*M_PI);
	}
	if(power==5){
		return (7*alphas_muR*TF*pow(-1 + x,4))/(24.*M_PI);
	}
	if(power==6){
		return (-11*alphas_muR*TF*pow(-1 + x,5))/(60.*M_PI);
	}
	if(power==7){
		return (2*alphas_muR*TF*pow(-1 + x,6))/(15.*M_PI);
	}
	if(power==8){
		return (-11*alphas_muR*TF*pow(-1 + x,7))/(105.*M_PI);
	}
	if(power==9){
		return (29*alphas_muR*TF*pow(-1 + x,8))/(336.*M_PI);
	}
	if(power==10){
		return (-37*alphas_muR*TF*pow(-1 + x,9))/(504.*M_PI);
	}
	if(power==-1){
		return DY_NLO_qg_full(x)-((alphas_muR*TF*(7759 - 22518*x + 62208*pow(x,2) - 105840*pow(x,3) + 111720*pow(x,4) - 87024*pow(x,5) + 47628*pow(x,6) - 17328*pow(x,7) + 3765*pow(x,8) - 370*pow(x,9) + 2520*(1 - 2*x + 2*pow(x,2))*log(Q2/muF2) + 5040*(1 - 2*x + 2*pow(x,2))*log(1 - x)))/(5040.*M_PI));
	}
	if(power==-2){
		return (alphas_muR*TF*(7759 - 22518*x + 62208*pow(x,2) - 105840*pow(x,3) + 111720*pow(x,4) - 87024*pow(x,5) + 47628*pow(x,6) - 17328*pow(x,7) + 3765*pow(x,8) - 370*pow(x,9) + 2520*(1 - 2*x + 2*pow(x,2))*log(Q2/muF2) + 5040*(1 - 2*x + 2*pow(x,2))*log(1 - x)))/(5040.*M_PI);
	}
}


