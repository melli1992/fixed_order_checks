#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "k_factors_nnlo_dy.h"
#include "k_factors_dy.h"
#include "polygamma.h"
using namespace std;

//////////////////////////////////////////////////////////
///
/// contains all K factors for drell yan  at NNLO
/// split up in LP, NLP, NNLP
/// base is mathematica code
///
//////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
///
/// qqbar channel 
///
////////////////////////////////////////////////////////////


/////////// integration routines

double vegas_DY_NNLO_qqbar_zdep(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; 
	double x = tau/z+k[1]*jacobian;
	double result = DY_LO_factor()*(DY_NNLO_qqbar_NS(z)*jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z
										+ DY_NNLO_BB_full(z)*jacobian*real(pdf_sum_qqbar_charge_unweighted(x,tau/z))/z
										+2.*DY_NNLO_BC_full(z)*jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z
										+DY_NNLO_qqbar_plus(z)*(jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z-(1.-tau)*real(pdf_sum_qqbar_charge_weighted(tau+k[1]*(1.-tau),tau)))
									);
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NNLO_qqbar_zdep_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double z = k[0];
	double jacobian = 1.-tau/z; 
	double x = tau/z+k[1]*jacobian;
	double result = DY_LO_factor()*(DY_NNLO_qqbar_NS(z)*jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z
										+ DY_NNLO_BB_full(z)*jacobian*real(fit_sum_qqbar_charge_unweighted(x,tau/z))/z
										+2.*DY_NNLO_BC_full(z)*jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z
										+DY_NNLO_qqbar_plus(z)*(jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z-(1.-tau)*real(fit_sum_qqbar_charge_weighted(tau+k[1]*(1.-tau),tau)))
									);
	
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NNLO_qqbar_zdepcorr(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double z = k[0]; // from 0 to tau
	double jacobian = 1.-tau; 
	double x = tau+k[1]*jacobian;
	double result =  DY_LO_factor()*(DY_NNLO_qqbar_plus(z))*jacobian*(-pdf_sum_qqbar_charge_weighted(x,tau));
	
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NNLO_qqbar_zdepcorr_fit(double *k, size_t dim, void *params){ /*correction from not integrating from 0 to 1 but tau to 1*/
	(void)(dim);
	(void)(params);

	double z = k[0]; // from 0 to tau
	double jacobian = 1.-tau; 
	double x = tau+k[1]*jacobian;
	
	double result =  DY_LO_factor()*(DY_NNLO_qqbar_plus(z))*jacobian*(-real(fit_sum_qqbar_charge_weighted(x,tau)));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NNLO_qqbar_delta(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double result = DY_LO_factor()*DY_NNLO_qqbar_delta()*real(pdf_sum_qqbar_charge_weighted(k[0],tau));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NNLO_qqbar_delta_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	
	double result = DY_LO_factor()*DY_NNLO_qqbar_delta()*real(fit_sum_qqbar_charge_weighted(k[0],tau));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NNLO_qqbar_LP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NNLO_qqbar_plus(z))*(jacobian*pdf_sum_qqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z - (1.-tau)*pdf_sum_qqbar_charge_weighted(tau+k[1]*(1.-tau),tau));
}
double vegas_DY_NNLO_qqbar_LP_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NNLO_qqbar_plus(z))*(jacobian*real(fit_sum_qqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z))/z - (1.-tau)*real(fit_sum_qqbar_charge_weighted(tau+k[1]*(1.-tau),tau)));
}
double vegas_DY_NNLO_qqbar_power(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	double result = DY_LO_factor()*(DY_NNLO_qqbar_NS_expansion(z,fp->power)*jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z
										+ DY_NNLO_BB_expansion(z,fp->power)*jacobian*real(pdf_sum_qqbar_charge_unweighted(x,tau/z))/z
										+2.*DY_NNLO_BC_expansion(z,fp->power)*jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z
									);
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NNLO_qqbar_power_fit(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	double result = DY_LO_factor()*(DY_NNLO_qqbar_NS_expansion(z,fp->power)*jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z
										+ DY_NNLO_BB_expansion(z,fp->power)*jacobian*real(fit_sum_qqbar_charge_unweighted(x,tau/z))/z
										+2.*DY_NNLO_BC_expansion(z,fp->power)*jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z
									);
	if (isnan(result)){return 0;}
	else{return result;}
}


double DY_NNLO_qqbar_plus(double x){
	return (pow(alphas_muR,2)*CF*(9*pow(log(Q2/muF2),2)*(11*CA - 2*(18*CF + nF) - 48*CF*log(1 - x)) + 6*log(Q2/muF2)*(-67*CA + 144*CF + 10*nF + 18*CA*zeta2 + 72*CF*zeta2 + 6*(11*CA - 2*(9*CF + nF))*log(1 - x) - 216*CF*pow(log(1 - x),2)) + 2*(202*CA - 28*nF - 198*CA*zeta2 + 36*nF*zeta2 - 189*CA*zeta3 - 864*CF*zeta3 + 6*(CA*(-67 + 18*zeta2) + 2*(5*nF + 36*CF*(2 + zeta2)))*log(1 - x) + 18*(11*CA - 2*nF)*pow(log(1 - x),2) - 432*CF*pow(log(1 - x),3))))/(108.*pow(M_PI,2)*(-1 + x));
}
double DY_NNLO_qqbar_delta(){
	return (pow(alphas_muR,2)*CF*(22995*CF + 3810*nF - 12600*CF*zeta2 - 2240*nF*zeta2 + 288*CF*pow(zeta2,2) - 10800*CF*zeta3 + 1440*nF*zeta3 + CA*(-23025 + 11840*zeta2 - 432*pow(zeta2,2) + 5040*zeta3) + 60*(-34*nF + CA*(193 - 72*zeta3) + 3*CF*(-93 + 24*zeta2 + 176*zeta3))*log(Q2/muF2) - 180*(11*CA - 2*(nF + CF*(9 - 16*zeta2)))*pow(log(Q2/muF2),2)))/(2880.*pow(M_PI,2));
}
double DY_NNLO_qqbar_NS(double x){
	return (pow(alphas_muR,2)*CF*(beta0*log(muR2/muF2)*(-4*(1 + x)*log(Q2/muF2) - 8*(1 + x)*log(1 - x) + (4*(1 + pow(x,2))*log(x))/(-1 + x)) + (2*nF*(94 - 206*x - 18*(1 + x)*pow(log(Q2/muF2),2) + 24*(-1 + 11*x)*log(1 - x) + 36*(2 - 3*x)*log(x) - (12*log(Q2/muF2)*(-1 + 12*x - 11*pow(x,2) + 6*(-1 + pow(x,2))*log(1 - x) - 6*(1 + pow(x,2))*log(x)))/(-1 + x) - (18*(1 + pow(x,2))*(log(x)*(5 - 8*log(1 - x) + 3*log(x)) - Li2(1 - x)))/(-1 + x) + 3*(1 + x)*(4*pow(M_PI,2) - 24*pow(log(1 - x),2) + 3*pow(log(x),2) + 6*Li2(1 - x))))/27. + 2*(-CA/2. + CF)*(94 - 78*x + 8*(-8 + 7*x)*log(1 - x) + 2*(22 - 9*x)*log(x) + (4*log(Q2/muF2)*(8 - 15*x + 7*pow(x,2) + (5 - 2*pow(x,2))*log(x) + (1 + pow(x,2))*pow(log(x),2) + 2*(1 + pow(x,2))*Li2(1 - x)))/(-1 + x) + ((1 + x)*(log(x)*(-168*log(1 - x) + log(x)*(69 + 4*log(x))) + 12*(-13 + 2*log(x))*Li2(1 - x) - 48*Li3(1 - x)))/6. - ((1 + pow(x,2))*(72*log(x) - 72*log(1 - x)*log(x) + 45*pow(log(x),2) - 156*log(1 - x)*pow(log(x),2) + 16*pow(log(x),3) - 12*(3 + 8*log(1 - x) + 6*log(x))*Li2(1 - x) - 216*log(x)*Li2(x) + 96*Li3(1 - x) + 216*Li3(x) - 216*zeta3))/(6.*(-1 + x))) + CF*(-72 + 48*x + 4*(64 + 3*x)*log(1 - x) - 64*(-1 + x)*(zeta2 - pow(log(1 - x),2)) - 16*log(x) + 8*(-4 + 13*x)*log(x) + 16*(7 - 6*x)*log(1 - x)*log(x) - 16*(3 + x)*log(1 - x)*log(x) + 16*(-2 + x)*pow(log(x),2) + 8*(3 + x)*pow(log(x),2) + pow(log(Q2/muF2),2)*(-8*(5 + x) + (16*(1 + pow(x,2))*log(x))/(-1 + x) + 8*(1 + x)*(-4*log(1 - x) + log(x))) + 8*(3 - 2*x)*Li2(1 - x) - 16*(3 + x)*Li2(1 - x) + (4*log(Q2/muF2)*(-30 + 26*x + 4*pow(x,2) - 8*zeta2 + 8*pow(x,2)*zeta2 - 24*(-1 + pow(x,2))*pow(log(1 - x),2) - 2*log(x) + 20*x*log(x) - 6*pow(x,2)*log(x) - 3*pow(log(x),2) - 9*pow(x,2)*pow(log(x),2) + 4*log(1 - x)*(7 - 8*x + pow(x,2) + (5 + 9*pow(x,2))*log(x)) + 4*(-3 + pow(x,2))*Li2(1 - x)))/(-1 + x) - (2*(1 + x)*(192*zeta3 - 96*zeta2*log(1 - x) + 96*pow(log(1 - x),3) + 48*zeta2*log(x) - 48*pow(log(1 - x),2)*log(x) + 24*log(1 - x)*pow(log(x),2) - 7*pow(log(x),3) - 72*log(1 - x)*Li2(1 - x) - 24*log(x)*Li2(x) + 60*Li3(1 - x) + 24*Li3(x) - 24*zeta3))/3. + (4*(1 + pow(x,2))*(-14*log(x) - 16*zeta2*log(x) + 31*pow(log(1 - x),2)*log(x) - 14*log(1 - x)*pow(log(x),2) + 3*pow(log(x),3) - 6*(log(1 - x) - log(x))*Li2(1 - x) + 8*log(x)*Li2(x) + 2*Li3(1 - x) - 8*Li3(x) + 8*zeta3))/(-1 + x)) + CA*(-16.51851851851852 + (2278*x)/27. - (4*(19 + 25*x)*zeta2)/3. + (22*(1 + x)*pow(log(Q2/muF2),2))/3. - (4*(38 + 239*x)*log(1 - x))/9. + (2*(-26 + 57*x)*log(x))/3. + 4*(-3 + x)*log(1 - x)*log(x) - 16*x*log(1 - x)*log(x) + ((23 - 25*x)*pow(log(x),2))/6. + 8*x*pow(log(x),2) + log(Q2/muF2)*((-4*(19 + 124*x))/9. + (1 + x)*(8*zeta2 + (88*log(1 - x))/3. - 6*log(x)) + ((1 + pow(x,2))*((70*log(x))/3. - 8*Li2(1 - x)))/(1 - x)) - 16*x*Li2(1 - x) - (4*(7 + x)*Li2(1 - x))/3. + ((1 + pow(x,2))*(208*log(x) - 48*zeta2*log(x) - 280*log(1 - x)*log(x) + 87*pow(log(x),2) + 12*log(1 - x)*pow(log(x),2) + 8*(-1 + 6*log(1 - x) - 6*log(x))*Li2(1 - x) + 24*log(x)*Li2(x) + 72*Li3(1 - x) - 24*Li3(x) + 24*zeta3))/(6.*(-1 + x)) + (1 + x)*(-28*zeta3 + 16*zeta2*log(1 - x) + (88*pow(log(1 - x),2))/3. + 8*log(1 - x)*Li2(1 - x) - 12*Li3(1 - x) + 8*(log(1 - x)*pow(log(x),2) + 2*log(x)*Li2(x) - 2*Li3(x) + 2*zeta3)))))/(16.*pow(M_PI,2));
}
double DY_NNLO_qqbar_NS_expansion(double x, int power){
	if(power==1){
		return (pow(alphas_muR,2)*CF*(CA*(701 - 504*zeta2 - 378*zeta3) - 4*(nF*(41 - 3*pow(M_PI,2)) + 216*CF*(1 + zeta2 + 2*zeta3)) - 54*beta0*log(muR2/muF2)*(-1 + log(Q2/muF2) + 2*log(1 - x)) + 3*(3*pow(log(Q2/muF2),2)*(11*CA - 2*(6*CF + nF) - 48*CF*log(1 - x)) + log(1 - x)*(-487*CA + 639*CF + 88*nF + 72*(CA + 4*CF)*zeta2 + 6*log(1 - x)*(22*CA + 93*CF - 4*nF - 48*CF*log(1 - x))) + 2*log(Q2/muF2)*(-133*CA + 225*CF + 22*nF + 18*(CA + 4*CF)*zeta2 + 6*log(1 - x)*(11*CA + 24*CF - 2*nF - 36*CF*log(1 - x))))))/(108.*pow(M_PI,2));
	}
	if(power==2){
		return (pow(alphas_muR,2)*CF*(1 - x)*(nF*(503 - 12*pow(M_PI,2)) + 27*CF*(-71 + 96*zeta2 + 64*zeta3) + CA*(-2921 + 558*zeta2 + 378*zeta3) + 54*beta0*log(muR2/muF2)*(-1 + log(Q2/muF2) + 2*log(1 - x)) + 9*pow(log(Q2/muF2),2)*(-11*CA - 36*CF + 2*nF + 48*CF*log(1 - x)) - 6*log(Q2/muF2)*(-235*CA - 135*CF + 34*nF + 18*(CA + 4*CF)*zeta2 + 6*log(1 - x)*(11*CA + 72*CF - 2*nF - 36*CF*log(1 - x))) + 6*log(1 - x)*(470*CA + 522*CF - 68*nF - 36*(CA + 4*CF)*zeta2 + 3*log(1 - x)*(-22*CA - 189*CF + 4*nF + 48*CF*log(1 - x)))))/(216.*pow(M_PI,2));
	}
	if(power==3){
		return (CF*pow(alphas_muR - alphas_muR*x,2)*(1125*CA + 1599*CF - 154*nF - 72*(CA + 8*CF)*zeta2 + 36*beta0*log(muR2/muF2) + 6*(24*CF*pow(log(Q2/muF2),2) + 2*log(Q2/muF2)*(-22*CA - 63*CF + 4*nF + 84*CF*log(1 - x)) + log(1 - x)*(-84*CA - 337*CF + 16*nF + 186*CF*log(1 - x)))))/(216.*pow(M_PI,2));
	}
	if(power==4){
		return (pow(alphas_muR,2)*CF*pow(1 - x,3)*(3937*CA + 2710*CF - 340*nF - 576*(CA + 6*CF)*zeta2 + 288*beta0*log(muR2/muF2) + 24*(36*CF*pow(log(Q2/muF2),2) + 2*log(Q2/muF2)*(-7*(5*CA + CF) + 8*nF + 144*CF*log(1 - x)) + log(1 - x)*(-131*CA - 23*CF + 32*nF + 324*CF*log(1 - x)))))/(3456.*pow(M_PI,2));
	}
	if(power==5){
		return (pow(alphas_muR,2)*CF*pow(-1 + x,4)*(102997*CA + 113974*CF + 3560*nF - 3600*(7*CA + 36*CF)*zeta2 + 12600*beta0*log(muR2/muF2) + 60*(540*CF*pow(log(Q2/muF2),2) + 10*log(Q2/muF2)*(-109*CA + 109*CF + 28*nF + 468*CF*log(1 - x)) + log(1 - x)*(-2036*CA + 3691*CF + 560*nF + 5310*CF*log(1 - x)))))/(216000.*pow(M_PI,2));
	}
	if(power==6){
		return (pow(alphas_muR,2)*CF*pow(1 - x,5)*(102013*CA + 246740*CF + 23560*nF - 3600*(11*CA + 52*CF)*zeta2 + 19800*beta0*log(muR2/muF2) + 60*(780*CF*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-1570*CA + 2662*CF + 440*nF + 7080*CF*log(1 - x)) + log(1 - x)*(-2940*CA + 9081*CF + 880*nF + 8070*CF*log(1 - x)))))/(432000.*pow(M_PI,2));
	}
	if(power==7){
		return (pow(alphas_muR,2)*CF*pow(-1 + x,6)*(4492397*CA + 22734990*CF + 2591540*nF - 1411200*(2*CA + 9*CF)*zeta2 + 1411200*beta0*log(muR2/muF2) + 840*(3780*CF*pow(log(Q2/muF2),2) + 7*log(Q2/muF2)*(-1060*CA + 2321*CF + 320*nF + 5040*CF*log(1 - x)) + log(1 - x)*(-13940*CA + 56661*CF + 4480*nF + 40320*CF*log(1 - x)))))/(3.7044e7*pow(M_PI,2));
	}
	if(power==8){
		return (pow(alphas_muR,2)*CF*pow(1 - x,7)*(14015705*CA + 150820218*CF + 18133640*nF - 1411200*(11*CA + 48*CF)*zeta2 + 7761600*beta0*log(muR2/muF2) + 1680*(10080*CF*pow(log(Q2/muF2),2) + 7*log(1 - x)*(-5137*CA + 25161*CF + 1760*nF + 15660*CF*log(1 - x)) + log(Q2/muF2)*(-19082*CA + 49432*CF + 6160*nF + 95760*CF*log(1 - x)))))/(2.370816e8*pow(M_PI,2));
	}
	if(power==9){
		return (pow(alphas_muR,2)*CF*pow(-1 + x,8)*(74056361*CA + 384*(5365123*CF + 657895*nF) - 6350400*(29*CA + 124*CF)*zeta2 + 92080800*beta0*log(muR2/muF2) + 2520*(78120*CF*pow(log(Q2/muF2),2) + 7*log(1 - x)*(-38284*CA + 216613*CF + 13920*nF + 122940*CF*log(1 - x)) + 6*log(Q2/muF2)*(-23639*CA + 69585*CF + 8120*nF + 125160*CF*log(1 - x)))))/(3.2006016e9*pow(M_PI,2));
	}
	if(power==10){
		return (pow(alphas_muR,2)*CF*pow(1 - x,9)*(6995601*CA + 2934902260*CF + 362705280*nF - 6350400*(37*CA + 156*CF)*zeta2 + 117482400*beta0*log(muR2/muF2) + 2520*(98280*CF*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-170682*CA + 558302*CF + 62160*nF + 952560*CF*log(1 - x)) + log(1 - x)*(-323220*CA + 2061983*CF + 124320*nF + 1092420*CF*log(1 - x)))))/(4.572288e9*pow(M_PI,2));
	}
	if(power==-1){
		return DY_NNLO_qqbar_NS(x)-((pow(alphas_muR,2)*CF*(148176000*(1 - x)*(nF*(503 - 12*pow(M_PI,2)) + 27*CF*(-71 + 96*zeta2 + 64*zeta3) + CA*(-2921 + 558*zeta2 + 378*zeta3) + 54*beta0*log(muR2/muF2)*(-1 + log(Q2/muF2) + 2*log(1 - x)) + 9*pow(log(Q2/muF2),2)*(-11*CA - 36*CF + 2*nF + 48*CF*log(1 - x)) - 6*log(Q2/muF2)*(-235*CA - 135*CF + 34*nF + 18*(CA + 4*CF)*zeta2 + 6*log(1 - x)*(11*CA + 72*CF - 2*nF - 36*CF*log(1 - x))) + 6*log(1 - x)*(470*CA + 522*CF - 68*nF - 36*(CA + 4*CF)*zeta2 + 3*log(1 - x)*(-22*CA - 189*CF + 4*nF + 48*CF*log(1 - x)))) + 148176000*pow(-1 + x,2)*(1125*CA + 1599*CF - 154*nF - 72*(CA + 8*CF)*zeta2 + 36*beta0*log(muR2/muF2) + 6*(24*CF*pow(log(Q2/muF2),2) + 2*log(Q2/muF2)*(-22*CA - 63*CF + 4*nF + 84*CF*log(1 - x)) + log(1 - x)*(-84*CA - 337*CF + 16*nF + 186*CF*log(1 - x)))) + 9261000*pow(1 - x,3)*(3937*CA + 2710*CF - 340*nF - 576*(CA + 6*CF)*zeta2 + 288*beta0*log(muR2/muF2) + 24*(36*CF*pow(log(Q2/muF2),2) + 2*log(Q2/muF2)*(-7*(5*CA + CF) + 8*nF + 144*CF*log(1 - x)) + log(1 - x)*(-131*CA - 23*CF + 32*nF + 324*CF*log(1 - x)))) + 148176*pow(-1 + x,4)*(102997*CA + 113974*CF + 3560*nF - 3600*(7*CA + 36*CF)*zeta2 + 12600*beta0*log(muR2/muF2) + 60*(540*CF*pow(log(Q2/muF2),2) + 10*log(Q2/muF2)*(-109*CA + 109*CF + 28*nF + 468*CF*log(1 - x)) + log(1 - x)*(-2036*CA + 3691*CF + 560*nF + 5310*CF*log(1 - x)))) + 74088*pow(1 - x,5)*(102013*CA + 246740*CF + 23560*nF - 3600*(11*CA + 52*CF)*zeta2 + 19800*beta0*log(muR2/muF2) + 60*(780*CF*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-1570*CA + 2662*CF + 440*nF + 7080*CF*log(1 - x)) + log(1 - x)*(-2940*CA + 9081*CF + 880*nF + 8070*CF*log(1 - x)))) + 864*pow(-1 + x,6)*(4492397*CA + 22734990*CF + 2591540*nF - 1411200*(2*CA + 9*CF)*zeta2 + 1411200*beta0*log(muR2/muF2) + 840*(3780*CF*pow(log(Q2/muF2),2) + 7*log(Q2/muF2)*(-1060*CA + 2321*CF + 320*nF + 5040*CF*log(1 - x)) + log(1 - x)*(-13940*CA + 56661*CF + 4480*nF + 40320*CF*log(1 - x)))) + 135*pow(1 - x,7)*(14015705*CA + 150820218*CF + 18133640*nF - 1411200*(11*CA + 48*CF)*zeta2 + 7761600*beta0*log(muR2/muF2) + 1680*(10080*CF*pow(log(Q2/muF2),2) + 7*log(1 - x)*(-5137*CA + 25161*CF + 1760*nF + 15660*CF*log(1 - x)) + log(Q2/muF2)*(-19082*CA + 49432*CF + 6160*nF + 95760*CF*log(1 - x)))) + 10*pow(-1 + x,8)*(74056361*CA + 384*(5365123*CF + 657895*nF) - 6350400*(29*CA + 124*CF)*zeta2 + 92080800*beta0*log(muR2/muF2) + 2520*(78120*CF*pow(log(Q2/muF2),2) + 7*log(1 - x)*(-38284*CA + 216613*CF + 13920*nF + 122940*CF*log(1 - x)) + 6*log(Q2/muF2)*(-23639*CA + 69585*CF + 8120*nF + 125160*CF*log(1 - x)))) + 7*pow(1 - x,9)*(6995601*CA + 2934902260*CF + 362705280*nF - 6350400*(37*CA + 156*CF)*zeta2 + 117482400*beta0*log(muR2/muF2) + 2520*(98280*CF*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-170682*CA + 558302*CF + 62160*nF + 952560*CF*log(1 - x)) + log(1 - x)*(-323220*CA + 2061983*CF + 124320*nF + 1092420*CF*log(1 - x)))) + 296352000*(CA*(701 - 504*zeta2 - 378*zeta3) - 4*(nF*(41 - 3*pow(M_PI,2)) + 216*CF*(1 + zeta2 + 2*zeta3)) - 54*beta0*log(muR2/muF2)*(-1 + log(Q2/muF2) + 2*log(1 - x)) + 3*(3*pow(log(Q2/muF2),2)*(11*CA - 2*(6*CF + nF) - 48*CF*log(1 - x)) + log(1 - x)*(-487*CA + 639*CF + 88*nF + 72*(CA + 4*CF)*zeta2 + 6*log(1 - x)*(22*CA + 93*CF - 4*nF - 48*CF*log(1 - x))) + 2*log(Q2/muF2)*(-133*CA + 225*CF + 22*nF + 18*(CA + 4*CF)*zeta2 + 6*log(1 - x)*(11*CA + 24*CF - 2*nF - 36*CF*log(1 - x)))))))/(3.2006016e10*pow(M_PI,2)));
	}
	if(power==-2){
		return (pow(alphas_muR,2)*CF*(148176000*(1 - x)*(nF*(503 - 12*pow(M_PI,2)) + 27*CF*(-71 + 96*zeta2 + 64*zeta3) + CA*(-2921 + 558*zeta2 + 378*zeta3) + 54*beta0*log(muR2/muF2)*(-1 + log(Q2/muF2) + 2*log(1 - x)) + 9*pow(log(Q2/muF2),2)*(-11*CA - 36*CF + 2*nF + 48*CF*log(1 - x)) - 6*log(Q2/muF2)*(-235*CA - 135*CF + 34*nF + 18*(CA + 4*CF)*zeta2 + 6*log(1 - x)*(11*CA + 72*CF - 2*nF - 36*CF*log(1 - x))) + 6*log(1 - x)*(470*CA + 522*CF - 68*nF - 36*(CA + 4*CF)*zeta2 + 3*log(1 - x)*(-22*CA - 189*CF + 4*nF + 48*CF*log(1 - x)))) + 148176000*pow(-1 + x,2)*(1125*CA + 1599*CF - 154*nF - 72*(CA + 8*CF)*zeta2 + 36*beta0*log(muR2/muF2) + 6*(24*CF*pow(log(Q2/muF2),2) + 2*log(Q2/muF2)*(-22*CA - 63*CF + 4*nF + 84*CF*log(1 - x)) + log(1 - x)*(-84*CA - 337*CF + 16*nF + 186*CF*log(1 - x)))) + 9261000*pow(1 - x,3)*(3937*CA + 2710*CF - 340*nF - 576*(CA + 6*CF)*zeta2 + 288*beta0*log(muR2/muF2) + 24*(36*CF*pow(log(Q2/muF2),2) + 2*log(Q2/muF2)*(-7*(5*CA + CF) + 8*nF + 144*CF*log(1 - x)) + log(1 - x)*(-131*CA - 23*CF + 32*nF + 324*CF*log(1 - x)))) + 148176*pow(-1 + x,4)*(102997*CA + 113974*CF + 3560*nF - 3600*(7*CA + 36*CF)*zeta2 + 12600*beta0*log(muR2/muF2) + 60*(540*CF*pow(log(Q2/muF2),2) + 10*log(Q2/muF2)*(-109*CA + 109*CF + 28*nF + 468*CF*log(1 - x)) + log(1 - x)*(-2036*CA + 3691*CF + 560*nF + 5310*CF*log(1 - x)))) + 74088*pow(1 - x,5)*(102013*CA + 246740*CF + 23560*nF - 3600*(11*CA + 52*CF)*zeta2 + 19800*beta0*log(muR2/muF2) + 60*(780*CF*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-1570*CA + 2662*CF + 440*nF + 7080*CF*log(1 - x)) + log(1 - x)*(-2940*CA + 9081*CF + 880*nF + 8070*CF*log(1 - x)))) + 864*pow(-1 + x,6)*(4492397*CA + 22734990*CF + 2591540*nF - 1411200*(2*CA + 9*CF)*zeta2 + 1411200*beta0*log(muR2/muF2) + 840*(3780*CF*pow(log(Q2/muF2),2) + 7*log(Q2/muF2)*(-1060*CA + 2321*CF + 320*nF + 5040*CF*log(1 - x)) + log(1 - x)*(-13940*CA + 56661*CF + 4480*nF + 40320*CF*log(1 - x)))) + 135*pow(1 - x,7)*(14015705*CA + 150820218*CF + 18133640*nF - 1411200*(11*CA + 48*CF)*zeta2 + 7761600*beta0*log(muR2/muF2) + 1680*(10080*CF*pow(log(Q2/muF2),2) + 7*log(1 - x)*(-5137*CA + 25161*CF + 1760*nF + 15660*CF*log(1 - x)) + log(Q2/muF2)*(-19082*CA + 49432*CF + 6160*nF + 95760*CF*log(1 - x)))) + 10*pow(-1 + x,8)*(74056361*CA + 384*(5365123*CF + 657895*nF) - 6350400*(29*CA + 124*CF)*zeta2 + 92080800*beta0*log(muR2/muF2) + 2520*(78120*CF*pow(log(Q2/muF2),2) + 7*log(1 - x)*(-38284*CA + 216613*CF + 13920*nF + 122940*CF*log(1 - x)) + 6*log(Q2/muF2)*(-23639*CA + 69585*CF + 8120*nF + 125160*CF*log(1 - x)))) + 7*pow(1 - x,9)*(6995601*CA + 2934902260*CF + 362705280*nF - 6350400*(37*CA + 156*CF)*zeta2 + 117482400*beta0*log(muR2/muF2) + 2520*(98280*CF*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-170682*CA + 558302*CF + 62160*nF + 952560*CF*log(1 - x)) + log(1 - x)*(-323220*CA + 2061983*CF + 124320*nF + 1092420*CF*log(1 - x)))) + 296352000*(CA*(701 - 504*zeta2 - 378*zeta3) - 4*(nF*(41 - 3*pow(M_PI,2)) + 216*CF*(1 + zeta2 + 2*zeta3)) - 54*beta0*log(muR2/muF2)*(-1 + log(Q2/muF2) + 2*log(1 - x)) + 3*(3*pow(log(Q2/muF2),2)*(11*CA - 2*(6*CF + nF) - 48*CF*log(1 - x)) + log(1 - x)*(-487*CA + 639*CF + 88*nF + 72*(CA + 4*CF)*zeta2 + 6*log(1 - x)*(22*CA + 93*CF - 4*nF - 48*CF*log(1 - x))) + 2*log(Q2/muF2)*(-133*CA + 225*CF + 22*nF + 18*(CA + 4*CF)*zeta2 + 6*log(1 - x)*(11*CA + 24*CF - 2*nF - 36*CF*log(1 - x)))))))/(3.2006016e10*pow(M_PI,2));
	}
}

double DY_NNLO_BB_full(double x){
	return -(pow(alphas_muR,2)*CF*((1 + x)*(15*(-1 + x) + pow(M_PI,2)*(1 + x)) - 3*pow(1 + x,2)*pow(log(x),2) + 3*log(x)*(-3 - 4*x - 3*pow(x,2) + 4*pow(1 + x,2)*log(1 + x)) + 12*pow(1 + x,2)*Li2(-x)))/(18.*pow(M_PI,2));
}

double DY_NNLO_BB_expansion(double x, int power){
	if(power==1){
		return 0;
	}
	if(power==2){
		return 0;
	}
	if(power==3){
		return 0;
	}
	if(power==4){
		return 0;
	}
	if(power==5){
		return 0;
	}
	if(power==6){
		return -(pow(alphas_muR,2)*CF*pow(-1 + x,5))/(60.*pow(M_PI,2));
	}
	if(power==7){
		return (pow(alphas_muR,2)*CF*pow(-1 + x,6))/(40.*pow(M_PI,2));
	}
	if(power==8){
		return (-37*pow(alphas_muR,2)*CF*pow(-1 + x,7))/(1260.*pow(M_PI,2));
	}
	if(power==9){
		return (2*pow(alphas_muR,2)*CF*pow(-1 + x,8))/(63.*pow(M_PI,2));
	}
	if(power==10){
		return (-107*pow(alphas_muR,2)*CF*pow(-1 + x,9))/(3240.*pow(M_PI,2));
	}
	if(power==-1){
		return DY_NNLO_BB_full(x)-((pow(alphas_muR,2)*(CA - 2*CF)*CF*(-72*pow(1 - x,7)*(2003715*pow(M_PI,2) + 4596480*log(Q2/muF2) + 9192960*log(1 - x) + 44100*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) - 2*(8859856 + 6011145*zeta2 + 66150*zeta3 - 66150*zeta3)) - 8467200*pow(1 - x,3)*(pow(M_PI,2)*(31 + log(64)) + 6*pow(M_PI,2)*log(Q2/muF2) - 36*zeta2*log((2*Q2)/muF2) + 12*(pow(M_PI,2) - 6*zeta2)*log(1 - x) - 6*(2 + 31*zeta2 + 3*zeta3 - 3*zeta3)) - 203212800*(12*zeta2 - 6*zeta3 + pow(M_PI,2)*(-2 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 6*zeta3) - 50803200*pow(-1 + x,2)*(-30*zeta2 - 6*zeta3 + pow(M_PI,2)*(5 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 6*zeta3) - 6350400*pow(-1 + x,4)*(-64 - 222*zeta2 - 12*zeta3 + pow(M_PI,2)*(37 + log(16)) + 4*pow(M_PI,2)*log(Q2/muF2) - 24*zeta2*log((2*Q2)/muF2) + 8*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 12*zeta3) + 101606400*(1 - x)*(-6 + 36*zeta2 - 30*zeta3 + pow(M_PI,2)*(-6 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 30*zeta3) - 42336*pow(1 - x,5)*(-22992 + 4745*pow(M_PI,2) - 28470*zeta2 - 900*zeta3 + 2880*log(Q2/muF2) + 5760*log(1 - x) + 300*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 900*zeta3) - 7056*pow(-1 + x,6)*(-178144 + 24075*pow(M_PI,2) - 144450*zeta2 - 2700*zeta3 + 34560*log(Q2/muF2) + 69120*log(1 - x) + 900*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 2700*zeta3) - 90*pow(-1 + x,8)*(-12677248 + 1376361*pow(M_PI,2) - 8258166*zeta2 - 52920*zeta3 + 4257792*log(Q2/muF2) + 8515584*log(1 - x) + 17640*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 52920*zeta3) - pow(1 - x,9)*(-940359296 + 107784495*pow(M_PI,2) - 646706970*zeta2 - 2381400*zeta3 + 409489920*log(Q2/muF2) + 818979840*log(1 - x) + 793800*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 2381400*zeta3)))/(9.7542144e9*pow(M_PI,2)));
	}
	if(power==-2){
		return (pow(alphas_muR,2)*(CA - 2*CF)*CF*(-72*pow(1 - x,7)*(2003715*pow(M_PI,2) + 4596480*log(Q2/muF2) + 9192960*log(1 - x) + 44100*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) - 2*(8859856 + 6011145*zeta2 + 66150*zeta3 - 66150*zeta3)) - 8467200*pow(1 - x,3)*(pow(M_PI,2)*(31 + log(64)) + 6*pow(M_PI,2)*log(Q2/muF2) - 36*zeta2*log((2*Q2)/muF2) + 12*(pow(M_PI,2) - 6*zeta2)*log(1 - x) - 6*(2 + 31*zeta2 + 3*zeta3 - 3*zeta3)) - 203212800*(12*zeta2 - 6*zeta3 + pow(M_PI,2)*(-2 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 6*zeta3) - 50803200*pow(-1 + x,2)*(-30*zeta2 - 6*zeta3 + pow(M_PI,2)*(5 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 6*zeta3) - 6350400*pow(-1 + x,4)*(-64 - 222*zeta2 - 12*zeta3 + pow(M_PI,2)*(37 + log(16)) + 4*pow(M_PI,2)*log(Q2/muF2) - 24*zeta2*log((2*Q2)/muF2) + 8*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 12*zeta3) + 101606400*(1 - x)*(-6 + 36*zeta2 - 30*zeta3 + pow(M_PI,2)*(-6 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 30*zeta3) - 42336*pow(1 - x,5)*(-22992 + 4745*pow(M_PI,2) - 28470*zeta2 - 900*zeta3 + 2880*log(Q2/muF2) + 5760*log(1 - x) + 300*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 900*zeta3) - 7056*pow(-1 + x,6)*(-178144 + 24075*pow(M_PI,2) - 144450*zeta2 - 2700*zeta3 + 34560*log(Q2/muF2) + 69120*log(1 - x) + 900*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 2700*zeta3) - 90*pow(-1 + x,8)*(-12677248 + 1376361*pow(M_PI,2) - 8258166*zeta2 - 52920*zeta3 + 4257792*log(Q2/muF2) + 8515584*log(1 - x) + 17640*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 52920*zeta3) - pow(1 - x,9)*(-940359296 + 107784495*pow(M_PI,2) - 646706970*zeta2 - 2381400*zeta3 + 409489920*log(Q2/muF2) + 818979840*log(1 - x) + 793800*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 2381400*zeta3)))/(9.7542144e9*pow(M_PI,2));
	}
}

double DY_NNLO_BC_full(double x){
	return (pow(alphas_muR,2)*CF*(-CA/2. + CF)*(81 - 42*x - 39*pow(x,2) + 6*(9 + 11*x)*log(x) - 3*(-6 + 8*x + 15*pow(x,2))*pow(log(x),2) + 2*(1 + 4*x + pow(x,2))*pow(log(x),3) - 54*(-1 + pow(x,2))*Li2(1 - x) + 24*(1 + 3*x + pow(x,2))*(log(x)*Li2(1 - x) + 2*S12(1 - x)) + pow(1 + x,2)*(3*pow(M_PI,2) + 2*pow(M_PI,2)*log(x) - 6*pow(M_PI,2)*log(1 + x) + 36*log(x)*log(1 + x) + 30*pow(log(x),2)*log(1 + x) - 36*log(x)*pow(log(1 + x),2) + 36*(1 + log(x) - 2*log(1 + x))*Li2(-x) - 12*Li3(-x) - 72*S12(-x))))/(24.*pow(M_PI,2));
}

double DY_NNLO_BC_expansion(double x, int power){
	if(power==1){
		return 0;
	}
	if(power==2){
		return 0;
	}
	if(power==3){
		return 0;
	}
	if(power==4){
		return -(pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,3))/(48.*pow(M_PI,2));
	}
	if(power==5){
		return (25*pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,4))/(576.*pow(M_PI,2));
	}
	if(power==6){
		return (-137*pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,5))/(2880.*pow(M_PI,2));
	}
	if(power==7){
		return (4159*pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,6))/(86400.*pow(M_PI,2));
	}
	if(power==8){
		return (-2059*pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,7))/(43200.*pow(M_PI,2));
	}
	if(power==9){
		return (198077*pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,8))/(4.2336e6*pow(M_PI,2));
	}
	if(power==10){
		return (-193621*pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,9))/(4.2336e6*pow(M_PI,2));
	}
	if(power==-1){
		return DY_NNLO_BC_full(x)-(-(pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,3)*(1270611 - 4157142*x + 6908540*pow(x,2) - 6864109*pow(x,3) + 4096482*pow(x,4) - 1359803*pow(x,5) + 193621*pow(x,6)))/(4.2336e6*pow(M_PI,2)));
	}
	if(power==-2){
		return -(pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,3)*(1270611 - 4157142*x + 6908540*pow(x,2) - 6864109*pow(x,3) + 4096482*pow(x,4) - 1359803*pow(x,5) + 193621*pow(x,6)))/(4.2336e6*pow(M_PI,2));
	}
}


////////////////////////////////////////////////////////////
///
/// qq + qbarqbar channel 
///
////////////////////////////////////////////////////////////

double vegas_DY_NNLO_qq_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	double result = DY_LO_factor()*jacobian/z*(DY_NNLO_CC_full(z)*real(pdf_sum_qq_charge_weighted_double(tau/z+k[1]*jacobian,tau/z))
												+ DY_NNLO_CD_full(z)*real(pdf_sum_qq_charge_weighted_double_vivj(tau/z+k[1]*jacobian,tau/z))
												+ DY_NNLO_CE_full(z)*real(pdf_sum_qq_charge_weighted_single(tau/z+k[1]*jacobian,tau/z))
												+ DY_NNLO_CF_full(z)*real(pdf_sum_qq_charge_weighted_single_vivi(tau/z+k[1]*jacobian,tau/z))  
												);
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NNLO_qq_full_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	double result = DY_LO_factor()*jacobian/z*(DY_NNLO_CC_full(z)*real(fit_sum_qq_charge_weighted_double(tau/z+k[1]*jacobian,tau/z))
												+ DY_NNLO_CD_full(z)*real(fit_sum_qq_charge_weighted_double_vivj(tau/z+k[1]*jacobian,tau/z))
												+ DY_NNLO_CE_full(z)*real(fit_sum_qq_charge_weighted_single(tau/z+k[1]*jacobian,tau/z))
												+ DY_NNLO_CF_full(z)*real(fit_sum_qq_charge_weighted_single_vivi(tau/z+k[1]*jacobian,tau/z))  
												);
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NNLO_qq_power(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	double result = DY_LO_factor()*jacobian/z*(DY_NNLO_CC_expansion(z, fp->power)*real(pdf_sum_qq_charge_weighted_double(tau/z+k[1]*jacobian,tau/z))
												+ DY_NNLO_CD_expansion(z, fp->power)*real(pdf_sum_qq_charge_weighted_double_vivj(tau/z+k[1]*jacobian,tau/z))
												+ DY_NNLO_CE_expansion(z, fp->power)*real(pdf_sum_qq_charge_weighted_single(tau/z+k[1]*jacobian,tau/z))
												+ DY_NNLO_CF_expansion(z, fp->power)*real(pdf_sum_qq_charge_weighted_single_vivi(tau/z+k[1]*jacobian,tau/z))  
												);
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_NNLO_qq_power_fit(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	double result = DY_LO_factor()*jacobian/z*(DY_NNLO_CC_expansion(z, fp->power)*real(fit_sum_qq_charge_weighted_double(tau/z+k[1]*jacobian,tau/z))
												+ DY_NNLO_CD_expansion(z, fp->power)*real(fit_sum_qq_charge_weighted_double_vivj(tau/z+k[1]*jacobian,tau/z))
												+ DY_NNLO_CE_expansion(z, fp->power)*real(fit_sum_qq_charge_weighted_single(tau/z+k[1]*jacobian,tau/z))
												+ DY_NNLO_CF_expansion(z, fp->power)*real(fit_sum_qq_charge_weighted_single_vivi(tau/z+k[1]*jacobian,tau/z))  
												);
	if (isnan(result)){return 0;}
	else{return result;}
}

//NNLO functions
double DY_NNLO_CC_full(double x){
	return (pow(alphas_muR,2)*CF*TF*(1779 + 116/x - 2598*x + 703*pow(x,2) - 24*(39 - 22/x - 39*x + 22*pow(x,2))*log(1 - x) + (72*(-4 - 3*x + 3*pow(x,2) + 4*pow(x,3))*(zeta2 - pow(log(1 - x),2)))/x + 6*(345 - 48*x + 20*pow(x,2))*log(x) + 216*(3 + 6*x + 4*pow(x,2))*log(1 - x)*log(x) - 45*(3 + 15*x + 8*pow(x,2))*pow(log(x),2) + 18*pow(log(Q2/muF2),2)*(3 + 4/x - 3*x - 4*pow(x,2) + 6*(1 + x)*log(x)) + 36*(39 + 16/x + 15*x + 8*pow(x,2))*Li2(1 - x) + 12*log(Q2/muF2)*(-39 + 22/x + 39*x - 22*pow(x,2) + 6*(3 + 4/x - 3*x - 4*pow(x,2))*log(1 - x) + 9*(3 + 6*x + 4*pow(x,2))*log(x) + 18*(1 + x)*((2*log(1 - x) - log(x))*log(x) + 2*Li2(1 - x))) - 18*(1 + x)*(4*pow(M_PI,2)*log(x) - 24*pow(log(1 - x),2)*log(x) - 12*log(1 - x)*pow(log(x),2) - 9*pow(log(x),3) - 12*(4*log(1 - x) + log(x))*Li2(1 - x) - 72*log(x)*Li2(x) + 48*Li3(1 - x) + 72*Li3(x) - 72*zeta3)))/(432.*pow(M_PI,2));
}

double DY_NNLO_CC_expansion(double x, int power){
	if(power==1){
		return 0;
	}
	if(power==2){
		return (pow(alphas_muR,2)*CF*TF*(1 - x)*(27 + 8*pow(M_PI,2) - 60*zeta2 + 3*pow(log(Q2/muF2),2) + 12*log(Q2/muF2)*(-1 + log(1 - x)) + 12*(-2 + log(1 - x))*log(1 - x)))/(24.*pow(M_PI,2));
	}
	if(power==3){
		return (3*pow(alphas_muR,2)*CF*TF*pow(-1 + x,2)*(-2 + log(Q2/muF2) + 2*log(1 - x)))/(8.*pow(M_PI,2));
	}
	if(power==4){
		return (pow(alphas_muR,2)*CF*TF*pow(1 - x,3)*(9*pow(log(Q2/muF2),2) + 18*log(Q2/muF2)*(-1 + 2*log(1 - x)) + 2*(27 + pow(M_PI,2) - 24*zeta2 + 18*(-1 + log(1 - x))*log(1 - x))))/(72.*pow(M_PI,2));
	}
	if(power==5){
		return (pow(alphas_muR,2)*CF*TF*pow(-1 + x,4)*(75 + 8*pow(M_PI,2) - 192*zeta2 + 36*pow(log(Q2/muF2),2) + 36*log(1 - x)*(1 + 4*log(1 - x)) + 18*log(Q2/muF2)*(1 + 8*log(1 - x))))/(288.*pow(M_PI,2));
	}
	if(power==6){
		return (pow(alphas_muR,2)*CF*TF*pow(1 - x,5)*(24137 + 1350*pow(M_PI,2) - 36000*zeta2 + 6975*pow(log(Q2/muF2),2) + 30*log(1 - x)*(523 + 930*log(1 - x)) + 15*log(Q2/muF2)*(523 + 1860*log(1 - x))))/(54000.*pow(M_PI,2));
	}
	if(power==7){
		return (pow(alphas_muR,2)*CF*TF*pow(-1 + x,6)*(5509 + 200*pow(M_PI,2) - 6000*zeta2 + 1200*pow(log(Q2/muF2),2) + 3430*log(1 - x) + 4800*pow(log(1 - x),2) + 5*log(Q2/muF2)*(343 + 960*log(1 - x))))/(9000.*pow(M_PI,2));
	}
	if(power==8){
		return (pow(alphas_muR,2)*CF*TF*pow(1 - x,7)*(3663971 + 98000*pow(M_PI,2) - 3292800*zeta2 + 676200*pow(log(Q2/muF2),2) + 560*log(1 - x)*(3947 + 4830*log(1 - x)) + 280*log(Q2/muF2)*(3947 + 9660*log(1 - x))))/(4.9392e6*pow(M_PI,2));
	}
	if(power==9){
		return (pow(alphas_muR,2)*CF*TF*pow(-1 + x,8)*(836373 + 17640*pow(M_PI,2) - 658560*zeta2 + 138180*pow(log(Q2/muF2),2) + 56*log(1 - x)*(8849 + 9870*log(1 - x)) + 28*log(Q2/muF2)*(8849 + 19740*log(1 - x))))/(987840.*pow(M_PI,2));
	}
	if(power==10){
		return (pow(alphas_muR,2)*CF*TF*pow(1 - x,9)*(21358847 + 370440*pow(M_PI,2) - 15240960*zeta2 + 3254580*pow(log(Q2/muF2),2) + 252*log(1 - x)*(49687 + 51660*log(1 - x)) + 126*log(Q2/muF2)*(49687 + 103320*log(1 - x))))/(2.286144e7*pow(M_PI,2));
	}
	if(power==-1){
		return DY_NNLO_CC_full(x)-((pow(alphas_muR,2)*CF*TF*(1 - x)*(-1500282000*(-1 + x)*(-2 + log(Q2/muF2) + 2*log(1 - x)) + 166698000*(27 + 8*pow(M_PI,2) - 60*zeta2 + 3*pow(log(Q2/muF2),2) + 12*log(Q2/muF2)*(-1 + log(1 - x)) + 12*(-2 + log(1 - x))*log(1 - x)) + 13891500*pow(1 - x,3)*(75 + 8*pow(M_PI,2) - 192*zeta2 + 36*pow(log(Q2/muF2),2) + 36*log(1 - x)*(1 + 4*log(1 - x)) + 18*log(Q2/muF2)*(1 + 8*log(1 - x))) + 444528*pow(1 - x,5)*(5509 + 200*pow(M_PI,2) - 6000*zeta2 + 1200*pow(log(Q2/muF2),2) + 10*log(1 - x)*(343 + 480*log(1 - x)) + 5*log(Q2/muF2)*(343 + 960*log(1 - x))) + 74088*pow(-1 + x,4)*(24137 + 1350*pow(M_PI,2) - 36000*zeta2 + 6975*pow(log(Q2/muF2),2) + 30*log(1 - x)*(523 + 930*log(1 - x)) + 15*log(Q2/muF2)*(523 + 1860*log(1 - x))) + 810*pow(-1 + x,6)*(3663971 + 98000*pow(M_PI,2) - 3292800*zeta2 + 676200*pow(log(Q2/muF2),2) + 560*log(1 - x)*(3947 + 4830*log(1 - x)) + 280*log(Q2/muF2)*(3947 + 9660*log(1 - x))) + 4050*pow(1 - x,7)*(836373 + 17640*pow(M_PI,2) - 658560*zeta2 + 138180*pow(log(Q2/muF2),2) + 56*log(1 - x)*(8849 + 9870*log(1 - x)) + 28*log(Q2/muF2)*(8849 + 19740*log(1 - x))) + 175*pow(-1 + x,8)*(21358847 + 370440*pow(M_PI,2) - 15240960*zeta2 + 3254580*pow(log(Q2/muF2),2) + 252*log(1 - x)*(49687 + 51660*log(1 - x)) + 126*log(Q2/muF2)*(49687 + 103320*log(1 - x))) + 55566000*pow(-1 + x,2)*(9*pow(log(Q2/muF2),2) + 18*log(Q2/muF2)*(-1 + 2*log(1 - x)) + 2*(27 + pow(M_PI,2) - 24*zeta2 + 18*(-1 + log(1 - x))*log(1 - x)))))/(4.000752e9*pow(M_PI,2)));
	}
	if(power==-2){
		return (pow(alphas_muR,2)*CF*TF*(1 - x)*(-1500282000*(-1 + x)*(-2 + log(Q2/muF2) + 2*log(1 - x)) + 166698000*(27 + 8*pow(M_PI,2) - 60*zeta2 + 3*pow(log(Q2/muF2),2) + 12*log(Q2/muF2)*(-1 + log(1 - x)) + 12*(-2 + log(1 - x))*log(1 - x)) + 13891500*pow(1 - x,3)*(75 + 8*pow(M_PI,2) - 192*zeta2 + 36*pow(log(Q2/muF2),2) + 36*log(1 - x)*(1 + 4*log(1 - x)) + 18*log(Q2/muF2)*(1 + 8*log(1 - x))) + 444528*pow(1 - x,5)*(5509 + 200*pow(M_PI,2) - 6000*zeta2 + 1200*pow(log(Q2/muF2),2) + 10*log(1 - x)*(343 + 480*log(1 - x)) + 5*log(Q2/muF2)*(343 + 960*log(1 - x))) + 74088*pow(-1 + x,4)*(24137 + 1350*pow(M_PI,2) - 36000*zeta2 + 6975*pow(log(Q2/muF2),2) + 30*log(1 - x)*(523 + 930*log(1 - x)) + 15*log(Q2/muF2)*(523 + 1860*log(1 - x))) + 810*pow(-1 + x,6)*(3663971 + 98000*pow(M_PI,2) - 3292800*zeta2 + 676200*pow(log(Q2/muF2),2) + 560*log(1 - x)*(3947 + 4830*log(1 - x)) + 280*log(Q2/muF2)*(3947 + 9660*log(1 - x))) + 4050*pow(1 - x,7)*(836373 + 17640*pow(M_PI,2) - 658560*zeta2 + 138180*pow(log(Q2/muF2),2) + 56*log(1 - x)*(8849 + 9870*log(1 - x)) + 28*log(Q2/muF2)*(8849 + 19740*log(1 - x))) + 175*pow(-1 + x,8)*(21358847 + 370440*pow(M_PI,2) - 15240960*zeta2 + 3254580*pow(log(Q2/muF2),2) + 252*log(1 - x)*(49687 + 51660*log(1 - x)) + 126*log(Q2/muF2)*(49687 + 103320*log(1 - x))) + 55566000*pow(-1 + x,2)*(9*pow(log(Q2/muF2),2) + 18*log(Q2/muF2)*(-1 + 2*log(1 - x)) + 2*(27 + pow(M_PI,2) - 24*zeta2 + 18*(-1 + log(1 - x))*log(1 - x)))))/(4.000752e9*pow(M_PI,2));
	}
}
double DY_NNLO_CD_full(double x){
	return (pow(alphas_muR,2)*CF*TF*(160*(-1 + x) - 24*(-6 + 4/x + x)*zeta3 - 16*(5 + 4*x)*log(x) + 8*(10 + x)*zeta2*log(x) - 52*x*pow(log(x),2) - (16*x*pow(log(x),3))/3. + 8*(5 - 4*x)*Li2(1 - x) - 8*(-10 + x)*log(x)*Li2(1 - x) + 32*(5/x + 2*x)*log(x)*Li2(-x) + 40*(1 + x)*(zeta2 + 2*log(x)*log(1 + x) + 2*Li2(-x)) + 8*(-6 + 4/x + 3*x)*Li3(1 - x) - 16*(-10 + 10/x + 3*x)*Li3(-x) - (8*(2 + 2*x + pow(x,2))*(log(1 + x)*(6*zeta2 - 5*pow(log(x),2) + 6*log(x)*log(1 + x) + 12*Li2(-x)) - 4*S12(1 - x) + 12*S12(-x)))/x))/(16.*pow(M_PI,2));
}

double DY_NNLO_CD_expansion(double x, int power){
	if(power==1){
		return (pow(alphas_muR,2)*CF*TF*(9*zeta3 + 5*pow(M_PI,2)*(-1 + log(8)) - 30*zeta2*(-1 + log(8)) - 9*zeta3))/(6.*pow(M_PI,2));
	}
	if(power==2){
		return -(pow(alphas_muR,2)*CF*TF*(-1 + x)*(pow(M_PI,2)*(1 + log(64)) - 6*(zeta2 + 9*zeta3 + zeta2*log(64) - 9*zeta3)))/(12.*pow(M_PI,2));
	}
	if(power==3){
		return (pow(alphas_muR,2)*CF*TF*pow(-1 + x,2)*(7 - 96*zeta3 + zeta2*(18 - 96*log(2)) + pow(M_PI,2)*(-3 + 16*log(2)) + 96*zeta3))/(16.*pow(M_PI,2));
	}
	if(power==4){
		return (pow(alphas_muR,2)*CF*TF*pow(1 - x,3)*(401 - 2592*zeta3 + 36*zeta2*(29 - 72*log(2)) + 6*pow(M_PI,2)*(-29 + 72*log(2)) + 2592*zeta3))/(432.*pow(M_PI,2));
	}
	if(power==5){
		return (pow(alphas_muR,2)*CF*TF*pow(-1 + x,4)*(3991 - 1671*pow(M_PI,2) + 10026*zeta2 - 20736*zeta3 + 3456*(pow(M_PI,2) - 6*zeta2)*log(2) + 20736*zeta3))/(3456.*pow(M_PI,2));
	}
	if(power==6){
		return (pow(alphas_muR,2)*CF*TF*pow(1 - x,5)*(533843 - 227925*pow(M_PI,2) + 1367550*zeta2 - 2592000*zeta3 + 432000*(pow(M_PI,2) - 6*zeta2)*log(2) + 2592000*zeta3))/(432000.*pow(M_PI,2));
	}
	if(power==7){
		return (pow(alphas_muR,2)*CF*TF*pow(-1 + x,6)*(1075556 - 480225*pow(M_PI,2) + 2881350*zeta2 - 5184000*zeta3 + 864000*(pow(M_PI,2) - 6*zeta2)*log(2) + 5184000*zeta3))/(864000.*pow(M_PI,2));
	}
	if(power==8){
		return (pow(alphas_muR,2)*CF*TF*pow(1 - x,7)*(60262802 + 170571450*zeta2 - 296352000*zeta3 + 1225*(-23207*pow(M_PI,2) + 40320*(pow(M_PI,2) - 6*zeta2)*log(2)) + 296352000*zeta3))/(4.9392e7*pow(M_PI,2));
	}
	if(power==9){
		return (pow(alphas_muR,2)*CF*TF*pow(-1 + x,8)*(5598119972 + 16792375950*zeta2 - 28449792000*zeta3 + 11025*(-253853*pow(M_PI,2) + 430080*(pow(M_PI,2) - 6*zeta2)*log(2)) + 28449792000*zeta3))/(4.741632e9*pow(M_PI,2));
	}
	if(power==10){
		return (pow(alphas_muR,2)*CF*TF*pow(1 - x,9)*(145435007884 + 462125950650*zeta2 - 768144384000*zeta3 + 33075*(-2328677*pow(M_PI,2) + 3870720*(pow(M_PI,2) - 6*zeta2)*log(2)) + 768144384000*zeta3))/(1.28024064e11*pow(M_PI,2));
	}
	if(power==-1){
		return DY_NNLO_CD_full(x)-((pow(alphas_muR,2)*CF*TF*(10668672000*(1 - x)*(pow(M_PI,2)*(1 + log(64)) - 6*(zeta2 + 9*zeta3 + zeta2*log(64) - 9*zeta3)) + 21337344000*(9*zeta3 + 5*pow(M_PI,2)*(-1 + log(8)) - 30*zeta2*(-1 + log(8)) - 9*zeta3) + 8001504000*pow(-1 + x,2)*(7 - 96*zeta3 + zeta2*(18 - 96*log(2)) + pow(M_PI,2)*(-3 + 16*log(2)) + 96*zeta3) + 296352000*pow(1 - x,3)*(401 - 2592*zeta3 + 36*zeta2*(29 - 72*log(2)) + 6*pow(M_PI,2)*(-29 + 72*log(2)) + 2592*zeta3) + 37044000*pow(-1 + x,4)*(3991 - 1671*pow(M_PI,2) + 10026*zeta2 - 20736*zeta3 + 3456*(pow(M_PI,2) - 6*zeta2)*log(2) + 20736*zeta3) + 296352*pow(1 - x,5)*(533843 - 227925*pow(M_PI,2) + 1367550*zeta2 - 2592000*zeta3 + 432000*(pow(M_PI,2) - 6*zeta2)*log(2) + 2592000*zeta3) + 148176*pow(-1 + x,6)*(1075556 - 480225*pow(M_PI,2) + 2881350*zeta2 - 5184000*zeta3 + 864000*(pow(M_PI,2) - 6*zeta2)*log(2) + 5184000*zeta3) + 2592*pow(1 - x,7)*(60262802 + 170571450*zeta2 - 296352000*zeta3 + 1225*(-23207*pow(M_PI,2) + 40320*(pow(M_PI,2) - 6*zeta2)*log(2)) + 296352000*zeta3) + 27*pow(-1 + x,8)*(5598119972 + 16792375950*zeta2 - 28449792000*zeta3 + 11025*(-253853*pow(M_PI,2) + 430080*(pow(M_PI,2) - 6*zeta2)*log(2)) + 28449792000*zeta3) + pow(1 - x,9)*(145435007884 + 462125950650*zeta2 - 768144384000*zeta3 + 33075*(-2328677*pow(M_PI,2) + 3870720*(pow(M_PI,2) - 6*zeta2)*log(2)) + 768144384000*zeta3)))/(1.28024064e11*pow(M_PI,2)));
	}
	if(power==-2){
		return (pow(alphas_muR,2)*CF*TF*(10668672000*(1 - x)*(pow(M_PI,2)*(1 + log(64)) - 6*(zeta2 + 9*zeta3 + zeta2*log(64) - 9*zeta3)) + 21337344000*(9*zeta3 + 5*pow(M_PI,2)*(-1 + log(8)) - 30*zeta2*(-1 + log(8)) - 9*zeta3) + 8001504000*pow(-1 + x,2)*(7 - 96*zeta3 + zeta2*(18 - 96*log(2)) + pow(M_PI,2)*(-3 + 16*log(2)) + 96*zeta3) + 296352000*pow(1 - x,3)*(401 - 2592*zeta3 + 36*zeta2*(29 - 72*log(2)) + 6*pow(M_PI,2)*(-29 + 72*log(2)) + 2592*zeta3) + 37044000*pow(-1 + x,4)*(3991 - 1671*pow(M_PI,2) + 10026*zeta2 - 20736*zeta3 + 3456*(pow(M_PI,2) - 6*zeta2)*log(2) + 20736*zeta3) + 296352*pow(1 - x,5)*(533843 - 227925*pow(M_PI,2) + 1367550*zeta2 - 2592000*zeta3 + 432000*(pow(M_PI,2) - 6*zeta2)*log(2) + 2592000*zeta3) + 148176*pow(-1 + x,6)*(1075556 - 480225*pow(M_PI,2) + 2881350*zeta2 - 5184000*zeta3 + 864000*(pow(M_PI,2) - 6*zeta2)*log(2) + 5184000*zeta3) + 2592*pow(1 - x,7)*(60262802 + 170571450*zeta2 - 296352000*zeta3 + 1225*(-23207*pow(M_PI,2) + 40320*(pow(M_PI,2) - 6*zeta2)*log(2)) + 296352000*zeta3) + 27*pow(-1 + x,8)*(5598119972 + 16792375950*zeta2 - 28449792000*zeta3 + 11025*(-253853*pow(M_PI,2) + 430080*(pow(M_PI,2) - 6*zeta2)*log(2)) + 28449792000*zeta3) + pow(1 - x,9)*(145435007884 + 462125950650*zeta2 - 768144384000*zeta3 + 33075*(-2328677*pow(M_PI,2) + 3870720*(pow(M_PI,2) - 6*zeta2)*log(2)) + 768144384000*zeta3)))/(1.28024064e11*pow(M_PI,2));
	}
}
double DY_NNLO_CE_full(double x){
	return (pow(alphas_muR,2)*CF*(-CA/2. + CF)*(2*(-9 + 7*x)*log(x) - 4*(1 + 3*x)*pow(log(x),2) + 8*(3 + x)*Li2(1 - x) + 4*(1 + x)*(zeta2 + 4*log(1 - x)*log(x) + 2*log(x)*log(1 + x) + 2*Li2(-x)) + log(Q2/muF2)*(-16*(-1 + x) + 8*(1 + x)*log(x) - (4*(1 + pow(x,2))*(2*zeta2 - pow(log(x),2) + 4*log(x)*log(1 + x) + 4*Li2(-x)))/(1 + x)) + (1 - x)*(-34 + 8*zeta3 + 32*log(1 - x) + 4*zeta2*log(x) - (2*pow(log(x),3))/3. - (4*pow(M_PI,2)*log(1 + x))/3. + 4*pow(log(x),2)*log(1 + x) - 8*log(x)*pow(log(1 + x),2) - 16*log(1 + x)*Li2(-x) + 8*Li3(-x) - 16*S12(-x)) - (4*(1 + pow(x,2))*(3*zeta3 + 12*zeta2*log(1 - x) - 9*zeta2*log(x) - 6*log(1 - x)*pow(log(x),2) + 2*pow(log(x),3) + 6*zeta2*log(1 + x) + 24*log(1 - x)*log(x)*log(1 + x) - 21*pow(log(x),2)*log(1 + x) + 6*log(x)*pow(log(1 + x),2) - 18*log(x)*Li2(1 - x) + 24*log(1 - x)*Li2(-x) - 24*log(x)*Li2(-x) + 12*log(1 + x)*Li2(-x) + 24*Li3(1 - x) + 6*Li3(-x) - 24*Li3((1 - x)/(1 + x)) + 24*Li3((-1 + x)/(1 + x)) - 24*S12(1 - x) + 12*S12(-x)))/(3.*(1 + x))))/(16.*pow(M_PI,2));
}

double DY_NNLO_CE_expansion(double x, int power){
	if(power==1){
		return -(pow(alphas_muR,2)*(CA - 2*CF)*CF*(12*zeta2 - 6*zeta3 + pow(M_PI,2)*(-2 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 6*zeta3))/(48.*pow(M_PI,2));
	}
	if(power==2){
		return (pow(alphas_muR,2)*(CA - 2*CF)*CF*(1 - x)*(-6 + 36*zeta2 - 30*zeta3 + pow(M_PI,2)*(-6 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 30*zeta3))/(96.*pow(M_PI,2));
	}
	if(power==3){
		return -((CA - 2*CF)*CF*pow(alphas_muR - alphas_muR*x,2)*(-30*zeta2 - 6*zeta3 + pow(M_PI,2)*(5 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 6*zeta3))/(192.*pow(M_PI,2));
	}
	if(power==4){
		return -(pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(1 - x,3)*(pow(M_PI,2)*(31 + log(64)) + 6*pow(M_PI,2)*log(Q2/muF2) - 36*zeta2*log((2*Q2)/muF2) + 12*(pow(M_PI,2) - 6*zeta2)*log(1 - x) - 6*(2 + 31*zeta2 + 3*zeta3 - 3*zeta3)))/(1152.*pow(M_PI,2));
	}
	if(power==5){
		return -(pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,4)*(-64 - 222*zeta2 - 12*zeta3 + pow(M_PI,2)*(37 + log(16)) + 4*pow(M_PI,2)*log(Q2/muF2) - 24*zeta2*log((2*Q2)/muF2) + 8*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 12*zeta3))/(1536.*pow(M_PI,2));
	}
	if(power==6){
		return -(pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(1 - x,5)*(-22992 + 4745*pow(M_PI,2) - 28470*zeta2 - 900*zeta3 + 2880*log(Q2/muF2) + 5760*log(1 - x) + 300*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 900*zeta3))/(230400.*pow(M_PI,2));
	}
	if(power==7){
		return -(pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,6)*(-178144 + 24075*pow(M_PI,2) - 144450*zeta2 - 2700*zeta3 + 34560*log(Q2/muF2) + 69120*log(1 - x) + 900*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 2700*zeta3))/(1.3824e6*pow(M_PI,2));
	}
	if(power==8){
		return -(pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(1 - x,7)*(2003715*pow(M_PI,2) + 4596480*log(Q2/muF2) + 9192960*log(1 - x) + 44100*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) - 2*(8859856 + 6011145*zeta2 + 66150*zeta3 - 66150*zeta3)))/(1.354752e8*pow(M_PI,2));
	}
	if(power==9){
		return -(pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,8)*(-12677248 + 1376361*pow(M_PI,2) - 8258166*zeta2 - 52920*zeta3 + 4257792*log(Q2/muF2) + 8515584*log(1 - x) + 17640*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 52920*zeta3))/(1.0838016e8*pow(M_PI,2));
	}
	if(power==10){
		return -(pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(1 - x,9)*(-940359296 + 107784495*pow(M_PI,2) - 646706970*zeta2 - 2381400*zeta3 + 409489920*log(Q2/muF2) + 818979840*log(1 - x) + 793800*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 2381400*zeta3))/(9.7542144e9*pow(M_PI,2));
	}
	if(power==-1){
		return DY_NNLO_CE_full(x)-((pow(alphas_muR,2)*(CA - 2*CF)*CF*(-72*pow(1 - x,7)*(2003715*pow(M_PI,2) + 4596480*log(Q2/muF2) + 9192960*log(1 - x) + 44100*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) - 2*(8859856 + 6011145*zeta2 + 66150*zeta3 - 66150*zeta3)) - 8467200*pow(1 - x,3)*(pow(M_PI,2)*(31 + log(64)) + 6*pow(M_PI,2)*log(Q2/muF2) - 36*zeta2*log((2*Q2)/muF2) + 12*(pow(M_PI,2) - 6*zeta2)*log(1 - x) - 6*(2 + 31*zeta2 + 3*zeta3 - 3*zeta3)) - 203212800*(12*zeta2 - 6*zeta3 + pow(M_PI,2)*(-2 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 6*zeta3) - 50803200*pow(-1 + x,2)*(-30*zeta2 - 6*zeta3 + pow(M_PI,2)*(5 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 6*zeta3) - 6350400*pow(-1 + x,4)*(-64 - 222*zeta2 - 12*zeta3 + pow(M_PI,2)*(37 + log(16)) + 4*pow(M_PI,2)*log(Q2/muF2) - 24*zeta2*log((2*Q2)/muF2) + 8*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 12*zeta3) + 101606400*(1 - x)*(-6 + 36*zeta2 - 30*zeta3 + pow(M_PI,2)*(-6 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 30*zeta3) - 42336*pow(1 - x,5)*(-22992 + 4745*pow(M_PI,2) - 28470*zeta2 - 900*zeta3 + 2880*log(Q2/muF2) + 5760*log(1 - x) + 300*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 900*zeta3) - 7056*pow(-1 + x,6)*(-178144 + 24075*pow(M_PI,2) - 144450*zeta2 - 2700*zeta3 + 34560*log(Q2/muF2) + 69120*log(1 - x) + 900*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 2700*zeta3) - 90*pow(-1 + x,8)*(-12677248 + 1376361*pow(M_PI,2) - 8258166*zeta2 - 52920*zeta3 + 4257792*log(Q2/muF2) + 8515584*log(1 - x) + 17640*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 52920*zeta3) - pow(1 - x,9)*(-940359296 + 107784495*pow(M_PI,2) - 646706970*zeta2 - 2381400*zeta3 + 409489920*log(Q2/muF2) + 818979840*log(1 - x) + 793800*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 2381400*zeta3)))/(9.7542144e9*pow(M_PI,2)));
	}
	if(power==-2){
		return (pow(alphas_muR,2)*(CA - 2*CF)*CF*(-72*pow(1 - x,7)*(2003715*pow(M_PI,2) + 4596480*log(Q2/muF2) + 9192960*log(1 - x) + 44100*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) - 2*(8859856 + 6011145*zeta2 + 66150*zeta3 - 66150*zeta3)) - 8467200*pow(1 - x,3)*(pow(M_PI,2)*(31 + log(64)) + 6*pow(M_PI,2)*log(Q2/muF2) - 36*zeta2*log((2*Q2)/muF2) + 12*(pow(M_PI,2) - 6*zeta2)*log(1 - x) - 6*(2 + 31*zeta2 + 3*zeta3 - 3*zeta3)) - 203212800*(12*zeta2 - 6*zeta3 + pow(M_PI,2)*(-2 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 6*zeta3) - 50803200*pow(-1 + x,2)*(-30*zeta2 - 6*zeta3 + pow(M_PI,2)*(5 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 6*zeta3) - 6350400*pow(-1 + x,4)*(-64 - 222*zeta2 - 12*zeta3 + pow(M_PI,2)*(37 + log(16)) + 4*pow(M_PI,2)*log(Q2/muF2) - 24*zeta2*log((2*Q2)/muF2) + 8*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 12*zeta3) + 101606400*(1 - x)*(-6 + 36*zeta2 - 30*zeta3 + pow(M_PI,2)*(-6 + log(4)) + 2*pow(M_PI,2)*log(Q2/muF2) - 12*zeta2*log((2*Q2)/muF2) + 4*(pow(M_PI,2) - 6*zeta2)*log(1 - x) + 30*zeta3) - 42336*pow(1 - x,5)*(-22992 + 4745*pow(M_PI,2) - 28470*zeta2 - 900*zeta3 + 2880*log(Q2/muF2) + 5760*log(1 - x) + 300*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 900*zeta3) - 7056*pow(-1 + x,6)*(-178144 + 24075*pow(M_PI,2) - 144450*zeta2 - 2700*zeta3 + 34560*log(Q2/muF2) + 69120*log(1 - x) + 900*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 2700*zeta3) - 90*pow(-1 + x,8)*(-12677248 + 1376361*pow(M_PI,2) - 8258166*zeta2 - 52920*zeta3 + 4257792*log(Q2/muF2) + 8515584*log(1 - x) + 17640*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 52920*zeta3) - pow(1 - x,9)*(-940359296 + 107784495*pow(M_PI,2) - 646706970*zeta2 - 2381400*zeta3 + 409489920*log(Q2/muF2) + 818979840*log(1 - x) + 793800*(pow(M_PI,2) - 6*zeta2)*(log((2*Q2)/muF2) + 2*log(1 - x)) + 2381400*zeta3)))/(9.7542144e9*pow(M_PI,2));
	}
}

// extra factor of 2 everywhere because I made a mistake in the matematica code
double DY_NNLO_CF_full(double x){
	return 2.*(pow(alphas_muR,2)*CF*(-CA/2. + CF)*(-30 + 56*x - 26*pow(x,2) + 4*(-7 + 6*x)*log(x) - (4*pow(-1 + x,2)*(9*pow(log(x),2) + 6*log(1 - x)*pow(log(x),2) + 2*pow(log(x),3) + 6*(3 + 2*log(x))*Li2(1 - x) + 12*log(x)*Li2(x) - 12*Li3(1 - x) - 12*Li3(x) + 12*zeta3))/3.))/(32.*pow(M_PI,2));
}

double DY_NNLO_CF_expansion(double x, int power){
	if(power==1){
		return 0;
	}
	if(power==2){
		return 0;
	}
	if(power==3){
		return 0;
	}
	if(power==4){
		return 2.*(pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,3))/(12.*pow(M_PI,2));
	}
	if(power==5){
		return 2.*(-5*pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,4))/(64.*pow(M_PI,2));
	}
	if(power==6){
		return 2.*(319*pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,5))/(4320.*pow(M_PI,2));
	}
	if(power==7){
		return 2.*(-809*pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,6))/(11520.*pow(M_PI,2));
	}
	if(power==8){
		return 2.*(33823*pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,7))/(504000.*pow(M_PI,2));
	}
	if(power==9){
		return 2.*(-19469*pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,8))/(302400.*pow(M_PI,2));
	}
	if(power==10){
		return 2.*(38257*pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,9))/(617400.*pow(M_PI,2));
	}
	if(power==-1){
		return DY_NNLO_CF_full(x)-2.*((pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,3)*(147874329 - 414483831*x + 669892119*pow(x,2) - 658426621*pow(x,3) + 390736424*pow(x,4) - 129259780*pow(x,5) + 18363360*pow(x,6)))/(2.96352e8*pow(M_PI,2)));
	}
	if(power==-2){
		return 2.*((pow(alphas_muR,2)*(CA - 2*CF)*CF*pow(-1 + x,3)*(147874329 - 414483831*x + 669892119*pow(x,2) - 658426621*pow(x,3) + 390736424*pow(x,4) - 129259780*pow(x,5) + 18363360*pow(x,6)))/(2.96352e8*pow(M_PI,2)));
	}
}


////////////////////////////////////////////////////////////
///
/// gg channel 
///
////////////////////////////////////////////////////////////

double vegas_DY_NNLO_gg_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NNLO_gg_full(z))*(jacobian*pdf_sum_gg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_DY_NNLO_gg_full_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NNLO_gg_full(z))*(jacobian*real(fit_sum_gg_charge_weighted(tau/z+k[1]*jacobian,tau/z))/z);
}
double vegas_DY_NNLO_gg_power(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NNLO_gg_expansion(z, fp->power))*(jacobian*pdf_sum_gg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_DY_NNLO_gg_power_fit(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NNLO_gg_expansion(z, fp->power))*(jacobian*real(fit_sum_gg_charge_weighted(tau/z+k[1]*jacobian,tau/z))/z);
}


///// NNLO functions
double DY_NNLO_gg_full(double x){
	return (pow(alphas_muR,2)*(-32 - 66*x + 98*pow(x,2) + 4*(5 + 9*x - 12*pow(x,2))*zeta2 + 4*(-1 - 2*x + 2*pow(x,2))*zeta3 + 2*(7 + 60*x - 67*pow(x,2))*log(1 - x) + 16*(-1 - 2*x + 3*pow(x,2))*pow(log(1 - x),2) + (-23 - 64*x + 105*pow(x,2))*log(x) + 4*(3 + 10*x + 10*pow(x,2))*zeta2*log(x) + 4*(1 + 8*x - 4*pow(x,2))*log(1 - x)*log(x) - 2*(3 + 7*x + 4*pow(x,2))*pow(log(x),2) - (2*(3 + 8*x + 8*pow(x,2))*pow(log(x),3))/3. - 2*pow(log(Q2/muF2),2)*(2 + 4*x - 6*pow(x,2) + pow(1 + 2*x,2)*log(x)) + log(Q2/muF2)*(7 + 60*x - 67*pow(x,2) + 16*(-1 - 2*x + 3*pow(x,2))*log(1 - x) + 2*(1 + 8*x - 4*pow(x,2))*log(x) + 2*pow(1 + 2*x,2)*(log(x)*(-4*log(1 - x) + log(x)) - 4*Li2(1 - x))) + 4*(-5 - 4*x + 14*pow(x,2))*Li2(1 - x) + 8*(2 + 4*x + pow(x,2))*log(x)*Li2(-x) + 8*(1 + x)*(log(x)*log(1 + x) + Li2(-x)) - 4*pow(1 + 2*x,2)*(log(1 - x)*(2*log(1 - x) - log(x))*log(x) + (4*log(1 - x) + log(x))*Li2(1 - x) - 4*Li3(1 - x)) + 8*(-1 - 2*x + pow(x,2))*Li3(-x) - 8*(1 + 10*x + 7*pow(x,2))*S12(1 - x) - 4*pow(1 + x,2)*(log(1 + x)*(2*zeta2 - 3*pow(log(x),2) + 2*log(x)*log(1 + x) + 4*Li2(-x)) + 4*S12(-x)) + (pow(CA,2)*(-47 - 144*x + 191*pow(x,2) - 2*(6 + 38*x + 75*pow(x,2))*log(x) + 2*(-2 + 2*x + 25*pow(x,2))*pow(log(x),2) - 24*pow(-1 + x,2)*S12(1 - x) + 4*pow(1 + x,2)*(2*zeta2 + 12*zeta3 + 6*zeta2*log(1 + x) + 4*log(x)*log(1 + x) - 9*pow(log(x),2)*log(1 + x) + 6*log(x)*pow(log(1 + x),2) + (4 - 18*log(x) + 12*log(1 + x))*Li2(-x) + 18*Li3(-x) + 12*S12(-x))))/(3.*(-1 + pow(CA,2)))))/(16.*pow(M_PI,2));
}

double DY_NNLO_gg_expansion(double x, int power){		
	if(power==1){
		return 0;
	}
	if(power==2){
		return (pow(alphas_muR,2)*(1 - x)*(3*pow(log(Q2/muF2),2) + 12*log(Q2/muF2)*(-1 + log(1 - x)) - 2*(-12 + pow(M_PI,2) - 6*(-2 + log(1 - x))*log(1 - x))))/(24.*pow(M_PI,2));
	}
	if(power==3){
		return (pow(alphas_muR - alphas_muR*x,2)*(-3*pow(log(Q2/muF2),2) - 12*log(Q2/muF2)*(-2 + log(1 - x)) + 2*(-27 + pow(M_PI,2) - 6*(-4 + log(1 - x))*log(1 - x))))/(16.*pow(M_PI,2));
	}
	if(power==4){
		return (pow(alphas_muR,2)*pow(1 - x,3)*(-461 + 455*pow(CA,2) - 12*(-1 + pow(CA,2))*pow(M_PI,2) + 18*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 24*(-1 + pow(CA,2))*log(1 - x)*(-13 + 3*log(1 - x)) + 12*(-1 + pow(CA,2))*log(Q2/muF2)*(-13 + 6*log(1 - x))))/(144.*(-1 + pow(CA,2))*pow(M_PI,2));
	}
	if(power==5){
		return (pow(alphas_muR,2)*pow(-1 + x,4)*(9*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 3*(-1 + pow(CA,2))*log(Q2/muF2)*(1 + 12*log(1 - x)) + 2*(4 + 3*pow(M_PI,2) - pow(CA,2)*(7 + 3*pow(M_PI,2)) + 3*(-1 + pow(CA,2))*log(1 - x)*(1 + 6*log(1 - x)))))/(288.*(-1 + pow(CA,2))*pow(M_PI,2));
	}
	if(power==6){
		return (pow(alphas_muR,2)*pow(1 - x,5)*(8623 - 11323*pow(CA,2) - 2400*(-1 + pow(CA,2))*pow(M_PI,2) + 3600*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 120*(-1 + pow(CA,2))*log(1 - x)*(277 + 120*log(1 - x)) + 60*(-1 + pow(CA,2))*log(Q2/muF2)*(277 + 240*log(1 - x))))/(216000.*(-1 + pow(CA,2))*pow(M_PI,2));
	}
	if(power==7){
		return (pow(alphas_muR,2)*pow(-1 + x,6)*(-2951 + 1201*pow(CA,2) - 1200*(-1 + pow(CA,2))*pow(M_PI,2) + 1800*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 120*(-1 + pow(CA,2))*log(1 - x)*(181 + 60*log(1 - x)) + 60*(-1 + pow(CA,2))*log(Q2/muF2)*(181 + 120*log(1 - x))))/(144000.*(-1 + pow(CA,2))*pow(M_PI,2));
	}
	if(power==8){
		return (pow(alphas_muR,2)*pow(1 - x,7)*(-1484437 + 1170837*pow(CA,2) - 176400*(-1 + pow(CA,2))*pow(M_PI,2) + 264600*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 840*(-1 + pow(CA,2))*log(1 - x)*(4001 + 1260*log(1 - x)) + 420*(-1 + pow(CA,2))*log(Q2/muF2)*(4001 + 2520*log(1 - x))))/(2.4696e7*(-1 + pow(CA,2))*pow(M_PI,2));
	}
	if(power==9){
		return (pow(alphas_muR,2)*pow(-1 + x,8)*(-9969413 + 8406313*pow(CA,2) - 764400*(-1 + pow(CA,2))*pow(M_PI,2) + 1146600*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 5040*(-1 + pow(CA,2))*log(1 - x)*(2897 + 910*log(1 - x)) + 2520*(-1 + pow(CA,2))*log(Q2/muF2)*(2897 + 1820*log(1 - x))))/(1.185408e8*(-1 + pow(CA,2))*pow(M_PI,2));
	}
	if(power==10){
		return (pow(alphas_muR,2)*pow(1 - x,9)*(-5863801 + 5062161*pow(CA,2) - 352800*(-1 + pow(CA,2))*pow(M_PI,2) + 529200*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 15120*(-1 + pow(CA,2))*log(1 - x)*(443 + 140*log(1 - x)) + 7560*(-1 + pow(CA,2))*log(Q2/muF2)*(443 + 280*log(1 - x))))/(5.92704e7*(-1 + pow(CA,2))*pow(M_PI,2));
	}
	if(power==-1){
		return DY_NNLO_gg_full(x)-((pow(alphas_muR,2)*(1 - x)*((4116000*pow(-1 + x,2)*(-461 + 455*pow(CA,2) - 12*(-1 + pow(CA,2))*pow(M_PI,2) + 18*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 24*(-1 + pow(CA,2))*log(1 - x)*(-13 + 3*log(1 - x)) + 12*(-1 + pow(CA,2))*log(Q2/muF2)*(-13 + 6*log(1 - x))))/(-1 + pow(CA,2)) + (4116*pow(1 - x,5)*(-2951 + 1201*pow(CA,2) - 1200*(-1 + pow(CA,2))*pow(M_PI,2) + 1800*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 120*(-1 + pow(CA,2))*log(1 - x)*(181 + 60*log(1 - x)) + 60*(-1 + pow(CA,2))*log(Q2/muF2)*(181 + 120*log(1 - x))))/(-1 + pow(CA,2)) + (2744*pow(-1 + x,4)*(8623 - 11323*pow(CA,2) - 2400*(-1 + pow(CA,2))*pow(M_PI,2) + 3600*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 120*(-1 + pow(CA,2))*log(1 - x)*(277 + 120*log(1 - x)) + 60*(-1 + pow(CA,2))*log(Q2/muF2)*(277 + 240*log(1 - x))))/(-1 + pow(CA,2)) + (10*pow(-1 + x,8)*(-5863801 + 5062161*pow(CA,2) - 352800*(-1 + pow(CA,2))*pow(M_PI,2) + 529200*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 15120*(-1 + pow(CA,2))*log(1 - x)*(443 + 140*log(1 - x)) + 7560*(-1 + pow(CA,2))*log(Q2/muF2)*(443 + 280*log(1 - x))))/(-1 + pow(CA,2)) + (5*pow(1 - x,7)*(-9969413 + 8406313*pow(CA,2) - 764400*(-1 + pow(CA,2))*pow(M_PI,2) + 1146600*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 5040*(-1 + pow(CA,2))*log(1 - x)*(2897 + 910*log(1 - x)) + 2520*(-1 + pow(CA,2))*log(Q2/muF2)*(2897 + 1820*log(1 - x))))/(-1 + pow(CA,2)) + (24*pow(-1 + x,6)*(-1484437 + 1170837*pow(CA,2) - 176400*(-1 + pow(CA,2))*pow(M_PI,2) + 264600*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 840*(-1 + pow(CA,2))*log(1 - x)*(4001 + 1260*log(1 - x)) + 420*(-1 + pow(CA,2))*log(Q2/muF2)*(4001 + 2520*log(1 - x))))/(-1 + pow(CA,2)) + 37044000*(1 - x)*(-3*pow(log(Q2/muF2),2) - 12*log(Q2/muF2)*(-2 + log(1 - x)) + 2*(-27 + pow(M_PI,2) - 6*(-4 + log(1 - x))*log(1 - x))) + 24696000*(3*pow(log(Q2/muF2),2) + 12*log(Q2/muF2)*(-1 + log(1 - x)) - 2*(-12 + pow(M_PI,2) - 6*(-2 + log(1 - x))*log(1 - x))) + (2058000*pow(1 - x,3)*(9*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 3*(-1 + pow(CA,2))*log(Q2/muF2)*(1 + 12*log(1 - x)) + 2*(4 + 3*pow(M_PI,2) - pow(CA,2)*(7 + 3*pow(M_PI,2)) + 3*(-1 + pow(CA,2))*log(1 - x)*(1 + 6*log(1 - x)))))/(-1 + pow(CA,2))))/(5.92704e8*pow(M_PI,2)));
	}
	if(power==-2){
		return (pow(alphas_muR,2)*(1 - x)*((4116000*pow(-1 + x,2)*(-461 + 455*pow(CA,2) - 12*(-1 + pow(CA,2))*pow(M_PI,2) + 18*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 24*(-1 + pow(CA,2))*log(1 - x)*(-13 + 3*log(1 - x)) + 12*(-1 + pow(CA,2))*log(Q2/muF2)*(-13 + 6*log(1 - x))))/(-1 + pow(CA,2)) + (4116*pow(1 - x,5)*(-2951 + 1201*pow(CA,2) - 1200*(-1 + pow(CA,2))*pow(M_PI,2) + 1800*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 120*(-1 + pow(CA,2))*log(1 - x)*(181 + 60*log(1 - x)) + 60*(-1 + pow(CA,2))*log(Q2/muF2)*(181 + 120*log(1 - x))))/(-1 + pow(CA,2)) + (2744*pow(-1 + x,4)*(8623 - 11323*pow(CA,2) - 2400*(-1 + pow(CA,2))*pow(M_PI,2) + 3600*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 120*(-1 + pow(CA,2))*log(1 - x)*(277 + 120*log(1 - x)) + 60*(-1 + pow(CA,2))*log(Q2/muF2)*(277 + 240*log(1 - x))))/(-1 + pow(CA,2)) + (10*pow(-1 + x,8)*(-5863801 + 5062161*pow(CA,2) - 352800*(-1 + pow(CA,2))*pow(M_PI,2) + 529200*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 15120*(-1 + pow(CA,2))*log(1 - x)*(443 + 140*log(1 - x)) + 7560*(-1 + pow(CA,2))*log(Q2/muF2)*(443 + 280*log(1 - x))))/(-1 + pow(CA,2)) + (5*pow(1 - x,7)*(-9969413 + 8406313*pow(CA,2) - 764400*(-1 + pow(CA,2))*pow(M_PI,2) + 1146600*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 5040*(-1 + pow(CA,2))*log(1 - x)*(2897 + 910*log(1 - x)) + 2520*(-1 + pow(CA,2))*log(Q2/muF2)*(2897 + 1820*log(1 - x))))/(-1 + pow(CA,2)) + (24*pow(-1 + x,6)*(-1484437 + 1170837*pow(CA,2) - 176400*(-1 + pow(CA,2))*pow(M_PI,2) + 264600*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 840*(-1 + pow(CA,2))*log(1 - x)*(4001 + 1260*log(1 - x)) + 420*(-1 + pow(CA,2))*log(Q2/muF2)*(4001 + 2520*log(1 - x))))/(-1 + pow(CA,2)) + 37044000*(1 - x)*(-3*pow(log(Q2/muF2),2) - 12*log(Q2/muF2)*(-2 + log(1 - x)) + 2*(-27 + pow(M_PI,2) - 6*(-4 + log(1 - x))*log(1 - x))) + 24696000*(3*pow(log(Q2/muF2),2) + 12*log(Q2/muF2)*(-1 + log(1 - x)) - 2*(-12 + pow(M_PI,2) - 6*(-2 + log(1 - x))*log(1 - x))) + (2058000*pow(1 - x,3)*(9*(-1 + pow(CA,2))*pow(log(Q2/muF2),2) + 3*(-1 + pow(CA,2))*log(Q2/muF2)*(1 + 12*log(1 - x)) + 2*(4 + 3*pow(M_PI,2) - pow(CA,2)*(7 + 3*pow(M_PI,2)) + 3*(-1 + pow(CA,2))*log(1 - x)*(1 + 6*log(1 - x)))))/(-1 + pow(CA,2))))/(5.92704e8*pow(M_PI,2));
	}
}


////////////////////////////////////////////////////////////
///
/// qg channel 
///
////////////////////////////////////////////////////////////

double vegas_DY_NNLO_qg_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NNLO_qg_full(z))*(jacobian*pdf_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_DY_NNLO_qg_full_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NNLO_qg_full(z))*(jacobian*real(fit_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z))/z);
}
double vegas_DY_NNLO_qg_power(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NNLO_qg_expansion(z, fp->power))*(jacobian*pdf_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_DY_NNLO_qg_power_fit(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return DY_LO_factor()*(DY_NNLO_qg_expansion(z, fp->power))*(jacobian*real(fit_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z))/z);
}

///// NNLO functions
double DY_NNLO_qg_full(double x){
	return (pow(alphas_muR,2)*TF*(beta0*log(muR2/muF2)*(1 + 6*x - 7*pow(x,2) + 2*(1 + 2*(-1 + x)*x)*(log(Q2/muF2) + 2*log(1 - x) - log(x))) + CA*(59.888888888888886 + 116/(27.*x) - (1226*x)/9. + (1837*pow(x,2))/27. + (2*pow(M_PI,2)*(15 - 8/x - 84*x + 107*pow(x,2)))/9. - (16*pow(M_PI,2)*(1 - x + 2*pow(x,2))*log(1 - x))/3. + (2*(-210 + 88/x + 75*x + 74*pow(x,2))*log(1 - x))/9. + (4*(6 + 8/x + 63*x - 77*pow(x,2))*pow(log(1 - x),2))/3. + (26*(1 - 2*x + 2*pow(x,2))*pow(log(1 - x),3))/3. - 8*x*log(x) + (8*pow(M_PI,2)*x*(-5 + 2*x)*log(x))/3. - (2*(-354 + 12*x + 457*pow(x,2))*log(x))/9. + 8*x*log(1 - x)*log(x) + 20*(1 - 2*x + 13*pow(x,2))*log(1 - x)*log(x) + 4*(1 + 22*x - 6*pow(x,2))*pow(log(1 - x),2)*log(x) - 4*x*pow(log(x),2) - (5 + (346*pow(x,2))/3.)*pow(log(x),2) + 4*(-3 - 14*x + 2*pow(x,2))*log(1 - x)*pow(log(x),2) + (2*(9 + 20*x)*pow(log(x),3))/3. + (2*pow(log(Q2/muF2),2)*(3 + 4/x + 24*x - 31*pow(x,2) + 6*(1 - 2*x + 2*pow(x,2))*log(1 - x) + 6*(1 + 4*x)*log(x)))/3. + 8*x*Li2(1 - x) + (4*(33 + 16/x + 90*x + 44*pow(x,2))*Li2(1 - x))/3. + 8*(7 + 10*x + 5*pow(x,2))*log(1 - x)*Li2(1 - x) + 8*(7 - 2*x)*x*log(x)*Li2(1 - x) + 8*(1 + 5*x + 4*pow(x,2))*(log(x)*log(1 + x) + Li2(-x)) + (2*log(Q2/muF2)*(-87 + 44/x - 12*x + 73*pow(x,2) - 12*pow(M_PI,2)*(1 - x + 2*pow(x,2)) + 6*(9 + 8/x + 54*x - 71*pow(x,2))*log(1 - x) + 54*(1 - 2*x + 2*pow(x,2))*pow(log(1 - x),2) + 18*(3 - 2*x + 28*pow(x,2))*log(x) + 36*(1 + 10*x - 2*pow(x,2))*log(1 - x)*log(x) - 36*(1 + 3*x)*pow(log(x),2) + 36*(3 + 6*x + 2*pow(x,2))*Li2(1 - x) - 36*(1 + 2*x + 2*pow(x,2))*(log(x)*log(1 + x) + Li2(-x))))/9. - 4*(15 + 34*x + 12*pow(x,2))*Li3(1 - x) - 4*(1 + 2*x + 2*pow(x,2))*(4*log(1 - x)*log(x)*log(1 + x) - 3*pow(log(x),2)*log(1 + x) + 4*(log(1 - x) - log(x))*Li2(-x) + 2*Li3(-x) - 4*Li3((1 - x)/(1 + x)) + 4*Li3((-1 + x)/(1 + x))) - 4*(1 + 4*x + 2*pow(x,2))*zeta3 + 4*(9 + 16*x + 4*pow(x,2))*(log(1 - x)*pow(log(x),2) + 2*log(x)*Li2(x) - 2*Li3(x) + 2*zeta3)) + CF*(-90.5 + 12*(-1 + x) + 233*x - (305*pow(x,2))/2. + (2*pow(M_PI,2)*(5 + 2*x - 12*pow(x,2)))/3. - 24*(-1 + x)*log(1 - x) + 2*(38 - 147*x + 88*pow(x,2))*log(1 - x) - 2*(23 - 80*x + 63*pow(x,2))*pow(log(1 - x),2) + 4*(-7 + 11*x)*log(x) + 4*pow(M_PI,2)*(1 - 2*x + 4*pow(x,2))*log(x) - (59 - 245*x + 174*pow(x,2))*log(x) + 8*(-3 + x)*log(1 - x)*log(x) + 4*(13 - 50*x + 48*pow(x,2))*log(1 - x)*log(x) - 6*(7 - 14*x + 22*pow(x,2))*pow(log(1 - x),2)*log(x) + 4*(-3 + x)*pow(log(x),2) - ((35 - 68*x + 4*pow(x,2))*pow(log(x),2))/2. + 8*(3 - 6*x + 10*pow(x,2))*log(1 - x)*pow(log(x),2) - ((17 - 34*x + 52*pow(x,2))*pow(log(x),3))/3. + 3*pow(log(Q2/muF2),2)*(-1 + 4*x + (4 - 8*x + 8*pow(x,2))*log(1 - x) + (-2 + 4*x - 8*pow(x,2))*log(x)) - 8*(-3 + x)*Li2(1 - x) + 2*(3 - 28*x + 40*pow(x,2))*Li2(1 - x) - 4*(3 - 6*x + 26*pow(x,2))*log(1 - x)*Li2(1 - x) + 4*(1 - 2*x)*log(x)*Li2(1 - x) + log(Q2/muF2)*(24 - 68*x + 22*pow(x,2) - 4*(8 - 34*x + 23*pow(x,2))*log(1 - x) - (4*(1 - 2*x + 2*pow(x,2))*(pow(M_PI,2) - 27*pow(log(1 - x),2)))/3. + 2*(5 - 40*x + 46*pow(x,2))*log(x) - 8*(5 - 10*x + 16*pow(x,2))*log(1 - x)*log(x) + 8*(1 - 2*x + 4*pow(x,2))*pow(log(x),2) - 48*pow(x,2)*Li2(1 - x)) - 16*(-1 + 2*x + 3*pow(x,2))*(log(x)*log(1 + x) + Li2(-x)) + 4*(-1 + 2*x + 18*pow(x,2))*Li3(1 - x) - 2*(11 - 22*x + 34*pow(x,2))*(log(1 - x)*pow(log(x),2) + 2*log(x)*Li2(x) - 2*Li3(x) + 2*zeta3) + (1 - 2*x + 2*pow(x,2))*((-8*pow(M_PI,2)*log(1 - x))/3. + (70*pow(log(1 - x),3))/3. - 16*log(x)*Li2(-x) + 32*Li3(-x) + 100*zeta3))))/(16.*pow(M_PI,2));
}

double DY_NNLO_qg_expansion(double x, int power){		
	if(power==1){
		return (pow(alphas_muR,2)*TF*(-12*CA + 6*CF*(-5 + pow(M_PI,2)) + 3*pow(log(Q2/muF2),2)*(9*CF + 4*(CA + 3*CF)*log(1 - x)) + 2*log(1 - x)*(9*(CA - 7*CF) - 2*(3*CA + 2*CF)*pow(M_PI,2) + 6*beta0*log(muR2/muF2) - 18*CF*log(1 - x) + (13*CA + 35*CF)*pow(log(1 - x),2)) + 2*log(Q2/muF2)*(6*CA - 33*CF - (3*CA + 2*CF)*pow(M_PI,2) + 3*beta0*log(muR2/muF2) + 18*log(1 - x)*(CF + (CA + 3*CF)*log(1 - x))) + 6*(CA + 38*CF)*zeta3))/(48.*pow(M_PI,2));
	}
	if(power==2){
		return (pow(alphas_muR,2)*TF*(1 - x)*(189*CA + 231*CF - 2*(8*CA + 5*CF)*pow(M_PI,2) - 3*beta0*log(muR2/muF2)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)) + 3*pow(log(Q2/muF2),2)*(4*CA + 3*CF - 4*(CA + 3*CF)*log(1 - x)) + log(Q2/muF2)*(-72*CA - 69*CF + 6*CA*pow(M_PI,2) + 4*CF*pow(M_PI,2) + 12*log(1 - x)*(7*CA + 17*CF - 3*(CA + 3*CF)*log(1 - x))) + log(1 - x)*(-3*(51*CA + 77*CF) + 4*(3*CA + 2*CF)*pow(M_PI,2) + log(1 - x)*(96*CA + 273*CF - 2*(13*CA + 35*CF)*log(1 - x))) - 6*(CA + 38*CF)*zeta3))/(24.*pow(M_PI,2));
	}
	if(power==3){
		return (TF*pow(alphas_muR - alphas_muR*x,2)*(-3*(453*CA + 962*CF) + 2*(55*CA + 36*CF)*pow(M_PI,2) + 12*beta0*log(muR2/muF2)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)) + 6*pow(log(Q2/muF2),2)*(-3*(4*CA + 9*CF) + 8*(CA + 3*CF)*log(1 - x)) - 2*log(Q2/muF2)*(-252*CA - 669*CF + 4*(3*CA + 2*CF)*pow(M_PI,2) + 24*log(1 - x)*(9*CA + 28*CF - 3*(CA + 3*CF)*log(1 - x))) + 2*log(1 - x)*(612*CA + 1713*CF - 8*(3*CA + 2*CF)*pow(M_PI,2) + log(1 - x)*(-3*(86*CA + 261*CF) + 4*(13*CA + 35*CF)*log(1 - x))) + 24*(CA + 38*CF)*zeta3))/(96.*pow(M_PI,2));
	}
	if(power==4){
		return (pow(alphas_muR,2)*TF*pow(1 - x,3)*(1195*CA + 3806*CF - 12*(6*CA + 7*CF)*pow(M_PI,2) + 36*beta0*log(muR2/muF2) + 54*(CA + 3*CF)*pow(log(Q2/muF2),2) + 3*log(1 - x)*(-310*CA - 1033*CF + 36*(4*CA + 9*CF)*log(1 - x)) + 3*log(Q2/muF2)*(-134*CA - 453*CF + 24*(5*CA + 13*CF)*log(1 - x))))/(216.*pow(M_PI,2));
	}
	if(power==5){
		return (pow(alphas_muR,2)*TF*pow(-1 + x,4)*(12305*CA + 11791*CF - 24*(57*CA + 40*CF)*pow(M_PI,2) + 504*beta0*log(muR2/muF2) + 12*(54*(2*CA + 3*CF)*pow(log(Q2/muF2),2) + log(1 - x)*(-5*(100*CA + 223*CF) + 342*(2*CA + 3*CF)*log(1 - x)) + log(Q2/muF2)*(-292*CA - 546*CF + 24*(25*CA + 41*CF)*log(1 - x)))))/(6912.*pow(M_PI,2));
	}
	if(power==6){
		return (pow(alphas_muR,2)*TF*pow(1 - x,5)*(139537*CA + 57107*CF - 200*(113*CA + 56*CF)*pow(M_PI,2) + 6600*beta0*log(muR2/muF2) + 600*(40*CA + 39*CF)*pow(log(Q2/muF2),2) + 80*log(Q2/muF2)*(-164*CA - 261*CF + 30*(51*CA + 61*CF)*log(1 - x)) + 20*log(1 - x)*(348*CA - 601*CF + 30*(226*CA + 255*CF)*log(1 - x))))/(144000.*pow(M_PI,2));
	}
	if(power==7){
		return (pow(alphas_muR,2)*TF*pow(-1 + x,6)*(182023*CA + 33254*CF - 600*(50*CA + 19*CF)*pow(M_PI,2) + 7200*beta0*log(muR2/muF2) + 900*(38*CA + 27*CF)*pow(log(Q2/muF2),2) + 30*log(1 - x)*(2845*CA + 2997*CF + 600*(10*CA + 9*CF)*log(1 - x)) + 60*log(Q2/muF2)*(290*CA + 279*CF + 60*(46*CA + 43*CF)*log(1 - x))))/(216000.*pow(M_PI,2));
	}
	if(power==8){
		return (pow(alphas_muR,2)*TF*pow(1 - x,7)*(64445971*CA + 7239998*CF - 58800*(163*CA + 50*CF)*pow(M_PI,2) + 1940400*beta0*log(muR2/muF2) + 840*(210*(65*CA + 36*CF)*pow(log(Q2/muF2),2) + log(Q2/muF2)*(14968*CA + 13731*CF + 1680*(38*CA + 29*CF)*log(1 - x)) + log(1 - x)*(49480*CA + 51243*CF + 210*(326*CA + 243*CF)*log(1 - x)))))/(7.4088e7*pow(M_PI,2));
	}
	if(power==9){
		return (pow(alphas_muR,2)*TF*pow(-1 + x,8)*(442321415*CA + 41879017*CF - 117600*(499*CA + 128*CF)*pow(M_PI,2) + 10231200*beta0*log(muR2/muF2) + 840*(420*(206*CA + 93*CF)*pow(log(Q2/muF2),2) + 8*log(Q2/muF2)*(15823*CA + 13064*CF + 210*(235*CA + 151*CF)*log(1 - x)) + 3*log(1 - x)*(123244*CA + 119319*CF + 140*(998*CA + 633*CF)*log(1 - x)))))/(4.741632e8*pow(M_PI,2));
	}
	if(power==10){
		return (pow(alphas_muR,2)*TF*pow(1 - x,9)*(6407234287*CA + 592700725*CF - 1058400*(727*CA + 160*CF)*pow(M_PI,2) + 117482400*beta0*log(muR2/muF2) + 2520*(1260*(308*CA + 117*CF)*pow(log(Q2/muF2),2) + 4*log(Q2/muF2)*(61*(2725*CA + 2031*CF) + 1260*(345*CA + 191*CF)*log(1 - x)) + log(1 - x)*(1817636*CA + 1635115*CF + 1260*(1454*CA + 801*CF)*log(1 - x)))))/(6.4012032e9*pow(M_PI,2));
	}
	if(power==-1){
		return DY_NNLO_qg_full(x)-((pow(alphas_muR,2)*TF*(296352000*pow(1 - x,3)*(1195*CA + 3806*CF - 12*(6*CA + 7*CF)*pow(M_PI,2) + 36*beta0*log(muR2/muF2) + 54*(CA + 3*CF)*pow(log(Q2/muF2),2) + 3*log(1 - x)*(-310*CA - 1033*CF + 36*(4*CA + 9*CF)*log(1 - x)) + 3*log(Q2/muF2)*(-134*CA - 453*CF + 24*(5*CA + 13*CF)*log(1 - x))) + 296352*pow(-1 + x,6)*(182023*CA + 33254*CF - 600*(50*CA + 19*CF)*pow(M_PI,2) + 7200*beta0*log(muR2/muF2) + 900*(38*CA + 27*CF)*pow(log(Q2/muF2),2) + 30*log(1 - x)*(2845*CA + 2997*CF + 600*(10*CA + 9*CF)*log(1 - x)) + 60*log(Q2/muF2)*(290*CA + 279*CF + 60*(46*CA + 43*CF)*log(1 - x))) + 444528*pow(1 - x,5)*(139537*CA + 57107*CF - 200*(113*CA + 56*CF)*pow(M_PI,2) + 6600*beta0*log(muR2/muF2) + 600*(40*CA + 39*CF)*pow(log(Q2/muF2),2) + 80*log(Q2/muF2)*(-164*CA - 261*CF + 30*(51*CA + 61*CF)*log(1 - x)) + 20*log(1 - x)*(348*CA - 601*CF + 30*(226*CA + 255*CF)*log(1 - x))) + 9261000*pow(-1 + x,4)*(12305*CA + 11791*CF - 24*(57*CA + 40*CF)*pow(M_PI,2) + 504*beta0*log(muR2/muF2) + 12*(54*(2*CA + 3*CF)*pow(log(Q2/muF2),2) + log(1 - x)*(-5*(100*CA + 223*CF) + 342*(2*CA + 3*CF)*log(1 - x)) + log(Q2/muF2)*(-292*CA - 546*CF + 24*(25*CA + 41*CF)*log(1 - x)))) + 864*pow(1 - x,7)*(64445971*CA + 7239998*CF - 58800*(163*CA + 50*CF)*pow(M_PI,2) + 1940400*beta0*log(muR2/muF2) + 840*(210*(65*CA + 36*CF)*pow(log(Q2/muF2),2) + log(Q2/muF2)*(14968*CA + 13731*CF + 1680*(38*CA + 29*CF)*log(1 - x)) + log(1 - x)*(49480*CA + 51243*CF + 210*(326*CA + 243*CF)*log(1 - x)))) + 135*pow(-1 + x,8)*(442321415*CA + 41879017*CF - 117600*(499*CA + 128*CF)*pow(M_PI,2) + 10231200*beta0*log(muR2/muF2) + 840*(420*(206*CA + 93*CF)*pow(log(Q2/muF2),2) + 8*log(Q2/muF2)*(15823*CA + 13064*CF + 210*(235*CA + 151*CF)*log(1 - x)) + 3*log(1 - x)*(123244*CA + 119319*CF + 140*(998*CA + 633*CF)*log(1 - x)))) + 10*pow(1 - x,9)*(6407234287*CA + 592700725*CF - 1058400*(727*CA + 160*CF)*pow(M_PI,2) + 117482400*beta0*log(muR2/muF2) + 2520*(1260*(308*CA + 117*CF)*pow(log(Q2/muF2),2) + 4*log(Q2/muF2)*(61*(2725*CA + 2031*CF) + 1260*(345*CA + 191*CF)*log(1 - x)) + log(1 - x)*(1817636*CA + 1635115*CF + 1260*(1454*CA + 801*CF)*log(1 - x)))) + 2667168000*(1 - x)*(189*CA + 231*CF - 2*(8*CA + 5*CF)*pow(M_PI,2) - 3*beta0*log(muR2/muF2)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)) + 3*pow(log(Q2/muF2),2)*(4*CA + 3*CF - 4*(CA + 3*CF)*log(1 - x)) + log(Q2/muF2)*(-72*CA - 69*CF + 6*CA*pow(M_PI,2) + 4*CF*pow(M_PI,2) + 12*log(1 - x)*(7*CA + 17*CF - 3*(CA + 3*CF)*log(1 - x))) + log(1 - x)*(-3*(51*CA + 77*CF) + 4*(3*CA + 2*CF)*pow(M_PI,2) + log(1 - x)*(96*CA + 273*CF - 2*(13*CA + 35*CF)*log(1 - x))) - 6*(CA + 38*CF)*zeta3) + 1333584000*(-12*CA + 6*CF*(-5 + pow(M_PI,2)) + 3*pow(log(Q2/muF2),2)*(9*CF + 4*(CA + 3*CF)*log(1 - x)) + 2*log(1 - x)*(9*(CA - 7*CF) - 2*(3*CA + 2*CF)*pow(M_PI,2) + 6*beta0*log(muR2/muF2) - 18*CF*log(1 - x) + (13*CA + 35*CF)*pow(log(1 - x),2)) + 2*log(Q2/muF2)*(6*CA - 33*CF - (3*CA + 2*CF)*pow(M_PI,2) + 3*beta0*log(muR2/muF2) + 18*log(1 - x)*(CF + (CA + 3*CF)*log(1 - x))) + 6*(CA + 38*CF)*zeta3) + 666792000*pow(-1 + x,2)*(-3*(453*CA + 962*CF) + 2*(55*CA + 36*CF)*pow(M_PI,2) + 12*beta0*log(muR2/muF2)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)) + 6*pow(log(Q2/muF2),2)*(-3*(4*CA + 9*CF) + 8*(CA + 3*CF)*log(1 - x)) - 2*log(Q2/muF2)*(-252*CA - 669*CF + 4*(3*CA + 2*CF)*pow(M_PI,2) + 24*log(1 - x)*(9*CA + 28*CF - 3*(CA + 3*CF)*log(1 - x))) + 2*log(1 - x)*(612*CA + 1713*CF - 8*(3*CA + 2*CF)*pow(M_PI,2) + log(1 - x)*(-3*(86*CA + 261*CF) + 4*(13*CA + 35*CF)*log(1 - x))) + 24*(CA + 38*CF)*zeta3)))/(6.4012032e10*pow(M_PI,2)));
	}
	if(power==-2){
		return (pow(alphas_muR,2)*TF*(296352000*pow(1 - x,3)*(1195*CA + 3806*CF - 12*(6*CA + 7*CF)*pow(M_PI,2) + 36*beta0*log(muR2/muF2) + 54*(CA + 3*CF)*pow(log(Q2/muF2),2) + 3*log(1 - x)*(-310*CA - 1033*CF + 36*(4*CA + 9*CF)*log(1 - x)) + 3*log(Q2/muF2)*(-134*CA - 453*CF + 24*(5*CA + 13*CF)*log(1 - x))) + 296352*pow(-1 + x,6)*(182023*CA + 33254*CF - 600*(50*CA + 19*CF)*pow(M_PI,2) + 7200*beta0*log(muR2/muF2) + 900*(38*CA + 27*CF)*pow(log(Q2/muF2),2) + 30*log(1 - x)*(2845*CA + 2997*CF + 600*(10*CA + 9*CF)*log(1 - x)) + 60*log(Q2/muF2)*(290*CA + 279*CF + 60*(46*CA + 43*CF)*log(1 - x))) + 444528*pow(1 - x,5)*(139537*CA + 57107*CF - 200*(113*CA + 56*CF)*pow(M_PI,2) + 6600*beta0*log(muR2/muF2) + 600*(40*CA + 39*CF)*pow(log(Q2/muF2),2) + 80*log(Q2/muF2)*(-164*CA - 261*CF + 30*(51*CA + 61*CF)*log(1 - x)) + 20*log(1 - x)*(348*CA - 601*CF + 30*(226*CA + 255*CF)*log(1 - x))) + 9261000*pow(-1 + x,4)*(12305*CA + 11791*CF - 24*(57*CA + 40*CF)*pow(M_PI,2) + 504*beta0*log(muR2/muF2) + 12*(54*(2*CA + 3*CF)*pow(log(Q2/muF2),2) + log(1 - x)*(-5*(100*CA + 223*CF) + 342*(2*CA + 3*CF)*log(1 - x)) + log(Q2/muF2)*(-292*CA - 546*CF + 24*(25*CA + 41*CF)*log(1 - x)))) + 864*pow(1 - x,7)*(64445971*CA + 7239998*CF - 58800*(163*CA + 50*CF)*pow(M_PI,2) + 1940400*beta0*log(muR2/muF2) + 840*(210*(65*CA + 36*CF)*pow(log(Q2/muF2),2) + log(Q2/muF2)*(14968*CA + 13731*CF + 1680*(38*CA + 29*CF)*log(1 - x)) + log(1 - x)*(49480*CA + 51243*CF + 210*(326*CA + 243*CF)*log(1 - x)))) + 135*pow(-1 + x,8)*(442321415*CA + 41879017*CF - 117600*(499*CA + 128*CF)*pow(M_PI,2) + 10231200*beta0*log(muR2/muF2) + 840*(420*(206*CA + 93*CF)*pow(log(Q2/muF2),2) + 8*log(Q2/muF2)*(15823*CA + 13064*CF + 210*(235*CA + 151*CF)*log(1 - x)) + 3*log(1 - x)*(123244*CA + 119319*CF + 140*(998*CA + 633*CF)*log(1 - x)))) + 10*pow(1 - x,9)*(6407234287*CA + 592700725*CF - 1058400*(727*CA + 160*CF)*pow(M_PI,2) + 117482400*beta0*log(muR2/muF2) + 2520*(1260*(308*CA + 117*CF)*pow(log(Q2/muF2),2) + 4*log(Q2/muF2)*(61*(2725*CA + 2031*CF) + 1260*(345*CA + 191*CF)*log(1 - x)) + log(1 - x)*(1817636*CA + 1635115*CF + 1260*(1454*CA + 801*CF)*log(1 - x)))) + 2667168000*(1 - x)*(189*CA + 231*CF - 2*(8*CA + 5*CF)*pow(M_PI,2) - 3*beta0*log(muR2/muF2)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)) + 3*pow(log(Q2/muF2),2)*(4*CA + 3*CF - 4*(CA + 3*CF)*log(1 - x)) + log(Q2/muF2)*(-72*CA - 69*CF + 6*CA*pow(M_PI,2) + 4*CF*pow(M_PI,2) + 12*log(1 - x)*(7*CA + 17*CF - 3*(CA + 3*CF)*log(1 - x))) + log(1 - x)*(-3*(51*CA + 77*CF) + 4*(3*CA + 2*CF)*pow(M_PI,2) + log(1 - x)*(96*CA + 273*CF - 2*(13*CA + 35*CF)*log(1 - x))) - 6*(CA + 38*CF)*zeta3) + 1333584000*(-12*CA + 6*CF*(-5 + pow(M_PI,2)) + 3*pow(log(Q2/muF2),2)*(9*CF + 4*(CA + 3*CF)*log(1 - x)) + 2*log(1 - x)*(9*(CA - 7*CF) - 2*(3*CA + 2*CF)*pow(M_PI,2) + 6*beta0*log(muR2/muF2) - 18*CF*log(1 - x) + (13*CA + 35*CF)*pow(log(1 - x),2)) + 2*log(Q2/muF2)*(6*CA - 33*CF - (3*CA + 2*CF)*pow(M_PI,2) + 3*beta0*log(muR2/muF2) + 18*log(1 - x)*(CF + (CA + 3*CF)*log(1 - x))) + 6*(CA + 38*CF)*zeta3) + 666792000*pow(-1 + x,2)*(-3*(453*CA + 962*CF) + 2*(55*CA + 36*CF)*pow(M_PI,2) + 12*beta0*log(muR2/muF2)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)) + 6*pow(log(Q2/muF2),2)*(-3*(4*CA + 9*CF) + 8*(CA + 3*CF)*log(1 - x)) - 2*log(Q2/muF2)*(-252*CA - 669*CF + 4*(3*CA + 2*CF)*pow(M_PI,2) + 24*log(1 - x)*(9*CA + 28*CF - 3*(CA + 3*CF)*log(1 - x))) + 2*log(1 - x)*(612*CA + 1713*CF - 8*(3*CA + 2*CF)*pow(M_PI,2) + log(1 - x)*(-3*(86*CA + 261*CF) + 4*(13*CA + 35*CF)*log(1 - x))) + 24*(CA + 38*CF)*zeta3)))/(6.4012032e10*pow(M_PI,2));
	}
}
