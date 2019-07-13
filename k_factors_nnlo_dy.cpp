#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "deriv_pdf.h"
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


//////////////////////////////////////////////
/// LP contribution for NNLO DY qq to gg 
//////////////////////////////////////////////
double vegas_NNLO_qqbar_LP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqbar_LP(z)+NNLO_qqbar_b0_LP(z))*(jacobian*pdf_sum_qqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z - (1.-tau)*pdf_sum_qqbar_charge_weighted(tau+k[1]*(1.-tau),tau));
}
//////////////////////////////////////////////////////////////
/// correction term stemming from not integrating from 0 to 1
/// for NNLO qq to gg 
//////////////////////////////////////////////////////////////
double vegas_NNLO_qqbar_LP_correction(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0]; // from 0 to tau
	double jacobian = 1.-tau; //needed to transform the boundary dependent terms
	double x = tau+k[1]*jacobian;
	return LO_factor()*(NNLO_qqbar_LP(z)+NNLO_qqbar_b0_LP(z))*jacobian*(-pdf_sum_qqbar_charge_weighted(x,tau));
}

////////////////////////////////////////
/// NLP contribution DY qq to gg nnlo
////////////////////////////////////////
double vegas_NNLO_qqbar_NLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqbar_NLP(z))*(jacobian*pdf_sum_qqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}


////////////////////////////////////////
/// NNLP contribution DY qq to gg nnlo
////////////////////////////////////////
double vegas_NNLO_qqbar_NNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqbar_NNLP(z))*(jacobian*pdf_sum_qqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}




////////////////////////////////////////
/// NNNLP contribution DY qq to gg nnlo
////////////////////////////////////////
double vegas_NNLO_qqbar_NNNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqbar_NNNLP(z))*(jacobian*pdf_sum_qqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}

////////////////////////////////////////////////////////////////
/// constant NNLO qq to gg  contribution 
/// (contains  alphas^2*delta(1-z))
/// already integrated over z
/// so make sure here no z integral 
////////////////////////////////////////////////////////////////
double vegas_NNLO_qqbar_delta(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double jacobian = 1.;//-tau; //needed to transform the boundary dependent terms
	double x = /*tau+*/k[0]*jacobian;
	return LO_factor()*(NNLO_qqbar_delta()+NNLO_qqbar_b0_const())*(jacobian*pdf_sum_qqbar_charge_weighted(/*tau+*/k[0]*jacobian,tau));
}


////////////////////////////////////////
/// Full contribution DY qq nnlo
////////////////////////////////////////
double vegas_NNLO_qqbar_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	double charge_factor = 0;
	for(int i = 1; i <=5; i++){
		charge_factor += eq[i-1]*eq[i-1]; //still need the charge factor
			}
	double B2term = LO_factor()*(NNLO_qqbar_B2(z))*charge_factor*(jacobian*pdf_sum_qqbar(tau/z+k[1]*jacobian,tau/z)/z);
	return LO_factor()*(NNLO_qqbar_HCF(z)+NNLO_qqbar_HCA(z)+NNLO_qqbar_A2(z)+B2term+2.*(NNLO_qqbar_AC(z)+NNLO_qqbar_BC(z))+NNLO_qqbar_b0_nonconst(z))*(jacobian*pdf_sum_qqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}


double NNLO_qqbar_LP(double z){
	return (pow(alphas_Q,2)*CF*(9*pow(log(Q2/muF2),2)*(11*CA - 2*(18*CF + nF) - 48*CF*log(1 - z)) + 
       6*log(Q2/muF2)*(-67*CA + 144*CF + 10*nF + 18*CA*zeta2 + 72*CF*zeta2 + 6*(11*CA - 2*(9*CF + nF))*log(1 - z) - 
          216*CF*pow(log(1 - z),2)) + 2*(202*CA - 28*nF - 198*CA*zeta2 + 36*nF*zeta2 - 189*CA*zeta3 - 864*CF*zeta3 + 
          6*(CA*(-67 + 18*zeta2) + 2*(5*nF + 36*CF*(2 + zeta2)))*log(1 - z) + 18*(11*CA - 2*nF)*pow(log(1 - z),2) - 432*CF*pow(log(1 - z),3)
          )))/(108.*pow(M_PI,2)*(-1 + z));
}

double NNLO_qqbar_b0_LP(double x){
	return -(pow(alphas_Q,2)*CF*(11*CA - 2*nF)*log(muR2/muF2)*
      (log(Q2/muF2) + 2*log(1 - x)))/
   (6.*pow(M_PI,2)*(-1 + x));
	
}

double NNLO_qqbar_NLP(double x){
	return (pow(alphas_Q,2)*CF*(6256*CA - 8208*CF - 1312*nF + 
       96*nF*pow(M_PI,2) - 4032*CA*zeta2 - 6912*CF*zeta2 - 
       3024*CA*zeta3 - 13824*CF*zeta3 - 
       6816*CA*log(Q2/muF2) + 11664*CF*log(Q2/muF2) + 
       1056*nF*log(Q2/muF2) + 864*CA*zeta2*log(Q2/muF2) + 
       3456*CF*zeta2*log(Q2/muF2) + 
       792*CA*pow(log(Q2/muF2),2) - 
       864*CF*pow(log(Q2/muF2),2) - 
       144*nF*pow(log(Q2/muF2),2) - 12552*CA*log(1 - x) + 
       17064*CF*log(1 - x) + 2112*nF*log(1 - x) + 
       1728*CA*zeta2*log(1 - x) + 6912*CF*zeta2*log(1 - x) + 
       3168*CA*log(Q2/muF2)*log(1 - x) + 
       6912*CF*log(Q2/muF2)*log(1 - x) - 
       576*nF*log(Q2/muF2)*log(1 - x) - 
       3456*CF*pow(log(Q2/muF2),2)*log(1 - x) + 
       3168*CA*pow(log(1 - x),2) + 
       13392*CF*pow(log(1 - x),2) - 
       576*nF*pow(log(1 - x),2) - 
       10368*CF*log(Q2/muF2)*pow(log(1 - x),2) - 
       6912*CF*pow(log(1 - x),3) - 
       144*(11*CA - 2*nF)*log(muR2/muF2)*
        (-1 + log(Q2/muF2) + 2*log(1 - x)) + 
       216*(CA - 2*CF)*(-3 + 2*log(Q2/muF2) + 4*log(1 - x))))
    /(864.*pow(M_PI,2));
}



double NNLO_qqbar_NNLP(double x){
	return (pow(alphas_Q,2)*CF*(1 - x)*
     (-11792*CA - 7452*CF + 2012*nF - 48*nF*pow(M_PI,2) + 
       2232*CA*zeta2 + 10368*CF*zeta2 + 1512*CA*zeta3 + 
       6912*CF*zeta3 + 5532*CA*log(Q2/muF2) + 
       3456*CF*log(Q2/muF2) - 816*nF*log(Q2/muF2) - 
       432*CA*zeta2*log(Q2/muF2) - 
       1728*CF*zeta2*log(Q2/muF2) - 
       396*CA*pow(log(Q2/muF2),2) - 
       1296*CF*pow(log(Q2/muF2),2) + 
       72*nF*pow(log(Q2/muF2),2) + 11064*CA*log(1 - x) + 
       12960*CF*log(1 - x) - 1632*nF*log(1 - x) - 
       864*CA*zeta2*log(1 - x) - 3456*CF*zeta2*log(1 - x) - 
       1584*CA*log(Q2/muF2)*log(1 - x) - 
       10368*CF*log(Q2/muF2)*log(1 - x) + 
       288*nF*log(Q2/muF2)*log(1 - x) + 
       1728*CF*pow(log(Q2/muF2),2)*log(1 - x) - 
       1584*CA*pow(log(1 - x),2) - 
       13608*CF*pow(log(1 - x),2) + 
       288*nF*pow(log(1 - x),2) + 
       5184*CF*log(Q2/muF2)*pow(log(1 - x),2) + 
       3456*CF*pow(log(1 - x),3) + 
       72*(11*CA - 2*nF)*log(muR2/muF2)*
        (-1 + log(Q2/muF2) + 2*log(1 - x)) + 
       108*(CA - 2*CF)*(1 + log(Q2/muF2) + 2*log(1 - x))))/
   (864.*pow(M_PI,2));
}


double NNLO_qqbar_NNNLP(double x){
	return (pow(alphas_Q,2)*CF*pow(1 - x,2)*
     (1125*CA + 1599*CF - 154*nF - 72*CA*zeta2 - 
       576*CF*zeta2 + 132*CA*log(muR2/muF2) - 
       24*nF*log(muR2/muF2) - 264*CA*log(Q2/muF2) - 
       756*CF*log(Q2/muF2) + 48*nF*log(Q2/muF2) + 
       144*CF*pow(log(Q2/muF2),2) - 504*CA*log(1 - x) - 
       2022*CF*log(1 - x) + 96*nF*log(1 - x) + 
       1008*CF*log(Q2/muF2)*log(1 - x) + 
       1116*CF*pow(log(1 - x),2)))/(216.*pow(M_PI,2));
}

double NNLO_qqbar_delta(){
	return (pow(alphas_Q,2)*(CA*CF*(-127.91666666666667 + (592*zeta2)/9. - (12*pow(zeta2,2))/5. + 28*zeta3 + 
          (64.33333333333333 - 24*zeta3)*log(Q2/muF2) - 11*pow(log(Q2/muF2),2)) + 
       CF*nF*(21.166666666666668 - (112*zeta2)/9. + 8*zeta3 - (34*log(Q2/muF2))/3. + 2*pow(log(Q2/muF2),2)) + 
       pow(CF,2)*(127.75 - 70*zeta2 + (8*pow(zeta2,2))/5. - 60*zeta3 + (-93 + 24*zeta2 + 176*zeta3)*log(Q2/muF2) + 
          (18 - 32*zeta2)*pow(log(Q2/muF2),2))))/(16.*pow(M_PI,2));
}

double NNLO_qqbar_HCF(double x){
	return (pow(alphas_Q,2)*pow(CF,2)*
     (-24*(3 - 2*x) + 8*(3 - 2*x)*Li2(1 - x) - 
       (48 + 16*x)*Li2(1 - x) + 4*(64 + 3*x)*log(1 - x) + 
       (1 - x)*(64*zeta2 - 64*pow(log(1 - x),2)) - 
       16*log(x) - 8*(4 - 13*x)*log(x) + 
       16*(7 - 6*x)*log(1 - x)*log(x) - 
       (48 + 16*x)*log(1 - x)*log(x) - 
       16*(2 - x)*pow(log(x),2) + 
       (24 + 8*x)*pow(log(x),2) + 
       pow(log(Q2/muF2),2)*
        (-8*(5 + x) - (16*(1 + pow(x,2))*log(x))/(1 - x) + 
          8*(1 + x)*(-4*log(1 - x) + log(x))) + 
       log(Q2/muF2)*(8*(15 + 2*x) - 16*(7 - x)*log(1 - x) + 
          16*(2 - 3*x)*log(x) + 
          (1 + x)*(32*zeta2 + 32*Li2(1 - x) - 
             96*pow(log(1 - x),2) + 32*log(1 - x)*log(x) - 
             12*pow(log(x),2)) + 
          ((1 + pow(x,2))*
             (16*Li2(1 - x) - 24*log(x) - 
               112*log(1 - x)*log(x) + 24*pow(log(x),2)))/
           (1 - x)) + ((1 + pow(x,2))*
          (-8*Li3(1 - x) + 24*Li2(1 - x)*log(1 - x) + 
            56*log(x) + 64*zeta2*log(x) - 
            24*Li2(1 - x)*log(x) - 
            124*pow(log(1 - x),2)*log(x) + 
            72*log(1 - x)*pow(log(x),2) - 
            12*pow(log(x),3) - 32*S12(1 - x)))/(1 - x) + 
       (1 + x)*(-128*zeta3 - 40*Li3(1 - x) + 
          64*zeta2*log(1 - x) + 48*Li2(1 - x)*log(1 - x) - 
          64*pow(log(1 - x),3) - 32*zeta2*log(x) + 
          32*pow(log(1 - x),2)*log(x) - 
          24*log(1 - x)*pow(log(x),2) + 
          (14*pow(log(x),3))/3. + 16*S12(1 - x))))/
   (16.*pow(M_PI,2));
}

double NNLO_qqbar_HCA(double x){
	return (pow(alphas_Q,2)*CA*CF*(-16.51851851851852 + (2278*x)/27. - (4*(19 + 25*x)*zeta2)/3. - 16*x*Li2(1 - x) - (4*(7 + x)*Li2(1 - x))/3. + (22*(1 + x)*pow(log(Q2/muF2),2))/3. - (4*(38 + 239*x)*log(1 - x))/9. - (2*(26 - 57*x)*log(x))/3. - 4*(3 - x)*log(1 - x)*log(x) - 16*x*log(1 - x)*log(x) + ((23 - 25*x)*pow(log(x),2))/6. + 8*x*pow(log(x),2) + log(Q2/muF2)*((-4*(19 + 124*x))/9. + (1 + x)*(8*zeta2 + (88*log(1 - x))/3. - 6*log(x)) + ((1 + pow(x,2))*(-8*Li2(1 - x) + (70*log(x))/3.))/(1 - x)) + ((1 + pow(x,2))*((4*Li2(1 - x))/3. - 12*Li3(1 - x) - 8*Li2(1 - x)*log(1 - x) - (104*log(x))/3. + 8*zeta2*log(x) + 8*Li2(1 - x)*log(x) + (140*log(1 - x)*log(x))/3. - (29*pow(log(x),2))/2. - 4*S12(1 - x)))/(1 - x) + (1 + x)*(-28*zeta3 - 12*Li3(1 - x) + 16*zeta2*log(1 - x) + 8*Li2(1 - x)*log(1 - x) + (88*pow(log(1 - x),2))/3. + 16*S12(1 - x))))/(16.*pow(M_PI,2));
}
   


double NNLO_qqbar_A2(double x){
	return (pow(alphas_Q,2)*CF*nF*((4*(47 - 103*x))/27. - 
       (4*(1 + x)*pow(log(Q2/muF2),2))/3. - 
       (16*(1 - 11*x)*log(1 - x))/9. + 
       (8*(2 - 3*x)*log(x))/3. + 
       log(Q2/muF2)*((-8*(1 - 11*x))/9. - 
          (16*(1 + x)*log(1 - x))/3. - 
          (16*(1 + pow(x,2))*log(x))/(3.*(1 - x))) + 
       (1 + x)*((8*pow(M_PI,2))/9. + (4*Li2(1 - x))/3. - 
          (16*pow(log(1 - x),2))/3. + (2*pow(log(x),2))/3.
          ) + ((1 + pow(x,2))*
          ((-4*Li2(1 - x))/3. + (20*log(x))/3. - 
            (32*log(1 - x)*log(x))/3. + 4*pow(log(x),2)))/
        (1 - x)))/(16.*pow(M_PI,2));

}


double NNLO_qqbar_AC(double x){
	return (pow(alphas_Q,2)*CF*(-CA/2. + CF)*
     (2*(47 - 39*x) - 8*(8 - 7*x)*log(1 - x) + 
       2*(22 - 9*x)*log(x) + 
       (1 + x)*(-26*Li2(1 - x) - 8*Li3(1 - x) + 
          4*Li2(1 - x)*log(x) - 28*log(1 - x)*log(x) + 
          (23*pow(log(x),2))/2. + (2*pow(log(x),3))/3.) + 
       log(Q2/muF2)*(-4*(8 - 7*x) - 14*(1 + x)*log(x) + 
          ((1 + pow(x,2))*
             (-8*Li2(1 - x) - 6*log(x) - 4*pow(log(x),2)))/
           (1 - x)) + ((1 + pow(x,2))*
          (-6*Li2(1 - x) + 16*Li3(1 - x) - 
            16*Li2(1 - x)*log(1 - x) + 12*log(x) - 
            12*Li2(1 - x)*log(x) - 12*log(1 - x)*log(x) + 
            (15*pow(log(x),2))/2. - 
            8*log(1 - x)*pow(log(x),2) + 
            (8*pow(log(x),3))/3. - 36*S12(1 - x)))/(1 - x)))
    /(16.*pow(M_PI,2));
}


double NNLO_qqbar_B2(double x){
	return (pow(alphas_Q,2)*CF*((40*(1 - pow(x,2)))/3. + 
       (8*(3 + 4*x + 3*pow(x,2))*log(x))/3. + 
       pow(1 + x,2)*((-8*pow(M_PI,2))/9. - 
          (32*Li2(-x))/3. + (8*pow(log(x),2))/3. - 
          (32*log(x)*log(1 + x))/3.)))/(16.*pow(M_PI,2));
}

double NNLO_qqbar_BC(double x){
	return (pow(alphas_Q,2)*CF*(-CA/2. + CF)*
     (-2*(-27 + 14*x + 13*pow(x,2)) + 
       36*(1 - pow(x,2))*Li2(1 - x) + 4*(9 + 11*x)*log(x) - 
       2*(-6 + 8*x + 15*pow(x,2))*pow(log(x),2) + 
       (4*(1 + 4*x + pow(x,2))*pow(log(x),3))/3. + 
       (1 + 3*x + pow(x,2))*
        (16*Li2(1 - x)*log(x) + 32*S12(1 - x)) + 
       pow(1 + x,2)*(2*pow(M_PI,2) + 24*Li2(-x) - 
          8*Li3(-x) + (4*pow(M_PI,2)*log(x))/3. + 
          24*Li2(-x)*log(x) - 4*pow(M_PI,2)*log(1 + x) - 
          48*Li2(-x)*log(1 + x) + 24*log(x)*log(1 + x) + 
          20*pow(log(x),2)*log(1 + x) - 
          24*log(x)*pow(log(1 + x),2) - 48*S12(-x))))/
   (16.*pow(M_PI,2));
}

double NNLO_qqbar_b0_nonconst(double x){
	return -(pow(alphas_Q,2)*CF*(11*CA - 2*nF)*(1 + pow(x,2))*
      log(muR2/muF2)*(log(Q2/muF2) + 2*log(1 - x) - log(x)))/
   (12.*pow(M_PI,2)*(-1 + x));
}


double NNLO_qqbar_b0_const(){
	return (pow(alphas_Q,2)*CF*(11*CA - 2*nF)*
     log(muR2/muF2)*(-24 + 2*pow(M_PI,2) + 9*log(Q2/muF2)))/
   (72.*pow(M_PI,2));
}

//note the extra minus sign because it is LP!
double vegas_NNLO_LP_int(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0]; // from tau to 1
	double x = k[1]; // from 0 to 1
	return -LO_factor()*(-(pow(alphas_Q,2)*CF*log(1 - z)*(9*pow(log(Q2/muF2),2)*(-11*CA + 36*CF + 2*nF + 24*CF*log(1 - z)) - 
        6*log(Q2/muF2)*(-67*CA + 144*CF + 10*nF + 18*CA*zeta2 + 72*CF*zeta2 + (33*CA - 6*(9*CF + nF))*log(1 - z) - 72*CF*pow(log(1 - z),2)) + 
        2*(-202*CA + 28*nF + 198*CA*zeta2 - 36*nF*zeta2 + 189*CA*zeta3 + 864*CF*zeta3 - 
           3*(CA*(-67 + 18*zeta2) + 2*(5*nF + 36*CF*(2 + zeta2)))*log(1 - z) + (-66*CA + 12*nF)*pow(log(1 - z),2) + 108*CF*pow(log(1 - z),3)
           )))/(108.*pow(M_PI,2)))*derivative_qq_pdf_jac(x, z, tau); // jacobian needs to be in the derivative
}

///////////////////////////////////////
///
/// the qg channel
///
///////////////////////////////////////

///////integration routines
double vegas_NNLO_qg_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qg_full(z))*(jacobian*pdf_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qg_NLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qg_NLP(z))*(jacobian*pdf_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qg_NNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qg_NNLP(z))*(jacobian*pdf_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qg_NNNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qg_NNNLP(z))*(jacobian*pdf_sum_qg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}

//////// NNLO functions
double NNLO_qg_NLP(double x){
	//return (pow(alphas_Q,2)*TF*(6*(-2*CA + CF*(-5 + pow(M_PI,2))) + 3*pow(log(Q2/muF2),2)*(9*CF + 4*(CA + 3*CF)*log(1 - x)) + 2*log(Q2/muF2)*(6*CA - 33*CF - (3*CA + 2*CF)*pow(M_PI,2) + (11*CA - 2*nF)*log(muR2/muF2) + 18*log(1 - x)*(CF + (CA + 3*CF)*log(1 - x))) + 2*log(1 - x)*(9*(CA - 7*CF) - 2*(3*CA + 2*CF)*pow(M_PI,2) + (22*CA - 4*nF)*log(muR2/muF2) + log(1 - x)*(-18*CF + (13*CA + 35*CF)*log(1 - x))) + 6*(CA + 38*CF)*zeta3))/(48.*pow(M_PI,2));
	return (pow(alphas_Q,2)*TF*(-12*CA + 6*CF*(-5 + pow(M_PI,2)) + 3*pow(log(Q2/muF2),2)*(9*CF + 4*(CA + 3*CF)*log(1 - x)) + 2*log(1 - x)*(9*(CA - 7*CF) - 2*(3*CA + 2*CF)*pow(M_PI,2) + 6*b0*log(muR2/muF2) - 18*CF*log(1 - x) + (13*CA + 35*CF)*pow(log(1 - x),2)) + 2*log(Q2/muF2)*(6*CA - 33*CF - (3*CA + 2*CF)*pow(M_PI,2) + 3*b0*log(muR2/muF2) + 18*log(1 - x)*(CF + (CA + 3*CF)*log(1 - x))) + 6*(CA + 38*CF)*zeta3))/(48.*pow(M_PI,2));
}
double NNLO_qg_NNLP(double x){
	//return -(pow(alphas_Q,2)*TF*(-1 + x)*(189*CA + 267*CF - 2*(8*CA + 5*CF)*pow(M_PI,2) - (11*CA - 2*nF)*log(muR2/muF2)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)) + 3*pow(log(Q2/muF2),2)*(4*CA + 3*CF - 4*(CA + 3*CF)*log(1 - x)) + log(Q2/muF2)*(CF*(-69 + 4*pow(M_PI,2)) + 6*CA*(pow(M_PI,2) - 2*(6 + log(32))) + 12*log(1 - x)*(12*CA + 17*CF - 3*(CA + 3*CF)*log(1 - x))) + log(1 - x)*(-153*CA - 303*CF + 4*(3*CA + 2*CF)*pow(M_PI,2) + log(1 - x)*(96*CA + 273*CF - 2*(13*CA + 35*CF)*log(1 - x))) - 6*(CA + 38*CF)*zeta3))/(24.*pow(M_PI,2));
	return -(pow(alphas_Q,2)*TF*(-1 + x)*(189*CA + 231*CF - 2*(8*CA + 5*CF)*pow(M_PI,2) - 3*b0*log(muR2/muF2)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)) + 3*pow(log(Q2/muF2),2)*(4*CA + 3*CF - 4*(CA + 3*CF)*log(1 - x)) + log(Q2/muF2)*(6*CA*(-12 + pow(M_PI,2)) + CF*(-69 + 4*pow(M_PI,2)) + 12*log(1 - x)*(7*CA + 17*CF - 3*(CA + 3*CF)*log(1 - x))) + log(1 - x)*(-3*(51*CA + 77*CF) + 4*(3*CA + 2*CF)*pow(M_PI,2) + log(1 - x)*(96*CA + 273*CF - 2*(13*CA + 35*CF)*log(1 - x))) - 6*(CA + 38*CF)*zeta3))/(24.*pow(M_PI,2));
}
double NNLO_qg_NNNLP(double x){
	//return -(pow(alphas_Q,2)*TF*pow(-1 + x,2)*(1359*CA + 3366*CF - 2*(55*CA + 36*CF)*pow(M_PI,2) - 4*(11*CA - 2*nF)*log(muR2/muF2)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)) + 6*pow(log(Q2/muF2),2)*(3*(4*CA + 9*CF) - 8*(CA + 3*CF)*log(1 - x)) + 2*log(Q2/muF2)*(CF*(-669 + 8*pow(M_PI,2)) + 12*CA*(-26 + pow(M_PI,2) - 7*log(2)) + 12*log(1 - x)*(25*CA + 56*CF - 6*(CA + 3*CF)*log(1 - x))) + 2*log(1 - x)*(-612*CA - 1713*CF + 8*(3*CA + 2*CF)*pow(M_PI,2) + log(1 - x)*(258*CA + 783*CF - 4*(13*CA + 35*CF)*log(1 - x))) - 24*(CA + 38*CF)*zeta3))/(96.*pow(M_PI,2));
	return -(pow(alphas_Q,2)*TF*pow(-1 + x,2)*(1359*CA + 2886*CF - 2*(55*CA + 36*CF)*pow(M_PI,2) - 12*b0*log(muR2/muF2)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)) + 6*pow(log(Q2/muF2),2)*(3*(4*CA + 9*CF) - 8*(CA + 3*CF)*log(1 - x)) + 2*log(Q2/muF2)*(12*CA*(-21 + pow(M_PI,2)) + CF*(-669 + 8*pow(M_PI,2)) + 24*log(1 - x)*(9*CA + 28*CF - 3*(CA + 3*CF)*log(1 - x))) + 2*log(1 - x)*(-612*CA - 1713*CF + 8*(3*CA + 2*CF)*pow(M_PI,2) + log(1 - x)*(258*CA + 783*CF - 4*(13*CA + 35*CF)*log(1 - x))) - 24*(CA + 38*CF)*zeta3))/(96.*pow(M_PI,2));
}
double NNLO_qg_full(double x){
	//return (pow(alphas_Q,2)*TF*((11*CA - 2*nF)*log(muR2/muF2)*(1 + 6*x - 7*pow(x,2) + 2*(1 + 2*(-1 + x)*x)*(log(Q2/muF2) + 2*log(1 - x) - log(x))) + 3*CA*(59.888888888888886 + 116/(27.*x) - (1226*x)/9. + (1837*pow(x,2))/27. + (2*pow(M_PI,2)*(15 - 8/x - 84*x + 107*pow(x,2)))/9. + 8*x*Li2(1 - x) + (4*(33 + 16/x + 90*x + 44*pow(x,2))*Li2(1 - x))/3. - 4*(15 + 34*x + 12*pow(x,2))*Li3(1 - x) - (16*pow(M_PI,2)*(1 - x + 2*pow(x,2))*log(1 - x))/3. + (2*(-210 + 88/x + 75*x + 74*pow(x,2))*log(1 - x))/9. + 8*(7 + 10*x + 5*pow(x,2))*Li2(1 - x)*log(1 - x) + (4*(6 + 8/x + 63*x - 77*pow(x,2))*pow(log(1 - x),2))/3. + (26*(1 - 2*x + 2*pow(x,2))*pow(log(1 - x),3))/3. - 8*x*log(x) + (8*pow(M_PI,2)*x*(-5 + 2*x)*log(x))/3. - (2*(-354 + 12*x + 457*pow(x,2))*log(x))/9. + 8*(7 - 2*x)*x*Li2(1 - x)*log(x) + 8*x*log(1 - x)*log(x) + 20*(1 - 2*x + 13*pow(x,2))*log(1 - x)*log(x) + 4*(1 + 22*x - 6*pow(x,2))*pow(log(1 - x),2)*log(x) - 4*x*pow(log(x),2) - (5 + (346*pow(x,2))/3.)*pow(log(x),2) + 4*(-3 - 14*x + 2*pow(x,2))*log(1 - x)*pow(log(x),2) + (2*(9 + 20*x)*pow(log(x),3))/3. + (2*pow(log(Q2/muF2),2)*(3 + 4/x + 24*x - 31*pow(x,2) + 6*(1 - 2*x + 2*pow(x,2))*log(1 - x) + 6*(1 + 4*x)*log(x)))/3. + (2*log(Q2/muF2)*(-87 + 44/x - 12*x + 73*pow(x,2) - 12*pow(M_PI,2)*(1 - x + 2*pow(x,2)) + 36*(3 + 6*x + 2*pow(x,2))*Li2(1 - x) + 6*(9 + 8/x + 54*x - 71*pow(x,2))*log(1 - x) + 54*(1 - 2*x + 2*pow(x,2))*pow(log(1 - x),2) + 18*(3 - 2*x + 28*pow(x,2))*log(x) + 36*(1 + 10*x - 2*pow(x,2))*log(1 - x)*log(x) - 36*(1 + 3*x)*pow(log(x),2) - 36*(1 + 2*x + 2*pow(x,2))*(Li2(-x) + log(1 - x)*log(x))))/9. + 8*(1 + 5*x + 4*pow(x,2))*(Li2(-x) + log(x)*log(1 + x)) - 4*(1 + 2*x + 2*pow(x,2))*(2*Li3(-x) - 4*Li3((1 - x)/(1 + x)) + 4*Li3((-1 + x)/(1 + x)) + 4*Li2(-x)*log(1 - x) - 4*Li2(-x)*log(x) + 4*log(1 - x)*log(x)*log(1 + x) - 3*pow(log(x),2)*log(1 + x)) + 8*(9 + 16*x + 4*pow(x,2))*S12(1 - x) - 4*(1 + 4*x + 2*pow(x,2))*zeta3) + 3*CF*(-90.5 - 12*(-1 + x) + 233*x - (305*pow(x,2))/2. + (2*pow(M_PI,2)*(5 + 2*x - 12*pow(x,2)))/3. + 8*(-3 + x)*Li2(1 - x) + 2*(3 - 28*x + 40*pow(x,2))*Li2(1 - x) + 4*(-1 + 2*x + 18*pow(x,2))*Li3(1 - x) + 24*(-1 + x)*log(1 - x) + 2*(38 - 147*x + 88*pow(x,2))*log(1 - x) - 4*(3 - 6*x + 26*pow(x,2))*Li2(1 - x)*log(1 - x) - 2*(23 - 80*x + 63*pow(x,2))*pow(log(1 - x),2) + (28 - 44*x)*log(x) + 4*pow(M_PI,2)*(1 - 2*x + 4*pow(x,2))*log(x) - (59 - 245*x + 174*pow(x,2))*log(x) + 4*(1 - 2*x)*Li2(1 - x)*log(x) + 8*(-3 + x)*log(1 - x)*log(x) + 4*(13 - 50*x + 48*pow(x,2))*log(1 - x)*log(x) - 6*(7 - 14*x + 22*pow(x,2))*pow(log(1 - x),2)*log(x) - 4*(-3 + x)*pow(log(x),2) - ((35 - 68*x + 4*pow(x,2))*pow(log(x),2))/2. + 8*(3 - 6*x + 10*pow(x,2))*log(1 - x)*pow(log(x),2) - ((17 - 34*x + 52*pow(x,2))*pow(log(x),3))/3. + 3*pow(log(Q2/muF2),2)*(-1 + 4*x + (4 - 8*x + 8*pow(x,2))*log(1 - x) + (-2 + 4*x - 8*pow(x,2))*log(x)) + log(Q2/muF2)*(24 - 68*x + 22*pow(x,2) - 48*pow(x,2)*Li2(1 - x) - 4*(8 - 34*x + 23*pow(x,2))*log(1 - x) - (4*(1 - 2*x + 2*pow(x,2))*(pow(M_PI,2) - 27*pow(log(1 - x),2)))/3. + 2*(5 - 40*x + 46*pow(x,2))*log(x) - 8*(5 - 10*x + 16*pow(x,2))*log(1 - x)*log(x) + 8*(1 - 2*x + 4*pow(x,2))*pow(log(x),2)) - 16*(-1 + 2*x + 3*pow(x,2))*(Li2(-x) + log(x)*log(1 + x)) - 4*(11 - 22*x + 34*pow(x,2))*S12(1 - x) + (2*(1 - 2*x + 2*pow(x,2))*(48*Li3(-x) - 4*pow(M_PI,2)*log(1 - x) + 35*pow(log(1 - x),3) - 24*Li2(-x)*log(x) + 150*zeta3))/3.)))/(48.*pow(M_PI,2));
	return (pow(alphas_Q,2)*TF*(b0*log(muR2/muF2)*(1 + 6*x - 7*pow(x,2) + 2*(1 + 2*(-1 + x)*x)*(log(Q2/muF2) + 2*log(1 - x) - log(x))) + CA*(59.888888888888886 + 116/(27.*x) - (1226*x)/9. + (1837*pow(x,2))/27. + (2*pow(M_PI,2)*(15 - 8/x - 84*x + 107*pow(x,2)))/9. + 8*x*Li2(1 - x) + (4*(33 + 16/x + 90*x + 44*pow(x,2))*Li2(1 - x))/3. - 4*(15 + 34*x + 12*pow(x,2))*Li3(1 - x) - (16*pow(M_PI,2)*(1 - x + 2*pow(x,2))*log(1 - x))/3. + (2*(-210 + 88/x + 75*x + 74*pow(x,2))*log(1 - x))/9. + 8*(7 + 10*x + 5*pow(x,2))*Li2(1 - x)*log(1 - x) + (4*(6 + 8/x + 63*x - 77*pow(x,2))*pow(log(1 - x),2))/3. + (26*(1 - 2*x + 2*pow(x,2))*pow(log(1 - x),3))/3. - 8*x*log(x) + (8*pow(M_PI,2)*x*(-5 + 2*x)*log(x))/3. - (2*(-354 + 12*x + 457*pow(x,2))*log(x))/9. + 8*(7 - 2*x)*x*Li2(1 - x)*log(x) + 8*x*log(1 - x)*log(x) + 20*(1 - 2*x + 13*pow(x,2))*log(1 - x)*log(x) + 4*(1 + 22*x - 6*pow(x,2))*pow(log(1 - x),2)*log(x) - 4*x*pow(log(x),2) - (5 + (346*pow(x,2))/3.)*pow(log(x),2) + 4*(-3 - 14*x + 2*pow(x,2))*log(1 - x)*pow(log(x),2) + (2*(9 + 20*x)*pow(log(x),3))/3. + (2*pow(log(Q2/muF2),2)*(3 + 4/x + 24*x - 31*pow(x,2) + 6*(1 - 2*x + 2*pow(x,2))*log(1 - x) + 6*(1 + 4*x)*log(x)))/3. + 8*(1 + 5*x + 4*pow(x,2))*(Li2(-x) + log(x)*log(1 + x)) - 4*(1 + 2*x + 2*pow(x,2))*(2*Li3(-x) - 4*Li3((1 - x)/(1 + x)) + 4*Li3((-1 + x)/(1 + x)) + 4*Li2(-x)*log(1 - x) - 4*Li2(-x)*log(x) + 4*log(1 - x)*log(x)*log(1 + x) - 3*pow(log(x),2)*log(1 + x)) + (2*log(Q2/muF2)*(-87 + 44/x - 12*x + 73*pow(x,2) - 12*pow(M_PI,2)*(1 - x + 2*pow(x,2)) + 36*(3 + 6*x + 2*pow(x,2))*Li2(1 - x) + 6*(9 + 8/x + 54*x - 71*pow(x,2))*log(1 - x) + 54*(1 - 2*x + 2*pow(x,2))*pow(log(1 - x),2) + 18*(3 - 2*x + 28*pow(x,2))*log(x) + 36*(1 + 10*x - 2*pow(x,2))*log(1 - x)*log(x) - 36*(1 + 3*x)*pow(log(x),2) - 36*(1 + 2*x + 2*pow(x,2))*(Li2(-x) + log(x)*log(1 + x))))/9. + 8*(9 + 16*x + 4*pow(x,2))*S12(1 - x) - 4*(1 + 4*x + 2*pow(x,2))*zeta3) + CF*(-90.5 + 12*(-1 + x) + 233*x - (305*pow(x,2))/2. + (2*pow(M_PI,2)*(5 + 2*x - 12*pow(x,2)))/3. - 8*(-3 + x)*Li2(1 - x) + 2*(3 - 28*x + 40*pow(x,2))*Li2(1 - x) + 4*(-1 + 2*x + 18*pow(x,2))*Li3(1 - x) - 24*(-1 + x)*log(1 - x) + 2*(38 - 147*x + 88*pow(x,2))*log(1 - x) - 4*(3 - 6*x + 26*pow(x,2))*Li2(1 - x)*log(1 - x) - 2*(23 - 80*x + 63*pow(x,2))*pow(log(1 - x),2) + 4*(-7 + 11*x)*log(x) + 4*pow(M_PI,2)*(1 - 2*x + 4*pow(x,2))*log(x) - (59 - 245*x + 174*pow(x,2))*log(x) + 4*(1 - 2*x)*Li2(1 - x)*log(x) + 8*(-3 + x)*log(1 - x)*log(x) + 4*(13 - 50*x + 48*pow(x,2))*log(1 - x)*log(x) - 6*(7 - 14*x + 22*pow(x,2))*pow(log(1 - x),2)*log(x) + 4*(-3 + x)*pow(log(x),2) - ((35 - 68*x + 4*pow(x,2))*pow(log(x),2))/2. + 8*(3 - 6*x + 10*pow(x,2))*log(1 - x)*pow(log(x),2) - ((17 - 34*x + 52*pow(x,2))*pow(log(x),3))/3. + 3*pow(log(Q2/muF2),2)*(-1 + 4*x + (4 - 8*x + 8*pow(x,2))*log(1 - x) + (-2 + 4*x - 8*pow(x,2))*log(x)) + log(Q2/muF2)*(24 - 68*x + 22*pow(x,2) - 48*pow(x,2)*Li2(1 - x) - 4*(8 - 34*x + 23*pow(x,2))*log(1 - x) - (4*(1 - 2*x + 2*pow(x,2))*(pow(M_PI,2) - 27*pow(log(1 - x),2)))/3. + 2*(5 - 40*x + 46*pow(x,2))*log(x) - 8*(5 - 10*x + 16*pow(x,2))*log(1 - x)*log(x) + 8*(1 - 2*x + 4*pow(x,2))*pow(log(x),2)) - 16*(-1 + 2*x + 3*pow(x,2))*(Li2(-x) + log(x)*log(1 + x)) - 4*(11 - 22*x + 34*pow(x,2))*S12(1 - x) + (2*(1 - 2*x + 2*pow(x,2))*(48*Li3(-x) - 4*pow(M_PI,2)*log(1 - x) + 35*pow(log(1 - x),3) - 24*Li2(-x)*log(x) + 150*zeta3))/3.)))/(16.*pow(M_PI,2));
}


///////////////
///
/// gg channel
///
///////////////

/////////// integration routines
double vegas_NNLO_gg_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_gg_full(z))*(jacobian*pdf_sum_gg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_gg_NLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_gg_NLP(z))*(jacobian*pdf_sum_gg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_gg_NNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_gg_NNLP(z))*(jacobian*pdf_sum_gg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_gg_NNNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_gg_NNNLP(z))*(jacobian*pdf_sum_gg_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}

//////// NNLO functions
double NNLO_gg_full(double x){
	return (pow(alphas_Q,2)*(-32 - 66*x + 98*pow(x,2) + 4*(5 + 9*x - 12*pow(x,2))*zeta2 + 4*(-1 - 2*x + 2*pow(x,2))*zeta3 + 4*(-5 - 4*x + 14*pow(x,2))*Li2(1 - x) + 8*(-1 - 2*x + pow(x,2))*Li3(-x) + 2*(7 + 60*x - 67*pow(x,2))*log(1 - x) + 16*(-1 - 2*x + 3*pow(x,2))*pow(log(1 - x),2) + (-23 - 64*x + 105*pow(x,2))*log(x) + 4*(3 + 10*x + 10*pow(x,2))*zeta2*log(x) + 8*(2 + 4*x + pow(x,2))*Li2(-x)*log(x) + 4*(1 + 8*x - 4*pow(x,2))*log(1 - x)*log(x) - 2*(3 + 7*x + 4*pow(x,2))*pow(log(x),2) - (2*(3 + 8*x + 8*pow(x,2))*pow(log(x),3))/3. - 2*pow(log(Q2/muF2),2)*(2 + 4*x - 6*pow(x,2) + pow(1 + 2*x,2)*log(x)) + 4*pow(1 + 2*x,2)*(4*Li3(1 - x) + log(1 - x)*log(x)*(-2*log(1 - x) + log(x)) - Li2(1 - x)*(4*log(1 - x) + log(x))) + log(Q2/muF2)*(7 + 60*x - 67*pow(x,2) + 16*(-1 - 2*x + 3*pow(x,2))*log(1 - x) + 2*(1 + 8*x - 4*pow(x,2))*log(x) + pow(1 + 2*x,2)*(-8*Li2(1 - x) + 2*log(x)*(-4*log(1 - x) + log(x)))) + 8*(1 + x)*(Li2(-x) + log(x)*log(1 + x)) - 8*(1 + 10*x + 7*pow(x,2))*S12(1 - x) - 4*pow(1 + x,2)*(4*Li2(-x)*log(1 + x) + log(1 + x)*(2*zeta2 - 3*pow(log(x),2) + 2*log(x)*log(1 + x)) + 4*S12(-x)) + (pow(CA,2)*(-47 - 144*x + 191*pow(x,2) - 2*(6 + 38*x + 75*pow(x,2))*log(x) + 2*(-2 + 2*x + 25*pow(x,2))*pow(log(x),2) - 24*pow(-1 + x,2)*S12(1 - x) + 4*pow(1 + x,2)*(2*zeta2 + 12*zeta3 + 18*Li3(-x) + 6*zeta2*log(1 + x) + 4*log(x)*log(1 + x) - 9*pow(log(x),2)*log(1 + x) + 6*log(x)*pow(log(1 + x),2) + Li2(-x)*(4 - 18*log(x) + 12*log(1 + x)) + 12*S12(-x))))/(3.*(-1 + pow(CA,2)))))/(16.*pow(M_PI,2));
}
double NNLO_gg_NLP(double x){
	return 0;
}
double NNLO_gg_NNLP(double x){
	return (pow(alphas_Q,2)*(1 - x)*(3*pow(log(Q2/muF2),2) + 12*log(Q2/muF2)*(-1 + log(1 - x)) - 2*(-12 + pow(M_PI,2) - 6*(-2 + log(1 - x))*log(1 - x))))/(24.*pow(M_PI,2));
}
double NNLO_gg_NNNLP(double x){
	return (pow(alphas_Q - alphas_Q*x,2)*(-3*pow(log(Q2/muF2),2) - 12*log(Q2/muF2)*(-2 + log(1 - x)) + 2*(-27 + pow(M_PI,2) - 6*(-4 + log(1 - x))*log(1 - x))))/(16.*pow(M_PI,2));
}

////////////////////////////////////////////////////////////
///
/// qq channel and qbarqbar channel (identical quarks)
///
////////////////////////////////////////////////////////////

/////////// integration routines (qq)
double vegas_NNLO_qq_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qq_full(z))*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qq_NLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qq_NLP(z))*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qq_NNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qq_NNLP(z))*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qq_NNNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qq_NNNLP(z))*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}

/////////// integration routines (qbarqbar)
double vegas_NNLO_qbarqbar_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qq_full(z))*(jacobian*pdf_sum_qbarqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qbarqbar_NLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qq_NLP(z))*(jacobian*pdf_sum_qbarqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qbarqbar_NNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qq_NNLP(z))*(jacobian*pdf_sum_qbarqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qbarqbar_NNNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qq_NNNLP(z))*(jacobian*pdf_sum_qbarqbar_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}

//////// NNLO functions
double NNLO_qq_full(double x){
	return (pow(alphas_Q,2)*CF*((2*TF*(1779 + 116/x - 2598*x + 703*pow(x,2) + 36*(39 + 16/x + 15*x + 8*pow(x,2))*Li2(1 - x) - 24*(39 - 22/x - 39*x + 22*pow(x,2))*log(1 - x) + (72*(-4 - 3*x + 3*pow(x,2) + 4*pow(x,3))*(zeta2 - pow(log(1 - x),2)))/x + 6*(345 - 48*x + 20*pow(x,2))*log(x) + 216*(3 + 6*x + 4*pow(x,2))*log(1 - x)*log(x) - 45*(3 + 15*x + 8*pow(x,2))*pow(log(x),2) + 18*pow(log(Q2/muF2),2)*(3 + 4/x - 3*x - 4*pow(x,2) + 6*(1 + x)*log(x)) + 12*log(Q2/muF2)*(-39 + 22/x + 39*x - 22*pow(x,2) + 6*(3 + 4/x - 3*x - 4*pow(x,2))*log(1 - x) + 9*(3 + 6*x + 4*pow(x,2))*log(x) + 18*(1 + x)*(2*Li2(1 - x) + (2*log(1 - x) - log(x))*log(x))) - 18*(1 + x)*(48*Li3(1 - x) + 4*pow(M_PI,2)*log(x) - 24*pow(log(1 - x),2)*log(x) + 24*log(1 - x)*pow(log(x),2) - 9*pow(log(x),3) - 12*Li2(1 - x)*(4*log(1 - x) + log(x)) - 72*S12(1 - x))))/27. + 2*(-CA/2. + CF)*(-30 + 56*x - 26*pow(x,2) + 4*(-7 + 6*x)*log(x) - (4*pow(-1 + x,2)*(-12*Li3(1 - x) + 9*pow(log(x),2) + 2*pow(log(x),3) + 6*Li2(1 - x)*(3 + 2*log(x)) + 12*S12(1 - x)))/3.) + TF*(160*(-1 + x) - 24*(-6 + 4/x + x)*zeta3 + 8*(5 - 4*x)*Li2(1 - x) + 8*(-6 + 4/x + 3*x)*Li3(1 - x) - 16*(-10 + 10/x + 3*x)*Li3(-x) - 16*(5 + 4*x)*log(x) + 8*(10 + x)*zeta2*log(x) - 8*(-10 + x)*Li2(1 - x)*log(x) + 32*(5/x + 2*x)*Li2(-x)*log(x) - 52*x*pow(log(x),2) - (16*x*pow(log(x),3))/3. + 40*(1 + x)*(zeta2 + 2*Li2(-x) + 2*log(x)*log(1 + x)) - (8*(2 + 2*x + pow(x,2))*(6*zeta2*log(1 + x) + 12*Li2(-x)*log(1 + x) - 5*pow(log(x),2)*log(1 + x) + 6*log(x)*pow(log(1 + x),2) - 4*S12(1 - x) + 12*S12(-x)))/x) + 4*(-CA/2. + CF)*(8*(3 + x)*Li2(1 - x) + 2*(-9 + 7*x)*log(x) - 4*(1 + 3*x)*pow(log(x),2) + 4*(1 + x)*(zeta2 + 2*Li2(-x) + 4*log(1 - x)*log(x) + 2*log(x)*log(1 + x)) + log(Q2/muF2)*(-16*(-1 + x) + 8*(1 + x)*log(x) - (4*(1 + pow(x,2))*(2*zeta2 + 4*Li2(-x) - pow(log(x),2) + 4*log(x)*log(1 + x)))/(1 + x)) - (4*(1 + pow(x,2))*(3*zeta3 + 24*Li3(1 - x) + 6*Li3(-x) - 24*Li3((1 - x)/(1 + x)) + 24*Li3((-1 + x)/(1 + x)) + 12*zeta2*log(1 - x) + 24*Li2(-x)*log(1 - x) - 9*zeta2*log(x) - 18*Li2(1 - x)*log(x) - 24*Li2(-x)*log(x) - 6*log(1 - x)*pow(log(x),2) + 2*pow(log(x),3) + 6*zeta2*log(1 + x) + 12*Li2(-x)*log(1 + x) + 24*log(1 - x)*log(x)*log(1 + x) - 21*pow(log(x),2)*log(1 + x) + 6*log(x)*pow(log(1 + x),2) - 24*S12(1 - x) + 12*S12(-x)))/(3.*(1 + x)) + (1 - x)*(8*Li3(-x) - (2*(51 - 12*zeta3 - 48*log(1 - x) - 6*zeta2*log(x) + pow(log(x),3) + 2*pow(M_PI,2)*log(1 + x) + 24*Li2(-x)*log(1 + x) - 6*pow(log(x),2)*log(1 + x) + 12*log(x)*pow(log(1 + x),2) + 24*S12(-x)))/3.))))/(32.*pow(M_PI,2));
}
double NNLO_qq_NLP(double x){
	return 0;
}
double NNLO_qq_NNLP(double x){
	return (pow(alphas_Q,2)*CF*(1 - x)*(-3*CA + 6*CF + (27 - 2*pow(M_PI,2))*TF + 3*TF*(-4 + log(Q2/muF2) + 2*log(1 - x))*(log(Q2/muF2) + 2*log(1 - x))))/(24.*pow(M_PI,2));
}
double NNLO_qq_NNNLP(double x){
	return (pow(alphas_Q,2)*CF*TF*pow(-1 + x,2)*(-17 + 12*log(Q2/muF2) + 24*log(1 - x)))/(32.*pow(M_PI,2));
}



////////////////////////////////////////////////////////////
///
/// qq channel and qbarqbar channel (non-identical quarks)
///
////////////////////////////////////////////////////////////

/////////// integration routines (qq)
double vegas_NNLO_qqNI_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqNI_full(z))*(jacobian*pdf_sum_qqNI_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qqNI_NLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqNI_NLP(z))*(jacobian*pdf_sum_qqNI_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qqNI_NNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqNI_NNLP(z))*(jacobian*pdf_sum_qqNI_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qqNI_NNNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqNI_NNNLP(z))*(jacobian*pdf_sum_qqNI_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}

/////////// integration routines (qbarqbar)
double vegas_NNLO_qbarqbarNI_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqNI_full(z))*(jacobian*pdf_sum_qbarqbarNI_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qbarqbarNI_NLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqNI_NLP(z))*(jacobian*pdf_sum_qbarqbarNI_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qbarqbarNI_NNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqNI_NNLP(z))*(jacobian*pdf_sum_qbarqbarNI_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qbarqbarNI_NNNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqNI_NNNLP(z))*(jacobian*pdf_sum_qbarqbarNI_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}

//////// NNLO functions
double NNLO_qqNI_full(double x){
	return (pow(alphas_Q,2)*CF*TF*(160*(-1 + x) - 24*(-6 + 4/x + x)*zeta3 + 8*(5 - 4*x)*Li2(1 - x) + 8*(-6 + 4/x + 3*x)*Li3(1 - x) - 16*(-10 + 10/x + 3*x)*Li3(-x) - 16*(5 + 4*x)*log(x) + 8*(10 + x)*zeta2*log(x) - 8*(-10 + x)*Li2(1 - x)*log(x) + 32*(5/x + 2*x)*Li2(-x)*log(x) - 52*x*pow(log(x),2) - (16*x*pow(log(x),3))/3. + 40*(1 + x)*(zeta2 + 2*Li2(-x) + 2*log(x)*log(1 + x)) + (2*(1779 + 116/x - 2598*x + 703*pow(x,2) + 36*(39 + 16/x + 15*x + 8*pow(x,2))*Li2(1 - x) - 24*(39 - 22/x - 39*x + 22*pow(x,2))*log(1 - x) + (72*(-4 - 3*x + 3*pow(x,2) + 4*pow(x,3))*(zeta2 - pow(log(1 - x),2)))/x + 6*(345 - 48*x + 20*pow(x,2))*log(x) + 216*(3 + 6*x + 4*pow(x,2))*log(1 - x)*log(x) - 45*(3 + 15*x + 8*pow(x,2))*pow(log(x),2) + 18*pow(log(Q2/muF2),2)*(3 + 4/x - 3*x - 4*pow(x,2) + 6*(1 + x)*log(x)) + 12*log(Q2/muF2)*(-39 + 22/x + 39*x - 22*pow(x,2) + 6*(3 + 4/x - 3*x - 4*pow(x,2))*log(1 - x) + 9*(3 + 6*x + 4*pow(x,2))*log(x) + 18*(1 + x)*(2*Li2(1 - x) + (2*log(1 - x) - log(x))*log(x))) - 18*(1 + x)*(48*Li3(1 - x) + 4*pow(M_PI,2)*log(x) - 24*pow(log(1 - x),2)*log(x) + 24*log(1 - x)*pow(log(x),2) - 9*pow(log(x),3) - 12*Li2(1 - x)*(4*log(1 - x) + log(x)) - 72*S12(1 - x))))/27. - (8*(2 + 2*x + pow(x,2))*(6*zeta2*log(1 + x) + 12*Li2(-x)*log(1 + x) - 5*pow(log(x),2)*log(1 + x) + 6*log(x)*pow(log(1 + x),2) - 4*S12(1 - x) + 12*S12(-x)))/x))/(16.*pow(M_PI,2));
}
double NNLO_qqNI_NLP(double x){
	return 0;
}
double NNLO_qqNI_NNLP(double x){
	return (pow(alphas_Q,2)*CF*TF*(1 - x)*(27 - 2*pow(M_PI,2) + 3*pow(log(Q2/muF2),2) + 12*log(Q2/muF2)*(-1 + log(1 - x)) + 12*(-2 + log(1 - x))*log(1 - x)))/(12.*pow(M_PI,2));
}
double NNLO_qqNI_NNNLP(double x){
	return (pow(alphas_Q,2)*CF*TF*pow(-1 + x,2)*(-17 + 12*log(Q2/muF2) + 24*log(1 - x)))/(16.*pow(M_PI,2));
}



////////////////////////////////////////////////////////////
///
/// qqbar channel (non-identical quarks)
///
////////////////////////////////////////////////////////////

/////////// integration routines (qqbar)
double vegas_NNLO_qqbarNI_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqbarNI_full(z))*(jacobian*pdf_sum_qqbarNI_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qqbarNI_NLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqbarNI_NLP(z))*(jacobian*pdf_sum_qqbarNI_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qqbarNI_NNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqbarNI_NNLP(z))*(jacobian*pdf_sum_qqbarNI_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_NNLO_qqbarNI_NNNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*(NNLO_qqbarNI_NNNLP(z))*(jacobian*pdf_sum_qqbarNI_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}

//////// NNLO functions
double NNLO_qqbarNI_full(double x){
	return (pow(alphas_Q,2)*CF*TF*(-160*(-1 + x) + 24*(-6 + 4/x + x)*zeta3 + 8*(-5 + 4*x)*Li2(1 - x) - 8*(-6 + 4/x + 3*x)*Li3(1 - x) + 16*(-10 + 10/x + 3*x)*Li3(-x) + 16*(5 + 4*x)*log(x) - 8*(10 + x)*zeta2*log(x) + 8*(-10 + x)*Li2(1 - x)*log(x) - 32*(5/x + 2*x)*Li2(-x)*log(x) + 52*x*pow(log(x),2) + (16*x*pow(log(x),3))/3. - 40*(1 + x)*(zeta2 + 2*Li2(-x) + 2*log(x)*log(1 + x)) + (2*(1779 + 116/x - 2598*x + 703*pow(x,2) + 36*(39 + 16/x + 15*x + 8*pow(x,2))*Li2(1 - x) - 24*(39 - 22/x - 39*x + 22*pow(x,2))*log(1 - x) + (72*(-4 - 3*x + 3*pow(x,2) + 4*pow(x,3))*(zeta2 - pow(log(1 - x),2)))/x + 6*(345 - 48*x + 20*pow(x,2))*log(x) + 216*(3 + 6*x + 4*pow(x,2))*log(1 - x)*log(x) - 45*(3 + 15*x + 8*pow(x,2))*pow(log(x),2) + 18*pow(log(Q2/muF2),2)*(3 + 4/x - 3*x - 4*pow(x,2) + 6*(1 + x)*log(x)) + 12*log(Q2/muF2)*(-39 + 22/x + 39*x - 22*pow(x,2) + 6*(3 + 4/x - 3*x - 4*pow(x,2))*log(1 - x) + 9*(3 + 6*x + 4*pow(x,2))*log(x) + 18*(1 + x)*(2*Li2(1 - x) + (2*log(1 - x) - log(x))*log(x))) - 18*(1 + x)*(48*Li3(1 - x) + 4*pow(M_PI,2)*log(x) - 24*pow(log(1 - x),2)*log(x) + 24*log(1 - x)*pow(log(x),2) - 9*pow(log(x),3) - 12*Li2(1 - x)*(4*log(1 - x) + log(x)) - 72*S12(1 - x))))/27. + (8*(2 + 2*x + pow(x,2))*(6*zeta2*log(1 + x) + 12*Li2(-x)*log(1 + x) - 5*pow(log(x),2)*log(1 + x) + 6*log(x)*pow(log(1 + x),2) - 4*S12(1 - x) + 12*S12(-x)))/x))/(16.*pow(M_PI,2));
}
double NNLO_qqbarNI_NLP(double x){
	return 0;
}
double NNLO_qqbarNI_NNLP(double x){
	return (pow(alphas_Q,2)*CF*TF*(1 - x)*(27 - 2*pow(M_PI,2) + 3*pow(log(Q2/muF2),2) + 12*log(Q2/muF2)*(-1 + log(1 - x)) + 12*(-2 + log(1 - x))*log(1 - x)))/(12.*pow(M_PI,2));
}
double NNLO_qqbarNI_NNNLP(double x){
	return (pow(alphas_Q,2)*CF*TF*pow(-1 + x,2)*(-31 + 12*log(Q2/muF2) + 24*log(1 - x)))/(16.*pow(M_PI,2));
}
