#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "deriv_pdf.h"
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
double vegas_NNLO_qqtogg_LP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*((pow(alphas_Q,2)*CF*(9*pow(log(Q2/muF2),2)*(11*CA - 2*(18*CF + nF) - 48*CF*log(1 - z)) + 
       6*log(Q2/muF2)*(-67*CA + 144*CF + 10*nF + 18*CA*zeta2 + 72*CF*zeta2 + 6*(11*CA - 2*(9*CF + nF))*log(1 - z) - 
          216*CF*pow(log(1 - z),2)) + 2*(202*CA - 28*nF - 198*CA*zeta2 + 36*nF*zeta2 - 189*CA*zeta3 - 864*CF*zeta3 + 
          6*(CA*(-67 + 18*zeta2) + 2*(5*nF + 36*CF*(2 + zeta2)))*log(1 - z) + 18*(11*CA - 2*nF)*pow(log(1 - z),2) - 432*CF*pow(log(1 - z),3)
          )))/(108.*pow(M_PI,2)*(-1 + z)))*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z - (1.-tau)*pdf_sum_qq_charge_weighted(tau+k[1]*(1.-tau),tau));
}
//////////////////////////////////////////////////////////////
/// correction term stemming from not integrating from 0 to 1
/// for NNLO qq to gg 
//////////////////////////////////////////////////////////////
double vegas_NNLO_qqtogg_LP_correction(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0]; // from 0 to tau
	double jacobian = 1.-tau; //needed to transform the boundary dependent terms
	double x = tau+k[1]*jacobian;
	return LO_factor()*((pow(alphas_Q,2)*CF*(9*pow(log(Q2/muF2),2)*(11*CA - 2*(18*CF + nF) - 48*CF*log(1 - z)) + 
       6*log(Q2/muF2)*(-67*CA + 144*CF + 10*nF + 18*CA*zeta2 + 72*CF*zeta2 + 6*(11*CA - 2*(9*CF + nF))*log(1 - z) - 
          216*CF*pow(log(1 - z),2)) + 2*(202*CA - 28*nF - 198*CA*zeta2 + 36*nF*zeta2 - 189*CA*zeta3 - 864*CF*zeta3 + 
          6*(CA*(-67 + 18*zeta2) + 2*(5*nF + 36*CF*(2 + zeta2)))*log(1 - z) + 18*(11*CA - 2*nF)*pow(log(1 - z),2) - 432*CF*pow(log(1 - z),3)
          )))/(108.*pow(M_PI,2)*(-1 + z)))*jacobian*(-pdf_sum_qq_charge_weighted(x,tau));
	
	//return log(1.-z)/(1.-z)*jacobian*(x*(z-tau) - x*(1.-tau)); // test function
}

////////////////////////////////////////
/// NLP contribution DY qq to gg nnlo
////////////////////////////////////////
double vegas_NNLO_qqtogg_NLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*((pow(alphas_Q,2)*CF*(-2*(27*CF*(19 + 16*zeta2 + 32*zeta3) + CA*(-391 + 252*zeta2 + 189*zeta3)) + 
       3*(9*CF*(16*(7 + zeta2) + z*(-33 + 16*zeta2)) + CA*(-194 + 36*zeta2 + z*(-329 + 36*zeta2)))*log(1 - z) + 
       18*(22*CA + 45*CF + 48*CF*z)*pow(log(1 - z),2) - 864*CF*pow(log(1 - z),3) + 
       9*pow(log(Q2/muF2),2)*(11*CA - 12*CF - 48*CF*log(1 - z)) + 
       6*log(Q2/muF2)*(9*CF*(27 + 8*zeta2) + 2*CA*(-71 + 9*zeta2) + 6*(11*CA + 24*CF)*log(1 - z) - 216*CF*pow(log(1 - z),2))))/
   (108.*pow(M_PI,2)))*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}


////////////////////////////////////////
/// NNLP contribution DY qq to gg nnlo
////////////////////////////////////////
double vegas_NNLO_qqtogg_NNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*((pow(alphas_Q,2)*CF*(-1 + z)*(2948*CA + 1863*CF - 144*CF*pow(M_PI,2) - 558*CA*zeta2 - 1728*CF*zeta2 - 378*CA*zeta3 - 1728*CF*zeta3 + 
       18*(27*CF*(1 - 4*z) + CA*(-53 + 9*z))*log(1 - z) + 18*(22*CA + 93*CF)*pow(log(1 - z),2) - 864*CF*pow(log(1 - z),3) + 
       9*pow(log(Q2/muF2),2)*(11*CA + 36*CF - 48*CF*log(1 - z)) + 
       3*log(Q2/muF2)*(144*CF*(-2 + zeta2) + CA*(-461 + 36*zeta2) + 12*(11*CA + 72*CF)*log(1 - z) - 432*CF*pow(log(1 - z),2))))/
   (216.*pow(M_PI,2)))*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}




////////////////////////////////////////
/// NNNLP contribution DY qq to gg nnlo
////////////////////////////////////////
double vegas_NNLO_qqtogg_NNNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*((pow(alphas_Q,2)*CF*(-1 + z)*(-(CA*(-1 + z)*(-1105 + 72*zeta2 + 288*log(Q2/muF2) + 714*log(1 - z))) + 
       CF*(-1639 + 144*pow(M_PI,2) + 1639*z - 288*zeta2 - 576*z*zeta2 + 144*(-1 + z)*pow(log(Q2/muF2),2) + 18*(-1 + z)*log(1 - z) + 
          1116*(-1 + z)*pow(log(1 - z),2) + 12*(-1 + z)*log(Q2/muF2)*(-59 + 84*log(1 - z)))))/(216.*pow(M_PI,2)))*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
}

////////////////////////////////////////////////////////////////
/// constant NNLO qq to gg  contribution 
/// (contains  alphas^2*delta(1-z))
/// already integrated over z
/// so make sure here no z integral 
////////////////////////////////////////////////////////////////
double vegas_NNLO_qqtogg_delta(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double jacobian = 1.;//-tau; //needed to transform the boundary dependent terms
	double x = /*tau+*/k[0]*jacobian;
	return LO_factor()*((pow(alphas_Q,2)*(CA*CF*(-127.91666666666667 + (592*zeta2)/9. - (12*pow(zeta2,2))/5. + 28*zeta3 + 
          (64.33333333333333 - 24*zeta3)*log(Q2/muF2) - 11*pow(log(Q2/muF2),2)) + 
       CF*nF*(21.166666666666668 - (112*zeta2)/9. + 8*zeta3 - (34*log(Q2/muF2))/3. + 2*pow(log(Q2/muF2),2)) + 
       pow(CF,2)*(127.75 - 70*zeta2 + (8*pow(zeta2,2))/5. - 60*zeta3 + (-93 + 24*zeta2 + 176*zeta3)*log(Q2/muF2) + 
          (18 - 32*zeta2)*pow(log(Q2/muF2),2))))/(16.*pow(M_PI,2)))*(jacobian*pdf_sum_qq_charge_weighted(/*tau+*/k[0]*jacobian,tau));
}




////////////////////////////////////////
/// Full contribution DY qq to gg nnlo
////////////////////////////////////////
double vegas_NNLO_qqtogg_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return LO_factor()*((pow(alphas_Q,2)*CF*(CF*(-72 + 48*z + 8*(3 - 2*z)*Li2(1 - z) - 16*(3 + z)*Li2(1 - z) + 4*(64 + 3*z)*log(1 - z) - 
          64*(-1 + z)*(zeta2 - pow(log(1 - z),2)) - 16*log(z) + 8*(-4 + 13*z)*log(z) + 16*(7 - 6*z)*log(1 - z)*log(z) - 
          16*(3 + z)*log(1 - z)*log(z) + 16*(-2 + z)*pow(log(z),2) + 8*(3 + z)*pow(log(z),2) + 
          pow(log(Q2/muF2),2)*(-8*(5 + z) + (16*(1 + pow(z,2))*log(z))/(-1 + z) + 8*(1 + z)*(-4*log(1 - z) + log(z))) + 
          (4*log(Q2/muF2)*(-30 + 26*z + 4*pow(z,2) - 8*zeta2 + 8*pow(z,2)*zeta2 + 4*(-3 + pow(z,2))*Li2(1 - z) - 
               24*(-1 + pow(z,2))*pow(log(1 - z),2) - 2*log(z) + 20*z*log(z) - 6*pow(z,2)*log(z) - 3*pow(log(z),2) - 
               9*pow(z,2)*pow(log(z),2) + 4*log(1 - z)*(7 - 8*z + pow(z,2) + (5 + 9*pow(z,2))*log(z))))/(-1 + z) + 
          (4*(1 + pow(z,2))*(2*Li3(1 - z) - 6*Li2(1 - z)*(log(1 - z) - log(z)) - 14*log(z) - 16*zeta2*log(z) + 
               31*pow(log(1 - z),2)*log(z) - 18*log(1 - z)*pow(log(z),2) + 3*pow(log(z),3) + 8*S12(1 - z)))/(-1 + z) + 
          (1 + z)*(-128*zeta3 - 40*Li3(1 - z) + 64*zeta2*log(1 - z) + 48*Li2(1 - z)*log(1 - z) - 64*pow(log(1 - z),3) - 32*zeta2*log(z) + 
             32*pow(log(1 - z),2)*log(z) - 24*log(1 - z)*pow(log(z),2) + (14*pow(log(z),3))/3. + 16*S12(1 - z))) + 
       CA*(-16.51851851851852 + (2278*z)/27. - (4*(19 + 25*z)*zeta2)/3. - 16*z*Li2(1 - z) - (4*(7 + z)*Li2(1 - z))/3. + 
          (22*(1 + z)*pow(log(Q2/muF2),2))/3. - (4*(38 + 239*z)*log(1 - z))/9. + 
          log(Q2/muF2)*((-4*(19 + 124*z))/9. + (2*(1 + pow(z,2))*(12*Li2(1 - z) - 35*log(z)))/(3.*(-1 + z)) + 
             (1 + z)*(8*zeta2 + (88*log(1 - z))/3. - 6*log(z))) + (2*(-26 + 57*z)*log(z))/3. + 4*(-3 + z)*log(1 - z)*log(z) - 
          16*z*log(1 - z)*log(z) + ((23 - 25*z)*pow(log(z),2))/6. + 8*z*pow(log(z),2) + 
          (1 + z)*(-28*zeta3 - 12*Li3(1 - z) + 16*zeta2*log(1 - z) + 8*Li2(1 - z)*log(1 - z) + (88*pow(log(1 - z),2))/3. + 16*S12(1 - z)) + 
          ((1 + pow(z,2))*(72*Li3(1 - z) + 8*Li2(1 - z)*(-1 + 6*log(1 - z) - 6*log(z)) + 208*log(z) - 48*zeta2*log(z) - 
               280*log(1 - z)*log(z) + 87*pow(log(z),2) + 24*S12(1 - z)))/(6.*(-1 + z)))))/(16.*pow(M_PI,2)))*(jacobian*pdf_sum_qq_charge_weighted(tau/z+k[1]*jacobian,tau/z)/z);
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
