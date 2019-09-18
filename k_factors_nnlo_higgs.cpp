#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "polygamma.h"
#include "gsl/gsl_sf_zeta.h"
#include "k_factors_higgs.h"
#include "k_factors_nnlo_higgs.h"
using namespace std;

//////////////////////////////////////////////////////////
///
/// contains all K factors for higgs
/// split up in NNLO and LP, NLP, NNLP, NNNLP, full
///
/// extracted log dependent parts from ihixs
///
//////////////////////////////////////////////////////////


///////////////////////////
/// gg channel
///////////////////////////
// the NNLO functions for the gg channel
// https://arxiv.org/pdf/hep-ph/0207004.pdf eqn. 47 and 48 (note that there is no scale dependence now) also we multiply with 1/x and take into account the extra 1/x = (1-x)/x + 1 term from the plus distribution
// checked with mathematica
double higgs_NNLO_gg_reg(double x){
	return 0.*logdep_gg(x)+pow(alphas_muR/M_PI,2)*(
	(1./x*((133.-90.*zeta2)*log(1.-x)+(-101./3.+33.*zeta2+351./2.*zeta3)-33.*pow(log(1.-x),2)+72.*pow(log(1.-x),3)) // this is the same
	+1./x*(9.*(38.*pow(x,2)-20.*pow(x,3)+18.*x-39.*pow(x,4)+14.+7.*pow(x,5))/(1.-pow(x,2))*Li3(x)-18.*pow((pow(x,2)+x+1),2)/(1.+x)*S12(pow(x,2))
			+9.*(4.*pow(x,4)+8.*pow(x,3)+21.*pow(x,2)+14.*x+7.)/(1.+x)*S12(-x)-9./2.*(5.*pow(x,5)-51.*pow(x,4)-57.*pow(x,3)+53.*pow(x,2)+59.*x-11.)/(1.-pow(x,2))*S12(x)
			-9./2.*(8.*pow(x,4)+8.*pow(x,3)-3.*pow(x,2)-2.*x-1.)/(1.+x)*Li3(-x)-9./2.*(16.+13.*pow(x,5)-40.*pow(x,3)-67.*pow(x,4)+64.*pow(x,2)+36.*x)/(1.-pow(x,2))*li2(x)*log(x)
			+9./2.*(2.*pow(x,4)-15.*pow(x,2)-10.*x-5.)/(1.+x)*li2(-x)*log(x)-9./4.*(59.+177.*pow(x,2)-116.*pow(x,3)+59.*pow(x,4)-118.*x)/(1.-x)*log(x)*pow(log(1.-x),2)
			+27.*(3.*pow(x,2)+2.*x+1)/(1.+x)*li2(-x)*log(1.+x)+9.*(6.-11.*pow(x,3)+18.*pow(x,2)-12.*x+6*pow(x,4))/(1.-x)*pow(log(x),2)*log(1.-x)
			+9./2.*(3.-8.*pow(x,3)+3.*pow(x,4)-6.*x+9.*pow(x,2))/(1.-x)*li2(x)*log(1.-x)-3./2.*(7.*x-7.*pow(x,3)+4.+18.*pow(x,2)-17.*pow(x,4)+9.*pow(x,5))/(1.-pow(x,2))*pow(log(x),3)
			+9./2.*(8.*pow(x,4)+16.*pow(x,3)+33.*pow(x,2)+22.*x+11.)/(1.+x)*zeta2*log(1.+x)-36.*pow(pow(x,2)+x+1.,2)/(1.+x)*li2(x)*log(1.+x)
			-9./4.*(4.*pow(x,4)+8.*pow(x,3)+27.*pow(x,2)+18.*x+9.)/(1.+x)*log(1.+x)*pow(log(x),2)+(-21.+63/2.*pow(x,2)-18.*x+33./2.*pow(x,3))*log(1.+x)*log(x)
			+27./2.*(3.*pow(x,2)+2.*x+1.)/(1.+x)*pow(log(1.+x),2)*log(x)-3./4.*(-280.*pow(x,3)+143.*pow(x,4)+394.*x-289+21.*pow(x,2))/(1.-x)*li2(x)
			+(-21.+63./2.*pow(x,2)-18.*x+33./2.*pow(x,3))*li2(-x)+(-2559./4.*pow(x,3)+1079./2.*pow(x,2)-2687./4.*x+2027./4.)*log(1.-x)
			-3./8.*(374.*pow(x,4)-389.*x+154.+699.*pow(x,2)-827.*pow(x,3))/(1.-x)*pow(log(x),2)+(330.*pow(x,3)-348.*pow(x,2)+381.*x-297.)*pow(log(1.-x),2)
			+3./4.*(-1180.*pow(x,3)+641.-1238.*x+1227.*pow(x,2)+605.*pow(x,4))/(1.-x)*log(x)*log(1.-x)-72.*(2.-x+pow(x,2))*x*pow(log(1.-x),3)
			-1./8.*(4318.*pow(x,4)-6955.*pow(x,3)+6447.*pow(x,2)-5611.*x+2333.)/(1.-x)*log(x)+3./4.*(495.*pow(x,4)-886.*pow(x,3)+564.*pow(x,2)-200.*x+16.)/(1.-x)*zeta2
			+9.*(6.*x+18.*pow(x,2)+2.+10.*pow(x,5)-6.*pow(x,3)-19.*pow(x,4))/(1.-pow(x,2))*zeta2*log(x)-9./2.*(-48.*pow(x,3)+23.*pow(x,4)-46.*x+3.+69.*pow(x,2))/(1.-x)*zeta2*log(1.-x)
			+9./2.*(-36.-15.*pow(x,4)-52.*x+19.*pow(x,2)+13.*pow(x,3)+33.*pow(x,5))/(1.-pow(x,2))*zeta3+7539./16.*pow(x,3)-24107./48.*pow(x,2)+22879./48.*x-18157./48.))

	+ nF*(1./x*(-10./3.*log(1.-x)+(14./9.-2.*zeta2)+2.*pow(log(1.-x),2)) // this piece same as mathematica
		+1./x*((31./6.*x+1./6.+65./12.*pow(x,2))*S12(x)+(-31./12.*pow(x,2)+1./6.-17./6.*x)*Li3(x)
				+(47./12.*pow(x,2)+25./6.*x-1./6.)*li2(x)*log(x)+(-1./12.*pow(x,2)+1./6.*x-1./6.)*zeta2*log(1.-x)-4.*x*(1.+x)*zeta2*log(x)
				+(-1./6.*x+1./6.+1./12.*pow(x,2))*li2(x)*log(1.-x)+(1./12.-1./12.*x+1./24.*pow(x,2))*log(1.-x)*log(x)*log((1.-x)/x)
				+5./9.*x*(1.+x)*pow(log(x),3)+(-17./6.*pow(x,2)-7./3.*x-1./3.)*zeta3+(-34./9.*pow(x,3)+2./3.*pow(x,2)-8./3.*x+16./9.)*(pow(log(1.-x),2)-zeta2)
				-2./9.*(21.*pow(x,2)+7.*x+25.*pow(x,4)+17.-61.*pow(x,3))/(1.-x)*log(x)*log(1.-x)+(785./54.*pow(x,3)-83./36.*pow(x,2)+49./18.*x-461./54.)*log(1.-x)
				+1./72.*(-351.*pow(x,3)+117.*pow(x,2)+68.+132.*pow(x,4)+52.*x)/(1.-x)*pow(log(x),2)+1./36.*(227.*pow(x,3)+68.+4.*pow(x,4)-302.*x+21.*pow(x,2))/(1.-x)*li2(1.-x)
				+1./216.*(333.*pow(x,2)+2384.*pow(x,4)-598.*x-3041.*pow(x,3)+1282.)/(1.-x)*log(x)-8887./648.*pow(x,3)+1267./432.*pow(x,2)-497./216.*x+12923./1296.)));
}

double logdep_gg(double x){
	return pow(alphas_muR/M_PI,2)*(log(muF2/Q2)*((pow(log(x),2)*(pow(CA,2)*(2. + 3.*x + 8.*pow(x,2) - 3.*pow(x,3) -
							8.*pow(x,4) + 	3.*pow(x,5))/(x*(-1. + pow(x,2)))+nF*(-1. + pow(CA,2))*(1. + x)/(2.*CA)))
						+ log(x)*(2.*log(muF2/Q2)*pow(CA,2)*(1. + 3.*pow(x,2) - 4.*pow(x,3) + pow(x,4))/((-1. + x)*x)
                      +(pow(CA,2)*(77. - 179.*x + 279.*pow(x,2) - 353.*pow(x,3) + 187.*pow(x,4)))/(6.*(-1. + x)*x)
                      + nF*((log(muF2/Q2)*(-1. + pow(CA,2))*(1. + x))/(4.*CA) + ((4. + 8.*x - 3.*pow(x,2) - 17.*pow(x,3) + 8*pow(x,4) - pow(CA,2)*(8. + 4.*x + 9.*pow(x,2) - 29.*pow(x,3) + 12.*pow(x,4))))/(12.*CA*(-1. + x)*x)))
            + ((log(muF2/Q2)*pow(CA,2)*(-99. + 94.*x - 83.*pow(x,2) + 99.*pow(x,3)))/(12.*x)
            + zeta2*(-((pow(CA,2)*(-4. - 9.*x - 24.*pow(x,2) - 14.*pow(x,3) + 6.*pow(x,4.)))/(x*(1. + x))))
            - (pow(CA,2)*(144.*(-1. + x)*pow(1. + x + pow(x,2),2)*(log(1.+x)*log(x)+li2(-x))+ (1 + x)*(-2161. + 3942.*x - 3294.*pow(x,2) + 3674.*pow(x,3) - 2161.*pow(x,4) + 288.*(3. - 2.*x + 9*pow(x,2) - 10.*pow(x,3) + 3.*pow(x,4))*li2(x) + 288.*(3. + 2.*x + 9.*pow(x,2) - 14.*pow(x,3) + 3.*pow(x,4))*(-log(1.-x)*log(x)-li2(x)))))/(72.*x*(-1. + pow(x,2)))
            + nF*((log(muF2/Q2)*(-4. - 3.*x + 3.*pow(x,2) + 4.*pow(x,3) + pow(CA,2)*(8. - 5.*x + pow(x,2) - 8.*pow(x,3))))/(24.*CA*x)
                  + zeta2*(-(((-1. + pow(CA,2))*(1. + x))/CA))
                  - ((16. - 73.*pow(CA,2) + 84.*x + 7.*pow(CA,2)*x - 66*pow(x,2) - 5*pow(CA,2)*pow(x,2) - 34*pow(x,3) + 91*pow(CA,2)*pow(x,3) + 36.*(-1. + pow(CA,2))*x*(1. + x)*li2(x) + 72.*(-1. + pow(CA,2))*x*(1. + x)*(-log(1.-x)*log(x)-li2(x))))/(36.*CA*x)
                )))
								+
								log(1.-x)*((-4*pow(log(muF2/Q2),2)*pow(CA,2)*(-1 + 2*x - pow(x,2) + pow(x,3)))/x  - (log(muF2/Q2)*pow(CA,2)*(110 - 237*x + 243*pow(x,2) - 226*pow(x,3) + 110*pow(x,4) + 6*(13 - 10*x + 39*pow(x,2) - 42*pow(x,3) + 13*pow(x,4))*log(x)))/(3.*(-1 + x)*x) + nF*((log(muF2/Q2)*(4 + 3*x - 3*pow(x,2) - 4*pow(x,3) + pow(CA,2)*(-8 + 5*x - pow(x,2) + 8*pow(x,3)) - 12*(-1 +pow(CA,2))*x*(1 + x)*log(x)))/(6.*CA*x)))
								+pow(log(1.-x),2)*((12*log(muF2/Q2)*pow(CA,2)*(-1 + 2*x - pow(x,2) + pow(x,3)))/x));
}

double logdep_gg_constant(){
  return pow(alphas_muR/M_PI,2)*((33./2. * zeta2  -171./2.*zeta3+27./2.  + (-zeta2 - 11./6.)*nF)*log(muF2/Q2) + (- 18. * zeta2)*pow(log(muF2/Q2),2.));
}

// https://arxiv.org/pdf/hep-ph/0207004.pdf eqn. 47 and 48 (note that there is no scale dependence now)
// same as in Mathematica code
double higgs_NNLO_gg_plus(double x){
	return pow(alphas_muR/M_PI,2)*((133.-90.*zeta2)*log(1.-x)/(1.-x)+(-101./3.+33.*zeta2+351./2.*zeta3)*1./(1.-x)-33.*pow(log(1.-x),2)/(1.-x)+72.*pow(log(1.-x),3)/(1.-x)+nF*(-10./3.*log(1.-x)/(1.-x)+(14./9.-2.*zeta2)*1./(1.-x)+2.*pow(log(1.-x),2)/(1.-x)));
}
// https://arxiv.org/pdf/hep-ph/0207004.pdf eqn. 47 and 48 (note that there is no scale dependence now)
// Lt = log(muR2/mt2)
// same as in Mathematica code
double higgs_NNLO_gg_delta(){
	return 0.*logdep_gg_constant()+ pow(alphas_muR/M_PI,2)*((11399./144.+133./2.*zeta2-165./4.*zeta3-9./20.*pow(zeta2,2)+19./8.*Lt)+nF*(-1189./144.+5./6.*zeta3-5./3.*zeta2+2./3.*Lt));
}
// power expansion, checked the coefficients with the mathematica code
double higgs_NNLO_gg_expansion(double x, int power){

	if(power==1){
		return pow(alphas_muR/M_PI,2)*(3147 - 684*pow(M_PI,2) + 2*nF*(-49 + 6*pow(M_PI,2)) + 3*log(1 - x)*(-1983 + 61*nF + 180*pow(M_PI,2) - 6*log(1 - x)*(-345 + 4*nF + 144*log(1 - x))) - 6318*zeta3)/36.;
	}
	if(power==2){
		return pow(alphas_muR/M_PI,2)*(1 - x)*((-55572 + 989*nF + 4734*pow(M_PI,2) - 64*nF*pow(M_PI,2) + 6*log(1 - x)*(12291 - 209*nF - 360*pow(M_PI,2) + log(1 - x)*(-5085 + 64*nF + 1728*log(1 - x))))/72. + 351*zeta3);
	}
	if(power==3){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,2)*(44828 - 965*nF - 2376*pow(M_PI,2) + 24*log(1 - x)*(-1682 + 28*nF + 1023*log(1 - x))))/96.;
	}
	if(power==4){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,3)*(-757.1067708333334 + 29*pow(M_PI,2) + nF*(11.765432098765432 - (5*pow(M_PI,2))/9.) + (798 - (28*nF)/3. - 15*pow(M_PI,2))*log(1 - x) + (5*(-855 + 16*nF)*pow(log(1 - x),2))/24. + 72*pow(log(1 - x),3) + (351*zeta3)/2.);
	}
	if(power==5){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,4)*(-327.39152083333335 + (391*pow(M_PI,2))/20. + nF*(5.709187885802469 - (5*pow(M_PI,2))/9.) + log(1 - x)*(485.4 - (272*nF)/45. - 15*pow(M_PI,2) + (-80.475 + (10*nF)/3.)*log(1 - x) + 72*pow(log(1 - x),2)) + (351*zeta3)/2.);
	}
	if(power==6){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,5)*(-233.03654166666666 + (609*pow(M_PI,2))/40. + nF*(4.4048391203703705 - (76*pow(M_PI,2))/135.) + (413.56 - (3197*nF)/675. - 15*pow(M_PI,2))*log(1 - x) + ((-6129 + 608*nF)*pow(log(1 - x),2))/180. + 72*pow(log(1 - x),3) + (351*zeta3)/2.);
	}
	if(power==7){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,6)*(-187.7184452137998 + (3457*pow(M_PI,2))/280. + nF*(3.8538812200806247 - (77*pow(M_PI,2))/135.) + (log(1 - x)*(3671379 - 36731*nF - 141750*pow(M_PI,2) + 60*log(1 - x)*(-405 + 539*nF + 11340*log(1 - x))))/9450. + (351*zeta3)/2.);
	}
	if(power==8){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,7)*(-157.21313551536687 + (1135*pow(M_PI,2))/112. + nF*(3.5441843393442034 - (109*pow(M_PI,2))/189.) + log(1 - x)*(379.50204081632654 - (85373*nF)/26460. - 15*pow(M_PI,2) + (21.776785714285715 + (218*nF)/63.)*log(1 - x) + 72*pow(log(1 - x),2)) + (351*zeta3)/2.);
	}
	if(power==9){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,8)*((-948743967014 + 23810115581*nF - 5174400*(-11421 + 800*nF)*pow(M_PI,2))/7.112448e9 - ((-9994419 + 70717*nF + 396900*pow(M_PI,2))*log(1 - x))/26460. + (41.87321428571428 + (220*nF)/63.)*pow(log(1 - x),2) + 72*pow(log(1 - x),3) + (351*zeta3)/2.);
	}
	if(power==10){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,9)*(-113.44591817050895 + (809*pow(M_PI,2))/120. + nF*(3.2171579041421463 - (95*pow(M_PI,2))/162.) + log(1 - x)*(379.74150793650796 - (7451*nF)/3402. - 15*pow(M_PI,2) + log(1 - x)*(59.1 + (95*nF)/27. + 72*log(1 - x))) + (351*zeta3)/2.);
	}
	if(power==-1){
		return higgs_NNLO_gg_reg(x)-pow(alphas_muR/M_PI,2)*((pow(-1 + x,2)*(44828 - 965*nF - 2376*pow(M_PI,2) + 24*log(1 - x)*(-1682 + 28*nF + 1023*log(1 - x))))/96. + (3147 - 684*pow(M_PI,2) + 2*nF*(-49 + 6*pow(M_PI,2)) + 3*log(1 - x)*(-1983 + 61*nF + 180*pow(M_PI,2) - 6*log(1 - x)*(-345 + 4*nF + 144*log(1 - x))) - 6318*zeta3)/36. + pow(-1 + x,8)*((-948743967014 + 23810115581*nF - 5174400*(-11421 + 800*nF)*pow(M_PI,2))/7.112448e9 - ((-9994419 + 70717*nF + 396900*pow(M_PI,2))*log(1 - x))/26460. + (41.87321428571428 + (220*nF)/63.)*pow(log(1 - x),2) + 72*pow(log(1 - x),3) + (351*zeta3)/2.) + pow(1 - x,3)*(-757.1067708333334 + 29*pow(M_PI,2) + nF*(11.765432098765432 - (5*pow(M_PI,2))/9.) + (798 - (28*nF)/3. - 15*pow(M_PI,2))*log(1 - x) + (5*(-855 + 16*nF)*pow(log(1 - x),2))/24. + 72*pow(log(1 - x),3) + (351*zeta3)/2.) + pow(1 - x,5)*(-233.03654166666666 + (609*pow(M_PI,2))/40. + nF*(4.4048391203703705 - (76*pow(M_PI,2))/135.) + (413.56 - (3197*nF)/675. - 15*pow(M_PI,2))*log(1 - x) + ((-6129 + 608*nF)*pow(log(1 - x),2))/180. + 72*pow(log(1 - x),3) + (351*zeta3)/2.) + pow(-1 + x,4)*(-327.39152083333335 + (391*pow(M_PI,2))/20. + nF*(5.709187885802469 - (5*pow(M_PI,2))/9.) + log(1 - x)*(485.4 - (272*nF)/45. - 15*pow(M_PI,2) + (-80.475 + (10*nF)/3.)*log(1 - x) + 72*pow(log(1 - x),2)) + (351*zeta3)/2.) + pow(1 - x,7)*(-157.21313551536687 + (1135*pow(M_PI,2))/112. + nF*(3.5441843393442034 - (109*pow(M_PI,2))/189.) + log(1 - x)*(379.50204081632654 - (85373*nF)/26460. - 15*pow(M_PI,2) + (21.776785714285715 + (218*nF)/63.)*log(1 - x) + 72*pow(log(1 - x),2)) + (351*zeta3)/2.) + pow(1 - x,9)*(-113.44591817050895 + (809*pow(M_PI,2))/120. + nF*(3.2171579041421463 - (95*pow(M_PI,2))/162.) + log(1 - x)*(379.74150793650796 - (7451*nF)/3402. - 15*pow(M_PI,2) + log(1 - x)*(59.1 + (95*nF)/27. + 72*log(1 - x))) + (351*zeta3)/2.) + pow(-1 + x,6)*(-187.7184452137998 + (3457*pow(M_PI,2))/280. + nF*(3.8538812200806247 - (77*pow(M_PI,2))/135.) + (log(1 - x)*(3671379 - 36731*nF - 141750*pow(M_PI,2) + 60*log(1 - x)*(-405 + 539*nF + 11340*log(1 - x))))/9450. + (351*zeta3)/2.) + (1 - x)*((-55572 + 989*nF + 4734*pow(M_PI,2) - 64*nF*pow(M_PI,2) + 6*log(1 - x)*(12291 - 209*nF - 360*pow(M_PI,2) + log(1 - x)*(-5085 + 64*nF + 1728*log(1 - x))))/72. + 351*zeta3));
	}
	if(power==-2){ // sanity check of above
		return pow(alphas_muR/M_PI,2)*((pow(-1 + x,2)*(44828 - 965*nF - 2376*pow(M_PI,2) + 24*log(1 - x)*(-1682 + 28*nF + 1023*log(1 - x))))/96. + (3147 - 684*pow(M_PI,2) + 2*nF*(-49 + 6*pow(M_PI,2)) + 3*log(1 - x)*(-1983 + 61*nF + 180*pow(M_PI,2) - 6*log(1 - x)*(-345 + 4*nF + 144*log(1 - x))) - 6318*zeta3)/36. + pow(-1 + x,8)*((-948743967014 + 23810115581*nF - 5174400*(-11421 + 800*nF)*pow(M_PI,2))/7.112448e9 - ((-9994419 + 70717*nF + 396900*pow(M_PI,2))*log(1 - x))/26460. + (41.87321428571428 + (220*nF)/63.)*pow(log(1 - x),2) + 72*pow(log(1 - x),3) + (351*zeta3)/2.) + pow(1 - x,3)*(-757.1067708333334 + 29*pow(M_PI,2) + nF*(11.765432098765432 - (5*pow(M_PI,2))/9.) + (798 - (28*nF)/3. - 15*pow(M_PI,2))*log(1 - x) + (5*(-855 + 16*nF)*pow(log(1 - x),2))/24. + 72*pow(log(1 - x),3) + (351*zeta3)/2.) + pow(1 - x,5)*(-233.03654166666666 + (609*pow(M_PI,2))/40. + nF*(4.4048391203703705 - (76*pow(M_PI,2))/135.) + (413.56 - (3197*nF)/675. - 15*pow(M_PI,2))*log(1 - x) + ((-6129 + 608*nF)*pow(log(1 - x),2))/180. + 72*pow(log(1 - x),3) + (351*zeta3)/2.) + pow(-1 + x,4)*(-327.39152083333335 + (391*pow(M_PI,2))/20. + nF*(5.709187885802469 - (5*pow(M_PI,2))/9.) + log(1 - x)*(485.4 - (272*nF)/45. - 15*pow(M_PI,2) + (-80.475 + (10*nF)/3.)*log(1 - x) + 72*pow(log(1 - x),2)) + (351*zeta3)/2.) + pow(1 - x,7)*(-157.21313551536687 + (1135*pow(M_PI,2))/112. + nF*(3.5441843393442034 - (109*pow(M_PI,2))/189.) + log(1 - x)*(379.50204081632654 - (85373*nF)/26460. - 15*pow(M_PI,2) + (21.776785714285715 + (218*nF)/63.)*log(1 - x) + 72*pow(log(1 - x),2)) + (351*zeta3)/2.) + pow(1 - x,9)*(-113.44591817050895 + (809*pow(M_PI,2))/120. + nF*(3.2171579041421463 - (95*pow(M_PI,2))/162.) + log(1 - x)*(379.74150793650796 - (7451*nF)/3402. - 15*pow(M_PI,2) + log(1 - x)*(59.1 + (95*nF)/27. + 72*log(1 - x))) + (351*zeta3)/2.) + pow(-1 + x,6)*(-187.7184452137998 + (3457*pow(M_PI,2))/280. + nF*(3.8538812200806247 - (77*pow(M_PI,2))/135.) + (log(1 - x)*(3671379 - 36731*nF - 141750*pow(M_PI,2) + 60*log(1 - x)*(-405 + 539*nF + 11340*log(1 - x))))/9450. + (351*zeta3)/2.) + (1 - x)*((-55572 + 989*nF + 4734*pow(M_PI,2) - 64*nF*pow(M_PI,2) + 6*log(1 - x)*(12291 - 209*nF - 360*pow(M_PI,2) + log(1 - x)*(-5085 + 64*nF + 1728*log(1 - x))))/72. + 351*zeta3));
	}

	else{return 0;}
}


///////////////////////////
/// qg channel
///////////////////////////
//NNLO coefficients  qg channel
// klopt met de code
double higgs_NNLO_qg_reg(double x){
	return  0.*logdep_qg(x)+pow(alphas_muR/M_PI,2)*
		(1./x*((170./3.*x+338./9.+119./3.*pow(x,2))*Li3(x)+(4.*x+4.+2.*pow(x,2))*Li3(-x)+(16.+8.*pow(x,2)+16.*x)*S12(-x)
				+(-614./9.*x-269./9.*pow(x,2)-74./9.)*S12(x)+(-2.*pow(x,2)-4.-4.*x)*S12(pow(x,2))+(367./27.+367./54.*pow(x,2)-367./27.*x)*pow(log(1.-x),3)
				+((2.+pow(x,2)-2.*x)*log(1.-x)-(446./9.*x+214./9.+281./9.*pow(x,2))*log(x)-(8.+4.*pow(x,2)+8.*x)*log(1.+x))*li2(x)
				+(8.+8.*x+4.*pow(x,2))*log((1.+x)/x)*li2(-x)+(-115./9.*pow(x,2)+230./9.*x-230./9.)*log(x)*pow(log(1.-x),2)
				+(107./9+107./18.*pow(x,2)-107./9.*x)*pow(log(x),2)*log(1.-x)+(-145./54.*pow(x,2)-71./27.*x-2.)*pow(log(x),3)
				+(-3.*pow(x,2)-6.-6.*x)*log(1.+x)*pow(log(x),2)+(4.*x+4.+2.*pow(x,2))*pow(log(1.+x),2)*log(x)
				+(-4./27.*pow(x,3)-74./9.*x-11./9.*pow(x,2)-166./27.)*li2(-x)+(2605./54.-146./9.*x+74./27.*pow(x,3)-79./6.*pow(x,2))*li2(x)
				+(1139./18.*x+37./12.*pow(x,2)+8.*pow(x,3)-72.)*pow(log(1.-x),2)+(-121./18.*pow(x,2)-326./27.*pow(x,3)-826./9.*x+5935./54.)*log(x)*log(1.-x)
				+(113./27.*pow(x,3)+244./9.*x-13./3.*pow(x,2)-31./2.)*pow(log(x),2)+(-4./27.*pow(x,3)-74./9.*x-11./9.*pow(x,2)-166./27.)*log(1.+x)*log(x)
				+zeta2*(-59./9.*pow(x,2)+118./9.*x-118./9)*log(1.-x)+zeta2*(140./9.*x+128./9.*pow(x,2)+52./9.)*log(x)
				+zeta2*(12.+12.*x+6.*pow(x,2))*log(1.+x)+(-392./81.*pow(x,3)-49./3.*pow(x,2)+23671./162.-106.*x)*log(1.-x)
				+(1985./108.*pow(x,2)+800./9.*x-12209./162.+616./81.*pow(x,3))*log(x)+(-292./27.*pow(x,3)-82./3.*x+16./3.*pow(x,2)+221./27.)*zeta2
				+(-18.*x+10.*pow(x,2)+92./9.)*zeta3-210115./1944.+1537./486.*pow(x,3)+16465./162.*x+2393./648.*pow(x,2))

		+nF*1./x*((1./18.*pow(x,2)-1./9.*x+1./9.)*pow(log(1.-x),2)+(-38./27.*x+19./27.*pow(x,2)+29./27.)*log(x)-209./81.*x+265./162. //klopt
				+((-4./9.+4./9.*x-2./9.*pow(x,2))*log(x)-pow(x,2)+16./9.*x-13./9.)*log(1.-x)+179./162.*pow(x,2)+(1./9.*pow(x,2)-2./9.*x+2./9.)*pow(log(x),2))
		);

}


double logdep_qg(double x){
return pow(alphas_muR/M_PI,2)*((-1 + pow(CA,2))*( - 1728*log(muF2/Q2)*pow(CA,2) - 1728*log(muF2/Q2)*pow(CA,2)*x - 864*log(muF2/Q2)*pow(CA,2)*pow(x,2))/(3456.*pow(CA,2)*x) +
(-1 + pow(CA,2))*(864*log(muF2/Q2)+ 4320*log(muF2/Q2)*pow(CA,2) + 15552*log(muF2/Q2)*pow(CA,2)*x + 7776*log(muF2/Q2)*pow(CA,2)*pow(x,2))/(3456.*pow(CA,2)*x)
 + ((-1 + pow(CA,2))*(636*log(muF2/Q2) - 29700*log(muF2/Q2)*pow(CA,2) - 8280*pow(log(muF2/Q2),2)*pow(CA,2) + 1392*log(muF2/Q2)*CA*nF + 288*pow(log(muF2/Q2),2)*CA*nF - 648*log(muF2/Q2)*x - 216*pow(log(muF2/Q2),2)*x + 21336*log(muF2/Q2)*pow(CA,2)*x + 6984*pow(log(muF2/Q2),2)*pow(CA,2)*x  - 1824*log(muF2/Q2)*CA*nF*x - 288*pow(log(muF2/Q2),2)*CA*nF*x - 216*log(muF2/Q2)*pow(x,2) + 54*pow(log(muF2/Q2),2)*pow(x,2)  + 4560*log(muF2/Q2)*pow(CA,2)*pow(x,2)  + 912*log(muF2/Q2)*CA*nF*pow(x,2) + 144*pow(log(muF2/Q2),2)*CA*nF*pow(x,2) - 744*log(muF2/Q2)*pow(x,3) + 1128*log(muF2/Q2)*pow(CA,2)*pow(x,3) + 864*pow(log(muF2/Q2),2)*pow(CA,2)*pow(x,3) - 864*log(muF2/Q2)*zeta2 +  6048*log(muF2/Q2)*pow(CA,2)*zeta2  + 1728*log(muF2/Q2)*x*zeta2 + 12096*log(muF2/Q2)*pow(CA,2)*x*zeta2 - 864*log(muF2/Q2)*pow(x,2)*zeta2 + 8640*log(muF2/Q2)*pow(CA,2)*pow(x,2)*zeta2))/(3456.*pow(CA,2)*x)
 +log(x)*(-((-1 + pow(CA,2))*(2232*log(muF2/Q2)*pow(CA,2) + 432*pow(log(muF2/Q2),2)*pow(CA,2) - 96*log(muF2/Q2)*CA*nF + 72*log(muF2/Q2)*x + 36*pow(log(muF2/Q2),2)*x - 2760*log(muF2/Q2)*pow(CA,2)*x + 396*pow(log(muF2/Q2),2)*pow(CA,2)*x  + 96*log(muF2/Q2)*CA*nF*x  - 342*log(muF2/Q2)*pow(x,2) - 18*pow(log(muF2/Q2),2)*pow(x,2)  + 390*log(muF2/Q2)*pow(CA,2)*pow(x,2) + 450*pow(log(muF2/Q2),2)*pow(CA,2)*pow(x,2) - 48*log(muF2/Q2)*CA*nF*pow(x,2) - 48*log(muF2/Q2)*pow(x,3)  - 528*log(muF2/Q2)*pow(CA,2)*pow(x,3)))/(576.*pow(CA,2)*x))
 +pow(log(x),2)*(((-1 + pow(CA,2))*(x*(24*log(muF2/Q2)*(-2 + x)) - pow(CA,2)*(24*log(muF2/Q2)*(12 + 14*x + 15*pow(x,2)))))/(384.*pow(CA,2)*x))
 +log(1-x)*(((-1 + pow(CA,2))*(- 162*log(muF2/Q2) - 36*pow(log(muF2/Q2),2) + 2766*log(muF2/Q2)*pow(CA,2) + 252*pow(log(muF2/Q2),2)*pow(CA,2)- 48*log(muF2/Q2)*CA*nF+ 288*log(muF2/Q2)*x + 36*pow(log(muF2/Q2),2)*x - 2496*log(muF2/Q2)*pow(CA,2)*x - 252*pow(log(muF2/Q2),2)*pow(CA,2)*x + 48*log(muF2/Q2)*CA*nF*x - 126*log(muF2/Q2)*pow(x,2) - 18*pow(log(muF2/Q2),2)*pow(x,2)  + 6*log(muF2/Q2)*pow(CA,2)*pow(x,2) + 126*pow(log(muF2/Q2),2)*pow(CA,2)*pow(x,2)  - 24*log(muF2/Q2)*CA*nF*pow(x,2)  - 288*log(muF2/Q2)*pow(CA,2)*pow(x,3)))/(288.*pow(CA,2)*x) + ((-1 + pow(CA,2))*(12*log(muF2/Q2)*(-pow(-1 + x,2) + pow(CA,2)*(17 + 6*x + 15*pow(x,2))))*(log(x)))/(48.*pow(CA,2)*x)
 +pow(log(1-x),2)*(((-1 + pow(CA,2))*(- 36*log(muF2/Q2)*(-1 + 7*pow(CA,2))*(2 - 2*x + pow(x,2))))/(192.*pow(CA,2)*x))));

}

// checked
double higgs_NNLO_qg_expansion(double x, int power){
	if(power==1){
		return pow(alphas_muR/M_PI,2)*(132 + 52*nF + 261*pow(M_PI,2) + 3*log(1 - x)*(2046 - 72*nF - 100*pow(M_PI,2) + log(1 - x)*(255 + 6*nF + 734*log(1 - x))) + 5598*zeta3)/324.;
	}
	if(power==2){
		return pow(alphas_muR/M_PI,2)*((1 - x)*(-17130 + 52*nF + 693*pow(M_PI,2) + 3*log(1 - x)*(7242 - 24*nF - 100*pow(M_PI,2) + log(1 - x)*(-1527 + 6*nF + 734*log(1 - x))) + 5598*zeta3))/324.;
	}
	if(power==3){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,2)*(-53805 + 772*nF + 2016*pow(M_PI,2) + 6*log(1 - x)*(-2*(-5451 + 60*nF + 100*pow(M_PI,2)) + log(1 - x)*(-2223 + 12*nF + 1468*log(1 - x))) + 22392*zeta3))/648.;
	}
	if(power==4){
		return pow(alphas_muR/M_PI,2)*(pow(1 - x,3)*(-31941 + 308*nF + 1696*pow(M_PI,2) + 6*log(1 - x)*(9018 - 88*nF - 200*pow(M_PI,2) + log(1 - x)*(-887 + 12*nF + 1468*log(1 - x))) + 22392*zeta3))/648.;
	}
	if(power==5){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,4)*(-140757 + 1264*nF + 8984*pow(M_PI,2) + 24*log(1 - x)*(15237 - 140*nF - 400*pow(M_PI,2) + log(1 - x)*(515 + 24*nF + 2936*log(1 - x))) + 179136*zeta3))/5184.;
	}
	if(power==6){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,5)*(-16.547468621399176 + (5*nF)/36. + (3527*pow(M_PI,2))/3240. + (log(1 - x)*(1087583 - 8580*nF - 30000*pow(M_PI,2) + 75*log(1 - x)*(2163 + 24*nF + 2936*log(1 - x))))/16200. + (311*zeta3)/9.);
	}
	if(power==7){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,6)*(-8.753997942386832 + (71*nF)/810. + (101*pow(M_PI,2))/180. + (log(1 - x)*(53977 - 354*nF - 1500*pow(M_PI,2) + 6*log(1 - x)*(2176 + 15*nF + 1835*log(1 - x))))/810. + (311*zeta3)/9.);
	}
	if(power==8){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,7)*(-2.3669754826092695 + (73*nF)/1134. + (409*pow(M_PI,2))/3780. + (log(1 - x)*(13392887 - 71610*nF - 367500*pow(M_PI,2) + 210*log(1 - x)*(20086 + 105*nF + 12845*log(1 - x))))/198450. + (311*zeta3)/9.);
	}
	if(power==9){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,8)*(3.1236286366288137 + (2603*nF)/45360. - (4399*pow(M_PI,2))/15120. + (log(1 - x)*(219331957 - 939960*nF - 5880000*pow(M_PI,2) + 210*log(1 - x)*(388701 + 1680*nF + 205520*log(1 - x))))/3.1752e6 + (311*zeta3)/9.);
	}
	if(power==10){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,9)*(7.980193194991841 + (1649*nF)/27216. - (88271*pow(M_PI,2))/136080. + (log(1 - x)*(225779997 - 760760*nF - 5880000*pow(M_PI,2) + 70*log(1 - x)*(1344743 + 5040*nF + 616560*log(1 - x))))/3.1752e6 + (311*zeta3)/9.);
	}
	if(power==-1){
		return higgs_NNLO_qqp_reg(x)-pow(alphas_muR/M_PI,2)*(pow(-1 + x,6)*(-8.753997942386832 + (71*nF)/810. + (101*pow(M_PI,2))/180. + (log(1 - x)*(53977 - 354*nF - 1500*pow(M_PI,2) + 6*log(1 - x)*(2176 + 15*nF + 1835*log(1 - x))))/810. + (311*zeta3)/9.) + pow(1 - x,5)*(-16.547468621399176 + (5*nF)/36. + (3527*pow(M_PI,2))/3240. + (log(1 - x)*(1087583 - 8580*nF - 30000*pow(M_PI,2) + 75*log(1 - x)*(2163 + 24*nF + 2936*log(1 - x))))/16200. + (311*zeta3)/9.) + pow(1 - x,7)*(-2.3669754826092695 + (73*nF)/1134. + (409*pow(M_PI,2))/3780. + (log(1 - x)*(13392887 - 71610*nF - 367500*pow(M_PI,2) + 210*log(1 - x)*(20086 + 105*nF + 12845*log(1 - x))))/198450. + (311*zeta3)/9.) + pow(-1 + x,8)*(3.1236286366288137 + (2603*nF)/45360. - (4399*pow(M_PI,2))/15120. + (log(1 - x)*(219331957 - 939960*nF - 5880000*pow(M_PI,2) + 210*log(1 - x)*(388701 + 1680*nF + 205520*log(1 - x))))/3.1752e6 + (311*zeta3)/9.) + pow(1 - x,9)*(7.980193194991841 + (1649*nF)/27216. - (88271*pow(M_PI,2))/136080. + (log(1 - x)*(225779997 - 760760*nF - 5880000*pow(M_PI,2) + 70*log(1 - x)*(1344743 + 5040*nF + 616560*log(1 - x))))/3.1752e6 + (311*zeta3)/9.) + ((1 - x)*(-17130 + 52*nF + 693*pow(M_PI,2) + 3*log(1 - x)*(7242 - 24*nF - 100*pow(M_PI,2) + log(1 - x)*(-1527 + 6*nF + 734*log(1 - x))) + 5598*zeta3))/324. + (132 + 52*nF + 261*pow(M_PI,2) + 3*log(1 - x)*(2046 - 72*nF - 100*pow(M_PI,2) + log(1 - x)*(255 + 6*nF + 734*log(1 - x))) + 5598*zeta3)/324. + (pow(-1 + x,2)*(-53805 + 772*nF + 2016*pow(M_PI,2) + 6*log(1 - x)*(-2*(-5451 + 60*nF + 100*pow(M_PI,2)) + log(1 - x)*(-2223 + 12*nF + 1468*log(1 - x))) + 22392*zeta3))/648. + (pow(1 - x,3)*(-31941 + 308*nF + 1696*pow(M_PI,2) + 6*log(1 - x)*(9018 - 88*nF - 200*pow(M_PI,2) + log(1 - x)*(-887 + 12*nF + 1468*log(1 - x))) + 22392*zeta3))/648. + (pow(-1 + x,4)*(-140757 + 1264*nF + 8984*pow(M_PI,2) + 24*log(1 - x)*(15237 - 140*nF - 400*pow(M_PI,2) + log(1 - x)*(515 + 24*nF + 2936*log(1 - x))) + 179136*zeta3))/5184.);
	}
	if(power==-2){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,6)*(-8.753997942386832 + (71*nF)/810. + (101*pow(M_PI,2))/180. + (log(1 - x)*(53977 - 354*nF - 1500*pow(M_PI,2) + 6*log(1 - x)*(2176 + 15*nF + 1835*log(1 - x))))/810. + (311*zeta3)/9.) + pow(1 - x,5)*(-16.547468621399176 + (5*nF)/36. + (3527*pow(M_PI,2))/3240. + (log(1 - x)*(1087583 - 8580*nF - 30000*pow(M_PI,2) + 75*log(1 - x)*(2163 + 24*nF + 2936*log(1 - x))))/16200. + (311*zeta3)/9.) + pow(1 - x,7)*(-2.3669754826092695 + (73*nF)/1134. + (409*pow(M_PI,2))/3780. + (log(1 - x)*(13392887 - 71610*nF - 367500*pow(M_PI,2) + 210*log(1 - x)*(20086 + 105*nF + 12845*log(1 - x))))/198450. + (311*zeta3)/9.) + pow(-1 + x,8)*(3.1236286366288137 + (2603*nF)/45360. - (4399*pow(M_PI,2))/15120. + (log(1 - x)*(219331957 - 939960*nF - 5880000*pow(M_PI,2) + 210*log(1 - x)*(388701 + 1680*nF + 205520*log(1 - x))))/3.1752e6 + (311*zeta3)/9.) + pow(1 - x,9)*(7.980193194991841 + (1649*nF)/27216. - (88271*pow(M_PI,2))/136080. + (log(1 - x)*(225779997 - 760760*nF - 5880000*pow(M_PI,2) + 70*log(1 - x)*(1344743 + 5040*nF + 616560*log(1 - x))))/3.1752e6 + (311*zeta3)/9.) + ((1 - x)*(-17130 + 52*nF + 693*pow(M_PI,2) + 3*log(1 - x)*(7242 - 24*nF - 100*pow(M_PI,2) + log(1 - x)*(-1527 + 6*nF + 734*log(1 - x))) + 5598*zeta3))/324. + (132 + 52*nF + 261*pow(M_PI,2) + 3*log(1 - x)*(2046 - 72*nF - 100*pow(M_PI,2) + log(1 - x)*(255 + 6*nF + 734*log(1 - x))) + 5598*zeta3)/324. + (pow(-1 + x,2)*(-53805 + 772*nF + 2016*pow(M_PI,2) + 6*log(1 - x)*(-2*(-5451 + 60*nF + 100*pow(M_PI,2)) + log(1 - x)*(-2223 + 12*nF + 1468*log(1 - x))) + 22392*zeta3))/648. + (pow(1 - x,3)*(-31941 + 308*nF + 1696*pow(M_PI,2) + 6*log(1 - x)*(9018 - 88*nF - 200*pow(M_PI,2) + log(1 - x)*(-887 + 12*nF + 1468*log(1 - x))) + 22392*zeta3))/648. + (pow(-1 + x,4)*(-140757 + 1264*nF + 8984*pow(M_PI,2) + 24*log(1 - x)*(15237 - 140*nF - 400*pow(M_PI,2) + log(1 - x)*(515 + 24*nF + 2936*log(1 - x))) + 179136*zeta3))/5184.);}
	else{return 0;}
}

///////////////////////////
/// qq channel (stricktly identical quarks!)
///////////////////////////
//NNLO coefficients  qq channel
//klopt met code
double higgs_NNLO_qq_reg(double x){
	return pow(alphas_muR/M_PI,2)*
		(1./x*((368./27.*x+104./27.*pow(x,2)+400./27.)*Li3(x)-32./9.*pow(x+2.,2)*S12(x)-4./27.*(2.+pow(x,2)-2.*x)*pow(log(x),2)*log(1.-x)
				-4./27.*pow(x+2.,2)*pow(log(x),3)-16./27.*(19.+5.*pow(x,2)+17.*x)*li2(x)*log(x)-32./9.*(x+3.)*(1.-x)*pow(log(1.-x),2)
				+16./3.*(x+3.)*(1.-x)*log(x)*log(1.-x)+4./27.*(26.*x-18.+9.*pow(x,2))*pow(log(x),2)-8./9.*(-6.+pow(x,2)+4.*x)*li2(x)
				+4./3.*(5.*x+17.)*(1.-x)*log(1.-x)+(8./9.*pow(x+2.,2)*zeta2-118./9.+248./27.*x+46./9.*pow(x,2))*log(x)
				+(-8./27.*pow(x,2)-16./27.+16./27.*x)*zeta3+(16./3.-32./9.*x-8./3.*pow(x,2))*zeta2-4./27.*(27.*x+160.)*(1.-x))
		+nF*1./x*(0.)
		);

}
// checked
double higgs_NNLO_qq_expansion(double x, int power){
	if(power==1){
		return 0;
	}
	if(power==2){
		return pow(alphas_muR/M_PI,2)*((1 - x)*(60 - 8*pow(M_PI,2) + 48*(-1 + log(1 - x))*log(1 - x)))/27.;
	}
	if(power==3){
		return pow(alphas_muR/M_PI,2)*(-4*pow(-1 + x,2)*(-31 + 3*pow(M_PI,2) + 18*log(1 - x) - 18*pow(log(1 - x),2)))/27.;
	}
	if(power==4){
		return pow(alphas_muR/M_PI,2)*(2*pow(-1 + x,3)*(-335 + 30*pow(M_PI,2) + 204*log(1 - x) - 180*pow(log(1 - x),2)))/81.;
	}
	if(power==5){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,4)*(3301 - 312*pow(M_PI,2) + 72*log(1 - x)*(-25 + 26*log(1 - x))))/324.;
	}
	if(power==6){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,5)*(-1438609 + 139800*pow(M_PI,2) + 676920*log(1 - x) - 838800*pow(log(1 - x),2)))/121500.;
	}
	if(power==7){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,6)*(4856549 - 478800*pow(M_PI,2) + 6840*log(1 - x)*(-283 + 420*log(1 - x))))/364500.;
	}
	if(power==8){
		return pow(alphas_muR/M_PI,2)*(pow(1 - x,7)*(1836307541 - 182221200*pow(M_PI,2) + 2520*log(1 - x)*(-242239 + 433860*log(1 - x))))/1.250235e8;
	}
	if(power==9){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,8)*(15.962012734605894 - (1499*pow(M_PI,2))/945. + (log(1 - x)*(-431723 + 944370*log(1 - x)))/99225.);
	}
	if(power==10){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,9)*(17.165953340438666 - (1609*pow(M_PI,2))/945. + (log(1 - x)*(-372563 + 1013670*log(1 - x)))/99225.);
	}
	if(power==-1){
		return higgs_NNLO_qq_reg(x)-pow(alphas_muR/M_PI,2)*(((-1 + x)*(-294819291614 + 1508892023605*x - 3838567340315*pow(x,2) + 6014649603745*pow(x,3) - 6178348500059*pow(x,4) + 4194625584295*pow(x,5) - 1821547831769*pow(x,6) + 459955373731*pow(x,7) - 51507541619*pow(x,8) + 1058400*pow(M_PI,2)*(27372 - 140205*x + 358675*pow(x,2) - 563365*pow(x,3) + 579147*pow(x,4) - 393265*pow(x,5) + 170767*pow(x,6) - 43113*pow(x,7) + 4827*pow(x,8))))/3.000564e9 + ((-3860578 + 20533963*x - 54178420*pow(x,2) + 91200410*pow(x,3) - 105339808*pow(x,4) + 85007958*pow(x,5) - 47301828*pow(x,6) + 17350530*pow(x,7) - 3784790*pow(x,8) + 372563*pow(x,9))*log(1 - x))/99225. - (2*(-27372 + 167577*x - 498880*pow(x,2) + 922040*pow(x,3) - 1142512*pow(x,4) + 972412*pow(x,5) - 564032*pow(x,6) + 213880*pow(x,7) - 47940*pow(x,8) + 4827*pow(x,9))*pow(log(1 - x),2))/945.);
	}
	if(power==-2){
		return pow(alphas_muR/M_PI,2)*(((-1 + x)*(-294819291614 + 1508892023605*x - 3838567340315*pow(x,2) + 6014649603745*pow(x,3) - 6178348500059*pow(x,4) + 4194625584295*pow(x,5) - 1821547831769*pow(x,6) + 459955373731*pow(x,7) - 51507541619*pow(x,8) + 1058400*pow(M_PI,2)*(27372 - 140205*x + 358675*pow(x,2) - 563365*pow(x,3) + 579147*pow(x,4) - 393265*pow(x,5) + 170767*pow(x,6) - 43113*pow(x,7) + 4827*pow(x,8))))/3.000564e9 + ((-3860578 + 20533963*x - 54178420*pow(x,2) + 91200410*pow(x,3) - 105339808*pow(x,4) + 85007958*pow(x,5) - 47301828*pow(x,6) + 17350530*pow(x,7) - 3784790*pow(x,8) + 372563*pow(x,9))*log(1 - x))/99225. - (2*(-27372 + 167577*x - 498880*pow(x,2) + 922040*pow(x,3) - 1142512*pow(x,4) + 972412*pow(x,5) - 564032*pow(x,6) + 213880*pow(x,7) - 47940*pow(x,8) + 4827*pow(x,9))*pow(log(1 - x),2))/945.);
	}
	else{return 0;}
}


double logdep_qq(double x){
	return pow(alphas_muR/M_PI,2)*((pow(-1 + pow(CA,2),2)*( 2*log(muF2/Q2)*pow(2 + x,2))*li2(x))/(8.*pow(CA,2)*x) + (pow(-1 + pow(CA,2),2)*( 192*log(muF2/Q2)*CA + 192*log(muF2/Q2)*CA*x + 48*log(muF2/Q2)*CA*pow(x,2))*(-log(1-x)*log(x)-li2(x))))/(96.*pow(CA,3)*x) + (pow(-1 + pow(CA,2),2)*(- 153*log(muF2/Q2)*CA - 36*pow(log(muF2/Q2),2)*CA + 108*log(muF2/Q2)*CA*x + 24*pow(log(muF2/Q2),2)*CA*x + 45*log(muF2/Q2)*CA*pow(x,2) + 12*pow(log(muF2/Q2),2)*CA*pow(x,2) + 96*log(muF2/Q2)*CA*zeta2  + 96*log(muF2/Q2)*CA*x*zeta2 + 24*log(muF2/Q2)*CA*pow(x,2)*zeta2 ))/(96.*pow(CA,3)*x)
	+ log(x)*(-(pow(-1 + pow(CA,2),2)*(CA*(6*pow(log(muF2/Q2),2)*pow(2 + x,2) - 6*log(muF2/Q2)*(-12 + 8*x + 5*pow(x,2)))))/(96.*pow(CA,3)*x))
	+pow(log(x),2)*(-(pow(-1 + pow(CA,2),2)*(CA*(log(muF2/Q2)*pow(2 + x,2))))/(16.*pow(CA,3)*x))
	+log(1.-x)*((pow(-1 + pow(CA,2),2)*( 72*log(muF2/Q2) - 48*log(muF2/Q2)*x  - 24*log(muF2/Q2)*pow(x,2)))/(48.*pow(CA,2)*x) + (pow(-1 + pow(CA,2),2)*(log(muF2/Q2)*pow(2 + x,2))*log(x))/(2.*pow(CA,2)*x));
}
///////////////////////////
/// qq' channel
///////////////////////////
//NNLO coefficients  qq' channel
// klopt met code
double higgs_NNLO_qqp_reg(double x){
	return pow(alphas_muR/M_PI,2)*
		(1./x*(32./9.*pow(x+2.,2)*(Li3(x)-S12(x))-8./3.*pow(x+2.,2)*log(x)*li2(x)-4./27.*pow(x+2.,2)*pow(log(x),3)
				-8./9.*(4.*x-6.+pow(x,2))*li2(x)-32./9.*(x+3.)*(1.-x)*pow(log(1.-x),2)+16./3.*(x+3.)*(1.-x)*log(x)*log(1.-x)
				+8./9.*(pow(x,2)+4.*x-3.)*pow(log(x),2)+8./9.*zeta2*pow(x+2.,2)*log(x)+4./3.*(5.*x+17.)*(1.-x)*log(1.-x)
				+2./9.*(29.*pow(x,2)+44.*x-59.)*log(x)+(16./3.-32./9.*x-8./3.*pow(x,2))*zeta2-2./9.*(11.*x+105.)*(1.-x))
		+nF*1./x*(0.)
		);

}
// checked
double higgs_NNLO_qqp_expansion(double x, int power){
	if(power==1){
		return 0;
	}
	if(power==2){
		return pow(alphas_muR/M_PI,2)*((1 - x)*(60 - 8*pow(M_PI,2) + 48*(-1 + log(1 - x))*log(1 - x)))/27.;
	}
	if(power==3){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,2)*(42 - 4*pow(M_PI,2) + 24*(-1 + log(1 - x))*log(1 - x)))/9.;
	}
	if(power==4){
		return pow(alphas_muR/M_PI,2)*(4*pow(-1 + x,3)*(-172 + 15*pow(M_PI,2) + 102*log(1 - x) - 90*pow(log(1 - x),2)))/81.;
	}
	if(power==5){
		return pow(alphas_muR/M_PI,2)*(-2*pow(-1 + x,4)*(-142 + 13*pow(M_PI,2) + 75*log(1 - x) - 78*pow(log(1 - x),2)))/27.;
	}
	if(power==6){
		return pow(alphas_muR/M_PI,2)*(pow(1 - x,5)*(744467 - 69900*pow(M_PI,2) + 60*log(1 - x)*(-5641 + 6990*log(1 - x))))/60750.;
	}
	if(power==7){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,6)*(1677283 - 159600*pow(M_PI,2) + 2280*log(1 - x)*(-283 + 420*log(1 - x))))/121500.;
	}
	if(power==8){
		return pow(alphas_muR/M_PI,2)*(pow(1 - x,7)*(634399147 - 60740400*pow(M_PI,2) + 840*log(1 - x)*(-242239 + 433860*log(1 - x))))/4.16745e7;
	}
	if(power==9){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,8)*(16.541672983479106 - (1499*pow(M_PI,2))/945. + (log(1 - x)*(-431723 + 944370*log(1 - x)))/99225.);
	}
	if(power==10){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,9)*(-2964401983 + 283827600*pow(M_PI,2) + 625905840*log(1 - x) - 1702965600*pow(log(1 - x),2)))/1.66698e8;
	}
	if(power==-1){
		return higgs_NNLO_qqp_reg(x)-pow(alphas_muR/M_PI,2)*(((-1 + x)*(-16921196098 + 86790453185*x - 220919323655*pow(x,2) + 346206684465*pow(x,3) - 355632299563*pow(x,4) + 241440062715*pow(x,5) - 104843098733*pow(x,6) + 26472679667*pow(x,7) - 2964401983*pow(x,8) + 58800*pow(M_PI,2)*(27372 - 140205*x + 358675*pow(x,2) - 563365*pow(x,3) + 579147*pow(x,4) - 393265*pow(x,5) + 170767*pow(x,6) - 43113*pow(x,7) + 4827*pow(x,8)) + 1680*(3860578 - 16673385*x + 37505035*pow(x,2) - 53695375*pow(x,3) + 51644433*pow(x,4) - 33363525*pow(x,5) + 13938303*pow(x,6) - 3412227*pow(x,7) + 372563*pow(x,8))*log(1 - x) - 352800*(27372 - 140205*x + 358675*pow(x,2) - 563365*pow(x,3) + 579147*pow(x,4) - 393265*pow(x,5) + 170767*pow(x,6) - 43113*pow(x,7) + 4827*pow(x,8))*pow(log(1 - x),2)))/1.66698e8);
	}
	if(power==-2){
		return pow(alphas_muR/M_PI,2)*(((-1 + x)*(-16921196098 + 86790453185*x - 220919323655*pow(x,2) + 346206684465*pow(x,3) - 355632299563*pow(x,4) + 241440062715*pow(x,5) - 104843098733*pow(x,6) + 26472679667*pow(x,7) - 2964401983*pow(x,8) + 58800*pow(M_PI,2)*(27372 - 140205*x + 358675*pow(x,2) - 563365*pow(x,3) + 579147*pow(x,4) - 393265*pow(x,5) + 170767*pow(x,6) - 43113*pow(x,7) + 4827*pow(x,8)) + 1680*(3860578 - 16673385*x + 37505035*pow(x,2) - 53695375*pow(x,3) + 51644433*pow(x,4) - 33363525*pow(x,5) + 13938303*pow(x,6) - 3412227*pow(x,7) + 372563*pow(x,8))*log(1 - x) - 352800*(27372 - 140205*x + 358675*pow(x,2) - 563365*pow(x,3) + 579147*pow(x,4) - 393265*pow(x,5) + 170767*pow(x,6) - 43113*pow(x,7) + 4827*pow(x,8))*pow(log(1 - x),2)))/1.66698e8);}
	else{return 0;}
}

///////////////////////////
/// qqbar channel
///////////////////////////
//NNLO coefficients  qqbar channel
//klopt met code
double higgs_NNLO_qqbar_reg(double x){
	return 0.*logdep_qqbar(x)+pow(alphas_muR/M_PI,2)*
		(1./x*((-16./9.-16./9.*x-8./9.*pow(x,2))*Li3(-x)+(-16./27*pow(x,2)-32./27.-32./27.*x)*S12(-x)
				+32./9.*pow(x+2.,2)*Li3(x)-32./9.*pow(x+2.,2)*S12(x)-4./27.*pow(x+2.,2)*pow(log(x),3)+4./9.*(2.+2.*x+pow(x,2))*log(1.+x)*pow(log(x),2)
				+(-8./27.*(2.+2.*x+pow(x,2))*pow(log(1.+x),2)-8./3.*pow(x+2.,2)*li2(x)+8./9.*(2.+2.*x+pow(x,2))*li2(-x))*log(x)
				-16./27.*(2.+2.*x+pow(x,2))*li2(-x)*log(1.+x)+32./81.*(1.-x)*(13.*pow(x,2)-35.*x-14.)*pow(log(1.-x),2)
				+-16./81.*(1.-x)*(37.*pow(x,2)-101.*x-44.)*log(x)*log(1.-x)-8./81.*(44.*pow(x,3)+39.*x-81.*pow(x,2)+27.)*pow(log(x),2)
				+16./27.*x*(x+6.*pow(x,2)+2.)*log(1.+x)*log(x)+8./81.*(42.*x-87.*pow(x,2)+12.+10.*pow(x,3))*li2(x)
				+16./27.*x*(x+6.*pow(x,2)+2.)*li2(-x)-4./81.*(1.-x)*(384.*pow(x,2)-967.*x-75.)*log(1.-x)+(-16./27.*pow(x,2)-32./27.-32./27.*x)*zeta3
				+(8./9.*pow(x+2.,2)*zeta2+4222./81*pow(x,2)-2896./81.*x-512./27.*pow(x,3)-10./3.)*log(x)-8./27.*(2.+2.*x+pow(x,2))*zeta2*log(1.+x)
				+(752./81.*pow(x,3)-544./27.*pow(x,2)+80./81.+400./27.*x)*zeta2+4./81.*(1.-x)*(783.*pow(x,2)-1925.*x+373.))
		+nF*1./x*(32./81.*pow(1.-x,3)*log(1.-x)+(-64./27.*pow(x,2)+64./81.*pow(x,3)-16./27.+80./27.*x)*log(x)-8./243.*(1.-x)*(41.*pow(x,2)-88.*x+23.))
		);

}


double logdep_qqbar(double x){
	return pow(alphas_muR/M_PI,2)*((pow(-1 + pow(CA,2),2)*(864*log(muF2/Q2)*CA+ 864*log(muF2/Q2)*CA*x + 216*log(muF2/Q2)*CA*pow(x,2))*(-log(1.-x)*log(x)-li2(x)))/(864.*pow(CA,3)*x)
  + (pow(-1 + pow(CA,2),2)*( - 1377*log(muF2/Q2)*CA - 324*pow(log(muF2/Q2),2)*CA  + 396*log(muF2/Q2)*pow(CA,2)  - 72*log(muF2/Q2)*CA*nF  + 144*log(muF2/Q2)*x  + 972*log(muF2/Q2)*CA*x + 216*pow(log(muF2/Q2),2)*CA*x  - 1332*log(muF2/Q2)*pow(CA,2)*x + 216*log(muF2/Q2)*CA*nF*x - 144*log(muF2/Q2)*pow(x,2) + 405*log(muF2/Q2)*CA*pow(x,2) + 108*pow(log(muF2/Q2),2)*CA*pow(x,2) + 1332*log(muF2/Q2)*pow(CA,2)*pow(x,2)  - 216*log(muF2/Q2)*CA*nF*pow(x,2) - 396*log(muF2/Q2)*pow(CA,2)*pow(x,3) + 72*log(muF2/Q2)*CA*nF*pow(x,3) + 864*log(muF2/Q2)*CA*zeta2 + 864*log(muF2/Q2)*CA*x*zeta2  + 216*log(muF2/Q2)*CA*pow(x,2)*zeta2))/(864.*pow(CA,3)*x)
 +log(x)*( -(pow(-1 + pow(CA,2),2)*(12*pow(CA,2)*( (6*log(muF2/Q2))*x - 2*(3*log(muF2/Q2))*pow(x,2) + (4*log(muF2/Q2))*pow(x,3)) - 6*x*(4*log(muF2/Q2)*(3 - 3*x + 2*pow(x,2))) + CA*(18*pow(log(muF2/Q2),2)*pow(2 + x,2) - 18*log(muF2/Q2)*(-12 + 8*x + 5*pow(x,2)))))/(288.*pow(CA,3)*x))
 +pow(log(x),2)*(-(pow(-1 + pow(CA,2),2)*( 3*CA*(log(muF2/Q2)*pow(2 + x,2))))/(48.*pow(CA,3)*x))
 +
 log(1.-x)*(
 -(pow(-1 + pow(CA,2),2)*((1 - x)*(12*pow(CA,2)*((2*log(muF2/Q2))*pow(-1 + x,2)) - 12*(2*log(muF2/Q2)*pow(-1 + x,2)) - CA*(72*log(muF2/Q2)*(3 + x)))))/(144.*pow(CA,3)*x)
 + (pow(-1 + pow(CA,2),2)*(2*log(muF2/Q2)*pow(2 + x,2)))*log(x)/(24.*pow(CA,3)*x)));
}

// checked
double higgs_NNLO_qqbar_expansion(double x, int power){
	if(power==1){
		return 0;
	}
	if(power==2){
		return pow(alphas_muR/M_PI,2)*((1 - x)*(60 - 8*pow(M_PI,2) + 48*(-1 + log(1 - x))*log(1 - x)))/27.;
	}
	if(power==3){
		return pow(alphas_muR/M_PI,2)*(-2*pow(-1 + x,2)*(-83 + 6*pow(M_PI,2) + 36*log(1 - x) - 36*pow(log(1 - x),2)))/27.;
	}
	if(power==4){
		return pow(alphas_muR/M_PI,2)*(4*pow(1 - x,3)*(10972 - 240*nF - 393*pow(M_PI,2) + 6*log(1 - x)*(-833 + 12*nF + 291*log(1 - x))))/729.;
	}
	if(power==5){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,4)*(33835 - 384*nF - 1734*pow(M_PI,2) + 6*log(1 - x)*(-2819 + 48*nF + 1326*log(1 - x))))/729.;
	}
	if(power==6){
		return pow(alphas_muR/M_PI,2)*(pow(1 - x,5)*(7689391 - 25800*nF - 467700*pow(M_PI,2) + 60*log(1 - x)*(-64043 + 1200*nF + 36570*log(1 - x))))/182250.;
	}
	if(power==7){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,6)*(4890253 + 13200*nF - 331600*pow(M_PI,2) + 40*log(1 - x)*(-59171 + 1200*nF + 39540*log(1 - x))))/121500.;
	}
	if(power==8){
		return pow(alphas_muR/M_PI,2)*(pow(1 - x,7)*(1634627677 + 12171600*nF - 119736400*pow(M_PI,2) + 280*log(1 - x)*(-2691197 + 58800*nF + 2065980*log(1 - x))))/4.16745e7;
	}
	if(power==9){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,8)*(6445572073 + 72676800*nF - 500407600*pow(M_PI,2) + 560*log(1 - x)*(-4999009 + 117600*nF + 4361910*log(1 - x))))/1.66698e8;
	}
	if(power==10){
		return pow(alphas_muR/M_PI,2)*(pow(1 - x,9)*(172752517031 + 2493237600*nF - 14034913200*pow(M_PI,2) + 5040*log(1 - x)*(-13912987 + 352800*nF + 13709430*log(1 - x))))/4.500846e9;
	}
	if(power==-1){
		return higgs_NNLO_qqbar_reg(x)-pow(alphas_muR/M_PI,2)*(((1 - x)*(166698000*(60 - 8*pow(M_PI,2) + 48*(-1 + log(1 - x))*log(1 - x)) + 333396000*(-1 + x)*(-83 + 6*pow(M_PI,2) + 36*log(1 - x) - 36*pow(log(1 - x),2)) + 24696000*pow(-1 + x,2)*(10972 - 240*nF - 393*pow(M_PI,2) + 6*log(1 - x)*(-833 + 12*nF + 291*log(1 - x))) + 6174000*pow(1 - x,3)*(33835 - 384*nF - 1734*pow(M_PI,2) + 6*log(1 - x)*(-2819 + 48*nF + 1326*log(1 - x))) + 24696*pow(-1 + x,4)*(7689391 - 25800*nF - 467700*pow(M_PI,2) + 60*log(1 - x)*(-64043 + 1200*nF + 36570*log(1 - x))) + 37044*pow(1 - x,5)*(4890253 + 13200*nF - 331600*pow(M_PI,2) + 40*log(1 - x)*(-59171 + 1200*nF + 39540*log(1 - x))) + 108*pow(-1 + x,6)*(1634627677 + 12171600*nF - 119736400*pow(M_PI,2) + 280*log(1 - x)*(-2691197 + 58800*nF + 2065980*log(1 - x))) + 27*pow(1 - x,7)*(6445572073 + 72676800*nF - 500407600*pow(M_PI,2) + 560*log(1 - x)*(-4999009 + 117600*nF + 4361910*log(1 - x))) + pow(-1 + x,8)*(172752517031 + 2493237600*nF - 14034913200*pow(M_PI,2) + 5040*log(1 - x)*(-13912987 + 352800*nF + 13709430*log(1 - x)))))/4.500846e9);
	}
	if(power==-2){
		return pow(alphas_muR/M_PI,2)*(((1 - x)*(166698000*(60 - 8*pow(M_PI,2) + 48*(-1 + log(1 - x))*log(1 - x)) + 333396000*(-1 + x)*(-83 + 6*pow(M_PI,2) + 36*log(1 - x) - 36*pow(log(1 - x),2)) + 24696000*pow(-1 + x,2)*(10972 - 240*nF - 393*pow(M_PI,2) + 6*log(1 - x)*(-833 + 12*nF + 291*log(1 - x))) + 6174000*pow(1 - x,3)*(33835 - 384*nF - 1734*pow(M_PI,2) + 6*log(1 - x)*(-2819 + 48*nF + 1326*log(1 - x))) + 24696*pow(-1 + x,4)*(7689391 - 25800*nF - 467700*pow(M_PI,2) + 60*log(1 - x)*(-64043 + 1200*nF + 36570*log(1 - x))) + 37044*pow(1 - x,5)*(4890253 + 13200*nF - 331600*pow(M_PI,2) + 40*log(1 - x)*(-59171 + 1200*nF + 39540*log(1 - x))) + 108*pow(-1 + x,6)*(1634627677 + 12171600*nF - 119736400*pow(M_PI,2) + 280*log(1 - x)*(-2691197 + 58800*nF + 2065980*log(1 - x))) + 27*pow(1 - x,7)*(6445572073 + 72676800*nF - 500407600*pow(M_PI,2) + 560*log(1 - x)*(-4999009 + 117600*nF + 4361910*log(1 - x))) + pow(-1 + x,8)*(172752517031 + 2493237600*nF - 14034913200*pow(M_PI,2) + 5040*log(1 - x)*(-13912987 + 352800*nF + 13709430*log(1 - x)))))/4.500846e9);
	}
	else{return 0;}
}
