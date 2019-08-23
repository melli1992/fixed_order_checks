#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "LHAPDF/LHAPDF.h"
#include "deriv_pdf.h"
#include "parameters.h"
using namespace std;

////////////////////////////////////////////////////////////////////////
///
/// this file contains all pdf related stuff, like sums of qqbar
/// and derivatives of pdfs
/// note that the xfxQ function of LHAPDF returns the x*f(x) value
///
////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////
/// this is the qqbar sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////

double pdf_sum_qqbar_charge_weighted(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		sum_pdf+= eq[i-1]*eq[i-1]*1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(i,x,muF)*pdfs[0]->xfxQ(-i,tau_over_z/x,muF)+pdfs[0]->xfxQ(i,tau_over_z/x,muF)*pdfs[0]->xfxQ(-i,x,muF));
		//cout << "in the sum " << sum_pdf << endl;
	}
	return 1./x*sum_pdf;
}
double pdf_sum_qqbar_charge_unweighted(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	double charge_factor = 0;
	for(int i = 1; i <=5; i++){
		charge_factor += eq[i-1]*eq[i-1];
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(i,x,muF)*pdfs[0]->xfxQ(-i,tau_over_z/x,muF)+pdfs[0]->xfxQ(i,tau_over_z/x,muF)*pdfs[0]->xfxQ(-i,x,muF));
	}
	return 1./x*charge_factor*sum_pdf;
}
//qq + qbarqbar +qqbar (identical+non-identical!), use it for DeltaNNLOC2
double pdf_sum_qq_charge_weighted_double(double x, double tau_over_z){ 
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = -5; i <=5; i++){
		for(int j = -5; j <=5; j++){
			if(j==0){continue;}
			if(i==0){continue;}
			sum_pdf+= eq[abs(i)-1]*eq[abs(i)-1]*1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(i,x,muF)*pdfs[0]->xfxQ(j,tau_over_z/x,muF)+pdfs[0]->xfxQ(i,tau_over_z/x,muF)*pdfs[0]->xfxQ(j,x,muF));
		}
	}
	return 1./x*sum_pdf;
}
//qq + qbarqbar +qqbar (identical+non-identical!), use it for DeltaNNLOCE
double pdf_sum_qq_charge_weighted_single(double x, double tau_over_z){ 
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = -5; i <=5; i++){
		if(i==0){continue;}
		sum_pdf+= eq[abs(i)-1]*eq[abs(i)-1]*1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(i,x,muF)*pdfs[0]->xfxQ(i,tau_over_z/x,muF)+pdfs[0]->xfxQ(i,tau_over_z/x,muF)*pdfs[0]->xfxQ(i,x,muF));
	}
	return 1./x*sum_pdf;
}


//qq + qbarqbar + qqbar (identical+nonidentical!), use it for DeltaNNLOCD
double pdf_sum_qq_charge_weighted_double_vivj(double x, double tau_over_z){ 
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = -5; i <=5; i++){
		for(int j = -5; j <=5; j++){
			if(j==0){continue;}
			if(i==0){continue;}
			sum_pdf+= eq[abs(i)-1]*eq[abs(j)-1]*1./x*1./(tau_over_z/x)*(-(float) i*j)/((float)abs(i*j))*(pdfs[0]->xfxQ(i,x,muF)*pdfs[0]->xfxQ(j,tau_over_z/x,muF));
		}
	}
	return 1./x*sum_pdf;
}


//qq + qbarqbar + qqbar (identical only), use it for DeltaNNLOCF
double pdf_sum_qq_charge_weighted_single_vivi(double x, double tau_over_z){ 
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = -5; i <=5; i++){
		if(i==0){continue;}
		sum_pdf+= eq[abs(i)-1]*eq[abs(i)-1]*1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(i,x,muF)*pdfs[0]->xfxQ(i,tau_over_z/x,muF));
	}
	return 1./x*sum_pdf;
}

/////////////////////////////////////////////////////////////////////////
/// this is the qg sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////

double pdf_sum_qg_charge_weighted(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		sum_pdf+= eq[i-1]*eq[i-1]*1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(i,x,muF)*pdfs[0]->xfxQ(0,tau_over_z/x,muF)+pdfs[0]->xfxQ(i,tau_over_z/x,muF)*pdfs[0]->xfxQ(0,x,muF));
		sum_pdf+= eq[i-1]*eq[i-1]*1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(-i,x,muF)*pdfs[0]->xfxQ(0,tau_over_z/x,muF)+pdfs[0]->xfxQ(-i,tau_over_z/x,muF)*pdfs[0]->xfxQ(0,x,muF));
	}
	return 1./x*sum_pdf;
}


/////////////////////////////////////////////////////////////////////////
/// this is the gg sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////

double pdf_sum_gg_charge_weighted(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	double charge_factor = 0;
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		charge_factor += eq[i-1]*eq[i-1]; //still need the charge factor
			}
	sum_pdf = 1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(0,x,muF)*pdfs[0]->xfxQ(0,tau_over_z/x,muF));
	return 1./x*sum_pdf*charge_factor;
}


///////////////////////////////////////
/// HIGGS STUFF
//////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/// this is the qqbar sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////
double pdf_sum_qqbar(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(i,x,muF)*pdfs[0]->xfxQ(-i,tau_over_z/x,muF)+pdfs[0]->xfxQ(i,tau_over_z/x,muF)*pdfs[0]->xfxQ(-i,x,muF));
		//cout << "in the sum " << sum_pdf << endl;
	}
	return 1./x*sum_pdf;
}
/////////////////////////////////////////////////////////////////////////
/// this is the qg sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////
double pdf_sum_qg(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(i,x,muF)*pdfs[0]->xfxQ(0,tau_over_z/x,muF)+pdfs[0]->xfxQ(i,tau_over_z/x,muF)*pdfs[0]->xfxQ(0,x,muF));
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(-i,x,muF)*pdfs[0]->xfxQ(0,tau_over_z/x,muF)+pdfs[0]->xfxQ(-i,tau_over_z/x,muF)*pdfs[0]->xfxQ(0,x,muF));
	}
	return 1./x*sum_pdf;
}
/////////////////////////////////////////////////////////////////////////
/// this is the gg sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////
double pdf_sum_gg(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	sum_pdf = 1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(0,x,muF)*pdfs[0]->xfxQ(0,tau_over_z/x,muF));
	return 1./x*sum_pdf;
}
/////////////////////////////////////////////////////////////////////////
/// this is the qq (+qbarqbar) sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////
double pdf_sum_qq(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(i,x,muF)*pdfs[0]->xfxQ(i,tau_over_z/x,muF)+pdfs[0]->xfxQ(-i,x,muF)*pdfs[0]->xfxQ(-i,tau_over_z/x,muF));
		//cout << "in the sum " << sum_pdf << endl;
	}
	return 1./x*sum_pdf;
}
/////////////////////////////////////////////////////////////////////////
/// this is the qq (+qbarqbar) sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////
double pdf_sum_qqNI(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		for(int j = 1; j <=5; j++){
			if(i>=j){continue;} //to avoid double counting!
			sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(i,x,muF)*pdfs[0]->xfxQ(j,tau_over_z/x,muF)+pdfs[0]->xfxQ(i,tau_over_z/x,muF)*pdfs[0]->xfxQ(j,x,muF));
			sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(-i,x,muF)*pdfs[0]->xfxQ(-j,tau_over_z/x,muF)+pdfs[0]->xfxQ(-i,tau_over_z/x,muF)*pdfs[0]->xfxQ(-j,x,muF));
		}
	}
	return 1./x*sum_pdf;
}

///////////////////////////////////////////////////////
/// numerical derivative of the weighted pdf sum for qg
///////////////////////////////////////////////////////
double derivative_qg_pdf(double x, double z, double tau, double eps){
	double tau_over_zpeps(tau/(z+eps)); //this is the argument of z + eps
	double tau_over_zmeps(tau/(z-eps)); //argument z - eps
	return (1./(z+eps)*pdf_sum_qg_charge_weighted(x,tau_over_zpeps)-1./(z-eps)*pdf_sum_qg_charge_weighted(x,tau_over_zmeps))/(2.*eps); //numerical derivative
	//return (jac_zpeps/(z+eps)*pdf_sum_charge_weighted(x,tau_over_zpeps)-jac_zmeps/(z)*pdf_sum_charge_weighted(x,tau/z))/(eps); //numerical derivative
}


////////////////////////////////////////////////////////////////////////////////////////
/// numerical derivative of the weighted pdf sum including a z-dependent jacobean for qg
////////////////////////////////////////////////////////////////////////////////////////
double derivative_qg_pdf_jac(double x, double z, double tau, double eps){

	double zpeps = z + eps;
	double zmeps = z - eps;
	double jacp = 1.-tau/zpeps, jacm = 1.-tau/zmeps;
	double xpeps = tau/zpeps+x*jacp, xmeps = tau/zmeps+x*jacm;

	if(xmeps < tau/(zpeps)){return 0;}
	if(xpeps < tau/(zmeps)){return 0;}
	return (jacp/(zpeps)*pdf_sum_qg_charge_weighted(xpeps,tau/(zpeps))-jacm/(zmeps)*pdf_sum_qg_charge_weighted(xmeps,tau/(zmeps)))/(2.*eps); //numerical derivative
}

///////////////////////////////////
/// checks for conservation charge
///////////////////////////////////
double vegas_pdf_up_minus_upbar(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return 1./k[0]*(pdfs[0]->xfxQ(2,k[0],muF)-pdfs[0]->xfxQ(-2,k[0],muF));

}

///////////////////////////////////
/// checks for momentum cons.
///////////////////////////////////
double vegas_pdf_mom_consv(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double sum_pdf(0);
	for(int i = 1; i <=5; i++){
		sum_pdf+= (pdfs[0]->xfxQ(i,k[0],muF)+pdfs[0]->xfxQ(-i,k[0],muF));
	}
	sum_pdf+=pdfs[0]->xfxQ(21,k[0],muF);
	return sum_pdf;
}
