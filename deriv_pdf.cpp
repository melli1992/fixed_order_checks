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
		sum_pdf+= eq[i-1]*eq[i-1]*1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
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
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
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
			sum_pdf+= eq[abs(i)-1]*eq[abs(i)-1]*1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(j,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(j,x,muF));
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
		sum_pdf+= eq[abs(i)-1]*eq[abs(i)-1]*1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(i,x,muF));
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
			sum_pdf+= eq[abs(i)-1]*eq[abs(j)-1]*1./x*1./(tau_over_z/x)*(-(float) i*j)/((float)abs(i*j))*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(j,tau_over_z/x,muF));
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
		sum_pdf+= eq[abs(i)-1]*eq[abs(i)-1]*1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(i,tau_over_z/x,muF));
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
		sum_pdf+= eq[i-1]*eq[i-1]*1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(0,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(0,x,muF));
		sum_pdf+= eq[i-1]*eq[i-1]*1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(-i,x,muF)*pdfs[use_member]->xfxQ(0,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(0,x,muF));
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
	sum_pdf = 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(0,x,muF)*pdfs[use_member]->xfxQ(0,tau_over_z/x,muF));
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
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
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
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(0,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(0,x,muF));
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(-i,x,muF)*pdfs[use_member]->xfxQ(0,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(0,x,muF));
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
	sum_pdf = 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(0,x,muF)*pdfs[use_member]->xfxQ(0,tau_over_z/x,muF));
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
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(-i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF));
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
	//int l = 0;
	for(int i = -5; i <=5; i++){
		for(int j = -5; j <=5; j++){
			if(abs(i)==abs(j)){continue;}
			if(i==0){continue;}
			if(j==0){continue;}
			//l += 1;
			//cout << "f_i(x1) i=" << i << ", f_j(x2), j=" << j << endl;//to avoid double counting!
			sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(j,tau_over_z/x,muF));
			//sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(-i,x,muF)*pdfs[use_member]->xfxQ(-j,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-j,x,muF));
			//sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(-i,x,muF)*pdfs[use_member]->xfxQ(j,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-j,x,muF));
		}
	}
	//cout << l << endl;
	return 1./x*sum_pdf;
}

//////////////////////////////////////////////////////
/// relevant for W+W-
//////////////////////////////////////////////////////
double pdf_sum_qqbarUP(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	int i = 2;
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
	i = 4;
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
	return 1./x*sum_pdf;
}
double pdf_sum_qqbarDOWN(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	int i = 1;
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
	i = 3;
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
	i = 5;
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
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



////////////////////////////////////////////////////////////////////////////////////////
/// numerical derivative of gg (for di-higgs)
////////////////////////////////////////////////////////////////////////////////////////
double derivative_gg_pdf(double x, double tau, double eps){
	eps = 1.E-6;
	double taupeps = tau + eps;
	double taumeps = tau - eps;
	return (pdf_sum_gg(x, taupeps) - pdf_sum_gg(x,taumeps))/(2.*eps);
}


///////////////////////////////////
/// checks for conservation charge
///////////////////////////////////
double vegas_pdf_up_minus_upbar(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return 1./k[0]*(pdfs[use_member]->xfxQ(2,k[0],muF)-pdfs[use_member]->xfxQ(-2,k[0],muF));

}

///////////////////////////////////
/// checks for momentum cons.
///////////////////////////////////
double vegas_pdf_mom_consv(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double sum_pdf(0);
	for(int i = 1; i <=5; i++){
		sum_pdf+= (pdfs[use_member]->xfxQ(i,k[0],muF)+pdfs[use_member]->xfxQ(-i,k[0],muF));
	}
	sum_pdf+=pdfs[use_member]->xfxQ(21,k[0],muF);
	return sum_pdf;
}

/////////////////////////////////////
/// lumi checks
/////////////////////////////////////
double vegas_lumi_gg(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return pdf_sum_gg(k[0],tau);
}
double vegas_lumi_qg(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return pdf_sum_qg(k[0],tau);
}
double vegas_lumi_qqbar(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return pdf_sum_qqbar(k[0],tau);
}
double vegas_lumi_qq(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return pdf_sum_qq(k[0],tau);
}
double vegas_lumi_qqNI(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return pdf_sum_qqNI(k[0],tau);
}
