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

double pdf_sum_qq_charge_weighted(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		sum_pdf+= eq[i-1]*eq[i-1]*1./x*1./(tau_over_z/x)*(pdfs[0]->xfxQ(i,x,muF)*pdfs[0]->xfxQ(-i,tau_over_z/x,muF)+pdfs[0]->xfxQ(i,tau_over_z/x,muF)*pdfs[0]->xfxQ(-i,x,muF));
		//cout << "in the sum " << sum_pdf << endl;
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

///////////////////////////////////////////////////////
/// numerical derivative of the weighted pdf sum for qq
///////////////////////////////////////////////////////
double derivative_qq_pdf(double x, double z, double tau, double eps){
	double tau_over_zpeps(tau/(z+eps)); //this is the argument of z + eps
	double tau_over_zmeps(tau/(z-eps)); //argument z - eps
	return (1./(z+eps)*pdf_sum_qq_charge_weighted(x,tau_over_zpeps)-1./(z-eps)*pdf_sum_qq_charge_weighted(x,tau_over_zmeps))/(2.*eps); //numerical derivative
	//return (jac_zpeps/(z+eps)*pdf_sum_charge_weighted(x,tau_over_zpeps)-jac_zmeps/(z)*pdf_sum_charge_weighted(x,tau/z))/(eps); //numerical derivative
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
/// numerical derivative of the weighted pdf sum including a z-dependent jacobean for qq
////////////////////////////////////////////////////////////////////////////////////////
double derivative_qq_pdf_jac(double x, double z, double tau, double eps){
	
	double zpeps = z + eps;
	double zmeps = z - eps;
	double jacp = 1.-tau/zpeps, jacm = 1.-tau/zmeps;
	double xpeps = tau/zpeps+x*jacp, xmeps = tau/zmeps+x*jacm;
	
	if(xmeps < tau/(zpeps)){return 0;}
	if(xpeps < tau/(zmeps)){return 0;}
	return (jacp/(zpeps)*pdf_sum_qq_charge_weighted(xpeps,tau/(zpeps))-jacm/(zmeps)*pdf_sum_qq_charge_weighted(xmeps,tau/(zmeps)))/(2.*eps); //numerical derivative
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
//////////////////////////////////////////////////////////////////
/// numerical derivative of the weighted pdf sum - usage for vegas
//////////////////////////////////////////////////////////////////
double vegas_sum_pdf_weigthed(double *k, size_t dim, void *params){
	(void)(dim);
	lumni_params * lp = (struct lumni_params *)params;
    // fail saves to catch errors of out of bounds pdfs -> commented because this makes it extremely slow
	//if((k[0] > xmax_pdfs) || (k[0] < xmin_pdfs)){return 0;}
	//if((tau/(lp->z)/k[0] < xmax_pdfs) || (tau/(lp->z)/k[0] < xmin_pdfs)){return 0;}
	return pdf_sum_qq_charge_weighted(k[0], tau/(lp->z));
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
