#include <bits/stdc++.h>
#include <stdlib.h>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "LHAPDF/LHAPDF.h"
#include "mellin_pdf.h"
#include "parameters.h"
#include "chebychev.h"
using namespace std;

//////////////////////////////////////////////////////////
/// contains the N space pdfs (2 methods)
///
/// we have the derivative method first
/// only first derivative is implemented
/// second derivative could be implemented as well
///
/// then we have the contour deformation to fitted PDFs
///
/// and we have the exact mellin space expression (need to implement still)
///
/// end of code contains chebychev fits, probably obsolete
///
//////////////////////////////////////////////////////////

//returns mellin space pdfs with 1/N^2 suppression via derivative method
complex<double> mellin_pdf_sum_qqbar_charge_weighted(double x1, double x2, complex<double> N){
	complex<double> sum_pdf(0,0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	for(int i = 1; i <=5; i++){
		sum_pdf+= eq[i-1]*eq[i-1]*1./pow(N,2)*pow(x1,N)*pow(x2,N)*(deriv_pdf(i,x1)*deriv_pdf(-i,x2)+deriv_pdf(i,x2)*deriv_pdf(-i,x1));
	}
	return sum_pdf;
}
double deriv_pdf(int i, double x, double eps){
	try{
		return ((pdfs[0]->xfxQ(i,x+eps,muF)/(x+eps))-(pdfs[0]->xfxQ(i,x-eps,muF)/(x-eps)))/(2.*eps);
	}
	catch (exception& e)
  {
    return 0;
  }
}

// test cases, these work, so numerical derivative is working
complex<double> mellin_test_pdf(double x1, double x2, complex<double> N){
	complex<double> sum_pdf(0,0);
	sum_pdf= 1./pow(N,2)*pow(x1,N)*pow(x2,N)*(test_deriv_pdf(x1)*test_deriv_pdf(x2));
	//sum_pdf= 1./pow(N,2)*pow(x1,N)*pow(x2,N)*(2.*x1*2.*x2);
	return sum_pdf;
}
double test_deriv_pdf(double x){
	double eps = 1.0E-5;
	if(((x+eps) > 1) or (x-eps < 0)){return 0;}
	else{return ((1.-pow(x+eps,2))-(1.-pow(x-eps,2)))/(2.*eps);}
}

/////////////////////////////////////////////////////
/// here is the code for the pdf fits
/// first we have the sum over all qqbar pairs
/// weigthed with the right normalization
/////////////////////////////////////////////////////
complex<double> fit_sum_qqbar_charge_weighted(complex<double> x1, complex<double> x2){
	complex<double> sum_pdf(0,0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	for(int i = 1; i <=5; i++){
		sum_pdf+= eq[i-1]*eq[i-1]*(fit_pdfs(5-i,x1)*fit_pdfs(5+i,x2)+fit_pdfs(5-i,x2)*fit_pdfs(5+i,x1));
	}
	return sum_pdf;
}
/// data structure that contains the fit coefficients
/// where i_coeff is a vector of elements, the elements are
/// data is stored like: {Qvalue, {bbar_coeff,cbar_coeff, sbar_coeff, ubar_coeff, dbar_coeff, g_coeff, d_coeff, u_coeff, s_coeff, c_coeff, b_coeff}}
std::unordered_map<double, vector<vector<double>>> fitcoeff = {
	{500,{{-1.2456100538,4.68897263286,-0.229793671739,0.843277587332,-0.158485801005,1.84621999137,5.692834127,-0.264099460625},{1.07402386079,6.76416560887,-0.266830011363,0.73563877183,-0.22857129335,-1.83997495695,7.78614109524,-0.225759159066},{0.251585742709,12.3621799749,-0.600174552326,-1.90306193949,0.904531006276,0.448216886868,16.6952673976,-0.293261776161},{1.41876010925,32.6996356528,-0.354479783283,-1.26371464848,0.397091748733,0.295343418204,9.06751475566,0.301948007994},{1.91020817555,6.35433712404,-0.323191939068,-0.756509180523,-0.0948132333114,-3.60943551453,6.39431343751,0.0825836279415},{-373.902422772,8.6481851032,-0.26428614758,0.843956145272,-0.204448010232,605.584363446,9.68471399425,-0.267290887874},{-0.546996119389,3.55165636983,-0.116292539393,-0.214583377237,0.530501004597,0.41704018547,3.12499332697,-0.302909501662},{0.334811915649,3.68827641672,-0.0247779747314,-1.73701373992,-1.99955578114,0.32830002086,15.2749119505,-0.312984968907},{-0.264407097187,2.95067243312,-0.0196211485701,1.45796491562,0.438374097732,0.272578776387,6.05118501114,-0.326817166142},{-0.0588000833014,5.39090752094,-0.618912875714,0.0885724491421,-1.08512810665,0.262670725013,6.40491840769,-0.335258374673},{-1.20691640506,4.69563014588,-0.229141056364,0.844776860091,-0.15708987512,1.78431928442,5.69975832033,-0.264587956825}}},
	{100,{{-1.2456100538,4.68897263286,-0.229793671739,0.843277587332,-0.158485801005,1.84621999137,5.692834127,-0.264099460625},{1.07402386079,6.76416560887,-0.266830011363,0.73563877183,-0.22857129335,-1.83997495695,7.78614109524,-0.225759159066},{0.251585742709,12.3621799749,-0.600174552326,-1.90306193949,0.904531006276,0.448216886868,16.6952673976,-0.293261776161},{1.41876010925,32.6996356528,-0.354479783283,-1.26371464848,0.397091748733,0.295343418204,9.06751475566,0.301948007994},{1.91020817555,6.35433712404,-0.323191939068,-0.756509180523,-0.0948132333114,-3.60943551453,6.39431343751,0.0825836279415},{-373.902422772,8.6481851032,-0.26428614758,0.843956145272,-0.204448010232,605.584363446,9.68471399425,-0.267290887874},{-0.546996119389,3.55165636983,-0.116292539393,-0.214583377237,0.530501004597,0.41704018547,3.12499332697,-0.302909501662},{0.334811915649,3.68827641672,-0.0247779747314,-1.73701373992,-1.99955578114,0.32830002086,15.2749119505,-0.312984968907},{-0.264407097187,2.95067243312,-0.0196211485701,1.45796491562,0.438374097732,0.272578776387,6.05118501114,-0.326817166142},{-0.0588000833014,5.39090752094,-0.618912875714,0.0885724491421,-1.08512810665,0.262670725013,6.40491840769,-0.335258374673},{-1.20691640506,4.69563014588,-0.229141056364,0.844776860091,-0.15708987512,1.78431928442,5.69975832033,-0.264587956825}}}
};
//this is xfx(x)
complex<double> xfit_pdfs(int i, complex<double> x){
	complex<double> y = 1.-2.*pow(x,0.5);
	return fitcoeff[Q][i][0]*pow(1.-x,fitcoeff[Q][i][1])*pow(x,fitcoeff[Q][i][2])*(1.+fitcoeff[Q][i][3]*y+fitcoeff[Q][i][4]*(2.*pow(y,2)-1.))+fitcoeff[Q][i][5]*pow(1.-x,fitcoeff[Q][i][6])*pow(x,fitcoeff[Q][i][7]);
}
//this is fx(x)
complex<double> fit_pdfs(int i, complex<double> x){
	complex<double> y = 1.-2.*pow(x,0.5);
	return 1./x*(fitcoeff[Q][i][0]*pow(1.-x,fitcoeff[Q][i][1])*pow(x,fitcoeff[Q][i][2])*(1.+fitcoeff[Q][i][3]*y+fitcoeff[Q][i][4]*(2.*pow(y,2)-1.))+fitcoeff[Q][i][5]*pow(1.-x,fitcoeff[Q][i][6])*pow(x,fitcoeff[Q][i][7]));
}

/////////////////////////////////////////////////////////////////////////////////
/// routines to calculate the chebychev polynomials (used only once for each muF2)
/// code became obsolete probably
/////////////////////////////////////////////////////////////////////////////////
using namespace Chebyshev ;
double vegas_coefficients_0(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * par =  (struct lumni_params *)params;
	if(k[0] <= 2.*1E-9-1.){return 0;}
	if(k[0] >= 1.){return 0;}
	double result = a0coeff(par->flavor,k[0]);
	if(isnan(result)){return 0;}
	else{return result;};
}
double vegas_coefficients_n(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * par =  (struct lumni_params *)params;
	if(k[0] <= 2.*1E-9-1.){return 0;}
	if(k[0] >= 1.){return 0;}
	double result = ancoeff(par->flavor,par->coefficient,k[0]);
	if(isnan(result)){return 0;}
	else{return result;};
}
double a0coeff(int i,double x){
	//return 1./M_PI*(pdfs[0]->xfxQ(i,1./2.*x+1./2.,muF))*Tn<double>(0, x)*pow(1.-pow(x,2),-0.5);
	return 1./M_PI*(pdfs[0]->xfxQ(i,pow(10.,4.5*x-4.5),muF))*Tn<double>(0, x)*pow(1.-pow(x,2),-0.5);
}
double ancoeff(int i,int n,double x){
	return 2./M_PI*(pdfs[0]->xfxQ(i,pow(10.,4.5*x-4.5),muF))*Tn<double>(n, x)*pow(1.-pow(x,2),-0.5);
}
