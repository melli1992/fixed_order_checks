#include <bits/stdc++.h>
#include <stdlib.h>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "LHAPDF/LHAPDF.h"
#include "mellin_pdf.h"
#include "polygamma.h"
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


//////////////////////////
/// the derivative method
/// returns mellin space pdfs with 1/N^2 suppression via derivative method
//////////////////////////
complex<double> mellin_pdf_sum_qqbar_charge_weighted(double x1, double x2, complex<double> N){
	complex<double> sum_pdf(0,0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; 
	for(int i = 1; i <=5; i++){
		sum_pdf+= eq[i-1]*eq[i-1]*1./pow(N,2)*pow(x1,N)*pow(x2,N)*(deriv_pdf(i,x1)*deriv_pdf(-i,x2)+deriv_pdf(i,x2)*deriv_pdf(-i,x1));
	}
	return sum_pdf;
}
// numerical derivative
double deriv_pdf(int i, double x, double eps){
	try{
		return ((pdfs[0]->xfxQ(i,x+eps,muF)/(x+eps))-(pdfs[0]->xfxQ(i,x-eps,muF)/(x-eps)))/(2.*eps);
	}
	catch (exception& e)
  {
    return 0;
  }
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
		sum_pdf+= eq[i-1]*eq[i-1]*1./x1*1./(x2/x1)*(xfit_pdfs(5-i,x1)*xfit_pdfs(5+i,x2/x1)+xfit_pdfs(5-i,x2/x1)*xfit_pdfs(5+i,x1));
	}
	return 1./x1*sum_pdf;
}
complex<double> fit_sum_qqbar_charge_unweighted(complex<double> x1, complex<double> x2){ //use for DeltaqqB2
	complex<double> sum_pdf(0,0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	double charge_factor = 0;
	for(int i = 1; i <=5; i++){
		charge_factor += eq[i-1]*eq[i-1];
		sum_pdf+= 1./x1*1./(x2/x1)*(xfit_pdfs(5-i,x1)*xfit_pdfs(5+i,x2/x1)+xfit_pdfs(5-i,x2/x1)*xfit_pdfs(5+i,x1));
	}
	return 1./x1*sum_pdf*charge_factor;
}
complex<double> fit_sum_qg_charge_weighted(complex<double> x1, complex<double> x2){
	complex<double> sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.};
	for(int i = 1; i <=5; i++){
		sum_pdf+= eq[i-1]*eq[i-1]*1./x1*1./(x2/x1)*(xfit_pdfs(5-i,x1)*xfit_pdfs(5,x2/x1)+xfit_pdfs(5,x1)*xfit_pdfs(5-i,x2/x1));
		sum_pdf+= eq[i-1]*eq[i-1]*1./x1*1./(x2/x1)*(xfit_pdfs(5+i,x1)*xfit_pdfs(5,x2/x1)+xfit_pdfs(5,x1)*xfit_pdfs(5+i,x2/x1));
	}
	return 1./x1*sum_pdf;
}

complex<double> fit_sum_gg_charge_weighted(complex<double> x1, complex<double> x2){
	complex<double> sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.};
	double charge_factor = 0;
	for(int i = 1; i <=5; i++){
		charge_factor += eq[i-1]*eq[i-1]; //still need the charge factor
			}
	sum_pdf = 1./x1*1./(x2/x1)*charge_factor*(xfit_pdfs(5,x1)*xfit_pdfs(5,x2/x1));
	return 1./x1*sum_pdf;
}
//qq + qbarqbar +qqbar (identical+non-identical!), use it for DeltaNNLOC2
complex<double> fit_sum_qq_charge_weighted_double(complex<double> x1, complex<double> x2){
	complex<double> sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	for(int i = -5; i <=5; i++){
		for(int j = -5; j <=5; j++){
			if(j==0){continue;}
			if(i==0){continue;}
			sum_pdf+= eq[abs(i)-1]*eq[abs(i)-1]*1./x1*1./(x2/x1)*(xfit_pdfs(5+i,x1)*xfit_pdfs(5+j,x2/x1)+xfit_pdfs(5+i,x2/x1)*xfit_pdfs(5+j,x1));
		}
	}
	return 1./x1*sum_pdf;
}
//qq + qbarqbar +qqbar (identical+non-identical!), use it for DeltaNNLOC2
complex<double> fit_sum_qq_charge_weighted_double_vivj(complex<double> x1, complex<double> x2){
	complex<double> sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	for(int i = -5; i <=5; i++){
		for(int j = -5; j <=5; j++){
			if(j==0){continue;}
			if(i==0){continue;}
			sum_pdf+= eq[abs(i)-1]*eq[abs(j)-1]*(-(float) i*j)/((float)abs(i*j))*1./x1*1./(x2/x1)*(xfit_pdfs(5+i,x1)*xfit_pdfs(5+j,x2/x1));
		}
	}
	return 1./x1*sum_pdf;
}
//qq + qbarqbar +qqbar (identical!), use it for DeltaNNLOC2
complex<double> fit_sum_qq_charge_weighted_single(complex<double> x1, complex<double> x2){
	complex<double> sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	for(int i = -5; i <=5; i++){
			if(i==0){continue;}
			sum_pdf+= eq[abs(i)-1]*eq[abs(i)-1]*1./x1*1./(x2/x1)*(xfit_pdfs(5+i,x1)*xfit_pdfs(5+i,x2/x1)+xfit_pdfs(5+i,x2/x1)*xfit_pdfs(5+i,x1));
	}
	return 1./x1*sum_pdf;
}
//qq + qbarqbar +qqbar (identical!), use it for DeltaNNLOC2
complex<double> fit_sum_qq_charge_weighted_single_vivi(complex<double> x1, complex<double> x2){
	complex<double> sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	for(int i = -5; i <=5; i++){
			if(i==0){continue;}
			sum_pdf+= eq[abs(i)-1]*eq[abs(i)-1]*1./x1*1./(x2/x1)*(xfit_pdfs(5+i,x1)*xfit_pdfs(5+i,x2/x1));
	}
	return 1./x1*sum_pdf;
}

// unweighted forms
complex<double> fit_sum_gg(complex<double> x1, complex<double> x2){
	complex<double> sum_pdf(0);
	sum_pdf = 1./x1*1./(x2/x1)*(xfit_pdfs(5,x1)*xfit_pdfs(5,x2/x1));
	return 1./x1*sum_pdf;
}
complex<double> fit_sum_qg(complex<double> x1, complex<double> x2){
	complex<double> sum_pdf(0);
	for(int i = 1; i <=5; i++){
		sum_pdf+= 1./x1*1./(x2/x1)*(xfit_pdfs(5-i,x1)*xfit_pdfs(5,x2/x1)+xfit_pdfs(5,x1)*xfit_pdfs(5-i,x2/x1));
		sum_pdf+= 1./x1*1./(x2/x1)*(xfit_pdfs(5+i,x1)*xfit_pdfs(5,x2/x1)+xfit_pdfs(5,x1)*xfit_pdfs(5+i,x2/x1));
	}
	return 1./x1*sum_pdf;
}
complex<double> fit_sum_qq(complex<double> x1, complex<double> x2){
	complex<double> sum_pdf(0);	
	for(int i = 1; i <=5; i++){
		sum_pdf+= 1./x1*1./(x2/x1)*(xfit_pdfs(5-i,x1)*xfit_pdfs(5-i,x2/x1)+xfit_pdfs(5+i,x1)*xfit_pdfs(5+i,x2/x1));
	}
	return 1./x1*sum_pdf;
}
complex<double> fit_sum_qqbar(complex<double> x1, complex<double> x2){
	complex<double> sum_pdf(0);	
	for(int i = 1; i <=5; i++){
		sum_pdf+= 1./x1*1./(x2/x1)*(xfit_pdfs(5-i,x1)*xfit_pdfs(5+i,x2/x1)+xfit_pdfs(5+i,x1)*xfit_pdfs(5-i,x2/x1));
	}
	return 1./x1*sum_pdf;
}
complex<double> fit_sum_qqNI(complex<double> x1, complex<double> x2){
	complex<double> sum_pdf(0);	
	for(int i = 1; i <=5; i++){
		for(int j = 1; j <=5; j++){
			if(i>=j){continue;} //to avoid double counting!
			sum_pdf+= 1./x1*1./(x2/x1)*(xfit_pdfs(5+i,x1)*xfit_pdfs(5+j,x2/x1)+xfit_pdfs(5+j,x1)*xfit_pdfs(5+i,x2/x1));
			sum_pdf+= 1./x1*1./(x2/x1)*(xfit_pdfs(5-i,x1)*xfit_pdfs(5-j,x2/x1)+xfit_pdfs(5-j,x1)*xfit_pdfs(5-i,x2/x1));
		}
	}
	return 1./x1*sum_pdf;
}

		
/////////////////////////////////////////////////////////////////////////////////
/// this is the convolution directly in N-space
/// note that this is the xfx(N) convoluted product of q and qbar. It works for CMP > 0.8 (if Nint not scaled!)
/////////////////////////////////////////////////////////////////////////////////
complex<double> fit_mellin_pdf_sum_qqbar_charge_weighted(complex<double> Nint){
	complex<double> sum_pdf(0,0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; 
	for(int i = 1; i <=5; i++){
		sum_pdf+= 2.*eq[i-1]*eq[i-1]*(xfit_Nspace_pdfs(5-i,Nint)*xfit_Nspace_pdfs(5+i,Nint));
	}
	return sum_pdf;
}

// gluon gluon convolution
complex<double> fit_mellin_pdf_sum_gg(complex<double> Nint){
	complex<double> sum_pdf(0,0);
	sum_pdf+= xfit_Nspace_pdfs(5,Nint)*xfit_Nspace_pdfs(5,Nint);
	return sum_pdf;
}
/// ---------------------------------------------------------------------------------------------



/////////////////////////////////////////////////////
/// and the functions that read from this structure
/////////////////////////////////////////////////////
// this is xfx(x)
complex<double> xfit_pdfs(int i, complex<double> x){
	complex<double> y = 1.-2.*pow(x,0.5);
	return fitcoeff[muF][i][0]*pow(1.-x,fitcoeff[muF][i][1])*pow(x,fitcoeff[muF][i][2])*(1.+fitcoeff[muF][i][3]*y+fitcoeff[muF][i][4]*(2.*pow(y,2)-1.))+fitcoeff[muF][i][5]*pow(1.-x,fitcoeff[muF][i][6])*pow(x,fitcoeff[muF][i][7]);
}
// this is fx(x)
complex<double> fit_pdfs(int i, complex<double> x){
	complex<double> y = 1.-2.*pow(x,0.5);
	return 1./x*(fitcoeff[muF][i][0]*pow(1.-x,fitcoeff[muF][i][1])*pow(x,fitcoeff[muF][i][2])*(1.+fitcoeff[muF][i][3]*y+fitcoeff[muF][i][4]*(2.*pow(y,2)-1.))+fitcoeff[muF][i][5]*pow(1.-x,fitcoeff[muF][i][6])*pow(x,fitcoeff[muF][i][7]));
}
// this is xfx(N) (mellin transform) checked that it gives the same answer when it is tranformed back to x-space
// note that if fx(N) is needed, then N-> N-1 works
complex<double> xfit_Nspace_pdfs(int i, complex<double> N){
	double A = fitcoeff[muF][i][0];
	double x3 = fitcoeff[muF][i][1];
	double x4 = fitcoeff[muF][i][2];
	double x5 = fitcoeff[muF][i][3];
	double x6 = fitcoeff[muF][i][4];
	double B = fitcoeff[muF][i][5];
	double x7 = fitcoeff[muF][i][6];
	double x8 = fitcoeff[muF][i][7];
	// these are optional, it looks like the bended contour now doesnt perform so well anymore
	// although it gets better for higher CMP (as expected)
	// checked that we can take these out, then it still works (it is quicker that way, as it doesn't have to evaluate the if statements)
	// also then we can bend the contour without making an error, as the function still exists there (it can be continued without introducing an error so it seems)
	//if(real(x4+N)<=0){return 0;}
	//if(real(x3)<=-1){return 0;}
	//if(real(x8+N)<=0){return 0;}
	//if(real(x7)<=-1){return 0;}
	return A*Gamma(1. + x3)*(-((2.*(x5 + 4.*x6)*Gamma(1./2. + N + x4))/Gamma(3./2. + N + x3 + x4)) + (((1. + N + x3 + x4)*(1. + x5) + (1. + 9.*N + x3 + 9.*x4)*x6)*Gamma(N + x4))/Gamma(2. + N + x3 + x4)) + (B*Gamma(1. + x7)*Gamma(N + x8))/Gamma(1. + N + x7 + x8);
}
