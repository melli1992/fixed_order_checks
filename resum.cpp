#include <cmath>
#include <iostream>
#include <fstream>
#include "cuba.h"
#include "cuba_integration.h"
#include <sstream>
#include "clooptools.h"
#include "deriv_pdf.h"
#include "monte_carlo.h"
#include "mellin_pdf.h"
#include <chrono>
#include "resum_functions.h"
#include "k_factors_dy.h"
#include "k_factors_higgs.h"
#include "k_factors_dihiggs.h"
#include "k_factors_diboson.h"
#include "k_factors_nnlo_dy.h"
#include "k_factors_nnlo_higgs.h"
#include "parameters.h"
#include "LHAPDF/LHAPDF.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <string.h>
#include <sstream>
#include "polygamma.h"
#include "mellin_pdf.h"

using namespace std;

string to_string(double Q){
	ostringstream q_to_str;
	q_to_str << Q;
	return q_to_str.str();
}

double polylog(int s, double z){
	s = s - 1;
	cout << "argument " << -log(z) << endl;
	return -gsl_sf_fermi_dirac_int(s,-log(z));
}


int main(int argc, char* argv[]){
	//////////////////////////////////////////////
	/// predefinition of everything, setting it up
	//////////////////////////////////////////////
	read_arguments(argc,argv);
	bool PDFmemuse = false;
    double scales[20] = {1000,900,800,750,700,650,600,550,500,480,450,420,400,390,370,350,330,300,270,260};
  // check whether the factorization scale is there
	ltini();
	bool TEST = false;

	// print of beginning of programme
	double z = 0.1;
	results lumi_pdfs;
	//int i = 0;
	//use_member = 50;
	//update_defaults();
	//for(int i = -5; i < 6; i++)
	/*{
		cout << "i=" << i << endl;
		for(int k = 1; k < 9; k++)
		{
			for(int j = 9; j > 0; j--)
			{
			complex<double> x = pow(10.,-k)*j;
			double xd = pow(10.,-k)*j;
			cout << xd << "," << real(xfit_pdfs(5+i, x)) << " " << real(Dxfit_pdfs(5+i, x)) << " " << real(NDxfit_pdfs(5+i, x, xd*1.E-6)) << ",";
			for(use_member = 0; use_member < 100; use_member++){ cout << pdfs[use_member]->xfxQ(i,xd,muF) << " "<< deriv_xpdf(i, xd, xd*1.E-6) << ",";}
			cout << endl;
			}
		}
	}*/
	/*
	ofstream output2;
	string q_str2 = "DihiggsSUSY3.txt";
	output2.open(q_str2.c_str()); //.c_str() needed to do a constant string conversion
	cout << q_str2 << endl;

	Q = 2.*sqrt(mH2);
	muF = Q;
	muR = Q;
	update_defaults();
	vector<results_c> diboson_LO_qqbar_full;
	set_SUSY(200,true);
	output2 << ", SUSY and mA2 = " << mA2 << endl;
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_hh","full",true);
	output2 << diboson_LO_qqbar_full[0].res << " " <<  diboson_LO_qqbar_full[0].err << endl;

	set_SUSY(950,true);
	output2 << ", SUSY and mA2 = " << mA2 << endl;
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_hh","full",true);
	output2 << diboson_LO_qqbar_full[0].res << " " <<  diboson_LO_qqbar_full[0].err << endl;
	set_SUSY(5000,true);
	output2 << ", SUSY and mA2 = " << mA2 << endl;
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_hh","full",true);
	output2 << diboson_LO_qqbar_full[0].res << " " <<  diboson_LO_qqbar_full[0].err << endl;
	output2.close();
	exit(0);

	for(int j = 0; j < 4; j++)
	{
		output2 << "j= " << j;
		cout << "j = " << j;
		//if(j==0) output2 << ", SM" << endl;
		//if(j==1) { set_SUSY(200,true);
		//			output2 << ", SUSY and mA2 = " << mA2 << endl;
		//			}
		if(j==0) { set_SUSY(950,true);
					output2 << ", SUSY and mA2 = " << mA2 << endl;
					}
		if(j==1) { set_SUSY(5000,true);
					output2 << ", SUSY and mA2 = " << mA2 << endl;
					}
		for(int i = 19; i > -1; i--){

		/*if(j == 0){ Q = scales[i];
					muF = Q;
					muR = Q;
					update_defaults(false);
					diboson_LO_qqbar_full = call_cuhre_dihiggs("SM","diff",true);
					output2 << "Q = " << Q << " GeV, " << diboson_LO_qqbar_full[0].res << " " <<  diboson_LO_qqbar_full[0].err << endl;
					cout << "Q = " << Q << " GeV, " << diboson_LO_qqbar_full[0].res << " " <<  diboson_LO_qqbar_full[0].err << endl; }
		if(j > -1)
		{ 			Q = scales[i];
					muF = Q;
					muR = Q;
					update_defaults(false);

					diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_hh","diff",true);
					output2 << "Q = " << Q << " GeV, " << diboson_LO_qqbar_full[0].res << " " <<  diboson_LO_qqbar_full[0].err << endl;
					cout << "Q = " << Q << " GeV, " << diboson_LO_qqbar_full[0].res << " " <<  diboson_LO_qqbar_full[0].err << endl;
					}
		}
	}

		output2.close();
	exit(0);
	//cout << fit_pdfs(5-i, x) << " " << Dfit_pdfs(5-i, x) << " " << NDfit_pdfs(5-i, x) << endl;
	//cout << pdfs[0]->xfxQ(i,xd,muF)/x << " " << deriv_pdf(i, xd) << endl;
	//exit(0);

	/*set_SUSY(950,true);
	vector<results_c> diboson_LO_qqbar_full;
	fitPDF = false;
	S = 13000;
	S2 = pow(S,2);
    Q = 2.*sqrt(mH2);
	muF = Q;
	muR = Q;
	update_defaults();

	diboson_LO_qqbar_full = call_cuhre_dihiggs("SM","full",false);
	cout << " SM hh =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_hh","full",false);
	cout << " hh =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	set_SUSY(5000,true);
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_hh","full",false);
	cout << " hh =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	set_SUSY(5900,true);
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_hh","full",false);
	cout << " hh =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;

	/*Q = sqrt(mH2)+sqrt(mHeavy2);
	muF = Q;
	muR = Q;
	update_defaults();
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_hH","full","delta",fitPDF);
	cout << " hH =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	Q = sqrt(mH2)+sqrt(mA2);
	muF = Q;
	muR = Q;
	update_defaults();
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_Ah","full","delta",fitPDF);
	cout << " Ah =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	Q = sqrt(mHeavy2)+sqrt(mHeavy2);
	muF = Q;
	muR = Q;
	update_defaults();
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_HH","full","delta",fitPDF);
	cout << " HH =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	Q = sqrt(mHeavy2)+sqrt(mA2);
	muF = Q;
	muR = Q;
	update_defaults();
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_AH","full","delta",fitPDF);
	cout << " AH =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	Q = sqrt(mA2)+sqrt(mA2);
	muF = Q;
	muR = Q;
	update_defaults();
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_AA","full","delta",fitPDF);
	cout << " AA =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;

	Q = 4.*sqrt(mH2);
	muF = 500.;
	muR = 500.;
	update_defaults();
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SM","diff","delta",fitPDF);
	cout << " hh =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SM","resumdiff","delta",fitPDF);
	cout << " hh =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_hh","diff","delta",fitPDF);
	cout << " hh =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_hh","resumdiff","delta",fitPDF);
	cout << " hh =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;

	Q = sqrt(mHeavy2)+2.*sqrt(mH2);
	muF = Q;
	muR = Q;
	update_defaults();
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_hh_diff","qqbar","delta",fitPDF);
	cout << " hh =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	diboson_LO_qqbar_full = call_cuhre_dihiggs("SUSY_hH_diff","qqbar","delta",fitPDF);
	cout << " hH =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	exit(0);

	//double s = 100*mW2;
	//partonic_up_wpwm(s);
	//exit(0);
	//diboson_LO_qqbar_full = call_cuhre_diboson("LOfull","ZZ",true);
	//cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	/*cout << "CMP " << CMP << " phiMP " << phiMP << endl;
	CMP = 1.3;
	diboson_LO_qqbar_full = call_cuhre_test("test","defor",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","nspace",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	CMP = 2.1-1;
	cout << "CMP " << CMP << " phiMP " << phiMP << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","defor",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","nspace",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	CMP = 1.9-1;
	cout << "CMP " << CMP << " phiMP " << phiMP << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","defor",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","nspace",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	CMP = 1.7-1;
	cout << "CMP " << CMP << " phiMP " << phiMP << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","defor",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","nspace",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	CMP = 1.5-1;
	cout << "CMP " << CMP << " phiMP " << phiMP << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","defor",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","nspace",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	CMP = 1.3-1;
	cout << "CMP " << CMP << " phiMP " << phiMP << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","defor",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","nspace",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	CMP = 1.1-1;
	cout << "CMP " << CMP << " phiMP " << phiMP << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","defor",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	diboson_LO_qqbar_full = call_cuhre_test("test","nspace",true);
	cout << " result =" << diboson_LO_qqbar_full[0].res << " " << diboson_LO_qqbar_full[0].err << " " << diboson_LO_qqbar_full[0].prob << endl;
	*/
	//////////////////////
	/// TEST FUNCTIONS
	//////////////////////
	bool MAKECOEF = false;
	if(MAKECOEF){
		Q2 = 1.;
		muF2 = 1.;
		muR2 = 1.;
		int Nsteps = 10000;
		ofstream output2;
		string q_str2 = "COEFFICIENTS_DY.txt";
		output2.open(q_str2.c_str()); //.c_str() needed to do a constant string conversion
		cout << q_str2 << endl;

		output2 << "x NLOqqbar NLOqg NNLOqqbarNS NNLOqg NNLOgg NNLOB2 NNLOBC NNLOC2 NNLOCD NNLOCE NNLOCF" << endl;
		for(int i=1;i<Nsteps;i++)
		{
			double x = double(float(i)/double(Nsteps));
			output2 << x << " ";
			output2 << /*M_PI/(alphas_muR)**/(DY_NLO_qqbar_reg(x)+DY_NLO_qqbar_plus(x)) << "," << /*M_PI/(alphas_muR)**/(DY_NLO_qqbar_plus(x)) << "," << /*M_PI/(alphas_muR)**/(DY_NLO_qqbar_plus(x)+DY_NLO_qqbar_expansion(x, 1)) << "," << /*M_PI/(alphas_muR)**/(DY_NLO_qqbar_plus(x)+DY_NLO_qqbar_expansion(x, 1)+DY_NLO_qqbar_expansion(x, 2)) << " ";
			output2 << /*M_PI/(alphas_muR)**/(DY_NLO_qg_full(x)) << "," << /*M_PI/(alphas_muR)**/(0) << "," << /*M_PI/(alphas_muR)**/(DY_NLO_qg_expansion(x, 1)) << "," << /*M_PI/(alphas_muR)**/(DY_NLO_qg_expansion(x, 1)+DY_NLO_qg_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qqbar_NS(x)+DY_NNLO_qqbar_plus(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qqbar_plus(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qqbar_plus(x)+DY_NNLO_qqbar_NS_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qqbar_plus(x)+DY_NNLO_qqbar_NS_expansion(x, 1)+DY_NNLO_qqbar_NS_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qg_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qg_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qg_expansion(x, 1)+DY_NNLO_qg_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_gg_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_gg_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_gg_expansion(x, 1)+DY_NNLO_gg_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_BB_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_BB_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_BB_expansion(x, 1)+DY_NNLO_BB_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_BC_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_BC_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_BC_expansion(x, 1)+DY_NNLO_BC_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CC_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CC_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CC_expansion(x, 1)+DY_NNLO_CC_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CD_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CD_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CD_expansion(x, 1)+DY_NNLO_CD_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CE_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CE_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CE_expansion(x, 1)+DY_NNLO_CE_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CF_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CF_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CF_expansion(x, 1)+DY_NNLO_CF_expansion(x, 2)) << " ";
			output2 << endl;
		}
		output2.close();

		q_str2 = "COEFFICIENTS_HIGGS.txt";
		output2.open(q_str2.c_str()); //.c_str() needed to do a constant string conversion
		cout << q_str2 << endl;

		output2 << "x NLOgg NLOqg NLOqqbar NNLOgg NNLOqg NNLOqq NNLOqqp NNLOqqbar" << endl;
		for(int i=1;i<Nsteps;i++)
		{
			double x = double(float(i)/double(Nsteps));
			output2 << x << " ";
			output2 << /*M_PI/(alphas_muR)**/(higgs_NLO_gg_reg(x)+higgs_NLO_gg_plus(x)) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_gg_plus(x)) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_gg_plus(x)+higgs_NLO_gg_expansion(x, 1)) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_gg_plus(x)+higgs_NLO_gg_expansion(x, 1)+higgs_NLO_gg_expansion(x, 2)) << " ";
			output2 << /*M_PI/(alphas_muR)**/(higgs_NLO_qg_full(x)) << "," << /*M_PI/(alphas_muR)**/(0) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_qg_expansion(x, 1)) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_qg_expansion(x, 1)+higgs_NLO_qg_expansion(x, 2)) << " ";
			output2 << /*M_PI/(alphas_muR)**/(higgs_NLO_qqbar_full(x)) << "," << /*M_PI/(alphas_muR)**/(0) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_qqbar_expansion(x, 1)) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_qqbar_expansion(x, 1)+higgs_NLO_qqbar_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_gg_reg(x)+higgs_NNLO_gg_plus(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_gg_plus(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_gg_plus(x)+higgs_NNLO_gg_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_gg_plus(x)+higgs_NNLO_gg_expansion(x, 1)+higgs_NNLO_gg_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qg_reg(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qg_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qg_expansion(x, 1)+higgs_NNLO_qg_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qq_reg(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qq_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qq_expansion(x, 1)+higgs_NNLO_qq_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qqp_reg(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qqp_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qqp_expansion(x, 1)+higgs_NNLO_qqp_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qqbar_reg(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qqbar_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qqbar_expansion(x, 1)+higgs_NNLO_qqbar_expansion(x, 2)) << " ";
			output2 << endl;
		}
		output2.close();
		exit(0);
	}
	if(TEST){
		z = 0.5;
		ofstream output2;
		string q_str2 = "COMPARISON_CODES.txt";
		output2.open(q_str2.c_str()); //.c_str() needed to do a constant string conversion
		cout << q_str2 << endl;
		output2 << "============================================================" << endl;
		output2 << "DY coefficients" << endl;
		output2 << "============================================================" << endl;
		output2 << "z=" << z <<", Q=" << Q << " GeV, muF=" << muF << " GeV, muR=" << muR << " GeV, tau=" << tau << ", alphas(muR2)=" << alphas_muR << endl;
		output2 << "LO factor (tau*sigma0): " << DY_LO_factor() << endl;
		output2 << "===============NLO-qqbar==========================" << endl;
		output2 << "NLO qqbar LP (*alphas/(4pi)): " << DY_LO_factor()*DY_NLO_qqbar_plus(z) << endl;
		output2 << "NLO qqbar z-regular: " << DY_LO_factor()*DY_NLO_qqbar_reg(z) << endl;
		output2 << "NLO qqbar delta: " << DY_LO_factor()*DY_NLO_qqbar_delta() << endl;

		output2 << "===============NLO-qg==========================" << endl;
		output2 << "NLO qg full: " << DY_LO_factor()*DY_NLO_qg_full(z) << endl;

		output2 << "===============NNLO-qqbar/qq/qbarqbar==========================" << endl;
		output2 << "NNLO qqbar LP (*alphas/(4pi))**2: " << DY_LO_factor()*DY_NNLO_qqbar_plus(z) << endl;
		output2 << "NNLO qqbar NS: " << DY_LO_factor()*DY_NNLO_qqbar_NS(z) << endl;
		output2 << "NNLO qqbar delta: " << DY_LO_factor()*DY_NNLO_qqbar_delta() << endl;
		output2 << "NNLO delta BB: " << DY_LO_factor()*DY_NNLO_BB_full(z) << endl;
		output2 << "NNLO delta BC: " << DY_LO_factor()*DY_NNLO_BC_full(z) << endl;
		output2 << "NNLO delta CC: " << DY_LO_factor()*DY_NNLO_CC_full(z) << endl;
		output2 << "NNLO delta CD: " << DY_LO_factor()*DY_NNLO_CD_full(z) << endl;
		output2 << "NNLO delta CE: " << DY_LO_factor()*DY_NNLO_CE_full(z) << endl;
		output2 << "NNLO delta CF: " << DY_LO_factor()*DY_NNLO_CF_full(z) << endl;
		output2 << "===============NNLO-qg==========================" << endl;
		output2 << "NNLO qg full: " << DY_LO_factor()*DY_NNLO_qg_full(z) << endl;
		output2 << "===============NNLO-gg==========================" << endl;
		output2 << "NNLO gg full: " << DY_LO_factor()*DY_NNLO_gg_full(z) << endl;

		output2 << "============================================================" << endl;
		output2 << "PDF convultions - weighted with charges" << endl;
		output2 << "============================================================" << endl;
		output2 << "Sum_i(Q_q[i]Q_q[i]*1/x*q_i(z)*qbar_i(tau/z) + <-> ), for DeltaNNLONS, DeltaNNLOBC: " << pdf_sum_qqbar_charge_weighted(z, tau) << endl;
		output2 << "Sum_i(Q_q[i]Q_q[i])*1/x*Sum_i(q_i(z)*qbar_i(tau/z) + <-> ), for DeltaNNLOB2: " << pdf_sum_qqbar_charge_unweighted(z, tau) << endl;
		output2 << "Sum_(i,j)(Q_q[i]Q_q[i]*1/x*q_i(z)*q_j(tau/z) + <-> + ...), for qq + qbarqbar +qqbar (identical+non-identical), for DeltaNNLOC2: " << pdf_sum_qq_charge_weighted_double(z, tau) << endl;
		output2 << "Sum_(i)(Q_q[i]Q_q[i]*1/x*q_i(z)*q_i(tau/z) + <-> + ...), for qq + qbarqbar +qqbar (identical), for DeltaNNLOCE: " << pdf_sum_qq_charge_weighted_single(z, tau) << endl;
		output2 << "Sum_(i,j)(Q_q[i]Q_q[j]*1/x*q_i(z)*q_j(tau/z) + ...), for qq + qbarqbar +qqbar (identical+NI), for DeltaNNLOCD: " << pdf_sum_qq_charge_weighted_double_vivj(z, tau) << endl;
		output2 << "Sum_(i)(Q_q[i]Q_q[i]*1/x*q_i(z)*q_i(tau/z) + ...), for qq + qbarqbar +qqbar (identical), for DeltaNNLOCF: " << pdf_sum_qq_charge_weighted_single_vivi(z, tau) << endl;
		output2 << "Sum_(i)(Q_q[i]Q_q[i]*1/x*q_i(z)*g(tau/z) + <-> + ...), for qg + qbarg: " << pdf_sum_qg_charge_weighted(z, tau) << endl;
		output2 << "Sum_(i)(Q_q[i]Q_q[i])*1/x*g(z)*g(tau/z)), for gg: " << pdf_sum_gg_charge_weighted(z, tau) << endl;


		output2 << "============================================================" << endl;
		output2 << "Higgs coefficients" << endl;
		output2 << "============================================================" << endl;
		output2 << "z=" << z <<", Q=" << Q << " GeV, muF=" << muF << " GeV, muR=" << muR << " GeV, tau=" << tau << ", alphas(muR2)=" << alphas_muR << endl;
		output2 << "LO factor (tau*sigma0): " << higgs_LO_factor() << endl;
		output2 << "===============NLO-gg==========================" << endl;
		output2 << "NLO gg LP (*alphas/(pi)): " << higgs_LO_factor()*higgs_NLO_gg_plus(z) << endl;
		output2 << "NLO gg z-regular: " << higgs_LO_factor()*higgs_NLO_gg_reg(z) << endl;
		output2 << "NLO gg delta: " << higgs_LO_factor()*higgs_NLO_gg_delta() << endl;

		output2 << "===============NLO-qg==========================" << endl;
		output2 << "NLO qg full: " << higgs_LO_factor()*higgs_NLO_qg_full(z) << endl;
		output2 << "===============NLO-qqbar==========================" << endl;
		output2 << "NLO qqbar full: " << higgs_LO_factor()*higgs_NLO_qqbar_full(z) << endl;


		output2 << "===============NNLO-gg==========================" << endl;
		output2 << "NNLO gg LP (*alphas/(pi))**2: " << higgs_LO_factor()*higgs_NNLO_gg_plus(z) << endl;
		output2 << "NNLO gg z-regular: " <<  higgs_LO_factor()*(higgs_NNLO_gg_reg(z) + 0*logdep_gg(z)) << endl;
		output2 << "NNLO gg delta: " << higgs_LO_factor()*(higgs_NNLO_gg_delta() + logdep_gg_constant()) << endl;
		nF = 0;
		output2 << "NNLO gg delta: " << higgs_LO_factor()*higgs_NNLO_gg_delta() << endl;
		nF = 5;
		output2 << "NNLO gg delta: " << higgs_LO_factor()*higgs_NNLO_gg_delta() << endl;
		output2 << "===============NNLO-qg==========================" << endl;
		output2 << "NNLO qg full: " << higgs_LO_factor()*(higgs_NNLO_qg_reg(z) +0*logdep_qg(z)) << endl;
		output2 << "===============NNLO-qqbar==========================" << endl;
		output2 << "NNLO qqbar full: " << higgs_LO_factor()*(higgs_NNLO_qqbar_reg(z)+0*logdep_qqbar(z)) << endl;
		output2 << "===============NNLO-qq==========================" << endl;
		output2 << "NNLO qq full: " << higgs_LO_factor()*(higgs_NNLO_qq_reg(z)+0*logdep_qq(z)) << endl;
		output2 << "===============NNLO-qq'==========================" << endl;
		output2 << "NNLO qq' full: " << higgs_LO_factor()*higgs_NNLO_qqp_reg(z) << endl;


		output2 << "============================================================" << endl;
		output2 << "PDF convultions" << endl;
		output2 << "============================================================" << endl;
		output2 << "Sum_i(1/x*q_i(z)*qbar_i(tau/z) + <-> ), for qqbar (identical): " << pdf_sum_qqbar(z, tau) << endl;
		output2 << "Sum_(i)(1/x*q_i(z)*g(tau/z) + <->), for qg (and qbarg): " << pdf_sum_qg(z, tau) << endl;
		output2 << "1/x*g(z)*g(tau/z)), for gg: " << pdf_sum_gg(z, tau) << endl;
		output2 << "Sum_i(1/x*q_i(z)*q_i(tau/z)+qbar_i(z)*qbar_i(tau/z)), for qq+qbarqbar (identical) : " << pdf_sum_qq(z, tau) << endl;
		output2 << "Sum_(i,j)(1/x*q_i(z)*q_j(tau/z)+<->+qbar_i(z)*qbar_j(tau/z)+<->), for qq+qbarqbar (non-identical) : " << pdf_sum_qqNI(z, tau) << endl;
	}
	//////////////////////
	/// LO declarations
	//////////////////////
	// DY
	vector<results_c> DY_LO_qqbar_full;
	// higgs
	vector<results_c> higgs_LO_gg_full;
	// dihiggs
	vector<results_c> dihiggs_LO_full;
	// wpwm
	vector<results_c> ww_LO_full;
	// zz
	vector<results_c> zz_LO_full;
	// wpwm
	vector<results_c> hh_LO_full;
	// wpwm
	vector<results_c> ww_LO_diff;
	// zz
	vector<results_c> zz_LO_diff;
	// zz
	vector<results_c> hh_LO_diff;
	////////////////////////////////////////////////////////////////////

	//////////////////////
	/// NLO declarations
	//////////////////////
	// DY
	vector<results_c> res_DY_NLO_qqbar_full, res_DY_NLO_qqbar_hard, res_DY_NLO_qqbar_LP, res_DY_NLO_qqbar_LP_part1, res_DY_NLO_qqbar_LP_cor,res_DY_NLO_qqbar_delta, res_DY_NLO_qqbar_accum, DY_NLO_qqbar_full;
	vector<results_c> res_DY_NLO_qg_full, res_DY_NLO_qg_accum;
	vector<results_c> DY_NLO_qg_powers, DY_NLO_qqbar_powers;
	// higgs
	vector<results_c> res_higgs_NLO_gg_full, res_higgs_NLO_gg_hard, res_higgs_NLO_gg_LP, res_higgs_NLO_gg_LP_part1, res_higgs_NLO_gg_LP_cor,res_higgs_NLO_gg_delta, res_higgs_NLO_gg_accum, higgs_NLO_gg_full;
	vector<results_c> res_higgs_NLO_qg_full, res_higgs_NLO_qg_accum;
	vector<results_c> res_higgs_NLO_qqbar_full,res_higgs_NLO_qqbar_accum;
	vector<results_c> higgs_NLO_qg_powers, higgs_NLO_gg_powers, higgs_NLO_qqbar_powers;
	// prompt photon (not complete yet)
	results res_pf_NLO_qqbar_hard, res_pf_NLO_qqbar_LP, res_pf_NLO_qqbar_LP_part1, res_pf_NLO_qqbar_LP_cor,res_pf_NLO_qqbar_delta, res_pf_NLO_qqbar_accum, res_pf_NLO_qqbar_full;
	results  pf_NLO_qqbar_powers;
	////////////////////////////////////////////////////////////////////

	//////////////////////
	/// NNLO declarations
	//////////////////////
	// DY
	vector<results_c> res_DY_NNLO_qqbar_full, DY_NNLO_qqbar_full, res_DY_NNLO_qqbar_hard, res_DY_NNLO_qqbar_LP, res_DY_NNLO_qqbar_LP_part1, res_DY_NNLO_qqbar_LP_cor,res_DY_NNLO_qqbar_delta, res_DY_NNLO_qqbar_accum;
	vector<results_c> res_DY_NNLO_qg_full, res_DY_NNLO_qg_accum;
	vector<results_c> res_DY_NNLO_qq_full, res_DY_NNLO_qq_accum;
	vector<results_c> res_DY_NNLO_gg_full, res_DY_NNLO_gg_accum;
	vector<results_c> DY_NNLO_gg_powers, DY_NNLO_qg_powers, DY_NNLO_qq_powers, DY_NNLO_qqbar_powers;
	// higgs
	vector<results_c> res_higgs_NNLO_gg_full, higgs_NNLO_gg_full, res_higgs_NNLO_gg_hard, res_higgs_NNLO_gg_LP, res_higgs_NNLO_gg_LP_part1, res_higgs_NNLO_gg_LP_cor,res_higgs_NNLO_gg_delta, res_higgs_NNLO_gg_accum;
	vector<results_c> res_higgs_NNLO_qg_full, res_higgs_NNLO_qg_accum;
	vector<results_c> res_higgs_NNLO_qq_full, res_higgs_NNLO_qq_accum;
	vector<results_c> res_higgs_NNLO_qqbar_full, res_higgs_NNLO_qqbar_accum;
	vector<results_c> res_higgs_NNLO_qqp_full, res_higgs_NNLO_qqp_accum;
	vector<results_c> higgs_NNLO_gg_powers, higgs_NNLO_qg_powers, higgs_NNLO_qq_powers, higgs_NNLO_qqp_powers, higgs_NNLO_qqbar_powers;
	////////////////////////////////////////////////////////////////////

	///////////////////////
	/// resummed declarations
	///////////////////////

	//DY
	vector<results_c> resummed_DY_LP_NNLL, resummed_DY_LP_NNLL_NLP_LL, resummed_DY_LP_NLL, resummed_DY_LP_NLL_NLP_LL, resummed_DY_LP_LL_NLP_LL, resummed_DY_LP_LL;
	vector<results_c> resummed_DY_LP_NNLL_exp_NNLO, resummed_DY_LP_NNLL_NLP_LL_exp_NNLO, resummed_DY_LP_NLL_exp_NNLO, resummed_DY_LP_NLL_NLP_LL_exp_NNLO, resummed_DY_LP_LL_NLP_LL_exp_NNLO, resummed_DY_LP_LL_exp_NNLO;
	vector<results_c> resummed_DY_LP_NNLL_exp_NLO, resummed_DY_LP_NNLL_NLP_LL_exp_NLO, resummed_DY_LP_NLL_exp_NLO, resummed_DY_LP_NLL_NLP_LL_exp_NLO, resummed_DY_LP_LL_NLP_LL_exp_NLO, resummed_DY_LP_LL_exp_NLO;

	// higgs
	vector<results_c> resummed_higgs_LP_NNLL, resummed_higgs_LP_NNLL_NLP_LL, resummed_higgs_LP_NLL, resummed_higgs_LP_NLL_NLP_LL, resummed_higgs_LP_LL_NLP_LL, resummed_higgs_LP_LL;
	vector<results_c> resummed_higgs_LP_NNLL_exp_NNLO, resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO, resummed_higgs_LP_NLL_exp_NNLO, resummed_higgs_LP_NLL_NLP_LL_exp_NNLO, resummed_higgs_LP_LL_NLP_LL_exp_NNLO, resummed_higgs_LP_LL_exp_NNLO;
	vector<results_c> resummed_higgs_LP_NNLL_exp_NLO, resummed_higgs_LP_NNLL_NLP_LL_exp_NLO, resummed_higgs_LP_NLL_exp_NLO, resummed_higgs_LP_NLL_NLP_LL_exp_NLO, resummed_higgs_LP_LL_NLP_LL_exp_NLO, resummed_higgs_LP_LL_exp_NLO;


	//wpwm
	vector<results_c> resummedM_ww_LP_NNLL, resummedM_ww_LP_NNLL_NLP_LL, resummedM_ww_LP_NLL, resummedM_ww_LP_NLL_NLP_LL, resummedM_ww_LP_LL_NLP_LL, resummedM_ww_LP_LL;
	vector<results_c> resummedQ_ww_LP_NNLL, resummedQ_ww_LP_NNLL_NLP_LL, resummedQ_ww_LP_NLL, resummedQ_ww_LP_NLL_NLP_LL, resummedQ_ww_LP_LL_NLP_LL, resummedQ_ww_LP_LL;
	vector<results_c> resummed_ww_diff_LP_NNLL, resummed_ww_diff_LP_NNLL_NLP_LL, resummed_ww_diff_LP_NLL, resummed_ww_diff_LP_NLL_NLP_LL, resummed_ww_diff_LP_LL_NLP_LL, resummed_ww_diff_LP_LL;

  //zz
	vector<results_c> resummedM_zz_LP_NNLL, resummedM_zz_LP_NNLL_NLP_LL, resummedM_zz_LP_NLL, resummedM_zz_LP_NLL_NLP_LL, resummedM_zz_LP_LL_NLP_LL, resummedM_zz_LP_LL;
	vector<results_c> resummedQ_zz_LP_NNLL, resummedQ_zz_LP_NNLL_NLP_LL, resummedQ_zz_LP_NLL, resummedQ_zz_LP_NLL_NLP_LL, resummedQ_zz_LP_LL_NLP_LL, resummedQ_zz_LP_LL;
	vector<results_c> resummed_zz_diff_LP_NNLL, resummed_zz_diff_LP_NNLL_NLP_LL, resummed_zz_diff_LP_NLL, resummed_zz_diff_LP_NLL_NLP_LL, resummed_zz_diff_LP_LL_NLP_LL, resummed_zz_diff_LP_LL;

  //hh
	vector<results_c> resummed_hh_diff_LP_NNLL, resummed_hh_diff_LP_NNLL_NLP_LL, resummed_hh_diff_LP_NLL, resummed_hh_diff_LP_NLL_NLP_LL, resummed_hh_diff_LP_LL_NLP_LL, resummed_hh_diff_LP_LL;

	////////////////////////////////////////////////////////////////////

	///////////////////
	/// output
	///////////////////
	ofstream output;
	ostringstream x_convert; // need this for the output
	x_convert << Q;
	string Qstring  = x_convert.str();
	ostringstream x_convert2; // need this for the output
	x_convert2 << alphas_Q;
	string asstring  = x_convert2.str();
	ostringstream x_convert3; // need this for the output
	x_convert3 << muF;
	string mufstring  = x_convert3.str();
	string q_str = "output_Q" + Qstring +"_as"+asstring+"_"+setname+"_muF"+mufstring;
	if(DY) q_str = q_str+"_DY";
	else if(hh && !SUSY) q_str = q_str+"_dihiggs";
	else if(hh && SUSY){ q_str = q_str+"_dihiggsSUSY";
							ostringstream x_convert4; // need this for the output
							x_convert4 << SUSYset;
							string SUS  = x_convert4.str();
							q_str = q_str+"_MA"+SUS;
							ostringstream x_convert5; // need this for the output
							x_convert5 << tb;
							string TA  = x_convert5.str();
							q_str = q_str+"_tb"+TA;
						}
	else if(higgs) q_str = q_str+"_Higgs";
	else if(WW) q_str = q_str+"_WW";
	else if(ZZ) q_str = q_str+"_ZZ";
	else {cout << "define the process!" << endl; exit(0);}
	if((DY) || (higgs)){
		if(LO) q_str = q_str+"_LO";
		if(NLO) q_str = q_str+"_NLO";
		if(NNLO) q_str = q_str+"_NNLO";
	}
	else{
		if(full) q_str = q_str+"_full";
		if(diff) q_str = q_str+"_diff";
	}
	if(PDFmemuse){q_str = q_str+"_members"; RES=false; fitPDF=false;}
	if(realPDF){q_str = q_str+"_real_pdfs"; RES=false; fitPDF=false;} //if real pdfs, the resummed cannot be used!
	if(fitPDF){q_str = q_str+"_new_fitted_pdfs";}
	if(RES) q_str = q_str+"_resummed";
	q_str = q_str + ".txt";
	cout << q_str.c_str() << endl;
	////////////////////////////////////////////////////////////////////

	/////////////////////////
	/// computation
	/////////////////////////
  if(higgs){
		cout << "working on higgs" << endl;
		Q = mH;
		Q2 = pow(Q,2);
		update_defaults();
		if((PDFmemuse == true) && (fitPDF == false))
		{
			for(int i = 0; i < 101;i++)
				{
					use_member = i;
					update_defaults();

					higgs_LO_gg_full.push_back(call_cuhre_higgs("LO","gg","LP",fitPDF)[0]);
					res_higgs_NLO_gg_hard.push_back(call_cuhre_higgs("NLO","gg","reg",fitPDF)[0]);
					res_higgs_NLO_gg_LP_part1.push_back(call_cuhre_higgs("NLO","gg","LP",fitPDF)[0]);
					res_higgs_NLO_gg_LP_cor.push_back(call_cuhre_higgs("NLO","gg","LP_corr",fitPDF)[0]);
					res_higgs_NLO_gg_delta.push_back(call_cuhre_higgs("NLO","gg","delta",fitPDF)[0]);
					res_higgs_NLO_qg_full.push_back(call_cuhre_higgs("NLO","qg","full",fitPDF)[0]);
					res_higgs_NLO_qqbar_full.push_back(call_cuhre_higgs("NLO","qqbar","full",fitPDF)[0]);
					res_higgs_NLO_gg_LP.push_back({res_higgs_NLO_gg_LP_part1[i].res + res_higgs_NLO_gg_LP_cor[i].res, res_higgs_NLO_gg_LP_part1[i].err + res_higgs_NLO_gg_LP_cor[i].err, res_higgs_NLO_gg_LP_part1[i].prob + res_higgs_NLO_gg_LP_cor[i].prob});
					higgs_NLO_gg_full.push_back({res_higgs_NLO_gg_hard[i].res+res_higgs_NLO_gg_LP_cor[i].res+res_higgs_NLO_gg_delta[i].res, res_higgs_NLO_gg_hard[i].err+res_higgs_NLO_gg_LP_cor[i].err+res_higgs_NLO_gg_delta[i].err, res_higgs_NLO_gg_hard[i].prob+res_higgs_NLO_gg_LP_cor[i].prob+res_higgs_NLO_gg_delta[i].prob});
					res_higgs_NNLO_gg_hard.push_back(call_cuhre_higgs("NNLO","gg","reg",fitPDF)[0]);
					res_higgs_NNLO_gg_LP_part1.push_back(call_cuhre_higgs("NNLO","gg","LP",fitPDF)[0]);
					res_higgs_NNLO_gg_LP_cor.push_back(call_cuhre_higgs("NNLO","gg","LP_corr",fitPDF)[0]);
					res_higgs_NNLO_gg_delta.push_back(call_cuhre_higgs("NNLO","gg","delta",fitPDF)[0]);
					res_higgs_NNLO_qg_full.push_back(call_cuhre_higgs("NNLO","qg","full",fitPDF)[0]);
					res_higgs_NNLO_qq_full.push_back(call_cuhre_higgs("NNLO","qq","full",fitPDF)[0]);
					res_higgs_NNLO_qqp_full.push_back(call_cuhre_higgs("NNLO","qqp","full",fitPDF)[0]);
					res_higgs_NNLO_qqbar_full.push_back(call_cuhre_higgs("NNLO","qqbar","full",fitPDF)[0]);
					res_higgs_NNLO_gg_LP.push_back({res_higgs_NNLO_gg_LP_part1[i].res + res_higgs_NNLO_gg_LP_cor[i].res, res_higgs_NNLO_gg_LP_part1[i].err + res_higgs_NNLO_gg_LP_cor[i].err, res_higgs_NNLO_gg_LP_part1[i].prob + res_higgs_NNLO_gg_LP_cor[i].prob});
					higgs_NNLO_gg_full.push_back({res_higgs_NNLO_gg_hard[i].res+res_higgs_NNLO_gg_LP_cor[i].res+res_higgs_NNLO_gg_delta[i].res, res_higgs_NNLO_gg_hard[i].err+res_higgs_NNLO_gg_LP_cor[i].err+res_higgs_NNLO_gg_delta[i].err, res_higgs_NNLO_gg_hard[i].prob+res_higgs_NNLO_gg_LP_cor[i].prob+res_higgs_NNLO_gg_delta[i].prob});
				}
		}
		else{
		if(LO&&higgs){
			cout << "working on LO" << endl;
			higgs_LO_gg_full = call_cuhre_higgs("LO","gg","LP",fitPDF);
		}
		if(NLO&&higgs){

			cout << "working on NLO" << endl;
			res_higgs_NLO_gg_hard = call_cuhre_higgs("NLO","gg","reg",fitPDF);
			res_higgs_NLO_gg_LP_part1 = call_cuhre_higgs("NLO","gg","LP",fitPDF);
			res_higgs_NLO_gg_LP_cor = call_cuhre_higgs("NLO","gg","LP_corr",fitPDF);
			res_higgs_NLO_gg_delta = call_cuhre_higgs("NLO","gg","delta",fitPDF);
			res_higgs_NLO_qg_full = call_cuhre_higgs("NLO","qg","full",fitPDF);
			res_higgs_NLO_qqbar_full = call_cuhre_higgs("NLO","qqbar","full",fitPDF);
			res_higgs_NLO_gg_LP.push_back({res_higgs_NLO_gg_LP_part1[0].res + res_higgs_NLO_gg_LP_cor[0].res, res_higgs_NLO_gg_LP_part1[0].err + res_higgs_NLO_gg_LP_cor[0].err, res_higgs_NLO_gg_LP_part1[0].prob + res_higgs_NLO_gg_LP_cor[0].prob});
			higgs_NLO_gg_powers = call_cuhre_higgs("NLO","gg","exp",fitPDF,10);
			higgs_NLO_qg_powers = call_cuhre_higgs("NLO","qg","exp",fitPDF,10);
			higgs_NLO_qqbar_powers = call_cuhre_higgs("NLO","qqbar","exp",fitPDF,10);
			higgs_NLO_gg_full.push_back({res_higgs_NLO_gg_hard[0].res+res_higgs_NLO_gg_LP_cor[0].res+res_higgs_NLO_gg_delta[0].res, res_higgs_NLO_gg_hard[0].err+res_higgs_NLO_gg_LP_cor[0].err+res_higgs_NLO_gg_delta[0].err, res_higgs_NLO_gg_hard[0].prob+res_higgs_NLO_gg_LP_cor[0].prob+res_higgs_NLO_gg_delta[0].prob});

		}
		if(NNLO&&higgs){
			cout << "working on NNLO" << endl;
			res_higgs_NNLO_gg_hard = call_cuhre_higgs("NNLO","gg","reg",fitPDF);
			res_higgs_NNLO_gg_LP_part1 = call_cuhre_higgs("NNLO","gg","LP",fitPDF);
			res_higgs_NNLO_gg_LP_cor = call_cuhre_higgs("NNLO","gg","LP_corr",fitPDF);
			res_higgs_NNLO_gg_delta = call_cuhre_higgs("NNLO","gg","delta",fitPDF);
			res_higgs_NNLO_qg_full = call_cuhre_higgs("NNLO","qg","full",fitPDF);
			res_higgs_NNLO_qq_full = call_cuhre_higgs("NNLO","qq","full",fitPDF);
			res_higgs_NNLO_qqp_full = call_cuhre_higgs("NNLO","qqp","full",fitPDF);
			res_higgs_NNLO_qqbar_full = call_cuhre_higgs("NNLO","qqbar","full",fitPDF);
			res_higgs_NNLO_gg_LP.push_back({res_higgs_NNLO_gg_LP_part1[0].res + res_higgs_NNLO_gg_LP_cor[0].res, res_higgs_NNLO_gg_LP_part1[0].err + res_higgs_NNLO_gg_LP_cor[0].err, res_higgs_NNLO_gg_LP_part1[0].prob + res_higgs_NNLO_gg_LP_cor[0].prob});
			higgs_NNLO_gg_powers = call_cuhre_higgs("NNLO","gg","exp",fitPDF,10);
			higgs_NNLO_qg_powers = call_cuhre_higgs("NNLO","qg","exp",fitPDF,10);
			higgs_NNLO_qq_powers = call_cuhre_higgs("NNLO","qq","exp",fitPDF,10);
			higgs_NNLO_qqp_powers = call_cuhre_higgs("NNLO","qqp","exp",fitPDF,10);
			higgs_NNLO_qqbar_powers = call_cuhre_higgs("NLO","qqbar","exp",fitPDF,10);
			higgs_NNLO_gg_full.push_back({res_higgs_NNLO_gg_hard[0].res+res_higgs_NNLO_gg_LP_cor[0].res+res_higgs_NNLO_gg_delta[0].res, res_higgs_NNLO_gg_hard[0].err+res_higgs_NNLO_gg_LP_cor[0].err+res_higgs_NNLO_gg_delta[0].err, res_higgs_NNLO_gg_hard[0].prob+res_higgs_NNLO_gg_LP_cor[0].prob+res_higgs_NNLO_gg_delta[0].prob});
			}
		if(RES&&higgs){
			cout << "computing the resummed results" << endl;
			cout << "LP NNLL + NLP LL" << endl;
			ISNNLL = 1;
			ISNLP = 1;
			ISLL = 1;
			ISNLL = 1;
			resummed_higgs_LP_NNLL_NLP_LL = call_cuhre_higgs("resum","gg","full",fitPDF);
			resummed_higgs_LP_NNLL_NLP_LL_exp_NLO =call_cuhre_higgs("resum","gg","NLO",fitPDF);
			resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO = call_cuhre_higgs("resum","gg","NNLO",fitPDF);

			cout << "LP NLL + NLP LL" << endl;
			ISNNLL = 0;
			ISNLP = 1;
			ISLL = 1;
			ISNLL = 1;
			resummed_higgs_LP_NLL_NLP_LL = call_cuhre_higgs("resum","gg","full",fitPDF);
			resummed_higgs_LP_NLL_NLP_LL_exp_NLO = call_cuhre_higgs("resum","gg","NLO",fitPDF);
			resummed_higgs_LP_NLL_NLP_LL_exp_NNLO = call_cuhre_higgs("resum","gg","NNLO",fitPDF);

			cout << "LP LL + NLP LL" << endl;
			ISNNLL = 0;
			ISNLP = 1;
			ISLL = 1;
			ISNLL = 0;
			resummed_higgs_LP_LL_NLP_LL = call_cuhre_higgs("resum","gg","full",fitPDF);
			resummed_higgs_LP_LL_NLP_LL_exp_NLO = call_cuhre_higgs("resum","gg","NLO",fitPDF);
			resummed_higgs_LP_LL_NLP_LL_exp_NNLO = call_cuhre_higgs("resum","gg","NNLO",fitPDF);

			cout << "LP NNLL" << endl;
			ISNNLL = 1;
			ISNLP = 0;
			ISLL = 1;
			ISNLL = 1;
			resummed_higgs_LP_NNLL = call_cuhre_higgs("resum","gg","full",fitPDF);
			resummed_higgs_LP_NNLL_exp_NLO = call_cuhre_higgs("resum","gg","NLO",fitPDF);
			resummed_higgs_LP_NNLL_exp_NNLO = call_cuhre_higgs("resum","gg","NNLO",fitPDF);

			cout << "LP NLL" << endl;
			ISNNLL = 0;
			ISNLP = 0;
			ISLL = 1;
			ISNLL = 1;
			resummed_higgs_LP_NLL = call_cuhre_higgs("resum","gg","full",fitPDF);
			resummed_higgs_LP_NLL_exp_NLO = call_cuhre_higgs("resum","gg","NLO",fitPDF);
			resummed_higgs_LP_NLL_exp_NNLO = call_cuhre_higgs("resum","gg","NNLO",fitPDF);

			cout << "LP LL" << endl;
			ISNNLL = 0;
			ISNLP = 0;
			ISLL = 1;
			ISNLL = 0;
			resummed_higgs_LP_LL = call_cuhre_higgs("resum","gg","full",fitPDF);
			resummed_higgs_LP_LL_exp_NLO = call_cuhre_higgs("resum","gg","NLO",fitPDF);
			resummed_higgs_LP_LL_exp_NNLO = call_cuhre_higgs("resum","gg","NNLO",fitPDF);


		}
	}
}
  if(DY){
		cout << "working on DY" << endl;
	    if((PDFmemuse == true) && (fitPDF == false))
		{
			for(int i = 0; i < 101;i++)
				{
					use_member = i;
					update_defaults();
					DY_LO_qqbar_full.push_back(call_cuhre_dy("LO","qqbar","LP",fitPDF)[0]);
					res_DY_NLO_qqbar_hard.push_back(call_cuhre_dy("NLO","qqbar","reg",fitPDF)[0]);
					res_DY_NLO_qqbar_LP_part1.push_back(call_cuhre_dy("NLO","qqbar","LP",fitPDF)[0]);
					res_DY_NLO_qqbar_LP_cor.push_back(call_cuhre_dy("NLO","qqbar","LP_corr",fitPDF)[0]);
					res_DY_NLO_qqbar_delta.push_back(call_cuhre_dy("NLO","qqbar","delta",fitPDF)[0]);
					res_DY_NLO_qg_full.push_back(call_cuhre_dy("NLO","qg","full",fitPDF)[0]);
					res_DY_NLO_qqbar_LP.push_back({res_DY_NLO_qqbar_LP_part1[i].res + res_DY_NLO_qqbar_LP_cor[i].res, res_DY_NLO_qqbar_LP_part1[i].err + res_DY_NLO_qqbar_LP_cor[i].err, res_DY_NLO_qqbar_LP_part1[i].prob + res_DY_NLO_qqbar_LP_cor[i].prob});
					DY_NLO_qqbar_full.push_back({res_DY_NLO_qqbar_hard[i].res+res_DY_NLO_qqbar_LP_cor[i].res+res_DY_NLO_qqbar_delta[i].res, res_DY_NLO_qqbar_hard[i].err+res_DY_NLO_qqbar_LP_cor[i].err+res_DY_NLO_qqbar_delta[i].err, res_DY_NLO_qqbar_hard[i].prob+res_DY_NLO_qqbar_LP_cor[i].prob+res_DY_NLO_qqbar_delta[i].prob});
					res_DY_NNLO_qqbar_hard.push_back(call_cuhre_dy("NNLO","qqbar","reg",fitPDF)[0]);
					res_DY_NNLO_qqbar_LP_part1.push_back(call_cuhre_dy("NNLO","qqbar","LP",fitPDF)[0]);
					res_DY_NNLO_qqbar_LP_cor.push_back(call_cuhre_dy("NNLO","qqbar","LP_corr",fitPDF)[0]);
					res_DY_NNLO_qqbar_delta.push_back(call_cuhre_dy("NNLO","qqbar","delta",fitPDF)[0]);
					res_DY_NNLO_qg_full.push_back(call_cuhre_dy("NNLO","qg","full",fitPDF)[0]);
					res_DY_NNLO_qq_full.push_back(call_cuhre_dy("NNLO","qq","full",fitPDF)[0]);
					res_DY_NNLO_gg_full.push_back(call_cuhre_dy("NNLO","gg","full",fitPDF)[0]);
					res_DY_NNLO_qqbar_LP.push_back({res_DY_NNLO_qqbar_LP_part1[i].res + res_DY_NNLO_qqbar_LP_cor[i].res, res_DY_NNLO_qqbar_LP_part1[i].err + res_DY_NNLO_qqbar_LP_cor[i].err, res_DY_NNLO_qqbar_LP_part1[i].prob + res_DY_NNLO_qqbar_LP_cor[i].prob});
					DY_NNLO_qqbar_full.push_back({res_DY_NNLO_qqbar_hard[i].res+res_DY_NNLO_qqbar_LP_cor[i].res+res_DY_NNLO_qqbar_delta[i].res, res_DY_NNLO_qqbar_hard[i].err+res_DY_NNLO_qqbar_LP_cor[i].err+res_DY_NNLO_qqbar_delta[i].err, res_DY_NNLO_qqbar_hard[i].prob+res_DY_NNLO_qqbar_LP_cor[i].prob+res_DY_NNLO_qqbar_delta[i].prob});
					cout << DY_NLO_qqbar_full[i].res << endl;
				}
		}
		else
		{
		update_defaults();

		if(LO&&DY){
			cout << "working on LO" << endl;
			DY_LO_qqbar_full = call_cuhre_dy("LO","qqbar","full",fitPDF);
		}
		if(NLO&&DY){
			cout << "working on NLO" << endl;
			res_DY_NLO_qqbar_hard = call_cuhre_dy("NLO","qqbar","reg",fitPDF);
			res_DY_NLO_qqbar_LP_part1 = call_cuhre_dy("NLO","qqbar","LP",fitPDF);
			res_DY_NLO_qqbar_LP_cor = call_cuhre_dy("NLO","qqbar","LP_corr",fitPDF);
			res_DY_NLO_qqbar_delta = call_cuhre_dy("NLO","qqbar","delta",fitPDF);
			res_DY_NLO_qg_full = call_cuhre_dy("NLO","qg","full",fitPDF);
			res_DY_NLO_qqbar_LP.push_back({res_DY_NLO_qqbar_LP_part1[0].res + res_DY_NLO_qqbar_LP_cor[0].res, res_DY_NLO_qqbar_LP_part1[0].err + res_DY_NLO_qqbar_LP_cor[0].err, res_DY_NLO_qqbar_LP_part1[0].prob + res_DY_NLO_qqbar_LP_cor[0].prob});
			DY_NLO_qqbar_powers = call_cuhre_dy("NLO","qqbar","exp",fitPDF,10);
			DY_NLO_qg_powers = call_cuhre_dy("NLO","qg","exp",fitPDF,10);
			DY_NLO_qqbar_full.push_back({res_DY_NLO_qqbar_hard[0].res+res_DY_NLO_qqbar_LP_cor[0].res+res_DY_NLO_qqbar_delta[0].res, res_DY_NLO_qqbar_hard[0].err+res_DY_NLO_qqbar_LP_cor[0].err+res_DY_NLO_qqbar_delta[0].err, res_DY_NLO_qqbar_hard[0].prob + res_DY_NLO_qqbar_LP_cor[0].prob+res_DY_NLO_qqbar_delta[0].prob});

		}
		if(NNLO&&DY){
			cout << "working on NNLO" << endl;
			res_DY_NNLO_qqbar_hard = call_cuhre_dy("NNLO","qqbar","reg",fitPDF);
			res_DY_NNLO_qqbar_LP_part1 = call_cuhre_dy("NNLO","qqbar","LP",fitPDF);
			res_DY_NNLO_qqbar_LP_cor = call_cuhre_dy("NNLO","qqbar","LP_corr",fitPDF);
			res_DY_NNLO_qqbar_delta = call_cuhre_dy("NNLO","qqbar","delta",fitPDF);
			res_DY_NNLO_qg_full = call_cuhre_dy("NNLO","qg","full",fitPDF);
			res_DY_NNLO_qq_full = call_cuhre_dy("NNLO","qq","full",fitPDF);
			res_DY_NNLO_gg_full = call_cuhre_dy("NNLO","gg","full",fitPDF);
			res_DY_NNLO_qqbar_LP.push_back({res_DY_NNLO_qqbar_LP_part1[0].res + res_DY_NNLO_qqbar_LP_cor[0].res, res_DY_NNLO_qqbar_LP_part1[0].err + res_DY_NNLO_qqbar_LP_cor[0].err, res_DY_NNLO_qqbar_LP_part1[0].prob + res_DY_NNLO_qqbar_LP_cor[0].prob});
			DY_NNLO_gg_powers = call_cuhre_dy("NNLO","gg","exp",fitPDF,10);
			DY_NNLO_qg_powers = call_cuhre_dy("NNLO","qg","exp",fitPDF,10);
			DY_NNLO_qq_powers = call_cuhre_dy("NNLO","qq","exp",fitPDF,10);
			DY_NNLO_qqbar_powers = call_cuhre_dy("NNLO","qqbar","exp",fitPDF,10);
			DY_NNLO_qqbar_full.push_back({res_DY_NNLO_qqbar_hard[0].res+res_DY_NNLO_qqbar_LP_cor[0].res+res_DY_NNLO_qqbar_delta[0].res, res_DY_NNLO_qqbar_hard[0].err+res_DY_NNLO_qqbar_LP_cor[0].err+res_DY_NNLO_qqbar_delta[0].err, res_DY_NNLO_qqbar_hard[0].prob+res_DY_NNLO_qqbar_LP_cor[0].prob+res_DY_NNLO_qqbar_delta[0].prob});
		}
		if(RES&&DY){
			cout << "computing the resummed results" << endl;
			cout << "LP NNLL + NLP LL" << endl;
			ISNNLL = 1;
			ISNLP = 1;
			ISLL = 1;
			ISNLL = 1;
			resummed_DY_LP_NNLL_NLP_LL = call_cuhre_dy("resum","gg","full",fitPDF);
			resummed_DY_LP_NNLL_NLP_LL_exp_NLO =call_cuhre_dy("resum","gg","NLO",fitPDF);
			resummed_DY_LP_NNLL_NLP_LL_exp_NNLO = call_cuhre_dy("resum","gg","NNLO",fitPDF);

			cout << "LP NLL + NLP LL" << endl;
			ISNNLL = 0;
			ISNLP = 1;
			ISLL = 1;
			ISNLL = 1;
			resummed_DY_LP_NLL_NLP_LL = call_cuhre_dy("resum","gg","full",fitPDF);
			resummed_DY_LP_NLL_NLP_LL_exp_NLO = call_cuhre_dy("resum","gg","NLO",fitPDF);
			resummed_DY_LP_NLL_NLP_LL_exp_NNLO = call_cuhre_dy("resum","gg","NNLO",fitPDF);

			cout << "LP LL + NLP LL" << endl;
			ISNNLL = 0;
			ISNLP = 1;
			ISLL = 1;
			ISNLL = 0;
			resummed_DY_LP_LL_NLP_LL = call_cuhre_dy("resum","gg","full",fitPDF);
			resummed_DY_LP_LL_NLP_LL_exp_NLO = call_cuhre_dy("resum","gg","NLO",fitPDF);
			resummed_DY_LP_LL_NLP_LL_exp_NNLO = call_cuhre_dy("resum","gg","NNLO",fitPDF);

			cout << "LP NNLL" << endl;
			ISNNLL = 1;
			ISNLP = 0;
			ISLL = 1;
			ISNLL = 1;
			resummed_DY_LP_NNLL = call_cuhre_dy("resum","gg","full",fitPDF);
			resummed_DY_LP_NNLL_exp_NLO = call_cuhre_dy("resum","gg","NLO",fitPDF);
			resummed_DY_LP_NNLL_exp_NNLO = call_cuhre_dy("resum","gg","NNLO",fitPDF);

			cout << "LP NLL" << endl;
			ISNNLL = 0;
			ISNLP = 0;
			ISLL = 1;
			ISNLL = 1;
			resummed_DY_LP_NLL = call_cuhre_dy("resum","gg","full",fitPDF);
			resummed_DY_LP_NLL_exp_NLO = call_cuhre_dy("resum","gg","NLO",fitPDF);
			resummed_DY_LP_NLL_exp_NNLO = call_cuhre_dy("resum","gg","NNLO",fitPDF);

			cout << "LP LL" << endl;
			ISNNLL = 0;
			ISNLP = 0;
			ISLL = 1;
			ISNLL = 0;
			resummed_DY_LP_LL = call_cuhre_dy("resum","gg","full",fitPDF);
			resummed_DY_LP_LL_exp_NLO = call_cuhre_dy("resum","gg","NLO",fitPDF);
			resummed_DY_LP_LL_exp_NNLO = call_cuhre_dy("resum","gg","NNLO",fitPDF);


		}
	}
	}
  if(ZZ){
		cout << "working on ZZ" << endl;
		if((PDFmemuse == true) && (fitPDF == false))
		{
			for(int i = 0; i < 101;i++)
			{
				use_member = i;
				if(full){
					Q = 2.*sqrt(mZ2);
					Q2 = pow(Q,2);
					update_defaults();
					zz_LO_full.push_back(call_cuhre_diboson("LOfull","ZZ",fitPDF)[0]);
				}
				if(diff){
					update_defaults();
					zz_LO_diff.push_back(call_cuhre_diboson("LOdiff","ZZ",fitPDF)[0]);
				}
			}
		}
		else
		{
		if(full){
				Q = 2.*sqrt(mZ2);
				Q2 = pow(Q,2);
				update_defaults();
				cout << "working on the total cross section" << endl;
				zz_LO_full = call_cuhre_diboson("LOfull","ZZ",fitPDF);
				if(RES){
					cout << "computing the resummed results (total cross-section)" << endl;
					cout << "LP NNLL + NLP LL" << endl;
					//ISNNLL = 1;
					//ISNLP = 1;
					//ISLL = 1;
					//ISNLL = 1;
					//resummedM_zz_LP_NNLL_NLP_LL = call_cuhre_diboson("resumfullM","ZZ",fitPDF);
					//resummedQ_zz_LP_NNLL_NLP_LL = call_cuhre_diboson("resumfullQ","ZZ",fitPDF);

					cout << "LP NLL + NLP LL" << endl;
					ISNNLL = 0;
					ISNLP = 1;
					ISLL = 1;
					ISNLL = 1;
					resummedM_zz_LP_NLL_NLP_LL = call_cuhre_diboson("resumfullM","ZZ",fitPDF);
					resummedQ_zz_LP_NLL_NLP_LL = call_cuhre_diboson("resumfullQ","ZZ",fitPDF);

					cout << "LP LL + NLP LL" << endl;
					ISNNLL = 0;
					ISNLP = 1;
					ISLL = 1;
					ISNLL = 0;

					resummedM_zz_LP_LL_NLP_LL = call_cuhre_diboson("resumfullM","ZZ",fitPDF);
					resummedQ_zz_LP_LL_NLP_LL = call_cuhre_diboson("resumfullQ","ZZ",fitPDF);

					//cout << "LP NNLL" << endl;
					//ISNNLL = 1;
					//ISNLP = 0;
					//ISLL = 1;
					//ISNLL = 1;
					//resummedM_zz_LP_NNLL = call_cuhre_diboson("resumfullM","ZZ",fitPDF);
					//resummedQ_zz_LP_NNLL = call_cuhre_diboson("resumfullQ","ZZ",fitPDF);
					cout << "LP NLL" << endl;
					ISNNLL = 0;
					ISNLP = 0;
					ISLL = 1;
					ISNLL = 1;
					resummedM_zz_LP_NLL = call_cuhre_diboson("resumfullM","ZZ",fitPDF);
					resummedQ_zz_LP_NLL = call_cuhre_diboson("resumfullQ","ZZ",fitPDF);

					cout << "LP LL" << endl;
					ISNNLL = 0;
					ISNLP = 0;
					ISLL = 1;
					ISNLL = 0;
					resummedM_zz_LP_LL = call_cuhre_diboson("resumfullM","ZZ",fitPDF);
					resummedQ_zz_LP_LL = call_cuhre_diboson("resumfullQ","ZZ",fitPDF);
			}
		}
		if(diff){

		  update_defaults();
			zz_LO_diff = call_cuhre_diboson("LOdiff","ZZ",fitPDF);
			if(RES){
				cout << "computing the resummed results (total cross-section)" << endl;
				//cout << "LP NNLL + NLP LL" << endl;
				//ISNNLL = 1;
				//ISNLP = 1;
				//ISLL = 1;
				//ISNLL = 1;
				//resummed_zz_diff_LP_NNLL_NLP_LL = call_cuhre_diboson("resumdiff","ZZ",fitPDF);

				cout << "LP NLL + NLP LL" << endl;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 1;
				resummed_zz_diff_LP_NLL_NLP_LL = call_cuhre_diboson("resumdiff","ZZ",fitPDF);

				cout << "LP LL + NLP LL" << endl;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 0;

				resummed_zz_diff_LP_LL_NLP_LL = call_cuhre_diboson("resumdiff","ZZ",fitPDF);

				//cout << "LP NNLL" << endl;
				//ISNNLL = 1;
				//ISNLP = 0;
				//ISLL = 1;
				//ISNLL = 1;
				//resummed_zz_diff_LP_NNLL = call_cuhre_diboson("resumdiff","ZZ",fitPDF);
				cout << "LP NLL" << endl;
				ISNNLL = 0;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 1;
				resummed_zz_diff_LP_NLL = call_cuhre_diboson("resumdiff","ZZ",fitPDF);

				cout << "LP LL" << endl;
				ISNNLL = 0;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 0;
				resummed_zz_diff_LP_LL = call_cuhre_diboson("resumdiff","ZZ",fitPDF);
			}
		}
	}
	}
	if(WW){
		cout << "working on WW" << endl;
		if((PDFmemuse == true) && (fitPDF == false))
		{
			for(int i = 0; i < 101;i++)
			{
				use_member = i;
				if(full){
					Q = 2.*sqrt(mW2);
					Q2 = pow(Q,2);
					update_defaults();
					ww_LO_full.push_back(call_cuhre_diboson("LOfull","WW",fitPDF)[0]);
				}
				if(diff){
					update_defaults();
					ww_LO_diff.push_back(call_cuhre_diboson("LOdiff","WW",fitPDF)[0]);
				}
			}
		}
		else{
			if(full){
				Q = 2.*sqrt(mW2);
				Q2 = pow(Q,2);
				update_defaults();
				cout << "working on the total cross section" << endl;
				ww_LO_full = call_cuhre_diboson("LOfull","WW",fitPDF);
				if(RES){
					cout << "computing the resummed results (total cross-section)" << endl;
					//cout << "LP NNLL + NLP LL" << endl;
					//ISNNLL = 1;
					//ISNLP = 1;
					//ISLL = 1;
					//ISNLL = 1;
					//resummedM_ww_LP_NNLL_NLP_LL = call_cuhre_diboson("resumfullM","WW",fitPDF);
					//resummedQ_ww_LP_NNLL_NLP_LL = call_cuhre_diboson("resumfullQ","WW",fitPDF);

					cout << "LP NLL + NLP LL" << endl;
					ISNNLL = 0;
					ISNLP = 1;
					ISLL = 1;
					ISNLL = 1;
					resummedM_ww_LP_NLL_NLP_LL = call_cuhre_diboson("resumfullM","WW",fitPDF);
					resummedQ_ww_LP_NLL_NLP_LL = call_cuhre_diboson("resumfullQ","WW",fitPDF);

					cout << "LP LL + NLP LL" << endl;
					ISNNLL = 0;
					ISNLP = 1;
					ISLL = 1;
					ISNLL = 0;

					resummedM_ww_LP_LL_NLP_LL = call_cuhre_diboson("resumfullM","WW",fitPDF);
					resummedQ_ww_LP_LL_NLP_LL = call_cuhre_diboson("resumfullQ","WW",fitPDF);

					//cout << "LP NNLL" << endl;
					//ISNNLL = 1;
					//ISNLP = 0;
					//ISLL = 1;
					//ISNLL = 1;
					//resummedM_ww_LP_NNLL = call_cuhre_diboson("resumfullM","WW",fitPDF);
					//resummedQ_ww_LP_NNLL = call_cuhre_diboson("resumfullQ","WW",fitPDF);
					cout << "LP NLL" << endl;
					ISNNLL = 0;
					ISNLP = 0;
					ISLL = 1;
					ISNLL = 1;
					resummedM_ww_LP_NLL = call_cuhre_diboson("resumfullM","WW",fitPDF);
					resummedQ_ww_LP_NLL = call_cuhre_diboson("resumfullQ","WW",fitPDF);

					cout << "LP LL" << endl;
					ISNNLL = 0;
					ISNLP = 0;
					ISLL = 1;
					ISNLL = 0;
					resummedM_ww_LP_LL = call_cuhre_diboson("resumfullM","WW",fitPDF);
					resummedQ_ww_LP_LL = call_cuhre_diboson("resumfullQ","WW",fitPDF);
			}
		}
		if(diff){

		  update_defaults();
			ww_LO_diff = call_cuhre_diboson("LOdiff","WW",fitPDF);
			if(RES){
				cout << "computing the resummed results (total cross-section)" << endl;
				//cout << "LP NNLL + NLP LL" << endl;
				//ISNNLL = 1;
				//ISNLP = 1;
				//ISLL = 1;
				//ISNLL = 1;
				//resummed_ww_diff_LP_NNLL_NLP_LL = call_cuhre_diboson("resumdiff","WW",fitPDF);

				cout << "LP NLL + NLP LL" << endl;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 1;
				resummed_ww_diff_LP_NLL_NLP_LL = call_cuhre_diboson("resumdiff","WW",fitPDF);

				cout << "LP LL + NLP LL" << endl;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 0;

				resummed_ww_diff_LP_LL_NLP_LL = call_cuhre_diboson("resumdiff","WW",fitPDF);

				//cout << "LP NNLL" << endl;
				//ISNNLL = 1;
				//ISNLP = 0;
				//ISLL = 1;
				//ISNLL = 1;
				//resummed_ww_diff_LP_NNLL = call_cuhre_diboson("resumdiff","WW",fitPDF);
				cout << "LP NLL" << endl;
				ISNNLL = 0;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 1;
				resummed_ww_diff_LP_NLL = call_cuhre_diboson("resumdiff","WW",fitPDF);

				cout << "LP LL" << endl;
				ISNNLL = 0;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 0;
				resummed_ww_diff_LP_LL = call_cuhre_diboson("resumdiff","WW",fitPDF);
			}
		}
	}
}
	if(hh && !SUSY){
		cout << "working on hh" << endl;
		if((PDFmemuse == true) && (fitPDF == false))
		{
			for(int i = 0; i < 101;i++)
			{
				use_member = i;
				if(full){
					Q = 2.*sqrt(mH2);
					Q2 = pow(Q,2);
					update_defaults();
					hh_LO_full.push_back(call_cuhre_dihiggs("SM","full", fitPDF)[0]);
				}
				if(diff){
					update_defaults();
					hh_LO_diff.push_back(call_cuhre_dihiggs("SM","diff",fitPDF)[0]);
				}
			}
		}
		else
		{
			if(full){
				Q = 2.*sqrt(mH2);
				Q2 = pow(Q,2);
				update_defaults();
				cout << "working on the total cross section" << endl;
				hh_LO_full = call_cuhre_dihiggs("SM","full",fitPDF);
		}
		if(diff){
	  	update_defaults();
			hh_LO_diff = call_cuhre_dihiggs("SM","diff",fitPDF);
			if(RES){
				cout << "computing the resummed results (total cross-section)" << endl;
				//cout << "LP NNLL + NLP LL" << endl;
				//ISNNLL = 1;
				//ISNLP = 1;
				//ISLL = 1;
				//ISNLL = 1;
				//resummed_hh_diff_LP_NNLL_NLP_LL = call_cuhre_dihiggs("SM","resumdiff",fitPDF);

				cout << "LP NLL + NLP LL" << endl;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 1;
				resummed_hh_diff_LP_NLL_NLP_LL = call_cuhre_dihiggs("SM","resumdiff",fitPDF);

				cout << "LP LL + NLP LL" << endl;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 0;

				resummed_hh_diff_LP_LL_NLP_LL = call_cuhre_dihiggs("SM","resumdiff",fitPDF);

				//cout << "LP NNLL" << endl;
				//ISNNLL = 1;
				//ISNLP = 0;
				//ISLL = 1;
				//ISNLL = 1;
				//resummed_hh_diff_LP_NNLL = call_cuhre_dihiggs("SM","resumdiff",fitPDF);
				cout << "LP NLL" << endl;
				ISNNLL = 0;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 1;
				resummed_hh_diff_LP_NLL = call_cuhre_dihiggs("SM","resumdiff",fitPDF);

				cout << "LP LL" << endl;
				ISNNLL = 0;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 0;
				resummed_hh_diff_LP_LL = call_cuhre_dihiggs("SM","resumdiff",fitPDF);
			}
		}
	}
}
	if(hh && SUSY){
		cout << tb << " " << SUSYset << endl;
		cout << "working on hh" << endl;
		set_SUSY(SUSYset,tb,true);
		if(full){
				Q = 2.*sqrt(mH2);
				Q2 = pow(Q,2);
				update_defaults();
				cout << "working on the total cross section" << endl;
				hh_LO_full = call_cuhre_dihiggs("SUSY_hh","full",fitPDF);
		}
		if(diff){
	  	update_defaults();
			hh_LO_diff = call_cuhre_dihiggs("SUSY_hh","diff",fitPDF);
			if(RES){
				cout << "computing the resummed results (total cross-section)" << endl;
				cout << "LP NLL + NLP LL" << endl;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 1;
				resummed_hh_diff_LP_NLL_NLP_LL = call_cuhre_dihiggs("SUSY_hh","resumdiff",fitPDF);

				cout << "LP LL + NLP LL" << endl;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 0;

				resummed_hh_diff_LP_LL_NLP_LL = call_cuhre_dihiggs("SUSY_hh","resumdiff",fitPDF);

				cout << "LP NLL" << endl;
				ISNNLL = 0;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 1;
				resummed_hh_diff_LP_NLL = call_cuhre_dihiggs("SUSY_hh","resumdiff",fitPDF);

				cout << "LP LL" << endl;
				ISNNLL = 0;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 0;
				resummed_hh_diff_LP_LL = call_cuhre_dihiggs("SUSY_hh","resumdiff",fitPDF);
			}
		}
	}


	/////////////////////////
	/// printouts
	/////////////////////////
	output.open(q_str.c_str()); //.c_str() needed to do a constant string conversion
		cout << "Making the printouts" << endl;
	if(higgs){
		output << "======================================================" << endl;
		output << "Higgs results" << endl;
		output << "======================================================" << endl;
		if(PDFmemuse)
		{
			output << "LO, NLOgg, NLOqg, NLOqqbar, NNLOgg, NNLOqg, NNLOqqbar, NNLOqq, NNLOqqp \n" << endl;
			for(int i = 0; i < higgs_LO_gg_full.size(); i++){
			output << "	member "<<i<<" : " <<  higgs_LO_gg_full[i].res << "," <<  higgs_NLO_gg_full[i].res << "," <<  res_higgs_NLO_qg_full[i].res << ","
			<<  res_higgs_NLO_qqbar_full[i].res << "," <<  higgs_NNLO_gg_full[i].res << ","  <<  res_higgs_NNLO_qg_full[i].res << ","
			 <<  res_higgs_NNLO_qqbar_full[i].res << ","  <<  res_higgs_NNLO_qq_full[i].res << ","  <<  res_higgs_NNLO_qqp_full[i].res << endl;

			}
		}
		else
		{
		if(LO&&higgs){

			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "LO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////

			output << "---------------------------------------" << endl;

			output << "Total xsec:					" << higgs_LO_gg_full[0].res << " pb +/- " << higgs_LO_gg_full[0].err << " " << higgs_LO_gg_full[0].prob << endl;
		}
		if(NLO&&higgs){
			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "NLO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////
			output << "gg channel" << endl;
			/////////////////////////////////////
			output << "Z-dependent, including LP (NLO):		" <<  res_higgs_NLO_gg_hard[0].res << " pb +/- " << res_higgs_NLO_gg_hard[0].err <<  endl;
			output << "Constant piece (NLO):				" <<  res_higgs_NLO_gg_delta[0].res << " pb +/- " << res_higgs_NLO_gg_delta[0].err <<  endl;
			output << "Fractional gg xsec (at NLO):			" <<  higgs_NLO_gg_full[0].res << " pb " <<  endl;
			res_higgs_NLO_gg_accum.push_back(res_higgs_NLO_gg_LP[0]);
			output << "	power "<<0<<" : " << res_higgs_NLO_gg_accum[0].res << " pb, fractional: "<< res_higgs_NLO_gg_accum[0].res/higgs_NLO_gg_full[0].res << endl;
			for(int i = 0; i < higgs_NLO_gg_powers.size(); i++){
						res_higgs_NLO_gg_accum[0].res = res_higgs_NLO_gg_accum[0].res+higgs_NLO_gg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NLO_gg_accum[0].res << " pb, fractional: "<< res_higgs_NLO_gg_accum[0].res/higgs_NLO_gg_full[0].res << "; increase is "  << higgs_NLO_gg_powers[i].res << " (+/-"<< higgs_NLO_gg_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qg channel" << endl;
			/////////////////////////////////////
			output << "Fractional qg xsec (at NLO):			" <<  res_higgs_NLO_qg_full[0].res << " pb +/-"  <<  res_higgs_NLO_qg_full[0].err  <<  endl;
			res_higgs_NLO_qg_accum.push_back({0,0,0});
			for(int i = 0; i < higgs_NLO_qg_powers.size(); i++){
						res_higgs_NLO_qg_accum[0].res = res_higgs_NLO_qg_accum[0].res+higgs_NLO_qg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NLO_qg_accum[0].res << " pb, fractional: "<< res_higgs_NLO_qg_accum[0].res/res_higgs_NLO_qg_full[0].res << "; increase is "  << higgs_NLO_qg_powers[i].res << " (+/-"<< higgs_NLO_qg_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qqbar channel" << endl;
			/////////////////////////////////////
			output << "Fractional qqbar xsec (at NLO):			" <<  res_higgs_NLO_qqbar_full[0].res << " pb +/-"  <<  res_higgs_NLO_qqbar_full[0].err  <<  endl;
			res_higgs_NLO_qqbar_accum.push_back({0,0,0});
			for(int i = 0; i < higgs_NLO_qqbar_powers.size(); i++){
						res_higgs_NLO_qqbar_accum[0].res = res_higgs_NLO_qqbar_accum[0].res+higgs_NLO_qqbar_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NLO_qqbar_accum[0].res << " pb, fractional: "<< res_higgs_NLO_qqbar_accum[0].res/res_higgs_NLO_qqbar_full[0].res << "; increase is "  << higgs_NLO_qqbar_powers[i].res << " (+/-"<< higgs_NLO_qqbar_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			output << "......................................................." << endl;
			output << "Total (LO+NLO): " << higgs_LO_gg_full[0].res + higgs_NLO_gg_full[0].res + res_higgs_NLO_qg_full[0].res + res_higgs_NLO_qqbar_full[0].res << endl;
			output << "Only gg channel (LO+NLO): " << higgs_LO_gg_full[0].res + higgs_NLO_gg_full[0].res << endl;
			output << "......................................................." << endl;
			////////////////////////////////////////////////////////////////////
		}
		if(NNLO&&higgs){
			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "NNLO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////
			output << "gg channel" << endl;
			/////////////////////////////////////
			output << "Z-dependent, including LP (NNLO):		" <<  res_higgs_NNLO_gg_hard[0].res << " pb +/- " << res_higgs_NNLO_gg_hard[0].err <<  endl;
			output << "Constant piece (NNLO):				" <<  res_higgs_NNLO_gg_delta[0].res << " pb +/- " << res_higgs_NNLO_gg_delta[0].err <<  endl;
			output << "Fractional gg xsec (at NNLO):			" <<  higgs_NNLO_gg_full[0].res << " pb "  <<  endl;
			res_higgs_NNLO_gg_accum.push_back(res_higgs_NNLO_gg_LP[0]);
			output << "	power "<<0<<" : " << res_higgs_NNLO_gg_accum[0].res << " pb, fractional: "<< res_higgs_NNLO_gg_accum[0].res/higgs_NNLO_gg_full[0].res << endl;
			for(int i = 0; i < higgs_NNLO_gg_powers.size(); i++){
						res_higgs_NNLO_gg_accum[0].res = res_higgs_NNLO_gg_accum[0].res+higgs_NNLO_gg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_gg_accum[0].res << " pb, fractional: "<< res_higgs_NNLO_gg_accum[0].res/higgs_NNLO_gg_full[0].res << "; increase is "  << higgs_NNLO_gg_powers[i].res << " (+/-"<< higgs_NNLO_gg_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qg channel" << endl;
			/////////////////////////////////////
			output << "Fractional qg xsec (at NNLO):			" <<  res_higgs_NNLO_qg_full[0].res << " pb +/-"  <<  res_higgs_NNLO_qg_full[0].err  <<  endl;
			res_higgs_NNLO_qg_accum.push_back({0,0,0});
			for(int i = 0; i < higgs_NNLO_qg_powers.size(); i++){
						res_higgs_NNLO_qg_accum[0].res = res_higgs_NNLO_qg_accum[0].res+higgs_NNLO_qg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_qg_accum[0].res << " pb, fractional: "<< res_higgs_NNLO_qg_accum[0].res/res_higgs_NNLO_qg_full[0].res << "; increase is "  << higgs_NNLO_qg_powers[i].res << " (+/-"<< higgs_NNLO_qg_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qqbar channel" << endl;
			/////////////////////////////////////
			output << "Fractional qqbar xsec (at NNLO):		" <<  res_higgs_NNLO_qqbar_full[0].res << " pb +/-"  <<  res_higgs_NNLO_qqbar_full[0].err  <<  endl;
			res_higgs_NNLO_qqbar_accum.push_back({0,0,0});
			for(int i = 0; i < higgs_NNLO_qqbar_powers.size(); i++){
						res_higgs_NNLO_qqbar_accum[0].res = res_higgs_NNLO_qqbar_accum[0].res+higgs_NNLO_qqbar_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_qqbar_accum[0].res << " pb, fractional: "<< res_higgs_NNLO_qqbar_accum[0].res/res_higgs_NNLO_qqbar_full[0].res << "; increase is "  << higgs_NNLO_qqbar_powers[i].res << " (+/-"<< higgs_NNLO_qqbar_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qq channel" << endl;
			/////////////////////////////////////
			output << "Fractional qq xsec (at NNLO):			" <<  res_higgs_NNLO_qq_full[0].res << " pb +/-"  <<  res_higgs_NNLO_qq_full[0].err  <<  endl;
			res_higgs_NNLO_qq_accum.push_back({0,0,0});
			for(int i = 0; i < higgs_NNLO_qq_powers.size(); i++){
						res_higgs_NNLO_qq_accum[0].res = res_higgs_NNLO_qq_accum[0].res+higgs_NNLO_qq_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_qq_accum[0].res << " pb, fractional: "<< res_higgs_NNLO_qq_accum[0].res/res_higgs_NNLO_qq_full[0].res << "; increase is "  << higgs_NNLO_qq_powers[i].res << " (+/-"<< higgs_NNLO_qq_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qq' channel" << endl;
			/////////////////////////////////////
			output << "Fractional qq' xsec (at NNLO):			" <<  res_higgs_NNLO_qqp_full[0].res << " pb +/-"  <<  res_higgs_NNLO_qqp_full[0].err  <<  endl;
			res_higgs_NNLO_qqp_accum.push_back({0,0,0});
			for(int i = 0; i < higgs_NNLO_qqp_powers.size(); i++){
						res_higgs_NNLO_qqp_accum[0].res = res_higgs_NNLO_qqp_accum[0].res+higgs_NNLO_qqp_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_qqp_accum[0].res << " pb, fractional: "<< res_higgs_NNLO_qqp_accum[0].res/res_higgs_NNLO_qqp_full[0].res << "; increase is "  << higgs_NNLO_qqp_powers[i].res << " (+/-"<< higgs_NNLO_qqp_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			output << "......................................................." << endl;
			output << "Total (LO+NLO+NNLO): " << higgs_LO_gg_full[0].res + higgs_NLO_gg_full[0].res + res_higgs_NLO_qg_full[0].res + res_higgs_NLO_qqbar_full[0].res + higgs_NNLO_gg_full[0].res + res_higgs_NNLO_qg_full[0].res + res_higgs_NNLO_qqbar_full[0].res+ res_higgs_NNLO_qq_full[0].res+ res_higgs_NNLO_qqp_full[0].res << " pb" << endl;
			output << "Only gg channel (LO+NLO+NNLO): " << higgs_LO_gg_full[0].res + higgs_NLO_gg_full[0].res+ higgs_NNLO_gg_full[0].res << " pb" << endl;
			output << "......................................................." << endl;
		}
		if(RES&&higgs){
			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "Resummed results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////

			output << "Resummed (LP NNLL + NLP LL): 		" << resummed_higgs_LP_NNLL_NLP_LL[0].res << " pb +/- " << resummed_higgs_LP_NNLL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NNLL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NNLL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NLL + NLP LL): 		" << resummed_higgs_LP_NLL_NLP_LL[0].res << " pb +/- " << resummed_higgs_LP_NLL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NLL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NLL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NLL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NLL_NLP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP LL + NLP LL): 		" << resummed_higgs_LP_LL_NLP_LL[0].res << " pb +/- " << resummed_higgs_LP_LL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_LL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_LL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_LL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_LL_NLP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NNLL): 		" << resummed_higgs_LP_NNLL[0].res << " pb +/- " << resummed_higgs_LP_NNLL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NNLL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NNLL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NNLL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NNLL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NLL): 		" << resummed_higgs_LP_NLL[0].res << " pb +/- " << resummed_higgs_LP_NLL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NLL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NLL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NLL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NLL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP LL): 		" << resummed_higgs_LP_LL[0].res << " pb +/- " << resummed_higgs_LP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_LL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "======================================================" << endl;
		}
	}
}
	if(DY){
		output << "======================================================" << endl;
		output << "DY results" << endl;
		output << "======================================================" << endl;
		if(PDFmemuse)
		{
			output << "LO, NLOqqbar, NLOqg, NNLOqqbar, NNLOqg, NNLOgg, NNLOqq \n" << endl;
			for(int i = 0; i < DY_LO_qqbar_full.size(); i++){
			output << "	member "<<i<<" : " <<  DY_LO_qqbar_full[i].res << "," <<  DY_NLO_qqbar_full[i].res << "," <<  res_DY_NLO_qg_full[i].res << ","
		     <<  DY_NNLO_qqbar_full[i].res << ","  <<  res_DY_NNLO_qg_full[i].res << ","
			 <<  res_DY_NNLO_gg_full[i].res << ","  <<  res_DY_NNLO_qq_full[i].res << endl;

			}
		}
		else{
		if(LO&&DY){

			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "LO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////

			output << "---------------------------------------" << endl;

			output << "Total xsec:					" << DY_LO_qqbar_full[0].res << " pb +/- " << DY_LO_qqbar_full[0].err <<  endl;
		}
		if(NLO&&DY){

			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "NLO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////
			output << "qqbar channel" << endl;
			/////////////////////////////////////
			output << "Z-dependent, including LP (NLO):		" <<  res_DY_NLO_qqbar_hard[0].res << " pb +/- " << res_DY_NLO_qqbar_hard[0].err <<  endl;
			output << "Constant piece (NLO):				" <<  res_DY_NLO_qqbar_delta[0].res << " pb +/- " << res_DY_NLO_qqbar_delta[0].err <<  endl;
			output << "Fractional qqbar xsec (at NLO):			" <<  DY_NLO_qqbar_full[0].res << " pb " <<  endl;
			res_DY_NLO_qqbar_accum.push_back(res_DY_NLO_qqbar_LP[0]);
			output << "	power "<<0<<" : " << res_DY_NLO_qqbar_accum[0].res << " pb, fractional: "<< res_DY_NLO_qqbar_accum[0].res/DY_NLO_qqbar_full[0].res << endl;
			for(int i = 0; i < DY_NLO_qqbar_powers.size(); i++){
						res_DY_NLO_qqbar_accum[0].res = res_DY_NLO_qqbar_accum[0].res+DY_NLO_qqbar_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NLO_qqbar_accum[0].res << " pb, fractional: "<< res_DY_NLO_qqbar_accum[0].res/DY_NLO_qqbar_full[0].res << "; increase is "  << DY_NLO_qqbar_powers[i].res << " (+/-"<< DY_NLO_qqbar_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
		/////////////////////////////////////
			output << "qg channel" << endl;
			/////////////////////////////////////
			output << "Fractional qg xsec (at NLO):			" <<  res_DY_NLO_qg_full[0].res << " pb +/-"  <<  res_DY_NLO_qg_full[0].err  <<  endl;
			res_DY_NLO_qg_accum.push_back({0,0,0});
			for(int i = 0; i < DY_NLO_qg_powers.size(); i++){
						res_DY_NLO_qg_accum[0].res = res_DY_NLO_qg_accum[0].res+DY_NLO_qg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NLO_qg_accum[0].res << " pb, fractional: "<< res_DY_NLO_qg_accum[0].res/res_DY_NLO_qg_full[0].res << "; increase is "  << DY_NLO_qg_powers[i].res << " (+/-"<< DY_NLO_qg_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			output << "......................................................." << endl;
			output << "Total (LO+NLO): " << DY_LO_qqbar_full[0].res + DY_NLO_qqbar_full[0].res + res_DY_NLO_qg_full[0].res << " pb" << endl;
			output << "Only qqbar channel (LO+NLO): " << DY_LO_qqbar_full[0].res + DY_NLO_qqbar_full[0].res << " pb" << endl;
			output << "......................................................." << endl;
			////////////////////////////////////////////////////////////////////
		}
		if(NNLO&&DY){
			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "NNLO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////
			output << "qqbar channel" << endl;
			/////////////////////////////////////
			output << "Z-dependent, including LP (NNLO):		" <<  res_DY_NNLO_qqbar_hard[0].res << " pb +/- " << res_DY_NNLO_qqbar_hard[0].err <<  endl;
			output << "Constant piece (NNLO):				" <<  res_DY_NNLO_qqbar_delta[0].res << " pb +/- " << res_DY_NNLO_qqbar_delta[0].err <<  endl;
			output << "Fractional qqbar xsec (at NNLO):			" <<  DY_NNLO_qqbar_full[0].res << " pb "  <<  endl;
			res_DY_NNLO_qqbar_accum.push_back(res_DY_NNLO_qqbar_LP[0]);
			output << "	power "<<0<<" : " << res_DY_NNLO_qqbar_accum[0].res << " pb, fractional: "<< res_DY_NNLO_qqbar_accum[0].res/DY_NNLO_qqbar_full[0].res << endl;
			for(int i = 0; i < DY_NNLO_qqbar_powers.size(); i++){
						res_DY_NNLO_qqbar_accum[0].res = res_DY_NNLO_qqbar_accum[0].res+DY_NNLO_qqbar_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NNLO_qqbar_accum[0].res << " pb, fractional: "<< res_DY_NNLO_qqbar_accum[0].res/DY_NNLO_qqbar_full[0].res << "; increase is "  << DY_NNLO_qqbar_powers[i].res << " (+/-"<< DY_NNLO_qqbar_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qg channel" << endl;
			/////////////////////////////////////
			output << "Fractional qg xsec (at NNLO):			" <<  res_DY_NNLO_qg_full[0].res << " pb +/-"  <<  res_DY_NNLO_qg_full[0].err  <<  endl;
			res_DY_NNLO_qg_accum.push_back({0,0,0});
			for(int i = 0; i < DY_NNLO_qg_powers.size(); i++){
						res_DY_NNLO_qg_accum[0].res = res_DY_NNLO_qg_accum[0].res+DY_NNLO_qg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NNLO_qg_accum[0].res << " pb, fractional: "<< res_DY_NNLO_qg_accum[0].res/res_DY_NNLO_qg_full[0].res << "; increase is "  << DY_NNLO_qg_powers[i].res << " (+/-"<< DY_NNLO_qg_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "gg channel" << endl;
			/////////////////////////////////////
			output << "Fractional gg xsec (at NNLO):		" <<  res_DY_NNLO_gg_full[0].res << " pb +/-"  <<  res_DY_NNLO_gg_full[0].err  <<  endl;
			res_DY_NNLO_gg_accum.push_back({0,0,0});
			for(int i = 0; i < DY_NNLO_gg_powers.size(); i++){
						res_DY_NNLO_gg_accum[0].res = res_DY_NNLO_gg_accum[0].res+DY_NNLO_gg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NNLO_gg_accum[0].res << " pb, fractional: "<< res_DY_NNLO_gg_accum[0].res/res_DY_NNLO_gg_full[0].res << "; increase is "  << DY_NNLO_gg_powers[i].res << " (+/-"<< DY_NNLO_gg_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qq channel" << endl;
			/////////////////////////////////////
			output << "Fractional qq xsec (at NNLO):			" <<  res_DY_NNLO_qq_full[0].res << " pb +/-"  <<  res_DY_NNLO_qq_full[0].err  <<  endl;
			res_DY_NNLO_qq_accum.push_back({0,0,0});
			for(int i = 0; i < DY_NNLO_qq_powers.size(); i++){
						res_DY_NNLO_qq_accum[0].res = res_DY_NNLO_qq_accum[0].res+DY_NNLO_qq_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NNLO_qq_accum[0].res << " pb, fractional: "<< res_DY_NNLO_qq_accum[0].res/res_DY_NNLO_qq_full[0].res << "; increase is "  << DY_NNLO_qq_powers[i].res << " (+/-"<< DY_NNLO_qq_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			output << "......................................................." << endl;
			output << "Total (LO+NLO+NNLO): " << DY_LO_qqbar_full[0].res + DY_NLO_qqbar_full[0].res + res_DY_NLO_qg_full[0].res + DY_NNLO_qqbar_full[0].res + + res_DY_NNLO_gg_full[0].res + res_DY_NNLO_qg_full[0].res+ res_DY_NNLO_qq_full[0].res << " pb" << endl;
			output << "Only qqbar channel (LO+NLO+NNLO): " << DY_LO_qqbar_full[0].res + DY_NLO_qqbar_full[0].res + DY_NNLO_qqbar_full[0].res << " pb" << endl;
			output << "......................................................." << endl;
		}
		if(RES&&DY){

			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "Resummed results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////

			output << "Resummed (LP NNLL + NLP LL): 		" << resummed_DY_LP_NNLL_NLP_LL[0].res << " pb +/- " << resummed_DY_LP_NNLL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NNLL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NNLL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NNLL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_NNLL_NLP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NLL + NLP LL): 		" << resummed_DY_LP_NLL_NLP_LL[0].res << " pb +/- " << resummed_DY_LP_NLL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NLL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NLL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NLL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_NLL_NLP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP LL + NLP LL): 		" << resummed_DY_LP_LL_NLP_LL[0].res << " pb +/- " << resummed_DY_LP_LL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_LL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_LL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_LL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_LL_NLP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NNLL): 		" << resummed_DY_LP_NNLL[0].res << " pb +/- " << resummed_DY_LP_NNLL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NNLL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NNLL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NNLL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_NNLL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NLL): 		" << resummed_DY_LP_NLL[0].res << " pb +/- " << resummed_DY_LP_NLL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NLL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NLL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NLL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_NLL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP LL): 		" << resummed_DY_LP_LL[0].res << " pb +/- " << resummed_DY_LP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_LL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_LL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "======================================================" << endl;
		}
	}
}
	if((hh && !SUSY)){
		output << "======================================================" << endl;
		output << "Di-higgs results" << endl;
		output << "======================================================" << endl;
		output << "Q = " << Q << ", muF = " << muF << ", muR = " << muR << endl;
		output << "======================================================" << endl;
		if(PDFmemuse)
		{
			if(full){
				for(int i = 0; i < hh_LO_full.size(); i++){
					output << "	member "<<i<<" : " <<  hh_LO_full[i].res << endl;
				}
			}
			if(diff){
				for(int i = 0; i < hh_LO_diff.size(); i++){
					output << "	member "<<i<<" : " <<  hh_LO_diff[i].res << endl;
				}
			 }
		}
		else
		{
			if(LO){

			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "LO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////

			output << "---------------------------------------" << endl;

			if(full) output << "Total xsec:					" << hh_LO_full[0].res << " pb +/- " << hh_LO_full[0].err << " " << hh_LO_full[0].prob << endl;
			if(diff) output << "Differential xsec:					" << hh_LO_diff[0].res << " pb/GeV +/- " << hh_LO_diff[0].err << " " << hh_LO_diff[0].prob << endl;
		}
		if(RES){
			output << "===================" << endl;
			output << "Resummed results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////

			if(diff){
			//output << "Resummed (LP NNLL + NLP LL): 		" << resummed_hh_diff_LP_NNLL_NLP_LL[0].res << " pb/GeV +/- " << resummed_hh_diff_LP_NNLL_NLP_LL[0].err <<  endl;
			output << "Resummed (LP NLL + NLP LL): 		" << resummed_hh_diff_LP_NLL_NLP_LL[0].res << " pb/GeV +/- " << resummed_hh_diff_LP_NLL_NLP_LL[0].err <<  endl;
			output << "Resummed (LP LL + NLP LL): 		" << resummed_hh_diff_LP_LL_NLP_LL[0].res << " pb/GeV +/- " << resummed_hh_diff_LP_LL_NLP_LL[0].err <<  endl;
			//output << "Resummed (LP NNLL): 		" << resummed_hh_diff_LP_NNLL[0].res << " pb/GeV +/- " << resummed_hh_diff_LP_NNLL[0].err <<  endl;
			output << "Resummed (LP NLL): 		" << resummed_hh_diff_LP_NLL[0].res << " pb/GeV +/- " << resummed_hh_diff_LP_NLL[0].err <<  endl;
			output << "Resummed (LP LL): 		" << resummed_hh_diff_LP_LL[0].res << " pb/GeV +/- " << resummed_hh_diff_LP_LL[0].err <<  endl;
			output << "======================================================" << endl;
			}
		}
	}
}
	if((hh && SUSY)){
		output << "======================================================" << endl;
		output << "Di-higgs SUSY results" << endl;
		output << "======================================================" << endl;
		output << "Q = " << Q << ", muF = " << muF << ", muR = " << muR << endl;
		output << "mH2 = " << sqrt(mH2) << ",mHeavy2 = " << sqrt(mHeavy2) << ",mA2 = " << sqrt(mA2) << endl;
		output << "======================================================" << endl;
		if(PDFmemuse)
		{
			if(full){
				for(int i = 0; i < hh_LO_full.size(); i++){
					output << "	member "<<i<<" : " <<  hh_LO_full[i].res << endl;
				}
			}
			if(diff){
				for(int i = 0; i < hh_LO_diff.size(); i++){
					output << "	member "<<i<<" : " <<  hh_LO_diff[i].res << endl;
				}
			 }
		}
		else
		{
			if(LO){

			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "LO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////

			output << "---------------------------------------" << endl;

			if(full) output << "Total xsec:					" << hh_LO_full[0].res << " pb +/- " << hh_LO_full[0].err << " " << hh_LO_full[0].prob << endl;
			if(diff) output << "Differential xsec:					" << hh_LO_diff[0].res << " pb/GeV +/- " << hh_LO_diff[0].err << " " << hh_LO_diff[0].prob << endl;
		}
		if(RES){
			output << "===================" << endl;
			output << "Resummed results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////

			if(diff){
			//output << "Resummed (LP NNLL + NLP LL): 		" << resummed_hh_diff_LP_NNLL_NLP_LL[0].res << " pb/GeV +/- " << resummed_hh_diff_LP_NNLL_NLP_LL[0].err <<  endl;
			output << "Resummed (LP NLL + NLP LL): 		" << resummed_hh_diff_LP_NLL_NLP_LL[0].res << " pb/GeV +/- " << resummed_hh_diff_LP_NLL_NLP_LL[0].err <<  endl;
			output << "Resummed (LP LL + NLP LL): 		" << resummed_hh_diff_LP_LL_NLP_LL[0].res << " pb/GeV +/- " << resummed_hh_diff_LP_LL_NLP_LL[0].err <<  endl;
			//output << "Resummed (LP NNLL): 		" << resummed_hh_diff_LP_NNLL[0].res << " pb/GeV +/- " << resummed_hh_diff_LP_NNLL[0].err <<  endl;
			output << "Resummed (LP NLL): 		" << resummed_hh_diff_LP_NLL[0].res << " pb/GeV +/- " << resummed_hh_diff_LP_NLL[0].err <<  endl;
			output << "Resummed (LP LL): 		" << resummed_hh_diff_LP_LL[0].res << " pb/GeV +/- " << resummed_hh_diff_LP_LL[0].err <<  endl;
			output << "======================================================" << endl;
			}
		}
	}
}
	if(ZZ){
		output << "======================================================" << endl;
		output << "ZZ results" << endl;
		output << "======================================================" << endl;
		output << "Q = " << Q << ", muF = " << muF << ", muR = " << muR << endl;
		output << "======================================================" << endl;

		if(PDFmemuse)
		{
			if(full){
				for(int i = 0; i < zz_LO_full.size(); i++){
					output << "	member "<<i<<" : " <<  zz_LO_full[i].res << endl;
				}
			}
			if(diff){
				for(int i = 0; i < zz_LO_diff.size(); i++){
					output << "	member "<<i<<" : " <<  zz_LO_diff[i].res << endl;
				}
			 }
		}
		else
		{
			if(LO){

				/////////////////////////////////////
				///
				output << "===================" << endl;
				output << "LO results" << endl;
				output << "===================" << endl;
				///
				/////////////////////////////////////

				output << "---------------------------------------" << endl;

				if(full) output << "Total xsec:					" << zz_LO_full[0].res << " pb +/- " << zz_LO_full[0].err << " " << zz_LO_full[0].prob << endl;
				if(diff) output << "Differential xsec:					" << zz_LO_diff[0].res << " pb/GeV +/- " << zz_LO_diff[0].err << " " << zz_LO_diff[0].prob << endl;
			}
			if(RES){
				output << "===================" << endl;
				output << "Resummed results" << endl;
				output << "===================" << endl;
				///
				/////////////////////////////////////

				if(full){
				//output << "ResummedQ (LP NNLL + NLP LL): 		" << resummedQ_zz_LP_NNLL_NLP_LL[0].res << " pb +/- " << resummedQ_zz_LP_NNLL_NLP_LL[0].err <<  endl;
				output << "ResummedQ (LP NLL + NLP LL): 		" << resummedQ_zz_LP_NLL_NLP_LL[0].res << " pb +/- " << resummedQ_zz_LP_NLL_NLP_LL[0].err <<  endl;
				output << "ResummedQ (LP LL + NLP LL): 		" << resummedQ_zz_LP_LL_NLP_LL[0].res << " pb +/- " << resummedQ_zz_LP_LL_NLP_LL[0].err <<  endl;
				//output << "ResummedQ (LP NNLL): 		" << resummedQ_zz_LP_NNLL[0].res << " pb +/- " << resummedQ_zz_LP_NNLL[0].err <<  endl;
				output << "ResummedQ (LP NLL): 		" << resummedQ_zz_LP_NLL[0].res << " pb +/- " << resummedQ_zz_LP_NLL[0].err <<  endl;
				output << "ResummedQ (LP LL): 		" << resummedQ_zz_LP_LL[0].res << " pb +/- " << resummedQ_zz_LP_LL[0].err <<  endl;
				output << "------------------------------------------------------" << endl;
				//output << "ResummedM (LP NNLL + NLP LL): 		" << resummedM_zz_LP_NNLL_NLP_LL[0].res << " pb +/- " << resummedM_zz_LP_NNLL_NLP_LL[0].err <<  endl;
				output << "ResummedM (LP NLL + NLP LL): 		" << resummedM_zz_LP_NLL_NLP_LL[0].res << " pb +/- " << resummedM_zz_LP_NLL_NLP_LL[0].err <<  endl;
				output << "ResummedM (LP LL + NLP LL): 		" << resummedM_zz_LP_LL_NLP_LL[0].res << " pb +/- " << resummedM_zz_LP_LL_NLP_LL[0].err <<  endl;
				//output << "ResummedM (LP NNLL): 		" << resummedM_zz_LP_NNLL[0].res << " pb +/- " << resummedM_zz_LP_NNLL[0].err <<  endl;
				output << "ResummedM (LP NLL): 		" << resummedM_zz_LP_NLL[0].res << " pb +/- " << resummedM_zz_LP_NLL[0].err <<  endl;
				output << "ResummedM (LP LL): 		" << resummedM_zz_LP_LL[0].res << " pb +/- " << resummedM_zz_LP_LL[0].err <<  endl;
				output << "======================================================" << endl;

				}
				if(diff){
				//output << "Resummed (LP NNLL + NLP LL): 		" << resummed_zz_diff_LP_NNLL_NLP_LL[0].res << " pb/GeV +/- " << resummed_zz_diff_LP_NNLL_NLP_LL[0].err <<  endl;
				output << "Resummed (LP NLL + NLP LL): 		" << resummed_zz_diff_LP_NLL_NLP_LL[0].res << " pb/GeV +/- " << resummed_zz_diff_LP_NLL_NLP_LL[0].err <<  endl;
				output << "Resummed (LP LL + NLP LL): 		" << resummed_zz_diff_LP_LL_NLP_LL[0].res << " pb/GeV +/- " << resummed_zz_diff_LP_LL_NLP_LL[0].err <<  endl;
				//output << "Resummed (LP NNLL): 		" << resummed_zz_diff_LP_NNLL[0].res << " pb/GeV +/- " << resummed_zz_diff_LP_NNLL[0].err <<  endl;
				output << "Resummed (LP NLL): 		" << resummed_zz_diff_LP_NLL[0].res << " pb/GeV +/- " << resummed_zz_diff_LP_NLL[0].err <<  endl;
				output << "Resummed (LP LL): 		" << resummed_zz_diff_LP_LL[0].res << " pb/GeV +/- " << resummed_zz_diff_LP_LL[0].err <<  endl;
				output << "======================================================" << endl;
				}
			}
		}
	}
	if(WW){
			output << "======================================================" << endl;
			output << "WW results" << endl;
			output << "======================================================" << endl;
			output << "Q = " << Q << ", muF = " << muF << ", muR = " << muR << endl;
			output << "======================================================" << endl;
		if(PDFmemuse)
		{
			if(full){
				for(int i = 0; i < ww_LO_full.size(); i++){
					output << "	member "<<i<<" : " <<  ww_LO_full[i].res << endl;
				}
			}
			if(diff){
				for(int i = 0; i < ww_LO_diff.size(); i++){
					output << "	member "<<i<<" : " <<  ww_LO_diff[i].res << endl;
				}
			 }
		}
		else
		{
			if(LO){

				/////////////////////////////////////
				///
				output << "===================" << endl;
				output << "LO results" << endl;
				output << "===================" << endl;
				///
				/////////////////////////////////////

				output << "---------------------------------------" << endl;

				if(full) output << "Total xsec:					" << ww_LO_full[0].res << " pb +/- " << ww_LO_full[0].err << " " << ww_LO_full[0].prob << endl;
				if(diff) output << "Differential xsec:					" << ww_LO_diff[0].res << " pb/GeV +/- " << ww_LO_diff[0].err << " " << ww_LO_diff[0].prob << endl;
			}
			if(RES){
				output << "===================" << endl;
				output << "Resummed results" << endl;
				output << "===================" << endl;
				///
				/////////////////////////////////////

				if(full){
				//output << "ResummedQ (LP NNLL + NLP LL): 		" << resummedQ_ww_LP_NNLL_NLP_LL[0].res << " pb +/- " << resummedQ_ww_LP_NNLL_NLP_LL[0].err <<  endl;
				output << "ResummedQ (LP NLL + NLP LL): 		" << resummedQ_ww_LP_NLL_NLP_LL[0].res << " pb +/- " << resummedQ_ww_LP_NLL_NLP_LL[0].err <<  endl;
				output << "ResummedQ (LP LL + NLP LL): 		" << resummedQ_ww_LP_LL_NLP_LL[0].res << " pb +/- " << resummedQ_ww_LP_LL_NLP_LL[0].err <<  endl;
				//output << "ResummedQ (LP NNLL): 		" << resummedQ_ww_LP_NNLL[0].res << " pb +/- " << resummedQ_ww_LP_NNLL[0].err <<  endl;
				output << "ResummedQ (LP NLL): 		" << resummedQ_ww_LP_NLL[0].res << " pb +/- " << resummedQ_ww_LP_NLL[0].err <<  endl;
				output << "ResummedQ (LP LL): 		" << resummedQ_ww_LP_LL[0].res << " pb +/- " << resummedQ_ww_LP_LL[0].err <<  endl;
				output << "------------------------------------------------------" << endl;
				//output << "ResummedM (LP NNLL + NLP LL): 		" << resummedM_ww_LP_NNLL_NLP_LL[0].res << " pb +/- " << resummedM_ww_LP_NNLL_NLP_LL[0].err <<  endl;
				output << "ResummedM (LP NLL + NLP LL): 		" << resummedM_ww_LP_NLL_NLP_LL[0].res << " pb +/- " << resummedM_ww_LP_NLL_NLP_LL[0].err <<  endl;
				output << "ResummedM (LP LL + NLP LL): 		" << resummedM_ww_LP_LL_NLP_LL[0].res << " pb +/- " << resummedM_ww_LP_LL_NLP_LL[0].err <<  endl;
				//output << "ResummedM (LP NNLL): 		" << resummedM_ww_LP_NNLL[0].res << " pb +/- " << resummedM_ww_LP_NNLL[0].err <<  endl;
				output << "ResummedM (LP NLL): 		" << resummedM_ww_LP_NLL[0].res << " pb +/- " << resummedM_ww_LP_NLL[0].err <<  endl;
				output << "ResummedM (LP LL): 		" << resummedM_ww_LP_LL[0].res << " pb +/- " << resummedM_ww_LP_LL[0].err <<  endl;
				output << "======================================================" << endl;

				}
				if(diff){
				//output << "Resummed (LP NNLL + NLP LL): 		" << resummed_ww_diff_LP_NNLL_NLP_LL[0].res << " pb/GeV +/- " << resummed_ww_diff_LP_NNLL_NLP_LL[0].err <<  endl;
				output << "Resummed (LP NLL + NLP LL): 		" << resummed_ww_diff_LP_NLL_NLP_LL[0].res << " pb/GeV +/- " << resummed_ww_diff_LP_NLL_NLP_LL[0].err <<  endl;
				output << "Resummed (LP LL + NLP LL): 		" << resummed_ww_diff_LP_LL_NLP_LL[0].res << " pb/GeV +/- " << resummed_ww_diff_LP_LL_NLP_LL[0].err <<  endl;
				//output << "Resummed (LP NNLL): 		" << resummed_ww_diff_LP_NNLL[0].res << " pb/GeV +/- " << resummed_ww_diff_LP_NNLL[0].err <<  endl;
				output << "Resummed (LP NLL): 		" << resummed_ww_diff_LP_NLL[0].res << " pb/GeV +/- " << resummed_ww_diff_LP_NLL[0].err <<  endl;
				output << "Resummed (LP LL): 		" << resummed_ww_diff_LP_LL[0].res << " pb/GeV +/- " << resummed_ww_diff_LP_LL[0].err <<  endl;
				output << "======================================================" << endl;
				}
			}
		}
	}
	output.close();
	return 0;
}
