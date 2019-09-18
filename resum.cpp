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

	bool TEST=false;
  double scales[17] = {5000,2500,2000,1500,1000,750,500,400,300,200,150,125,100,75,50,25,10};
  // check whether the factorization scale is there
	unordered_map<double, vector<vector<double>>>::const_iterator got = fitcoeff.find(muF);
	ltini();
	if ( got == fitcoeff.end() ){
		cout << "Setname " << setname <<" is not implemented with factorization scale ";
		cout << "mu_F="<<muF<<"GeV" << endl;
		cout << "scales that are here: ";
		for(int i = 0; i < 17; i++){cout << scales[i] << "GeV, ";}
		cout << endl;
		cout << "Exiting program" << endl;
		exit(0);
		return 0;
	}
	// print of beginning of programme
	update_defaults();
	double z = 0.1;
	results lumi_pdfs;

	double eta = 1.5;
	lumni_params params = {z, Q, 2*Q/S, exp(eta), exp(-eta), 0,0,0};
	params.z = 0.5;
	double j = 10.;

	/*
	ofstream output3;
	ostringstream x_convert3; // need this for the output

	double scaling = muR/Q;
	cout << scaling << endl;
	x_convert3 << scaling;
	string SCAL  = x_convert3.str();
	string q_str3 = "APROX_dihiggs_"+setname+"_"+SCAL+".txt";
	output3.open(q_str3.c_str()); //.c_str() needed to do a constant string conversion
	cout << q_str3 << endl;


	results dihiggs_LO_gg_full;
	vector<results_c> result;
	scaling = 0.5;

 	printf("\n-------------------- diHiggslets --------------------\n");
	double scales2[20] = {255.,300.,325.,350.,375.,400.,425.,450.,475.,500.,550.,600.,650.,700.,800.,900,1000.,1100.,1200.,1400.};
	for(int i = 0; i < 20; i ++){
		Q = scales2[i];
		cout << "Q= " << Q;
		Q2 = pow(Q,2);
		tau = Q2/S2;
		muR = scaling*Q;
		muR2 = pow(muR,2);
		muF = scaling*Q;
		muF2 = pow(muF,2);
		update_defaults(false);
		result = call_cuhre_dihiggs("LOdiff","gg","LP",false);
		cout << " result =" << result[0][0].res << " " << result[0][0].err << " " << result[0].prob << endl;

	}*/
/*
	scaling = 0.5;
	//Q = 2.*mH;
	Q = 400.;
	Q2 = pow(Q,2);
	tau = Q2/S2;
	muR = scaling*Q;
	muR2 = pow(muR,2);
	muF = scaling*Q;
	muF2 = pow(muF,2);
	update_defaults(false);

	result = call_cuhre_dihiggs("resum","gg","LOinversediff",false);
	cout << " TEST resum LO  =" << result[0][0].res << " " << result[0][0].err << " " << result[0].prob << endl;
	result = call_cuhre_dihiggs("resum","gg","diff",false);
	cout << " Resum LO  =" << result[0][0].res << " " << result[0][0].err << " " << result[0].prob << endl;

		result = call_cuhre_dihiggs("LOdiff","gg","test2",true);
		cout << " TEST resum LO  =" << result[0][0].res << " " << result[0][0].err << " " << result[0].prob << endl;

	dihiggs_LO_gg_full = call_vegas(init_vegas_dihiggs("LOapprox"),params);
  cout << "Total cross section " << dihiggs_LO_gg_full[0].res << " fb +/-" << dihiggs_LO_gg_full[0].err << " fb, as(Q)= " << alphas_Q << ", as(muR)= " << alphas_muR << endl;
	dihiggs_LO_gg_full = call_vegas(init_vegas_mellin("dihiggs","LOmellin"),params);
	cout << "Total cross section (fitted) " << dihiggs_LO_gg_full[0].res << " fb +/-" << dihiggs_LO_gg_full[0].err << " fb, as(Q)= " << alphas_Q << ", as(muR)= " << alphas_muR << endl;

	exit(0);
*/



	//////////////////////
	/// TEST FUNCTIONS
	//////////////////////
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
		output2 << "NLO qqbar LP (*alphas/(4pi)): " << DY_NLO_qqbar_plus(z) << endl;
		output2 << "NLO qqbar z-regular: " << DY_NLO_qqbar_reg(z) << endl;
		output2 << "NLO qqbar delta: " << DY_NLO_qqbar_delta() << endl;

		output2 << "===============NLO-qg==========================" << endl;
		output2 << "NLO qg full: " << DY_NLO_qg_full(z) << endl;

		output2 << "===============NNLO-qqbar/qq/qbarqbar==========================" << endl;
		output2 << "NNLO qqbar LP (*alphas/(4pi))**2: " << DY_NNLO_qqbar_plus(z) << endl;
		output2 << "NNLO qqbar NS: " << DY_NNLO_qqbar_NS(z) << endl;
		output2 << "NNLO qqbar delta: " << DY_NNLO_qqbar_delta() << endl;
		output2 << "NNLO delta BB: " << DY_NNLO_BB_full(z) << endl;
		output2 << "NNLO delta BC: " << DY_NNLO_BC_full(z) << endl;
		output2 << "NNLO delta CC: " << DY_NNLO_CC_full(z) << endl;
		output2 << "NNLO delta CD: " << DY_NNLO_CD_full(z) << endl;
		output2 << "NNLO delta CE: " << DY_NNLO_CE_full(z) << endl;
		output2 << "NNLO delta CF: " << DY_NNLO_CF_full(z) << endl;
		output2 << "===============NNLO-qg==========================" << endl;
		output2 << "NNLO qg full: " << DY_NNLO_qg_full(z) << endl;
		output2 << "===============NNLO-gg==========================" << endl;
		output2 << "NNLO gg full: " << DY_NNLO_gg_full(z) << endl;

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
		output2 << "NLO gg LP (*alphas/(pi)): " << higgs_NLO_gg_plus(z) << endl;
		output2 << "NLO gg z-regular: " << higgs_NLO_gg_reg(z) << endl;
		output2 << "NLO gg delta: " << higgs_NLO_gg_delta() << endl;

		output2 << "===============NLO-qg==========================" << endl;
		output2 << "NLO qg full: " << higgs_NLO_qg_full(z) << endl;
		output2 << "===============NLO-qqbar==========================" << endl;
		output2 << "NLO qqbar full: " << higgs_NLO_qqbar_full(z) << endl;


		output2 << "===============NNLO-gg==========================" << endl;
		output2 << "NNLO gg LP (*alphas/(pi))**2: " << higgs_NNLO_gg_plus(z) << endl;
		output2 << "NNLO gg z-regular: " <<  higgs_NNLO_gg_reg(z) + logdep_gg(z) << endl;
		output2 << "NNLO gg delta: " << higgs_NNLO_gg_delta() + logdep_gg_constant() << endl;
		nF = 0;
		output2 << "NNLO gg delta: " << higgs_NNLO_gg_delta() << endl;
		nF = 5;
		output2 << "NNLO gg delta: " << higgs_NNLO_gg_delta() << endl;
		output2 << "===============NNLO-qg==========================" << endl;
		output2 << "NNLO qg full: " << higgs_NNLO_qg_reg(z) +logdep_qg(z) << endl;
		output2 << "===============NNLO-qqbar==========================" << endl;
		output2 << "NNLO qqbar full: " << higgs_NNLO_qqbar_reg(z)+logdep_qqbar(z) << endl;
		output2 << "===============NNLO-qq==========================" << endl;
		output2 << "NNLO qq full: " << higgs_NNLO_qq_reg(z)+logdep_qq(z) << endl;
		output2 << "===============NNLO-qq'==========================" << endl;
		output2 << "NNLO qq' full: " << higgs_NNLO_qqp_reg(z) << endl;


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
	string q_str = "output_Q" + Qstring +"_as"+asstring+"_"+setname;
	if(DY) q_str = q_str+"_DY";
	if(higgs) q_str = q_str+"_Higgs";
	if(LO) q_str = q_str+"_LO";
	if(NLO) q_str = q_str+"_NLO";
	if(NNLO) q_str = q_str+"_NNLO";
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
		update_defaults(false,false);

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
if(DY){
		cout << "working on DY" << endl;
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
	/*if(PF){
			cout << "computing NLO qqbar (prompt photon)" << endl;
		    res_pf_NLO_qqbar_hard = call_vegas(init_vegas_pf("NLO"), params);
			res_pf_NLO_qqbar_LP_part1 = call_vegas(init_vegas_pf("NLO","LP"), params);
			res_pf_NLO_qqbar_LP_cor = call_vegas(init_vegas_pf("NLO","LP_corr"), params);
			res_pf_NLO_qqbar_delta = call_vegas(init_vegas_pf("NLO","delta"), params);
			res_pf_NLO_qqbar_LP[0].res = res_pf_NLO_qqbar_LP_part1[0].res + res_pf_NLO_qqbar_LP_cor[0].res;
			res_pf_NLO_qqbar_LP[0].err = res_pf_NLO_qqbar_LP_part1[0].err + res_pf_NLO_qqbar_LP_cor[0].err;
			res_pf_NLO_qqbar_full[0].res = res_pf_NLO_qqbar_LP[0].res + res_pf_NLO_qqbar_delta[0].res + res_pf_NLO_qqbar_hard[0].res;
			res_pf_NLO_qqbar_full[0].err = res_pf_NLO_qqbar_LP[0].err + res_pf_NLO_qqbar_delta[0].err + res_pf_NLO_qqbar_hard[0].err;
			res_pf_NLO_qqbar_accum[0].res = res_pf_NLO_qqbar_LP[0].res+res_pf_NLO_qqbar_delta[0].res;
			res_pf_NLO_qqbar_accum[0].err = res_pf_NLO_qqbar_LP[0].err+res_pf_NLO_qqbar_delta[0].err;

		}*/


	/////////////////////////
	/// printouts
	/////////////////////////
	output.open(q_str.c_str()); //.c_str() needed to do a constant string conversion
		cout << "Making the printouts" << endl;

	if(higgs){
		output << "======================================================" << endl;
		output << "Higgs results" << endl;
		output << "======================================================" << endl;
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
	if(DY){
		output << "======================================================" << endl;
		output << "DY results" << endl;
		output << "======================================================" << endl;
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





	/*if(PF){
		output << endl << "RESULTS (prompt photon)" << endl;
		if(NLO){
			output << "Total (NLO): " << res_pf_NLO_qqbar_full.res+res_pf_NLO_qqbar_full.res << " pb (" << res_pf_NLO_qqbar_full.err+res_pf_NLO_qqbar_full.err << ")" << endl;
			output << "Total (NLO, gg): " << res_pf_NLO_qqbar_full.res << " pb (" << res_pf_NLO_qqbar_full.err << ")" << endl;
			output << "Total (NLO, LP): " << res_pf_NLO_qqbar_LP.res << " pb (" << res_pf_NLO_qqbar_LP.err << ")" << endl;
	    	output << "Total (NLO, gg delta): " << res_pf_NLO_qqbar_delta.res << " pb (" << res_pf_NLO_qqbar_delta.err << ")" << endl;
			output << "Total (NLO, gg z): " << res_pf_NLO_qqbar_full.res-res_pf_NLO_qqbar_delta.res << " pb (" << -res_pf_NLO_qqbar_delta.err+res_pf_NLO_qqbar_full.err << ")" << endl;
		}
		output << "---------------------------------------------------------------------------" << endl;

	}*/
	output.close();
	return 0;
}
