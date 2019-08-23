#include <cmath>
#include <iostream>
#include <fstream>
#include "cuba.h"
#include <sstream>
#include "deriv_pdf.h"
#include "monte_carlo.h"
#include "mellin_pdf.h"
#include "resum_functions.h"
#include "k_factors_dy.h"
#include "k_factors_higgs.h"
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
    double scales[17] = {5000,2500,2000,1500,1000,750,500,400,300,200,150,125,100,75,50,25,10}; 
    
    // check whether the factorization scale is there
    unordered_map<double, vector<vector<double>>>::const_iterator got = fitcoeff.find(muF);
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
	
	double z;
	double eta = 1.5;
	lumni_params params = {z, Q, 2*Q/S, exp(eta), exp(-eta), 0,0,0};
	params.z = 0.5;
	////////////////////////////////////////////////////////////////////
	
	//////////////////////
	/// LO declarations
	//////////////////////
	// DY
	results DY_LO_qqbar_full;
	// higgs
	results higgs_LO_gg_full;
	////////////////////////////////////////////////////////////////////
	
	//////////////////////
	/// NLO declarations
	//////////////////////
	// DY
	results res_DY_NLO_qqbar_full, res_DY_NLO_qqbar_hard, res_DY_NLO_qqbar_LP, res_DY_NLO_qqbar_LP_part1, res_DY_NLO_qqbar_LP_cor,res_DY_NLO_qqbar_delta, res_DY_NLO_qqbar_accum, DY_NLO_qqbar_full;
	results res_DY_NLO_qg_full, res_DY_NLO_qg_accum;
	vector<results> DY_NLO_qg_powers, DY_NLO_qqbar_powers;
	// higgs
	results res_higgs_NLO_gg_full, res_higgs_NLO_gg_hard, res_higgs_NLO_gg_LP, res_higgs_NLO_gg_LP_part1, res_higgs_NLO_gg_LP_cor,res_higgs_NLO_gg_delta, res_higgs_NLO_gg_accum, higgs_NLO_gg_full;
	results res_higgs_NLO_qg_full, res_higgs_NLO_qg_accum;
	results res_higgs_NLO_qqbar_full,res_higgs_NLO_qqbar_accum;
	vector<results> higgs_NLO_qg_powers, higgs_NLO_gg_powers, higgs_NLO_qqbar_powers;
	// prompt photon (not complete yet)
	results res_pf_NLO_qqbar_hard, res_pf_NLO_qqbar_LP, res_pf_NLO_qqbar_LP_part1, res_pf_NLO_qqbar_LP_cor,res_pf_NLO_qqbar_delta, res_pf_NLO_qqbar_accum, res_pf_NLO_qqbar_full;
	vector<results> pf_NLO_qqbar_powers;
	////////////////////////////////////////////////////////////////////
	
	//////////////////////
	/// NNLO declarations
	//////////////////////
	// DY
	results res_DY_NNLO_qqbar_full, DY_NNLO_qqbar_full, res_DY_NNLO_qqbar_hard, res_DY_NNLO_qqbar_LP, res_DY_NNLO_qqbar_LP_part1, res_DY_NNLO_qqbar_LP_cor,res_DY_NNLO_qqbar_delta, res_DY_NNLO_qqbar_accum;
	results res_DY_NNLO_qg_full, res_DY_NNLO_qg_accum;
	results res_DY_NNLO_qq_full, res_DY_NNLO_qq_accum;
	results res_DY_NNLO_gg_full, res_DY_NNLO_gg_accum;
	vector<results> DY_NNLO_gg_powers, DY_NNLO_qg_powers, DY_NNLO_qq_powers, DY_NNLO_qqbar_powers;
	// higgs
	results res_higgs_NNLO_gg_full, higgs_NNLO_gg_full, res_higgs_NNLO_gg_hard, res_higgs_NNLO_gg_LP, res_higgs_NNLO_gg_LP_part1, res_higgs_NNLO_gg_LP_cor,res_higgs_NNLO_gg_delta, res_higgs_NNLO_gg_accum;
	results res_higgs_NNLO_qg_full, res_higgs_NNLO_qg_accum;
	results res_higgs_NNLO_qq_full, res_higgs_NNLO_qq_accum;
	results res_higgs_NNLO_qqbar_full, res_higgs_NNLO_qqbar_accum;
	results res_higgs_NNLO_qqp_full, res_higgs_NNLO_qqp_accum;
	vector<results> higgs_NNLO_gg_powers, higgs_NNLO_qg_powers, higgs_NNLO_qq_powers, higgs_NNLO_qqp_powers, higgs_NNLO_qqbar_powers;
	////////////////////////////////////////////////////////////////////
	
	///////////////////////
	/// resummed declarations
	///////////////////////
	
	//DY
	results resummed_DY_LP_NNLL, resummed_DY_LP_NNLL_NLP_LL, resummed_DY_LP_NLL, resummed_DY_LP_NLL_NLP_LL, resummed_DY_LP_LL_NLP_LL, resummed_DY_LP_LL;
	results resummed_DY_LP_NNLL_exp_NNLO, resummed_DY_LP_NNLL_NLP_LL_exp_NNLO, resummed_DY_LP_NLL_exp_NNLO, resummed_DY_LP_NLL_NLP_LL_exp_NNLO, resummed_DY_LP_LL_NLP_LL_exp_NNLO, resummed_DY_LP_LL_exp_NNLO;
	results resummed_DY_LP_NNLL_exp_NLO, resummed_DY_LP_NNLL_NLP_LL_exp_NLO, resummed_DY_LP_NLL_exp_NLO, resummed_DY_LP_NLL_NLP_LL_exp_NLO, resummed_DY_LP_LL_NLP_LL_exp_NLO, resummed_DY_LP_LL_exp_NLO;
	
	// higgs
	results resummed_higgs_LP_NNLL, resummed_higgs_LP_NNLL_NLP_LL, resummed_higgs_LP_NLL, resummed_higgs_LP_NLL_NLP_LL, resummed_higgs_LP_LL_NLP_LL, resummed_higgs_LP_LL;
	results resummed_higgs_LP_NNLL_exp_NNLO, resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO, resummed_higgs_LP_NLL_exp_NNLO, resummed_higgs_LP_NLL_NLP_LL_exp_NNLO, resummed_higgs_LP_LL_NLP_LL_exp_NNLO, resummed_higgs_LP_LL_exp_NNLO;
	results resummed_higgs_LP_NNLL_exp_NLO, resummed_higgs_LP_NNLL_NLP_LL_exp_NLO, resummed_higgs_LP_NLL_exp_NLO, resummed_higgs_LP_NLL_NLP_LL_exp_NLO, resummed_higgs_LP_LL_NLP_LL_exp_NLO, resummed_higgs_LP_LL_exp_NLO;
	
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
	if(fitPDF){q_str = q_str+"_fitted_pdfs";}
	if(RES) q_str = q_str+"_resummed";
	q_str = q_str + ".txt";
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
			higgs_LO_gg_full = call_vegas(init_vegas_higgs("LO","full","gg",0,fitPDF), params);
		}
		if(NLO&&higgs){
			cout << "working on NLO" << endl;
			res_higgs_NLO_gg_hard = call_vegas(init_vegas_higgs("NLO","reg","gg",0,fitPDF),params);
			res_higgs_NLO_gg_LP_part1 = call_vegas(init_vegas_higgs("NLO","LP","gg",0,fitPDF),params);
			res_higgs_NLO_gg_LP_cor = call_vegas(init_vegas_higgs("NLO","LP_corr","gg",0,fitPDF),params);
			res_higgs_NLO_gg_delta = call_vegas(init_vegas_higgs("NLO","delta","gg",0,fitPDF),params);
			res_higgs_NLO_qg_full = call_vegas(init_vegas_higgs("NLO","full","qg",0,fitPDF),params);
			res_higgs_NLO_qqbar_full = call_vegas(init_vegas_higgs("NLO","full","qqbar",0,fitPDF),params);
			
			res_higgs_NLO_gg_LP.res = res_higgs_NLO_gg_LP_part1.res + res_higgs_NLO_gg_LP_cor.res;
			res_higgs_NLO_gg_LP.err = res_higgs_NLO_gg_LP_part1.err + res_higgs_NLO_gg_LP_cor.err;
			
			cout << "doing the powers" << endl;
			for(int i =1; i<11;i++){
						params.power = i;
						higgs_NLO_gg_powers.push_back(call_vegas(init_vegas_higgs("NLO","powers","gg", params.power,fitPDF), params));
						higgs_NLO_qg_powers.push_back(call_vegas(init_vegas_higgs("NLO","powers","qg", params.power,fitPDF), params));
						higgs_NLO_qqbar_powers.push_back(call_vegas(init_vegas_higgs("NLO","powers","qqbar", params.power,fitPDF), params));
					}
					
			cout << "finishing powers" << endl;
			params.power = -1;
			higgs_NLO_gg_powers.push_back(call_vegas(init_vegas_higgs("NLO","powers","gg", params.power,fitPDF), params));
			higgs_NLO_qg_powers.push_back(call_vegas(init_vegas_higgs("NLO","powers","qg", params.power,fitPDF), params));
			higgs_NLO_qqbar_powers.push_back(call_vegas(init_vegas_higgs("NLO","powers","qqbar", params.power,fitPDF), params));
			params.power = -2;
			higgs_NLO_gg_powers.push_back(call_vegas(init_vegas_higgs("NLO","powers","gg", params.power,fitPDF), params));
			higgs_NLO_qg_powers.push_back(call_vegas(init_vegas_higgs("NLO","powers","qg", params.power,fitPDF), params));
			higgs_NLO_qqbar_powers.push_back(call_vegas(init_vegas_higgs("NLO","powers","qqbar", params.power,fitPDF), params));
			params.power = 0;
			higgs_NLO_gg_full.res = res_higgs_NLO_gg_hard.res+res_higgs_NLO_gg_LP_cor.res+res_higgs_NLO_gg_delta.res;
			higgs_NLO_gg_full.err = res_higgs_NLO_gg_hard.err+res_higgs_NLO_gg_LP_cor.err+res_higgs_NLO_gg_delta.err;
		}
		if(NNLO&&higgs){
			cout << "working on NNLO" << endl;
			res_higgs_NNLO_gg_hard = call_vegas(init_vegas_higgs("NNLO","reg","gg",0,fitPDF),params);
			res_higgs_NNLO_qg_full = call_vegas(init_vegas_higgs("NNLO","full","qg",0,fitPDF),params);
			res_higgs_NNLO_qq_full = call_vegas(init_vegas_higgs("NNLO","full","qq",0,fitPDF),params);
			res_higgs_NNLO_qqp_full = call_vegas(init_vegas_higgs("NNLO","full","qqp",0,fitPDF),params);
			res_higgs_NNLO_qqbar_full = call_vegas(init_vegas_higgs("NNLO","full","qqbar",0,fitPDF),params);
			res_higgs_NNLO_gg_LP_cor = call_vegas(init_vegas_higgs("NNLO","LP_corr","gg",0,fitPDF),params);
			res_higgs_NNLO_gg_LP_part1 = call_vegas(init_vegas_higgs("NNLO","LP","gg",0,fitPDF),params);
			res_higgs_NNLO_gg_delta = call_vegas(init_vegas_higgs("NNLO","delta","gg",0,fitPDF),params);
			
			res_higgs_NNLO_gg_LP.res = res_higgs_NNLO_gg_LP_part1.res + res_higgs_NNLO_gg_LP_cor.res;
			res_higgs_NNLO_gg_LP.err = res_higgs_NNLO_gg_LP_part1.err + res_higgs_NNLO_gg_LP_cor.err;
			cout << "doing the powers" << endl;
			for(int i =1; i<11;i++){
				
						params.power = i;
						higgs_NNLO_gg_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","gg", params.power,fitPDF), params));
						higgs_NNLO_qg_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","qg", params.power,fitPDF), params));
						higgs_NNLO_qq_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","qq", params.power,fitPDF), params));
						higgs_NNLO_qqbar_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","qqbar", params.power,fitPDF), params));
						higgs_NNLO_qqp_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","qqp", params.power,fitPDF), params));
					}
			cout << "finishing powers" << endl;
			params.power = -1;
			higgs_NNLO_gg_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","gg", params.power,fitPDF), params));
			higgs_NNLO_qg_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","qg", params.power,fitPDF), params));
			higgs_NNLO_qq_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","qq", params.power,fitPDF), params));
			higgs_NNLO_qqbar_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","qqbar", params.power,fitPDF), params));
			higgs_NNLO_qqp_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","qqp", params.power,fitPDF), params));
			params.power = -2;
			higgs_NNLO_gg_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","gg", params.power,fitPDF), params));
			higgs_NNLO_qg_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","qg", params.power,fitPDF), params));
			higgs_NNLO_qq_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","qq", params.power,fitPDF), params));
			higgs_NNLO_qqbar_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","qqbar", params.power,fitPDF), params));
			higgs_NNLO_qqp_powers.push_back(call_vegas(init_vegas_higgs("NNLO","powers","qqp", params.power,fitPDF), params));
			params.power = 0;
			higgs_NNLO_gg_full.res = res_higgs_NNLO_gg_hard.res+res_higgs_NNLO_gg_LP_cor.res+res_higgs_NNLO_gg_delta.res;
			higgs_NNLO_gg_full.err = res_higgs_NNLO_gg_hard.err+res_higgs_NNLO_gg_LP_cor.err+res_higgs_NNLO_gg_delta.err;
		}
		if(RES&&higgs){
			cout << "computing the resummed results" << endl;
			cout << "LP NNLL + NLP LL" << endl;
			ISNNLL = 1;
			ISNLP = 1;
			ISLL = 1;
			ISNLL = 1;
			resummed_higgs_LP_NNLL_NLP_LL = call_vegas(init_vegas_mellin("higgs","resum"),params);
			resummed_higgs_LP_NNLL_NLP_LL_exp_NLO = call_vegas(init_vegas_mellin("higgs","resumexpandedNLO"),params);
			resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO = call_vegas(init_vegas_mellin("higgs","resumexpandedNNLO"),params);
			
			cout << "LP NLL + NLP LL" << endl;
			ISNNLL = 0;
			ISNLP = 1;
			ISLL = 1;
			ISNLL = 1;
			resummed_higgs_LP_NLL_NLP_LL = call_vegas(init_vegas_mellin("higgs","resum"),params);
			resummed_higgs_LP_NLL_NLP_LL_exp_NLO = call_vegas(init_vegas_mellin("higgs","resumexpandedNLO"),params);
			resummed_higgs_LP_NLL_NLP_LL_exp_NNLO = call_vegas(init_vegas_mellin("higgs","resumexpandedNNLO"),params);
				
			cout << "LP LL + NLP LL" << endl;
			ISNNLL = 0;
			ISNLP = 1;
			ISLL = 1;
			ISNLL = 0;
			resummed_higgs_LP_LL_NLP_LL = call_vegas(init_vegas_mellin("higgs","resum"),params);
			resummed_higgs_LP_LL_NLP_LL_exp_NLO = call_vegas(init_vegas_mellin("higgs","resumexpandedNLO"),params);
			resummed_higgs_LP_LL_NLP_LL_exp_NNLO = call_vegas(init_vegas_mellin("higgs","resumexpandedNNLO"),params);
			
			cout << "LP NNLL" << endl;
			ISNNLL = 1;
			ISNLP = 0;
			ISLL = 1;
			ISNLL = 1;
			resummed_higgs_LP_NNLL = call_vegas(init_vegas_mellin("higgs","resum"),params);
			resummed_higgs_LP_NNLL_exp_NLO = call_vegas(init_vegas_mellin("higgs","resumexpandedNLO"),params);
			resummed_higgs_LP_NNLL_exp_NNLO = call_vegas(init_vegas_mellin("higgs","resumexpandedNNLO"),params);
			
			cout << "LP NLL" << endl;
			ISNNLL = 0;
			ISNLP = 0;
			ISLL = 1;
			ISNLL = 1;
			resummed_higgs_LP_NLL = call_vegas(init_vegas_mellin("higgs","resum"),params);
			resummed_higgs_LP_NLL_exp_NLO = call_vegas(init_vegas_mellin("higgs","resumexpandedNLO"),params);
			resummed_higgs_LP_NLL_exp_NNLO = call_vegas(init_vegas_mellin("higgs","resumexpandedNNLO"),params);
			
			cout << "LP LL" << endl;
			ISNNLL = 0;
			ISNLP = 0;
			ISLL = 1;
			ISNLL = 0;
			resummed_higgs_LP_LL = call_vegas(init_vegas_mellin("higgs","resum"),params);
			resummed_higgs_LP_LL_exp_NLO = call_vegas(init_vegas_mellin("higgs","resumexpandedNLO"),params);
			resummed_higgs_LP_LL_exp_NNLO = call_vegas(init_vegas_mellin("higgs","resumexpandedNNLO"),params);
			
		
		}
	}
	if(DY){
	
		cout << "working on DY" << endl;
		if(LO&&DY){
			cout << "working on LO" << endl;
			cout << fitPDF << endl;
			DY_LO_qqbar_full = call_vegas(init_vegas_DY("LO","full","qqbar",0,fitPDF), params);
		}
		if(NLO&&DY){
			cout << "working on NLO" << endl;
			res_DY_NLO_qqbar_hard = call_vegas(init_vegas_DY("NLO","reg","qqbar",0,fitPDF),params);
			res_DY_NLO_qqbar_LP_part1 = call_vegas(init_vegas_DY("NLO","LP","qqbar",0,fitPDF),params);
			res_DY_NLO_qqbar_LP_cor = call_vegas(init_vegas_DY("NLO","LP_corr","qqbar",0,fitPDF),params);
			res_DY_NLO_qqbar_delta = call_vegas(init_vegas_DY("NLO","delta","qqbar",0,fitPDF),params);
			res_DY_NLO_qg_full = call_vegas(init_vegas_DY("NLO","full","qg",0,fitPDF),params);
			
			res_DY_NLO_qqbar_LP.res = res_DY_NLO_qqbar_LP_part1.res + res_DY_NLO_qqbar_LP_cor.res;
			res_DY_NLO_qqbar_LP.err = res_DY_NLO_qqbar_LP_part1.err + res_DY_NLO_qqbar_LP_cor.err;
			
			cout << "doing the powers" << endl;
			for(int i =1; i<11;i++){
						params.power = i;
						DY_NLO_qg_powers.push_back(call_vegas(init_vegas_DY("NLO","powers","qg", params.power,fitPDF), params));
						DY_NLO_qqbar_powers.push_back(call_vegas(init_vegas_DY("NLO","powers","qqbar", params.power,fitPDF), params));
					}
					
			cout << "finishing powers" << endl;
			params.power = -1;
			DY_NLO_qg_powers.push_back(call_vegas(init_vegas_DY("NLO","powers","qg", params.power,fitPDF), params));
			DY_NLO_qqbar_powers.push_back(call_vegas(init_vegas_DY("NLO","powers","qqbar", params.power,fitPDF), params));
			params.power = -2;
			DY_NLO_qg_powers.push_back(call_vegas(init_vegas_DY("NLO","powers","qg", params.power,fitPDF), params));
			DY_NLO_qqbar_powers.push_back(call_vegas(init_vegas_DY("NLO","powers","qqbar", params.power,fitPDF), params));
			params.power = 0;
			DY_NLO_qqbar_full.res = res_DY_NLO_qqbar_hard.res+res_DY_NLO_qqbar_LP_cor.res+res_DY_NLO_qqbar_delta.res;
			DY_NLO_qqbar_full.err = res_DY_NLO_qqbar_hard.err+res_DY_NLO_qqbar_LP_cor.err+res_DY_NLO_qqbar_delta.err;
		}
		if(NNLO&&DY){
			cout << "working on NNLO" << endl;
			res_DY_NNLO_qqbar_hard = call_vegas(init_vegas_DY("NNLO","reg","qqbar",0,fitPDF),params);
			res_DY_NNLO_qg_full = call_vegas(init_vegas_DY("NNLO","full","qg",0,fitPDF),params);
			res_DY_NNLO_gg_full = call_vegas(init_vegas_DY("NNLO","full","gg",0,fitPDF),params);
			res_DY_NNLO_qq_full = call_vegas(init_vegas_DY("NNLO","full","qq",0,fitPDF),params);
			res_DY_NNLO_qqbar_LP_cor = call_vegas(init_vegas_DY("NNLO","LP_corr","qqbar",0,fitPDF),params);
			res_DY_NNLO_qqbar_LP_part1 = call_vegas(init_vegas_DY("NNLO","LP","qqbar",0,fitPDF),params);
			res_DY_NNLO_qqbar_delta = call_vegas(init_vegas_DY("NNLO","delta","qqbar",0,fitPDF),params);
			
			res_DY_NNLO_qqbar_LP.res = res_DY_NNLO_qqbar_LP_part1.res + res_DY_NNLO_qqbar_LP_cor.res;
			res_DY_NNLO_qqbar_LP.err = res_DY_NNLO_qqbar_LP_part1.err + res_DY_NNLO_qqbar_LP_cor.err;
			cout << "doing the powers" << endl;
			for(int i =1; i<11;i++){
				
						params.power = i;
						DY_NNLO_gg_powers.push_back(call_vegas(init_vegas_DY("NNLO","powers","gg", params.power,fitPDF), params));
						DY_NNLO_qg_powers.push_back(call_vegas(init_vegas_DY("NNLO","powers","qg", params.power,fitPDF), params));
						DY_NNLO_qq_powers.push_back(call_vegas(init_vegas_DY("NNLO","powers","qq", params.power,fitPDF), params));
						DY_NNLO_qqbar_powers.push_back(call_vegas(init_vegas_DY("NNLO","powers","qqbar", params.power,fitPDF), params));
					}
			cout << "finishing powers" << endl;
			params.power = -1;
			DY_NNLO_gg_powers.push_back(call_vegas(init_vegas_DY("NNLO","powers","gg", params.power,fitPDF), params));
			DY_NNLO_qg_powers.push_back(call_vegas(init_vegas_DY("NNLO","powers","qg", params.power,fitPDF), params));
			DY_NNLO_qq_powers.push_back(call_vegas(init_vegas_DY("NNLO","powers","qq", params.power,fitPDF), params));
			DY_NNLO_qqbar_powers.push_back(call_vegas(init_vegas_DY("NNLO","powers","qqbar", params.power,fitPDF), params));
			params.power = -2;
			DY_NNLO_gg_powers.push_back(call_vegas(init_vegas_DY("NNLO","powers","gg", params.power,fitPDF), params));
			DY_NNLO_qg_powers.push_back(call_vegas(init_vegas_DY("NNLO","powers","qg", params.power,fitPDF), params));
			DY_NNLO_qq_powers.push_back(call_vegas(init_vegas_DY("NNLO","powers","qq", params.power,fitPDF), params));
			DY_NNLO_qqbar_powers.push_back(call_vegas(init_vegas_DY("NNLO","powers","qqbar", params.power,fitPDF), params));
			params.power = 0;
			DY_NNLO_qqbar_full.res = res_DY_NNLO_qqbar_hard.res+res_DY_NNLO_qqbar_LP_cor.res+res_DY_NNLO_qqbar_delta.res;
			DY_NNLO_qqbar_full.err = res_DY_NNLO_qqbar_hard.err+res_DY_NNLO_qqbar_LP_cor.err+res_DY_NNLO_qqbar_delta.err;
		}
		if(RES&&DY){
			cout << "computing the resummed results" << endl;
			cout << "LP NNLL + NLP LL" << endl;
			ISNNLL = 1;
			ISNLP = 1;
			ISLL = 1;
			ISNLL = 1;
			resummed_DY_LP_NNLL_NLP_LL = call_vegas(init_vegas_mellin("DY","resum"),params);
			resummed_DY_LP_NNLL_NLP_LL_exp_NLO = call_vegas(init_vegas_mellin("DY","resumexpandedNLO"),params);
			resummed_DY_LP_NNLL_NLP_LL_exp_NNLO = call_vegas(init_vegas_mellin("DY","resumexpandedNNLO"),params);
			
			cout << "LP NLL + NLP LL" << endl;
			ISNNLL = 0;
			ISNLP = 1;
			ISLL = 1;
			ISNLL = 1;
			resummed_DY_LP_NLL_NLP_LL = call_vegas(init_vegas_mellin("DY","resum"),params);
			resummed_DY_LP_NLL_NLP_LL_exp_NLO = call_vegas(init_vegas_mellin("DY","resumexpandedNLO"),params);
			resummed_DY_LP_NLL_NLP_LL_exp_NNLO = call_vegas(init_vegas_mellin("DY","resumexpandedNNLO"),params);
				
			cout << "LP LL + NLP LL" << endl;
			ISNNLL = 0;
			ISNLP = 1;
			ISLL = 1;
			ISNLL = 0;
			resummed_DY_LP_LL_NLP_LL = call_vegas(init_vegas_mellin("DY","resum"),params);
			resummed_DY_LP_LL_NLP_LL_exp_NLO = call_vegas(init_vegas_mellin("DY","resumexpandedNLO"),params);
			resummed_DY_LP_LL_NLP_LL_exp_NNLO = call_vegas(init_vegas_mellin("DY","resumexpandedNNLO"),params);
			
			cout << "LP NNLL" << endl;
			ISNNLL = 1;
			ISNLP = 0;
			ISLL = 1;
			ISNLL = 1;
			resummed_DY_LP_NNLL = call_vegas(init_vegas_mellin("DY","resum"),params);
			resummed_DY_LP_NNLL_exp_NLO = call_vegas(init_vegas_mellin("DY","resumexpandedNLO"),params);
			resummed_DY_LP_NNLL_exp_NNLO = call_vegas(init_vegas_mellin("DY","resumexpandedNNLO"),params);
			
			cout << "LP NLL" << endl;
			ISNNLL = 0;
			ISNLP = 0;
			ISLL = 1;
			ISNLL = 1;
			resummed_DY_LP_NLL = call_vegas(init_vegas_mellin("DY","resum"),params);
			resummed_DY_LP_NLL_exp_NLO = call_vegas(init_vegas_mellin("DY","resumexpandedNLO"),params);
			resummed_DY_LP_NLL_exp_NNLO = call_vegas(init_vegas_mellin("DY","resumexpandedNNLO"),params);
			
			cout << "LP LL" << endl;
			ISNNLL = 0;
			ISNLP = 0;
			ISLL = 1;
			ISNLL = 0;
			resummed_DY_LP_LL = call_vegas(init_vegas_mellin("DY","resum"),params);
			resummed_DY_LP_LL_exp_NLO = call_vegas(init_vegas_mellin("DY","resumexpandedNLO"),params);
			resummed_DY_LP_LL_exp_NNLO = call_vegas(init_vegas_mellin("DY","resumexpandedNNLO"),params);
			
		
		}
	}
	if(PF){
			cout << "computing NLO qqbar (prompt photon)" << endl;
		    res_pf_NLO_qqbar_hard = call_vegas(init_vegas_pf("NLO"), params);
			res_pf_NLO_qqbar_LP_part1 = call_vegas(init_vegas_pf("NLO","LP"), params);
			res_pf_NLO_qqbar_LP_cor = call_vegas(init_vegas_pf("NLO","LP_corr"), params);
			res_pf_NLO_qqbar_delta = call_vegas(init_vegas_pf("NLO","delta"), params);
			res_pf_NLO_qqbar_LP.res = res_pf_NLO_qqbar_LP_part1.res + res_pf_NLO_qqbar_LP_cor.res;
			res_pf_NLO_qqbar_LP.err = res_pf_NLO_qqbar_LP_part1.err + res_pf_NLO_qqbar_LP_cor.err;
			res_pf_NLO_qqbar_full.res = res_pf_NLO_qqbar_LP.res + res_pf_NLO_qqbar_delta.res + res_pf_NLO_qqbar_hard.res;
			res_pf_NLO_qqbar_full.err = res_pf_NLO_qqbar_LP.err + res_pf_NLO_qqbar_delta.err + res_pf_NLO_qqbar_hard.err;
			res_pf_NLO_qqbar_accum.res = res_pf_NLO_qqbar_LP.res+res_pf_NLO_qqbar_delta.res;
			res_pf_NLO_qqbar_accum.err = res_pf_NLO_qqbar_LP.err+res_pf_NLO_qqbar_delta.err;

		}
	
	
	/////////////////////////
	/// printouts
	/////////////////////////
	
	output.open(q_str.c_str()); //.c_str() needed to do a constant string conversion
	
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
			
			output << "Total xsec:					" << higgs_LO_gg_full.res << " pb +/- " << higgs_LO_gg_full.err <<  endl;
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
			output << "Z-dependent, including LP (NLO):		" <<  res_higgs_NLO_gg_hard.res << " pb +/- " << res_higgs_NLO_gg_hard.err <<  endl;
			output << "Constant piece (NLO):				" <<  res_higgs_NLO_gg_delta.res << " pb +/- " << res_higgs_NLO_gg_delta.err <<  endl;
			output << "Fractional gg xsec (at NLO):			" <<  higgs_NLO_gg_full.res << " pb " <<  endl;
			res_higgs_NLO_gg_accum.res = res_higgs_NLO_gg_LP.res;
			res_higgs_NLO_gg_accum.err = res_higgs_NLO_gg_LP.err;
			output << "	power "<<0<<" : " << res_higgs_NLO_gg_accum.res << " pb, fractional: "<< res_higgs_NLO_gg_accum.res/higgs_NLO_gg_full.res << endl;		
			for(int i = 0; i < higgs_NLO_gg_powers.size()-2; i++){
						res_higgs_NLO_gg_accum.res = res_higgs_NLO_gg_accum.res+higgs_NLO_gg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NLO_gg_accum.res << " pb, fractional: "<< res_higgs_NLO_gg_accum.res/higgs_NLO_gg_full.res << "; increase is "  << higgs_NLO_gg_powers[i].res << " (+/-"<< higgs_NLO_gg_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << higgs_NLO_gg_powers[higgs_NLO_gg_powers.size()-1].res << " (+/-"<< higgs_NLO_gg_powers[higgs_NLO_gg_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << higgs_NLO_gg_powers[higgs_NLO_gg_powers.size()-2].res << " (+/-"<< higgs_NLO_gg_powers[higgs_NLO_gg_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qg channel" << endl;
			/////////////////////////////////////
			output << "Fractional qg xsec (at NLO):			" <<  res_higgs_NLO_qg_full.res << " pb +/-"  <<  res_higgs_NLO_qg_full.err  <<  endl;
			res_higgs_NLO_qg_accum.res = 0;
			res_higgs_NLO_qg_accum.err = 0;
			for(int i = 0; i < higgs_NLO_qg_powers.size()-2; i++){
						res_higgs_NLO_qg_accum.res = res_higgs_NLO_qg_accum.res+higgs_NLO_qg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NLO_qg_accum.res << " pb, fractional: "<< res_higgs_NLO_qg_accum.res/res_higgs_NLO_qg_full.res << "; increase is "  << higgs_NLO_qg_powers[i].res << " (+/-"<< higgs_NLO_qg_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << higgs_NLO_qg_powers[higgs_NLO_qg_powers.size()-1].res << " (+/-"<< higgs_NLO_qg_powers[higgs_NLO_qg_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << higgs_NLO_qg_powers[higgs_NLO_qg_powers.size()-2].res << " (+/-"<< higgs_NLO_qg_powers[higgs_NLO_qg_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qqbar channel" << endl;
			/////////////////////////////////////
			output << "Fractional qqbar xsec (at NLO):			" <<  res_higgs_NLO_qqbar_full.res << " pb +/-"  <<  res_higgs_NLO_qqbar_full.err  <<  endl;
			res_higgs_NLO_qqbar_accum.res = 0;
			res_higgs_NLO_qqbar_accum.err = 0;
			for(int i = 0; i < higgs_NLO_qqbar_powers.size()-2; i++){
						res_higgs_NLO_qqbar_accum.res = res_higgs_NLO_qqbar_accum.res+higgs_NLO_qqbar_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NLO_qqbar_accum.res << " pb, fractional: "<< res_higgs_NLO_qqbar_accum.res/res_higgs_NLO_qqbar_full.res << "; increase is "  << higgs_NLO_qqbar_powers[i].res << " (+/-"<< higgs_NLO_qqbar_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << higgs_NLO_qqbar_powers[higgs_NLO_qqbar_powers.size()-1].res << " (+/-"<< higgs_NLO_qqbar_powers[higgs_NLO_qqbar_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << higgs_NLO_qqbar_powers[higgs_NLO_qqbar_powers.size()-2].res << " (+/-"<< higgs_NLO_qqbar_powers[higgs_NLO_qqbar_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			output << "......................................................." << endl;
			output << "Total (LO+NLO): " << higgs_LO_gg_full.res + higgs_NLO_gg_full.res + res_higgs_NLO_qg_full.res + res_higgs_NLO_qqbar_full.res << endl;
			output << "Only gg channel (LO+NLO): " << higgs_LO_gg_full.res + higgs_NLO_gg_full.res << endl;
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
			output << "Z-dependent, including LP (NNLO):		" <<  res_higgs_NNLO_gg_hard.res << " pb +/- " << res_higgs_NNLO_gg_hard.err <<  endl;
			output << "Constant piece (NNLO):				" <<  res_higgs_NNLO_gg_delta.res << " pb +/- " << res_higgs_NNLO_gg_delta.err <<  endl;
			output << "Fractional gg xsec (at NNLO):			" <<  higgs_NNLO_gg_full.res << " pb "  <<  endl;
			res_higgs_NNLO_gg_accum.res = res_higgs_NNLO_gg_LP.res;
			res_higgs_NNLO_gg_accum.err = res_higgs_NNLO_gg_LP.err;
			output << "	power "<<0<<" : " << res_higgs_NNLO_gg_accum.res << " pb, fractional: "<< res_higgs_NNLO_gg_accum.res/higgs_NNLO_gg_full.res << endl;	
			for(int i = 0; i < higgs_NNLO_gg_powers.size()-2; i++){
						res_higgs_NNLO_gg_accum.res = res_higgs_NNLO_gg_accum.res+higgs_NNLO_gg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_gg_accum.res << " pb, fractional: "<< res_higgs_NNLO_gg_accum.res/higgs_NNLO_gg_full.res << "; increase is "  << higgs_NNLO_gg_powers[i].res << " (+/-"<< higgs_NNLO_gg_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << higgs_NNLO_gg_powers[higgs_NNLO_gg_powers.size()-1].res << " (+/-"<< higgs_NNLO_gg_powers[higgs_NNLO_gg_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << higgs_NNLO_gg_powers[higgs_NNLO_gg_powers.size()-2].res << " (+/-"<< higgs_NNLO_gg_powers[higgs_NNLO_gg_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qg channel" << endl;
			/////////////////////////////////////
			output << "Fractional qg xsec (at NNLO):			" <<  res_higgs_NNLO_qg_full.res << " pb +/-"  <<  res_higgs_NNLO_qg_full.err  <<  endl;
			res_higgs_NNLO_qg_accum.res = 0;
			res_higgs_NNLO_qg_accum.err = 0;
			for(int i = 0; i < higgs_NNLO_qg_powers.size()-2; i++){
						res_higgs_NNLO_qg_accum.res = res_higgs_NNLO_qg_accum.res+higgs_NNLO_qg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_qg_accum.res << " pb, fractional: "<< res_higgs_NNLO_qg_accum.res/res_higgs_NNLO_qg_full.res << "; increase is "  << higgs_NNLO_qg_powers[i].res << " (+/-"<< higgs_NNLO_qg_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << higgs_NNLO_qg_powers[higgs_NNLO_qg_powers.size()-1].res << " (+/-"<< higgs_NNLO_qg_powers[higgs_NNLO_qg_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << higgs_NNLO_qg_powers[higgs_NNLO_qg_powers.size()-2].res << " (+/-"<< higgs_NNLO_qg_powers[higgs_NNLO_qg_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qqbar channel" << endl;
			/////////////////////////////////////
			output << "Fractional qqbar xsec (at NNLO):		" <<  res_higgs_NNLO_qqbar_full.res << " pb +/-"  <<  res_higgs_NNLO_qqbar_full.err  <<  endl;
			res_higgs_NNLO_qqbar_accum.res = 0;
			res_higgs_NNLO_qqbar_accum.err = 0;
			for(int i = 0; i < higgs_NNLO_qqbar_powers.size()-2; i++){
						res_higgs_NNLO_qqbar_accum.res = res_higgs_NNLO_qqbar_accum.res+higgs_NNLO_qqbar_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_qqbar_accum.res << " pb, fractional: "<< res_higgs_NNLO_qqbar_accum.res/res_higgs_NNLO_qqbar_full.res << "; increase is "  << higgs_NNLO_qqbar_powers[i].res << " (+/-"<< higgs_NNLO_qqbar_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << higgs_NNLO_qqbar_powers[higgs_NNLO_qqbar_powers.size()-1].res << " (+/-"<< higgs_NNLO_qqbar_powers[higgs_NNLO_qqbar_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << higgs_NNLO_qqbar_powers[higgs_NNLO_qqbar_powers.size()-2].res << " (+/-"<< higgs_NNLO_qqbar_powers[higgs_NNLO_qqbar_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qq channel" << endl;
			/////////////////////////////////////
			output << "Fractional qq xsec (at NNLO):			" <<  res_higgs_NNLO_qq_full.res << " pb +/-"  <<  res_higgs_NNLO_qq_full.err  <<  endl;
			res_higgs_NNLO_qq_accum.res = 0;
			res_higgs_NNLO_qq_accum.err = 0;
			for(int i = 0; i < higgs_NNLO_qq_powers.size()-2; i++){
						res_higgs_NNLO_qq_accum.res = res_higgs_NNLO_qq_accum.res+higgs_NNLO_qq_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_qq_accum.res << " pb, fractional: "<< res_higgs_NNLO_qq_accum.res/res_higgs_NNLO_qq_full.res << "; increase is "  << higgs_NNLO_qq_powers[i].res << " (+/-"<< higgs_NNLO_qq_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << higgs_NNLO_qq_powers[higgs_NNLO_qq_powers.size()-1].res << " (+/-"<< higgs_NNLO_qq_powers[higgs_NNLO_qq_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << higgs_NNLO_qq_powers[higgs_NNLO_qq_powers.size()-2].res << " (+/-"<< higgs_NNLO_qq_powers[higgs_NNLO_qq_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qq' channel" << endl;
			/////////////////////////////////////
			output << "Fractional qq' xsec (at NNLO):			" <<  res_higgs_NNLO_qqp_full.res << " pb +/-"  <<  res_higgs_NNLO_qqp_full.err  <<  endl;
			res_higgs_NNLO_qqp_accum.res = 0;
			res_higgs_NNLO_qqp_accum.err = 0;
			for(int i = 0; i < higgs_NNLO_qqp_powers.size()-2; i++){
						res_higgs_NNLO_qqp_accum.res = res_higgs_NNLO_qqp_accum.res+higgs_NNLO_qqp_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_qqp_accum.res << " pb, fractional: "<< res_higgs_NNLO_qqp_accum.res/res_higgs_NNLO_qqp_full.res << "; increase is "  << higgs_NNLO_qqp_powers[i].res << " (+/-"<< higgs_NNLO_qqp_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << higgs_NNLO_qqp_powers[higgs_NNLO_qqp_powers.size()-1].res << " (+/-"<< higgs_NNLO_qqp_powers[higgs_NNLO_qqp_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << higgs_NNLO_qqp_powers[higgs_NNLO_qqp_powers.size()-2].res << " (+/-"<< higgs_NNLO_qqp_powers[higgs_NNLO_qqp_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			output << "......................................................." << endl;
			output << "Total (LO+NLO+NNLO): " << higgs_LO_gg_full.res + higgs_NLO_gg_full.res + res_higgs_NLO_qg_full.res + res_higgs_NLO_qqbar_full.res + higgs_NNLO_gg_full.res + res_higgs_NNLO_qg_full.res + res_higgs_NNLO_qqbar_full.res+ res_higgs_NNLO_qq_full.res+ res_higgs_NNLO_qqp_full.res << " pb" << endl;
			output << "Only gg channel (LO+NLO+NNLO): " << higgs_LO_gg_full.res + higgs_NLO_gg_full.res+ higgs_NNLO_gg_full.res << " pb" << endl;
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
			
			output << "Resummed (LP NNLL + NLP LL): 		" << resummed_higgs_LP_NNLL_NLP_LL.res << " pb +/- " << resummed_higgs_LP_NNLL_NLP_LL.err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NNLL_NLP_LL_exp_NLO.res << " pb +/- " << resummed_higgs_LP_NNLL_NLP_LL_exp_NLO.err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO.res << " pb +/- " << resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO.err <<  endl;
			output << "--------------------------------" << endl;
			
			output << "Resummed (LP NLL + NLP LL): 		" << resummed_higgs_LP_NLL_NLP_LL.res << " pb +/- " << resummed_higgs_LP_NLL_NLP_LL.err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NLL_NLP_LL_exp_NLO.res << " pb +/- " << resummed_higgs_LP_NLL_NLP_LL_exp_NLO.err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NLL_NLP_LL_exp_NNLO.res << " pb +/- " << resummed_higgs_LP_NLL_NLP_LL_exp_NNLO.err <<  endl;
			output << "--------------------------------" << endl;
			
			output << "Resummed (LP LL + NLP LL): 		" << resummed_higgs_LP_LL_NLP_LL.res << " pb +/- " << resummed_higgs_LP_LL_NLP_LL.err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_LL_NLP_LL_exp_NLO.res << " pb +/- " << resummed_higgs_LP_LL_NLP_LL_exp_NLO.err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_LL_NLP_LL_exp_NNLO.res << " pb +/- " << resummed_higgs_LP_LL_NLP_LL_exp_NNLO.err <<  endl;
			output << "--------------------------------" << endl;
			
			output << "Resummed (LP NNLL): 		" << resummed_higgs_LP_NNLL.res << " pb +/- " << resummed_higgs_LP_NNLL.err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NNLL_exp_NLO.res << " pb +/- " << resummed_higgs_LP_NNLL_exp_NLO.err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NNLL_exp_NNLO.res << " pb +/- " << resummed_higgs_LP_NNLL_exp_NNLO.err <<  endl;
			output << "--------------------------------" << endl;
			
			output << "Resummed (LP NLL): 		" << resummed_higgs_LP_NLL.res << " pb +/- " << resummed_higgs_LP_NLL.err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NLL_exp_NLO.res << " pb +/- " << resummed_higgs_LP_NLL_exp_NLO.err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NLL_exp_NNLO.res << " pb +/- " << resummed_higgs_LP_NLL_exp_NNLO.err <<  endl;
			output << "--------------------------------" << endl;
			
			output << "Resummed (LP LL): 		" << resummed_higgs_LP_LL.res << " pb +/- " << resummed_higgs_LP_LL.err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_LL_exp_NLO.res << " pb +/- " << resummed_higgs_LP_LL_exp_NLO.err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_LL_exp_NNLO.res << " pb +/- " << resummed_higgs_LP_LL_exp_NNLO.err <<  endl;
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
			
			output << "Total xsec:					" << DY_LO_qqbar_full.res << " pb +/- " << DY_LO_qqbar_full.err <<  endl;
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
			output << "Z-dependent, including LP (NLO):		" <<  res_DY_NLO_qqbar_hard.res << " pb +/- " << res_DY_NLO_qqbar_hard.err <<  endl;
			output << "Constant piece (NLO):				" <<  res_DY_NLO_qqbar_delta.res << " pb +/- " << res_DY_NLO_qqbar_delta.err <<  endl;
			output << "Fractional qqbar xsec (at NLO):			" <<  DY_NLO_qqbar_full.res << " pb " <<  endl;
			res_DY_NLO_qqbar_accum.res = res_DY_NLO_qqbar_LP.res;
			res_DY_NLO_qqbar_accum.err = res_DY_NLO_qqbar_LP.err;
			output << "	power "<<0<<" : " << res_DY_NLO_qqbar_accum.res << " pb, fractional: "<< res_DY_NLO_qqbar_accum.res/DY_NLO_qqbar_full.res << endl;		
			for(int i = 0; i < DY_NLO_qqbar_powers.size()-2; i++){
						res_DY_NLO_qqbar_accum.res = res_DY_NLO_qqbar_accum.res+DY_NLO_qqbar_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NLO_qqbar_accum.res << " pb, fractional: "<< res_DY_NLO_qqbar_accum.res/DY_NLO_qqbar_full.res << "; increase is "  << DY_NLO_qqbar_powers[i].res << " (+/-"<< DY_NLO_qqbar_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << DY_NLO_qqbar_powers[DY_NLO_qqbar_powers.size()-1].res << " (+/-"<< DY_NLO_qqbar_powers[DY_NLO_qqbar_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << DY_NLO_qqbar_powers[DY_NLO_qqbar_powers.size()-2].res << " (+/-"<< DY_NLO_qqbar_powers[DY_NLO_qqbar_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qg channel" << endl;
			/////////////////////////////////////
			output << "Fractional qg xsec (at NLO):			" <<  res_DY_NLO_qg_full.res << " pb +/-"  <<  res_DY_NLO_qg_full.err  <<  endl;
			res_DY_NLO_qg_accum.res = 0;
			res_DY_NLO_qg_accum.err = 0;
			for(int i = 0; i < DY_NLO_qg_powers.size()-2; i++){
						res_DY_NLO_qg_accum.res = res_DY_NLO_qg_accum.res+DY_NLO_qg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NLO_qg_accum.res << " pb, fractional: "<< res_DY_NLO_qg_accum.res/res_DY_NLO_qg_full.res << "; increase is "  << DY_NLO_qg_powers[i].res << " (+/-"<< DY_NLO_qg_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << DY_NLO_qg_powers[DY_NLO_qg_powers.size()-1].res << " (+/-"<< DY_NLO_qg_powers[DY_NLO_qg_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << DY_NLO_qg_powers[DY_NLO_qg_powers.size()-2].res << " (+/-"<< DY_NLO_qg_powers[DY_NLO_qg_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			output << "......................................................." << endl;
			output << "Total (LO+NLO): " << DY_LO_qqbar_full.res + DY_NLO_qqbar_full.res + res_DY_NLO_qg_full.res << " pb" << endl;
			output << "Only qqbar channel (LO+NLO): " << DY_LO_qqbar_full.res + DY_NLO_qqbar_full.res << " pb" << endl;
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
			output << "Z-dependent, including LP (NNLO):		" <<  res_DY_NNLO_qqbar_hard.res << " pb +/- " << res_DY_NNLO_qqbar_hard.err <<  endl;
			output << "Constant piece (NNLO):				" <<  res_DY_NNLO_qqbar_delta.res << " pb +/- " << res_DY_NNLO_qqbar_delta.err <<  endl;
			output << "Fractional qqbar xsec (at NNLO):			" <<  DY_NNLO_qqbar_full.res << " pb "  <<  endl;
			res_DY_NNLO_qqbar_accum.res = res_DY_NNLO_qqbar_LP.res;
			res_DY_NNLO_qqbar_accum.err = res_DY_NNLO_qqbar_LP.err;
			output << "	power "<<0<<" : " << res_DY_NNLO_qqbar_accum.res << " pb, fractional: "<< res_DY_NNLO_qqbar_accum.res/DY_NNLO_qqbar_full.res << endl;	
			for(int i = 0; i < DY_NNLO_qqbar_powers.size()-2; i++){
						res_DY_NNLO_qqbar_accum.res = res_DY_NNLO_qqbar_accum.res+DY_NNLO_qqbar_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NNLO_qqbar_accum.res << " pb, fractional: "<< res_DY_NNLO_qqbar_accum.res/DY_NNLO_qqbar_full.res << "; increase is "  << DY_NNLO_qqbar_powers[i].res << " (+/-"<< DY_NNLO_qqbar_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << DY_NNLO_qqbar_powers[DY_NNLO_qqbar_powers.size()-1].res << " (+/-"<< DY_NNLO_qqbar_powers[DY_NNLO_qqbar_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << DY_NNLO_qqbar_powers[DY_NNLO_qqbar_powers.size()-2].res << " (+/-"<< DY_NNLO_qqbar_powers[DY_NNLO_qqbar_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qg channel" << endl;
			/////////////////////////////////////
			output << "Fractional qg xsec (at NNLO):			" <<  res_DY_NNLO_qg_full.res << " pb +/-"  <<  res_DY_NNLO_qg_full.err  <<  endl;
			res_DY_NNLO_qg_accum.res = 0;
			res_DY_NNLO_qg_accum.err = 0;
			for(int i = 0; i < DY_NNLO_qg_powers.size()-2; i++){
						res_DY_NNLO_qg_accum.res = res_DY_NNLO_qg_accum.res+DY_NNLO_qg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NNLO_qg_accum.res << " pb, fractional: "<< res_DY_NNLO_qg_accum.res/res_DY_NNLO_qg_full.res << "; increase is "  << DY_NNLO_qg_powers[i].res << " (+/-"<< DY_NNLO_qg_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << DY_NNLO_qg_powers[DY_NNLO_qg_powers.size()-1].res << " (+/-"<< DY_NNLO_qg_powers[DY_NNLO_qg_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << DY_NNLO_qg_powers[DY_NNLO_qg_powers.size()-2].res << " (+/-"<< DY_NNLO_gg_powers[DY_NNLO_qg_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			/////////////////////////////////////
			output << "gg channel" << endl;
			/////////////////////////////////////
			output << "Fractional gg xsec (at NNLO):		" <<  res_DY_NNLO_gg_full.res << " pb +/-"  <<  res_DY_NNLO_gg_full.err  <<  endl;
			res_DY_NNLO_gg_accum.res = 0;
			res_DY_NNLO_gg_accum.err = 0;
			for(int i = 0; i < DY_NNLO_gg_powers.size()-2; i++){
						res_DY_NNLO_gg_accum.res = res_DY_NNLO_gg_accum.res+DY_NNLO_gg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NNLO_gg_accum.res << " pb, fractional: "<< res_DY_NNLO_gg_accum.res/res_DY_NNLO_gg_full.res << "; increase is "  << DY_NNLO_gg_powers[i].res << " (+/-"<< DY_NNLO_gg_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << DY_NNLO_gg_powers[DY_NNLO_gg_powers.size()-1].res << " (+/-"<< DY_NNLO_gg_powers[DY_NNLO_gg_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << DY_NNLO_gg_powers[DY_NNLO_gg_powers.size()-2].res << " (+/-"<< DY_NNLO_gg_powers[DY_NNLO_gg_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qq channel" << endl;
			/////////////////////////////////////
			output << "Fractional qq xsec (at NNLO):			" <<  res_DY_NNLO_qq_full.res << " pb +/-"  <<  res_DY_NNLO_qq_full.err  <<  endl;
			res_DY_NNLO_qq_accum.res = 0;
			res_DY_NNLO_qq_accum.err = 0;
			for(int i = 0; i < DY_NNLO_qq_powers.size()-2; i++){
						res_DY_NNLO_qq_accum.res = res_DY_NNLO_qq_accum.res+DY_NNLO_qq_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NNLO_qq_accum.res << " pb, fractional: "<< res_DY_NNLO_qq_accum.res/res_DY_NNLO_qq_full.res << "; increase is "  << DY_NNLO_qq_powers[i].res << " (+/-"<< DY_NNLO_qq_powers[i].err << ") pb" << endl;
					}
			output << "	Tot powers  : " << DY_NNLO_qq_powers[DY_NNLO_qq_powers.size()-1].res << " (+/-"<< DY_NNLO_qq_powers[DY_NNLO_qq_powers.size()-1].err << ") pb" << endl;
			output << "	Full - Tot powers  : " << DY_NNLO_qq_powers[DY_NNLO_qq_powers.size()-2].res << " (+/-"<< DY_NNLO_gg_powers[DY_NNLO_qq_powers.size()-2].err << ") pb" << endl;
			output << "---------" << endl;
			output << "......................................................." << endl;
			output << "Total (LO+NLO+NNLO): " << DY_LO_qqbar_full.res + DY_NLO_qqbar_full.res + res_DY_NLO_qg_full.res + DY_NNLO_qqbar_full.res + + res_DY_NNLO_gg_full.res + res_DY_NNLO_qg_full.res+ res_DY_NNLO_qq_full.res << " pb" << endl;
			output << "Only qqbar channel (LO+NLO+NNLO): " << DY_LO_qqbar_full.res + DY_NLO_qqbar_full.res + DY_NNLO_qqbar_full.res << " pb" << endl;
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
			
			output << "Resummed (LP NNLL + NLP LL): 		" << resummed_DY_LP_NNLL_NLP_LL.res << " pb +/- " << resummed_DY_LP_NNLL_NLP_LL.err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NNLL_NLP_LL_exp_NLO.res << " pb +/- " << resummed_DY_LP_NNLL_NLP_LL_exp_NLO.err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NNLL_NLP_LL_exp_NNLO.res << " pb +/- " << resummed_DY_LP_NNLL_NLP_LL_exp_NNLO.err <<  endl;
			output << "--------------------------------" << endl;
			
			output << "Resummed (LP NLL + NLP LL): 		" << resummed_DY_LP_NLL_NLP_LL.res << " pb +/- " << resummed_DY_LP_NLL_NLP_LL.err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NLL_NLP_LL_exp_NLO.res << " pb +/- " << resummed_DY_LP_NLL_NLP_LL_exp_NLO.err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NLL_NLP_LL_exp_NNLO.res << " pb +/- " << resummed_DY_LP_NLL_NLP_LL_exp_NNLO.err <<  endl;
			output << "--------------------------------" << endl;
			
			output << "Resummed (LP LL + NLP LL): 		" << resummed_DY_LP_LL_NLP_LL.res << " pb +/- " << resummed_DY_LP_LL_NLP_LL.err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_LL_NLP_LL_exp_NLO.res << " pb +/- " << resummed_DY_LP_LL_NLP_LL_exp_NLO.err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_LL_NLP_LL_exp_NNLO.res << " pb +/- " << resummed_DY_LP_LL_NLP_LL_exp_NNLO.err <<  endl;
			output << "--------------------------------" << endl;
			
			output << "Resummed (LP NNLL): 		" << resummed_DY_LP_NNLL.res << " pb +/- " << resummed_DY_LP_NNLL.err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NNLL_exp_NLO.res << " pb +/- " << resummed_DY_LP_NNLL_exp_NLO.err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NNLL_exp_NNLO.res << " pb +/- " << resummed_DY_LP_NNLL_exp_NNLO.err <<  endl;
			output << "--------------------------------" << endl;
			
			output << "Resummed (LP NLL): 		" << resummed_DY_LP_NLL.res << " pb +/- " << resummed_DY_LP_NLL.err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NLL_exp_NLO.res << " pb +/- " << resummed_DY_LP_NLL_exp_NLO.err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NLL_exp_NNLO.res << " pb +/- " << resummed_DY_LP_NLL_exp_NNLO.err <<  endl;
			output << "--------------------------------" << endl;
			
			output << "Resummed (LP LL): 		" << resummed_DY_LP_LL.res << " pb +/- " << resummed_DY_LP_LL.err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_LL_exp_NLO.res << " pb +/- " << resummed_DY_LP_LL_exp_NLO.err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_LL_exp_NNLO.res << " pb +/- " << resummed_DY_LP_LL_exp_NNLO.err <<  endl;
			output << "--------------------------------" << endl;
			
			output << "======================================================" << endl;
		}
	}
	
	
	


	if(PF){
		output << endl << "RESULTS (prompt photon)" << endl;
		if(NLO){
			output << "Total (NLO): " << res_pf_NLO_qqbar_full.res+res_pf_NLO_qqbar_full.res << " pb (" << res_pf_NLO_qqbar_full.err+res_pf_NLO_qqbar_full.err << ")" << endl;
			output << "Total (NLO, gg): " << res_pf_NLO_qqbar_full.res << " pb (" << res_pf_NLO_qqbar_full.err << ")" << endl;
			output << "Total (NLO, LP): " << res_pf_NLO_qqbar_LP.res << " pb (" << res_pf_NLO_qqbar_LP.err << ")" << endl;
	    	output << "Total (NLO, gg delta): " << res_pf_NLO_qqbar_delta.res << " pb (" << res_pf_NLO_qqbar_delta.err << ")" << endl;
			output << "Total (NLO, gg z): " << res_pf_NLO_qqbar_full.res-res_pf_NLO_qqbar_delta.res << " pb (" << -res_pf_NLO_qqbar_delta.err+res_pf_NLO_qqbar_full.err << ")" << endl;
		}
		output << "---------------------------------------------------------------------------" << endl;

		}
	output.close();
	return 0;
}
