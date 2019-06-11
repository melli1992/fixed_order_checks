#include <cmath>
#include <iostream>
#include <fstream>
#include "cuba.h"
#include "deriv_pdf.h"
#include "monte_carlo.h"
#include "k_factors_dy.h"
#include "k_factors_higgs.h"
#include "k_factors_nnlo_dy.h"
#include "parameters.h"
#include "LHAPDF/LHAPDF.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <string.h>
#include <sstream>
#include "polygamma.h"

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
	LHAPDF::PDFSet setk(setname);
	print_defaults();
	int nmem(0.); //number of members
	vector<int> pids; //number of flavors, span from -5 to 5 with 0 = 21 gluon
	nmem = setk.size()-1;
	pdfs = setk.mkPDFs();
	pids = pdfs[0]->flavors();
	xmin_pdfs = pdfs[0]->xMin();
	xmax_pdfs = pdfs[0]->xMax();
	alphas_muF = pdfs[0]->alphasQ(muF);
	alphas_Q = pdfs[0]->alphasQ(Q);
	alphas_muR = pdfs[0]->alphasQ(muR);
	cout << "alphas_Q " << alphas_Q << endl;
	cout << "alphas_muF " << alphas_muF << endl;
	cout << "alphas_muR " << alphas_muR << endl;
	S12(1.-0.1);
	double z;
	lumni_params params = {z, pdfs};
	params.z = 0.5;
	///////////////
	bool DY = true, higgs = true;
	bool TEST=true;
	//LO
	bool LO = true;
	results LO_qqbar_full;
	results higgs_LO_gg_full;
	//NLO
	bool NLO = true;
	results NLO_qqbar_hard, NLO_qqbar_LP_part1, NLO_qqbar_LP_cor, NLO_qqbar_LPint, NLO_qqbar_NLP, NLO_qqbar_NNLP,NLO_qqbar_delta, NLO_qqbar_LP, NLO_qqbar_full, res_NLO_qg_full,res_NLO_qg_NLP,res_NLO_qg_NNLP,res_NLO_qg_NNNLP;
	results res_higgs_NLO_gg_hard, res_higgs_NLO_gg_LP, res_higgs_NLO_gg_LP_part1,res_higgs_NLO_gg_LP_cor, res_higgs_NLO_gg_NLP, res_higgs_NLO_gg_NNLP,res_higgs_NLO_gg_NNNLP,res_higgs_NLO_gg_delta, res_higgs_NLO_gg_full, res_higgs_NLO_qg_full,res_higgs_NLO_qg_NLP,res_higgs_NLO_qg_NNLP,res_higgs_NLO_qg_NNNLP;
	//NNLO
	bool NNLO = false;
	results res_NNLO_qqbar_hard, res_NNLO_qqbar_LP_part1, res_NNLO_qqbar_LP_cor, res_NNLO_qqbar_LPint, res_NNLO_qqbar_NLP, res_NNLO_qqbar_NNLP,res_NNLO_qqbar_NNNLP,res_NNLO_qqbar_delta, res_NNLO_qqbar_LP, res_NNLO_qqbar_full, res_NNLO_qg_full,res_NNLO_qg_NLP,res_NNLO_qg_NNLP,res_NNLO_qg_NNNLP, res_NNLO_gg_full,res_NNLO_gg_NLP,res_NNLO_gg_NNLP,res_NNLO_gg_NNNLP, res_NNLO_qq_full,res_NNLO_qq_NLP,res_NNLO_qq_NNLP,res_NNLO_qq_NNNLP, res_NNLO_qqbarNI_full,res_NNLO_qqbarNI_NLP,res_NNLO_qqbarNI_NNLP,res_NNLO_qqbarNI_NNNLP,res_NNLO_qqNI_full,res_NNLO_qqNI_NLP,res_NNLO_qqNI_NNLP,res_NNLO_qqNI_NNNLP,res_NNLO_qbarqbarNI_full,res_NNLO_qbarqbarNI_NLP,res_NNLO_qbarqbarNI_NNLP,res_NNLO_qbarqbarNI_NNNLP,res_NNLO_qbarqbar_full,res_NNLO_qbarqbar_NLP,res_NNLO_qbarqbar_NNLP,res_NNLO_qbarqbar_NNNLP;
	ofstream output;
	string q_str = "output_Q" + to_string(Q) +"_as"+to_string(alphas_Q)+"_"+setname;
	if(DY) q_str = q_str+"_DY";
	if(higgs) q_str = q_str+"_Higgs";
	if(LO) q_str = q_str+"_LO";
	if(NLO) q_str = q_str+"_NLO";
	if(NNLO) q_str = q_str+"_NNLO";
	q_str = q_str + ".txt";
	output.open(q_str.c_str()); //.c_str() needed to do a constant string conversion
	//////////////////////////////////////////////
	
	
	if(TEST && DY){
		string q_str = "test_values_" + to_string(Q) +"_as"+to_string(alphas_Q)+"_"+setname+".txt";
		ofstream output2;
		output2.open(q_str.c_str()); //.c_str() needed to do a constant string conversion
		cout << "---------- TESTING PARAMETERS -----------" << endl;
		for(int i =1; i < 10; i++)
		{ 
			double z = (i)*0.1;
			output2 << "=====================================" << endl;
			output2 << "z = " << z << endl;
			output2 << "NNLO qqbar LP = " << NNLO_qqbar_LP(z) << endl;
			output2 << "NNLO qqbar b0 LP = " << NNLO_qqbar_b0_LP(z) << endl;
			output2 << "NNLO qqbar NLP = " << NNLO_qqbar_NLP(z) << endl;
			output2 << "NNLO qqbar NNLP = " << NNLO_qqbar_NNLP(z) << endl;
			output2 << "NNLO qqbar NNNLP = " << NNLO_qqbar_NNNLP(z) << endl;
			output2 << "NNLO qqbar HCA = " << NNLO_qqbar_HCA(z) << endl;
			output2 << "NNLO qqbar HCF = " << NNLO_qqbar_HCF(z) << endl;
			output2 << "NNLO qqbar A2 = " << NNLO_qqbar_A2(z) << endl;
			output2 << "NNLO qqbar B2 = " << NNLO_qqbar_B2(z) << endl;
			output2 << "NNLO qqbar AC = " << NNLO_qqbar_AC(z) << endl;
			output2 << "NNLO qqbar BC = " << NNLO_qqbar_BC(z) << endl;
			output2 << "NNLO qqbar delta = " << NNLO_qqbar_delta() << endl;
			output2 << "NNLO qqbar b0H = " << NNLO_qqbar_b0_nonconst(z) << endl;
			output2 << "NNLO qqbar b0 delta = " << NNLO_qqbar_b0_const() << endl;
			output2 << "NNLO qg full = " << NNLO_qg_full(z) << endl;
			output2 << "NNLO qg NLP = " << NNLO_qg_NLP(z) << endl;
			output2 << "NNLO qg NNLP = " << NNLO_qg_NNLP(z) << endl;
			output2 << "NNLO qg NNNLP = " << NNLO_qg_NNNLP(z) << endl;
			output2 << "NNLO gg full = " << NNLO_gg_full(z) << endl;
			output2 << "NNLO gg NLP = " << NNLO_gg_NLP(z) << endl;
			output2 << "NNLO gg NNLP = " << NNLO_gg_NNLP(z) << endl;
			output2 << "NNLO gg NNNLP = " << NNLO_gg_NNNLP(z) << endl;
			output2 << "NNLO qq (=qbarqbar) full = " << NNLO_qq_full(z) << endl;
			output2 << "NNLO qq (=qbarqbar)  NLP = " << NNLO_qq_NLP(z) << endl;
			output2 << "NNLO qq (=qbarqbar)  NNLP = " << NNLO_qq_NNLP(z) << endl;
			output2 << "NNLO qq (=qbarqbar)  NNNLP = " << NNLO_qq_NNNLP(z) << endl;
			output2 << "NNLO qqNI (=qbarqbarNI) full = " << NNLO_qqNI_full(z) << endl;
			output2 << "NNLO qqNI (=qbarqbarNI) NLP = " << NNLO_qqNI_NLP(z) << endl;
			output2 << "NNLO qqNI (=qbarqbarNI) NNLP = " << NNLO_qqNI_NNLP(z) << endl;
			output2 << "NNLO qqNI (=qbarqbarNI) NNNLP = " << NNLO_qqNI_NNNLP(z) << endl;
			output2 << "NNLO qqbarNI full = " << NNLO_qqbarNI_full(z) << endl;
			output2 << "NNLO qqbarNI NLP = " << NNLO_qqbarNI_NLP(z) << endl;
			output2 << "NNLO qqbarNI NNLP = " << NNLO_qqbarNI_NNLP(z) << endl;
			output2 << "NNLO qqbarNI NNNLP = " << NNLO_qqbarNI_NNNLP(z) << endl;
			
		}
		output2.close();
	}
	
	if(TEST && higgs){
		string q_str = "test_values_" + to_string(Q) +"_as"+to_string(alphas_Q)+"_"+setname+"_higgs.txt";
		ofstream output2;
		output2.open(q_str.c_str()); //.c_str() needed to do a constant string conversion
		cout << "---------- TESTING PARAMETERS -----------" << endl;
		output2 << alphas_Q << endl;
		for(int i =1; i < 10; i++)
		{ 
			double z = (i)*0.1;
			output2 << "=====================================" << endl;
			output2 << "z = " << z << endl;
			output2 << "NLO gg full = " << higgs_NLO_gg_full(z) << endl;
			output2 << "NLO gg LP = " << higgs_NLO_gg_LP(z) << endl;
			output2 << "NLO gg NLP = " << higgs_NLO_gg_NLP(z) << endl;
			output2 << "NLO gg NNLP = " << higgs_NLO_gg_NNLP(z) << endl;
			output2 << "NLO gg NNNLP = " << higgs_NLO_gg_NNNLP(z) << endl;
			output2 << "NLO gg delta = " << higgs_NLO_gg_delta() << endl;
			
		}
		output2.close();
	}
	
	/////////////////////////
	/// integration routines 
	/////////////////////////
	
	//LO
	if(LO){
		if(DY){
			/// LO qqbar
			cout << "computing LO (DY)" << endl;
			LO_qqbar_full = call_vegas(init_vegas_dy("LO"), params);
		}
		else if(higgs){
			/// LO gg
			cout << "computing LO (higgs)" << endl;
			higgs_LO_gg_full = call_vegas(init_vegas_higgs("LO"), params);
		}
	}
	
	//NLO
	if(NLO){
		if(DY){
			/// NLO qqbar
			cout << "computing NLO qqbar (DY)" << endl;
			NLO_qqbar_hard = call_vegas(init_vegas_dy("NLO"), params);
			NLO_qqbar_LP_part1 = call_vegas(init_vegas_dy("NLO","LP"), params);
			NLO_qqbar_LP_cor = call_vegas(init_vegas_dy("NLO","LP_corr"), params);
			NLO_qqbar_LPint = call_vegas(init_vegas_dy("NLO","LP","qqbar",true), params);
			NLO_qqbar_NLP = call_vegas(init_vegas_dy("NLO","NLP"), params);
			NLO_qqbar_NNLP = call_vegas(init_vegas_dy("NLO","NNLP"), params);
			NLO_qqbar_delta = call_vegas(init_vegas_dy("NLO","delta"), params);
			NLO_qqbar_LP.res = NLO_qqbar_LP_part1.res + NLO_qqbar_LP_cor.res;
			NLO_qqbar_LP.err = NLO_qqbar_LP_part1.err + NLO_qqbar_LP_cor.err;
			NLO_qqbar_full.res = NLO_qqbar_LP.res + NLO_qqbar_delta.res + NLO_qqbar_hard.res;
			NLO_qqbar_full.err = NLO_qqbar_LP.err + NLO_qqbar_delta.err + NLO_qqbar_hard.err;
			/// NLO qg
			cout << "computing NLO qg (DY)" << endl;
			res_NLO_qg_full = call_vegas(init_vegas_dy("NLO","full","qg"), params);
			res_NLO_qg_NLP = call_vegas(init_vegas_dy("NLO","NLP","qg"), params);
			res_NLO_qg_NNLP = call_vegas(init_vegas_dy("NLO","NNLP","qg"), params);
			res_NLO_qg_NNNLP = call_vegas(init_vegas_dy("NLO","NNNLP","qg"), params);
		}
		else if(higgs){
			/// NLO gg
			cout << "computing NLO gg (Higgs)" << endl;
			res_higgs_NLO_gg_hard = call_vegas(init_vegas_higgs("NLO"), params);
			res_higgs_NLO_gg_LP_part1 = call_vegas(init_vegas_higgs("NLO","LP"), params);
			res_higgs_NLO_gg_LP_cor = call_vegas(init_vegas_higgs("NLO","LP_corr"), params);
			cout << res_higgs_NLO_gg_LP_cor.res << endl;
			res_higgs_NLO_gg_NLP = call_vegas(init_vegas_higgs("NLO","NLP"), params);
			res_higgs_NLO_gg_NNLP = call_vegas(init_vegas_higgs("NLO","NNLP"), params);
			res_higgs_NLO_gg_NNNLP = call_vegas(init_vegas_higgs("NLO","NNNLP"), params);
			res_higgs_NLO_gg_delta = call_vegas(init_vegas_higgs("NLO","delta"), params);
			res_higgs_NLO_gg_LP.res = res_higgs_NLO_gg_LP_part1.res + res_higgs_NLO_gg_LP_cor.res;
			res_higgs_NLO_gg_LP.err = res_higgs_NLO_gg_LP_part1.err + res_higgs_NLO_gg_LP_cor.err;
			res_higgs_NLO_gg_full.res = res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res + res_higgs_NLO_gg_hard.res;
			res_higgs_NLO_gg_full.err = res_higgs_NLO_gg_LP.err + res_higgs_NLO_gg_delta.err + res_higgs_NLO_gg_hard.err;
			/// NLO qg
			//cout << "computing NLO qg" << endl;
			//res_NLO_qg_full = call_vegas(init_vegas_dy("NLO","full","qg"), params);
			//res_NLO_qg_NLP = call_vegas(init_vegas_dy("NLO","NLP","qg"), params);
			//res_NLO_qg_NNLP = call_vegas(init_vegas_dy("NLO","NNLP","qg"), params);
			//res_NLO_qg_NNNLP = call_vegas(init_vegas_dy("NLO","NNNLP","qg"), params);
		}
	}
	if(NNLO&&DY){
		/// NNLO
		cout << "computing NNLO qqbar" << endl;
		res_NNLO_qqbar_hard = call_vegas(init_vegas_dy("NNLO"), params);
		res_NNLO_qqbar_LP_part1 = call_vegas(init_vegas_dy("NNLO","LP"), params);
		res_NNLO_qqbar_LP_cor = call_vegas(init_vegas_dy("NNLO","LP_corr"), params);
		res_NNLO_qqbar_LPint = call_vegas(init_vegas_dy("NNLO","LP","qqbar",true), params);
		res_NNLO_qqbar_NLP = call_vegas(init_vegas_dy("NNLO","NLP"), params);
		res_NNLO_qqbar_NNLP = call_vegas(init_vegas_dy("NNLO","NNLP"), params);
		res_NNLO_qqbar_NNNLP = call_vegas(init_vegas_dy("NNLO","NNNLP"), params);
		res_NNLO_qqbar_delta = call_vegas(init_vegas_dy("NNLO","delta"), params);
		res_NNLO_qqbar_LP.res = res_NNLO_qqbar_LP_part1.res + res_NNLO_qqbar_LP_cor.res;
		res_NNLO_qqbar_LP.err = res_NNLO_qqbar_LP_part1.err + res_NNLO_qqbar_LP_cor.err;
		res_NNLO_qqbar_full.res = res_NNLO_qqbar_LP.res + res_NNLO_qqbar_delta.res + res_NNLO_qqbar_hard.res;
		res_NNLO_qqbar_full.err = res_NNLO_qqbar_LP.err + res_NNLO_qqbar_delta.err + res_NNLO_qqbar_hard.err;
		
		cout << "computing NNLO qg" << endl;
		res_NNLO_qg_full = call_vegas(init_vegas_dy("NNLO","full","qg"), params);
		res_NNLO_qg_NLP = call_vegas(init_vegas_dy("NNLO","NLP","qg"), params);
		res_NNLO_qg_NNLP = call_vegas(init_vegas_dy("NNLO","NNLP","qg"), params);
		res_NNLO_qg_NNNLP = call_vegas(init_vegas_dy("NNLO","NNNLP","qg"), params);
		cout << "computing NNLO gg" << endl;
		res_NNLO_gg_full = call_vegas(init_vegas_dy("NNLO","full","gg"), params);
		res_NNLO_gg_NLP = call_vegas(init_vegas_dy("NNLO","NLP","gg"), params);
		res_NNLO_gg_NNLP = call_vegas(init_vegas_dy("NNLO","NNLP","gg"), params);
		res_NNLO_gg_NNNLP = call_vegas(init_vegas_dy("NNLO","NNNLP","gg"), params);
		cout << "computing NNLO qq" << endl;
		res_NNLO_qq_full = call_vegas(init_vegas_dy("NNLO","full","qq"), params);
		res_NNLO_qq_NLP = call_vegas(init_vegas_dy("NNLO","NLP","qq"), params);
		res_NNLO_qq_NNLP = call_vegas(init_vegas_dy("NNLO","NNLP","qq"), params);
		res_NNLO_qq_NNNLP = call_vegas(init_vegas_dy("NNLO","NNNLP","qq"), params);
		cout << "computing NNLO qbarqbar" << endl;
		res_NNLO_qbarqbar_full = call_vegas(init_vegas_dy("NNLO","full","qbarqbar"), params);
		res_NNLO_qbarqbar_NLP = call_vegas(init_vegas_dy("NNLO","NLP","qbarqbar"), params);
		res_NNLO_qbarqbar_NNLP = call_vegas(init_vegas_dy("NNLO","NNLP","qbarqbar"), params);
		res_NNLO_qbarqbar_NNNLP = call_vegas(init_vegas_dy("NNLO","NNNLP","qbarqbar"), params);
		cout << "computing NNLO qqNI" << endl;
		res_NNLO_qqNI_full = call_vegas(init_vegas_dy("NNLO","full","qqNI"), params);
		res_NNLO_qqNI_NLP = call_vegas(init_vegas_dy("NNLO","NLP","qqNI"), params);
		res_NNLO_qqNI_NNLP = call_vegas(init_vegas_dy("NNLO","NNLP","qqNI"), params);
		res_NNLO_qqNI_NNNLP = call_vegas(init_vegas_dy("NNLO","NNNLP","qqNI"), params);
		cout << "computing NNLO qqbarNI" << endl;
		res_NNLO_qqbarNI_full = call_vegas(init_vegas_dy("NNLO","full","qqbarNI"), params);
		res_NNLO_qqbarNI_NLP = call_vegas(init_vegas_dy("NNLO","NLP","qqbarNI"), params);
		res_NNLO_qqbarNI_NNLP = call_vegas(init_vegas_dy("NNLO","NNLP","qqbarNI"), params);
		res_NNLO_qqbarNI_NNNLP = call_vegas(init_vegas_dy("NNLO","NNNLP","qqbarNI"), params);
		cout << "computing NNLO qbarqbarNI" << endl;
		res_NNLO_qbarqbarNI_full = call_vegas(init_vegas_dy("NNLO","full","qbarqbarNI"), params);
		res_NNLO_qbarqbarNI_NLP = call_vegas(init_vegas_dy("NNLO","NLP","qbarqbarNI"), params);
		res_NNLO_qbarqbarNI_NNLP = call_vegas(init_vegas_dy("NNLO","NNLP","qbarqbarNI"), params);
		res_NNLO_qbarqbarNI_NNNLP = call_vegas(init_vegas_dy("NNLO","NNNLP","qbarqbarNI"), params);
	}
	
	/////////////////////////////////////
	/// output of results to output file
	/////////////////////////////////////

    output << endl << "INPUT PARAMETERS" << endl ;
    output << "..........................................................................." << endl;
    output << "sqrt[S] = " << S << " GeV" << endl<< "Q = " << Q << " GeV" << endl << "tau = " << tau << endl << "tau*sigma_0 = "<< LO_factor() << endl;
    output << "alphas_Q " << alphas_Q << endl;
	output << "alphas_muF " << alphas_muF << endl;
	output << "alphas_muR " << alphas_muR << endl;
	
    output << "---------------------------------------------------------------------------" << endl;
    
    if(DY){
		output << endl << "RESULTS (DY)" << endl;
    
		output << "..........................................................................." << endl;
		if(LO){
			output << "Total (LO): " << LO_qqbar_full.res << " pb/GeV^2 (" <<  LO_qqbar_full.err << ")" << endl;
		}
		if(NLO){
			output << "Total (NLO): " << NLO_qqbar_full.res+res_NLO_qg_full.res << " pb/GeV^2 (" << NLO_qqbar_full.err+res_NLO_qg_full.err << ")" << endl;
			output << "Total (NLO, qqbar): " << NLO_qqbar_full.res << " pb/GeV^2 (" << NLO_qqbar_full.err << ")" << endl;
			output << "Total (NLO, qg): " << res_NLO_qg_full.res << " pb/GeV^2 (" << res_NLO_qg_full.err << ")" << endl;
		}
		if(NNLO){
			output << "Total (NNLO): " << res_NNLO_qqbar_full.res+res_NNLO_qg_full.res+res_NNLO_gg_full.res+res_NNLO_qq_full.res+res_NNLO_qbarqbar_full.res+res_NNLO_qqNI_full.res+res_NNLO_qqbarNI_full.res+res_NNLO_qbarqbarNI_full.res << " pb/GeV^2 (" << res_NNLO_qqbar_full.err+res_NNLO_qg_full.err+res_NNLO_gg_full.err+res_NNLO_qq_full.err+res_NNLO_qbarqbar_full.res+res_NNLO_qqNI_full.err+res_NNLO_qqbarNI_full.err+res_NNLO_qbarqbarNI_full.err << ")" << endl;
			output << "Total (NNLO, qqbar): " << res_NNLO_qqbar_full.res << " pb/GeV^2 (" << res_NNLO_qqbar_full.err << ")" << endl;
			output << "Total (NNLO, qg): " << res_NNLO_qg_full.res << " pb/GeV^2 (" << res_NNLO_qg_full.err << ")" << endl;
			output << "Total (NNLO, gg): " << res_NNLO_gg_full.res << " pb/GeV^2 (" << res_NNLO_gg_full.err << ")" << endl;
			output << "Total (NNLO, qq): " << res_NNLO_qq_full.res << " pb/GeV^2 (" << res_NNLO_qq_full.err << ")" << endl;
			output << "Total (NNLO, qbarqbar): " << res_NNLO_qbarqbar_full.res << " pb/GeV^2 (" << res_NNLO_qbarqbar_full.err << ")" << endl;
			output << "Total (NNLO, qq non-identical): " << res_NNLO_qqNI_full.res << " pb/GeV^2 (" << res_NNLO_qqNI_full.err << ")" << endl;
			output << "Total (NNLO, qqbar non-identical): " << res_NNLO_qqbarNI_full.res << " pb/GeV^2 (" << res_NNLO_qqbarNI_full.err << ")" << endl;
			output << "Total (NNLO, qbarqbar non-identical): " << res_NNLO_qbarqbarNI_full.res << " pb/GeV^2 (" << res_NNLO_qbarqbarNI_full.err << ")" << endl;
		}
		output << "---------------------------------------------------------------------------" << endl;

		if(NLO){
			output << endl << "POWER EXPANSION (NLO)" << endl;
			output << "==================================qqbar results ==============================" << endl;
			output << "LP         : " << NLO_qqbar_LP.res  << " pb/GeV^2 (" << NLO_qqbar_LP.err  << ")"  << endl;
			output << "LP (int)   : " << NLO_qqbar_LPint.res  << " pb/GeV^2 (" << NLO_qqbar_LPint.err  << ")"  << endl;
			output << "LP (+delta): " << NLO_qqbar_LP.res + NLO_qqbar_delta.res << " pb/GeV^2 (" << NLO_qqbar_LP.err + NLO_qqbar_delta.err << ")"  << "  ; fractional: "<<(NLO_qqbar_LP.res + NLO_qqbar_delta.res)/(NLO_qqbar_full.res) << endl;
			output << "NLP        : " << NLO_qqbar_NLP.res << " pb/GeV^2 (" << NLO_qqbar_NLP.err << ")"  << "  ; fractional: "<<(NLO_qqbar_NLP.res)/(NLO_qqbar_full.res) << endl;
			output << "NNLP       : " << NLO_qqbar_NNLP.res << " pb/GeV^2 (" << NLO_qqbar_NNLP.err << ")"  << " ; fractional: "<<(NLO_qqbar_NNLP.res)/(NLO_qqbar_full.res) << endl;
		   
			output << endl << "CUMULATIVE" << endl;
			output << "..........................................................................." << endl;
			output << "LP : " << NLO_qqbar_LP.res  << " pb/GeV^2 " << " ; fractional: "<<(NLO_qqbar_LP.res)/(NLO_qqbar_full.res)<< endl;
			output << "LP + delta : " << NLO_qqbar_LP.res + NLO_qqbar_delta.res  << " pb/GeV^2 " << " ; fractional: "<<(NLO_qqbar_LP.res + NLO_qqbar_delta.res)/(NLO_qqbar_full.res)<< endl;
			output << "LP + NLP + delta : " << NLO_qqbar_NLP.res+NLO_qqbar_LP.res + NLO_qqbar_delta.res << " pb/GeV^2 " << " ; fractional: "<<(NLO_qqbar_NLP.res+NLO_qqbar_LP.res + NLO_qqbar_delta.res )/(NLO_qqbar_full.res) << endl;
			output << "LP + NLP + NNLP + delta: " << NLO_qqbar_NNLP.res+NLO_qqbar_NLP.res+NLO_qqbar_LP.res + NLO_qqbar_delta.res << " pb/GeV^2 " << " ; fractional: "<<(NLO_qqbar_NNLP.res+NLO_qqbar_NLP.res+NLO_qqbar_LP.res + NLO_qqbar_delta.res )/(NLO_qqbar_full.res) << endl;
			output << "---------------------------------------------------------------------------" << endl;
		
			output << "==================================qg results ==============================" << endl;
			output << "NLP        : " << res_NLO_qg_NLP.res << " pb/GeV^2 (" << res_NLO_qg_NLP.err << ")"  << "  ; fractional: "<<(res_NLO_qg_NLP.res)/(res_NLO_qg_full.res) << endl;
			output << "NNLP       : " << res_NLO_qg_NNLP.res << " pb/GeV^2 (" << res_NLO_qg_NNLP.err << ")"  << " ; fractional: "<<(res_NLO_qg_NNLP.res)/(res_NLO_qg_full.res) << endl;
			output << "NNLP       : " << res_NLO_qg_NNNLP.res << " pb/GeV^2 (" << res_NLO_qg_NNNLP.err << ")"  << " ; fractional: "<<(res_NLO_qg_NNNLP.res)/(res_NLO_qg_full.res) << endl;
		   
			output << endl << "CUMULATIVE" << endl;
			output << "..........................................................................." << endl;
			output << "NLP : " << res_NLO_qg_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NLO_qg_NLP.res)/(res_NLO_qg_full.res) << endl;
			output << "NLP + NNLP: " << res_NLO_qg_NNLP.res+res_NLO_qg_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NLO_qg_NNLP.res+res_NLO_qg_NLP.res)/(res_NLO_qg_full.res) << endl;
			output << "NLP + NNLP + NNNLP: " << res_NLO_qg_NNNLP.res+res_NNLO_qg_NLP.res+res_NLO_qg_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NLO_qg_NNNLP.res+res_NLO_qg_NNLP.res+res_NLO_qg_NLP.res)/(res_NLO_qg_full.res) << endl;
			
		
		}
		if(NNLO){
			output << endl << "POWER EXPANSION (NNLO)" << endl;
			output << "..........................................................................." << endl;
			output << "===============================qqbar results ==============================" << endl;
			output << "LP         : " << res_NNLO_qqbar_LP.res  << " pb/GeV^2 (" << res_NNLO_qqbar_LP.err  << ")"  << endl;
			output << "LP (int)   : " << res_NNLO_qqbar_LPint.res  << " pb/GeV^2 (" << res_NNLO_qqbar_LPint.err  << ")"  << endl;
			output << "LP (+delta): " << res_NNLO_qqbar_LP.res + res_NNLO_qqbar_delta.res << " pb/GeV^2 (" << res_NNLO_qqbar_LP.err + res_NNLO_qqbar_delta.err << ")"  << "  ; fractional: "<<(res_NNLO_qqbar_LP.res + res_NNLO_qqbar_delta.res)/(res_NNLO_qqbar_full.res) << endl;
			output << "NLP        : " << res_NNLO_qqbar_NLP.res << " pb/GeV^2 (" << res_NNLO_qqbar_NLP.err << ")"  << "  ; fractional: "<<(res_NNLO_qqbar_NLP.res)/(res_NNLO_qqbar_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qqbar_NNLP.res << " pb/GeV^2 (" << res_NNLO_qqbar_NNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qqbar_NNLP.res)/(res_NNLO_qqbar_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qqbar_NNNLP.res << " pb/GeV^2 (" << res_NNLO_qqbar_NNNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qqbar_NNNLP.res)/(res_NNLO_qqbar_full.res) << endl;
		   
			output << endl << "CUMULATIVE" << endl;
			output << "..........................................................................." << endl;
			output << "LP : " << res_NNLO_qqbar_LP.res  << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qqbar_LP.res)/(res_NNLO_qqbar_full.res)<< endl;
			output << "LP + delta : " << res_NNLO_qqbar_LP.res + res_NNLO_qqbar_delta.res  << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qqbar_LP.res + res_NNLO_qqbar_delta.res)/(res_NNLO_qqbar_full.res)<< endl;
			output << "LP + NLP + delta : " << res_NNLO_qqbar_NLP.res+res_NNLO_qqbar_LP.res + res_NNLO_qqbar_delta.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qqbar_NLP.res+res_NNLO_qqbar_LP.res + res_NNLO_qqbar_delta.res )/(res_NNLO_qqbar_full.res) << endl;
			output << "LP + NLP + NNLP + delta: " << res_NNLO_qqbar_NNLP.res+res_NNLO_qqbar_NLP.res+res_NNLO_qqbar_LP.res + res_NNLO_qqbar_delta.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qqbar_NNLP.res+res_NNLO_qqbar_NLP.res+res_NNLO_qqbar_LP.res + res_NNLO_qqbar_delta.res )/(res_NNLO_qqbar_full.res) << endl;
			output << "LP + NLP + NNLP + NNNLP + delta: " << res_NNLO_qqbar_NNNLP.res+res_NNLO_qqbar_NNLP.res+res_NNLO_qqbar_NLP.res+res_NNLO_qqbar_LP.res + res_NNLO_qqbar_delta.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qqbar_NNNLP.res+res_NNLO_qqbar_NNLP.res+res_NNLO_qqbar_NLP.res+res_NNLO_qqbar_LP.res + res_NNLO_qqbar_delta.res )/(res_NNLO_qqbar_full.res) << endl;
			output << "..........................................................................." << endl;
			
			output << "==================================qg results ==============================" << endl;
			output << "NLP        : " << res_NNLO_qg_NLP.res << " pb/GeV^2 (" << res_NNLO_qg_NLP.err << ")"  << "  ; fractional: "<<(res_NNLO_qg_NLP.res)/(res_NNLO_qg_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qg_NNLP.res << " pb/GeV^2 (" << res_NNLO_qg_NNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qg_NNLP.res)/(res_NNLO_qg_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qg_NNNLP.res << " pb/GeV^2 (" << res_NNLO_qg_NNNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qg_NNNLP.res)/(res_NNLO_qg_full.res) << endl;
		   
			output << endl << "CUMULATIVE" << endl;
			output << "..........................................................................." << endl;
			output << "NLP : " << res_NNLO_qg_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qg_NLP.res)/(res_NNLO_qg_full.res) << endl;
			output << "NLP + NNLP: " << res_NNLO_qg_NNLP.res+res_NNLO_qg_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qg_NNLP.res+res_NNLO_qg_NLP.res)/(res_NNLO_qg_full.res) << endl;
			output << "NLP + NNLP + NNNLP: " << res_NNLO_qg_NNNLP.res+res_NNLO_qg_NNLP.res+res_NNLO_qg_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qg_NNNLP.res+res_NNLO_qg_NNLP.res+res_NNLO_qg_NLP.res)/(res_NNLO_qg_full.res) << endl;
				
			output << "==================================gg results ==============================" << endl;
			output << "NLP        : " << res_NNLO_gg_NLP.res << " pb/GeV^2 (" << res_NNLO_gg_NLP.err << ")"  << "  ; fractional: "<<(res_NNLO_gg_NLP.res)/(res_NNLO_gg_full.res) << endl;
			output << "NNLP       : " << res_NNLO_gg_NNLP.res << " pb/GeV^2 (" << res_NNLO_gg_NNLP.err << ")"  << " ; fractional: "<<(res_NNLO_gg_NNLP.res)/(res_NNLO_gg_full.res) << endl;
			output << "NNLP       : " << res_NNLO_gg_NNNLP.res << " pb/GeV^2 (" << res_NNLO_gg_NNNLP.err << ")"  << " ; fractional: "<<(res_NNLO_gg_NNNLP.res)/(res_NNLO_gg_full.res) << endl;
		   
			output << endl << "CUMULATIVE" << endl;
			output << "..........................................................................." << endl;
			output << "NLP : " << res_NNLO_gg_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_gg_NLP.res)/(res_NNLO_gg_full.res) << endl;
			output << "NLP + NNLP: " << res_NNLO_gg_NNLP.res+res_NNLO_gg_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_gg_NNLP.res+res_NNLO_gg_NLP.res)/(res_NNLO_gg_full.res) << endl;
			output << "NLP + NNLP + NNNLP: " << res_NNLO_gg_NNNLP.res+res_NNLO_gg_NNLP.res+res_NNLO_gg_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_gg_NNNLP.res+res_NNLO_gg_NNLP.res+res_NNLO_gg_NLP.res)/(res_NNLO_gg_full.res) << endl;
				
			output << "==================================qq results ==============================" << endl;
			output << "NLP        : " << res_NNLO_qq_NLP.res << " pb/GeV^2 (" << res_NNLO_qq_NLP.err << ")"  << "  ; fractional: "<<(res_NNLO_qq_NLP.res)/(res_NNLO_qq_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qq_NNLP.res << " pb/GeV^2 (" << res_NNLO_qq_NNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qq_NNLP.res)/(res_NNLO_qq_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qq_NNNLP.res << " pb/GeV^2 (" << res_NNLO_qq_NNNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qq_NNNLP.res)/(res_NNLO_qq_full.res) << endl;
		   
			output << endl << "CUMULATIVE" << endl;
			output << "..........................................................................." << endl;
			output << "NLP : " << res_NNLO_qq_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qq_NLP.res)/(res_NNLO_qq_full.res) << endl;
			output << "NLP + NNLP: " << res_NNLO_qq_NNLP.res+res_NNLO_qq_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qq_NNLP.res+res_NNLO_qq_NLP.res)/(res_NNLO_qq_full.res) << endl;
			output << "NLP + NNLP + NNNLP: " << res_NNLO_qq_NNNLP.res+res_NNLO_qq_NNLP.res+res_NNLO_qq_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qq_NNNLP.res+res_NNLO_qq_NNLP.res+res_NNLO_qq_NLP.res)/(res_NNLO_qq_full.res) << endl;
				
			output << "==================================qbarqbar results ==============================" << endl;
			output << "NLP        : " << res_NNLO_qbarqbar_NLP.res << " pb/GeV^2 (" << res_NNLO_qbarqbar_NLP.err << ")"  << "  ; fractional: "<<(res_NNLO_qbarqbar_NLP.res)/(res_NNLO_qbarqbar_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qbarqbar_NNLP.res << " pb/GeV^2 (" << res_NNLO_qbarqbar_NNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qbarqbar_NNLP.res)/(res_NNLO_qbarqbar_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qbarqbar_NNNLP.res << " pb/GeV^2 (" << res_NNLO_qbarqbar_NNNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qbarqbar_NNNLP.res)/(res_NNLO_qbarqbar_full.res) << endl;
		   
			output << endl << "CUMULATIVE" << endl;
			output << "..........................................................................." << endl;
			output << "NLP : " << res_NNLO_qbarqbar_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qbarqbar_NLP.res)/(res_NNLO_qbarqbar_full.res) << endl;
			output << "NLP + NNLP: " << res_NNLO_qbarqbar_NNLP.res+res_NNLO_qbarqbar_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qbarqbar_NNLP.res+res_NNLO_qbarqbar_NLP.res)/(res_NNLO_qbarqbar_full.res) << endl;
			output << "NLP + NNLP + NNNLP: " << res_NNLO_qbarqbar_NNNLP.res+res_NNLO_qbarqbar_NNLP.res+res_NNLO_qbarqbar_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qbarqbar_NNNLP.res+res_NNLO_qbarqbar_NNLP.res+res_NNLO_qbarqbar_NLP.res)/(res_NNLO_qbarqbar_full.res) << endl;
				
			output << "========================qq (non-identical) results ==============================" << endl;
			output << "NLP        : " << res_NNLO_qqNI_NLP.res << " pb/GeV^2 (" << res_NNLO_qqNI_NLP.err << ")"  << "  ; fractional: "<<(res_NNLO_qqNI_NLP.res)/(res_NNLO_qqNI_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qqNI_NNLP.res << " pb/GeV^2 (" << res_NNLO_qqNI_NNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qqNI_NNLP.res)/(res_NNLO_qqNI_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qqNI_NNNLP.res << " pb/GeV^2 (" << res_NNLO_qqNI_NNNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qqNI_NNNLP.res)/(res_NNLO_qqNI_full.res) << endl;
		   
			output << endl << "CUMULATIVE" << endl;
			output << "..........................................................................." << endl;
			output << "NLP : " << res_NNLO_qqNI_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qqNI_NLP.res)/(res_NNLO_qqNI_full.res) << endl;
			output << "NLP + NNLP: " << res_NNLO_qqNI_NNLP.res+res_NNLO_qqNI_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qqNI_NNLP.res+res_NNLO_qqNI_NLP.res)/(res_NNLO_qqNI_full.res) << endl;
			output << "NLP + NNLP + NNNLP: " << res_NNLO_qqNI_NNNLP.res+res_NNLO_qqNI_NNLP.res+res_NNLO_qqNI_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qqNI_NNNLP.res+res_NNLO_qqNI_NNLP.res+res_NNLO_qqNI_NLP.res)/(res_NNLO_qqNI_full.res) << endl;
				
			output << "======================qqbar (non-identical) results ==============================" << endl;
			output << "NLP        : " << res_NNLO_qqbarNI_NLP.res << " pb/GeV^2 (" << res_NNLO_qqbarNI_NLP.err << ")"  << "  ; fractional: "<<(res_NNLO_qqbarNI_NLP.res)/(res_NNLO_qqbarNI_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qqbarNI_NNLP.res << " pb/GeV^2 (" << res_NNLO_qqbarNI_NNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qqbarNI_NNLP.res)/(res_NNLO_qqbarNI_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qqbarNI_NNNLP.res << " pb/GeV^2 (" << res_NNLO_qqbarNI_NNNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qqbarNI_NNNLP.res)/(res_NNLO_qqbarNI_full.res) << endl;
		   
			output << endl << "CUMULATIVE" << endl;
			output << "..........................................................................." << endl;
			output << "NLP : " << res_NNLO_qqbarNI_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qqbarNI_NLP.res)/(res_NNLO_qqbarNI_full.res) << endl;
			output << "NLP + NNLP: " << res_NNLO_qqbarNI_NNLP.res+res_NNLO_qqbarNI_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qqbarNI_NNLP.res+res_NNLO_qqbarNI_NLP.res)/(res_NNLO_qqbarNI_full.res) << endl;
			output << "NLP + NNLP + NNNLP: " << res_NNLO_qqbarNI_NNNLP.res+res_NNLO_qqbarNI_NNLP.res+res_NNLO_qqbarNI_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qqbarNI_NNNLP.res+res_NNLO_qqbarNI_NNLP.res+res_NNLO_qqbarNI_NLP.res)/(res_NNLO_qqbarNI_full.res) << endl;
			
				
			output << "====================qbarqbar (non-identical) results ============================" << endl;
			output << "NLP        : " << res_NNLO_qbarqbarNI_NLP.res << " pb/GeV^2 (" << res_NNLO_qbarqbarNI_NLP.err << ")"  << "  ; fractional: "<<(res_NNLO_qbarqbarNI_NLP.res)/(res_NNLO_qbarqbarNI_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qbarqbarNI_NNLP.res << " pb/GeV^2 (" << res_NNLO_qbarqbarNI_NNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qbarqbarNI_NNLP.res)/(res_NNLO_qbarqbarNI_full.res) << endl;
			output << "NNLP       : " << res_NNLO_qbarqbarNI_NNNLP.res << " pb/GeV^2 (" << res_NNLO_qbarqbarNI_NNNLP.err << ")"  << " ; fractional: "<<(res_NNLO_qbarqbarNI_NNNLP.res)/(res_NNLO_qbarqbarNI_full.res) << endl;
		   
			output << endl << "CUMULATIVE" << endl;
			output << "..........................................................................." << endl;
			output << "NLP : " << res_NNLO_qbarqbarNI_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qbarqbarNI_NLP.res)/(res_NNLO_qbarqbarNI_full.res) << endl;
			output << "NLP + NNLP: " << res_NNLO_qbarqbarNI_NNLP.res+res_NNLO_qbarqbarNI_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qbarqbarNI_NNLP.res+res_NNLO_qbarqbarNI_NLP.res)/(res_NNLO_qbarqbarNI_full.res) << endl;
			output << "NLP + NNLP + NNNLP: " << res_NNLO_qbarqbarNI_NNNLP.res+res_NNLO_qbarqbarNI_NNLP.res+res_NNLO_qbarqbarNI_NLP.res << " pb/GeV^2 " << " ; fractional: "<<(res_NNLO_qbarqbarNI_NNNLP.res+res_NNLO_qbarqbarNI_NNLP.res+res_NNLO_qbarqbarNI_NLP.res)/(res_NNLO_qbarqbarNI_full.res) << endl;
			
			output << "---------------------------------------------------------------------------" << endl;
		}
	}
    if(higgs){
		output << endl << "RESULTS (higgs)" << endl;
    
		output << "..........................................................................." << endl;
		if(LO){
			output << "Total (LO): " << higgs_LO_gg_full.res << " pb (" <<  higgs_LO_gg_full.err << ")" << endl;
		}
		if(NLO){
			output << "Total (NLO, gg): " << res_higgs_NLO_gg_full.res << " pb (" << res_higgs_NLO_gg_full.err << ")" << endl;
			output << "Total (NLO, gg delta): " << higgs_LO_gg_full.res+res_higgs_NLO_gg_delta.res << " pb (" << res_higgs_NLO_gg_delta.err+higgs_LO_gg_full.err << ")" << endl;
			output << "Total (NLO, gg z): " << res_higgs_NLO_gg_full.res-res_higgs_NLO_gg_delta.res << " pb (" << -res_higgs_NLO_gg_delta.err+res_higgs_NLO_gg_full.err << ")" << endl;
			output << "Total (NLO, gg z): " << res_higgs_NLO_gg_full.res-res_higgs_NLO_gg_delta.res-res_higgs_NLO_gg_LP_cor.res << " pb (" << -res_higgs_NLO_gg_delta.err+res_higgs_NLO_gg_full.err << ")" << endl;
		}
		output << "---------------------------------------------------------------------------" << endl;

		if(NLO){
			output << endl << "POWER EXPANSION (NLO)" << endl;
			output << "==================================gg results ==============================" << endl;
			output << "LP         : " << res_higgs_NLO_gg_LP.res  << " pb (" << res_higgs_NLO_gg_LP.err  << ")"  << endl;
			output << "delta      : " << res_higgs_NLO_gg_delta.res << " pb (" << res_higgs_NLO_gg_delta.err << ")"  << "  ; fractional: "<<(res_higgs_NLO_gg_delta.res)/(res_higgs_NLO_gg_full.res) << endl;
			output << "LP +delta  : " << res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res << " pb (" << res_higgs_NLO_gg_LP.err + res_higgs_NLO_gg_delta.err << ")"  << "  ; fractional: "<<(res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res)/(res_higgs_NLO_gg_full.res) << endl;
			output << "NLP        : " << res_higgs_NLO_gg_NLP.res << " pb (" << res_higgs_NLO_gg_NLP.err << ")"  << "  ; fractional: "<<(res_higgs_NLO_gg_NLP.res)/(res_higgs_NLO_gg_full.res) << endl;
			output << "NNLP       : " << res_higgs_NLO_gg_NNLP.res << " pb (" << res_higgs_NLO_gg_NNLP.err << ")"  << " ; fractional: "<<(res_higgs_NLO_gg_NNLP.res)/(res_higgs_NLO_gg_full.res) << endl;
		    output << "NNNLP       : " << res_higgs_NLO_gg_NNNLP.res << " pb (" << res_higgs_NLO_gg_NNNLP.err << ")"  << " ; fractional: "<<(res_higgs_NLO_gg_NNNLP.res)/(res_higgs_NLO_gg_full.res) << endl;
		   
			output << endl << "CUMULATIVE" << endl;
			output << "..........................................................................." << endl;
			output << "LP : " << res_higgs_NLO_gg_LP.res  << " pb " << " ; fractional: "<<(res_higgs_NLO_gg_LP.res)/(res_higgs_NLO_gg_full.res)<< endl;
			output << "LP + delta : " << res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res  << " pb " << " ; fractional: "<<(res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res)/(res_higgs_NLO_gg_full.res)<< endl;
			output << "LP + NLP + delta : " << res_higgs_NLO_gg_NLP.res+res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res << " pb " << " ; fractional: "<<(res_higgs_NLO_gg_NLP.res+res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res )/(res_higgs_NLO_gg_full.res) << endl;
			output << "LP + NLP + NNLP + delta: " << res_higgs_NLO_gg_NNLP.res+res_higgs_NLO_gg_NLP.res+res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res << " pb " << " ; fractional: "<<(res_higgs_NLO_gg_NNLP.res+res_higgs_NLO_gg_NLP.res+res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res )/(res_higgs_NLO_gg_full.res) << endl;
			output << "LP + NLP + NNLP + NNNLP + delta: " << res_higgs_NLO_gg_NNNLP.res+res_higgs_NLO_gg_NNLP.res+res_higgs_NLO_gg_NLP.res+res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res << " pb " << " ; fractional: "<<(res_higgs_NLO_gg_NNNLP.res+res_higgs_NLO_gg_NNLP.res+res_higgs_NLO_gg_NLP.res+res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res )/(res_higgs_NLO_gg_full.res) << endl;
			output << "---------------------------------------------------------------------------" << endl;
			
			output << endl << "COMPARISON" << endl;
			output << "..........................................................................." << endl;
			output << "power 0 : " << res_higgs_NLO_gg_LP.res + higgs_LO_gg_full.res+res_higgs_NLO_gg_delta.res << " pb " << " ; fractional: "<<(res_higgs_NLO_gg_LP.res)/(res_higgs_NLO_gg_full.res)<< endl;
			output << "power 1 : " << res_higgs_NLO_gg_LP.res + higgs_LO_gg_full.res+res_higgs_NLO_gg_delta.res  << " pb " << " ; fractional: "<<(res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res)/(res_higgs_NLO_gg_full.res)<< endl;
			output << "power 2 : " << res_higgs_NLO_gg_NLP.res+res_higgs_NLO_gg_LP.res + higgs_LO_gg_full.res+res_higgs_NLO_gg_delta.res << " pb " << " ; fractional: "<<(res_higgs_NLO_gg_NLP.res+res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res )/(res_higgs_NLO_gg_full.res) << endl;
			output << "power 3 : " << res_higgs_NLO_gg_NNLP.res+res_higgs_NLO_gg_NLP.res+res_higgs_NLO_gg_LP.res + higgs_LO_gg_full.res+res_higgs_NLO_gg_delta.res << " pb " << " ; fractional: "<<(res_higgs_NLO_gg_NNLP.res+res_higgs_NLO_gg_NLP.res+res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res )/(res_higgs_NLO_gg_full.res) << endl;
			output << "power 4 : " << res_higgs_NLO_gg_NNNLP.res+res_higgs_NLO_gg_NNLP.res+res_higgs_NLO_gg_NLP.res+res_higgs_NLO_gg_LP.res + higgs_LO_gg_full.res+res_higgs_NLO_gg_delta.res << " pb " << " ; fractional: "<<(res_higgs_NLO_gg_NNNLP.res+res_higgs_NLO_gg_NNLP.res+res_higgs_NLO_gg_NLP.res+res_higgs_NLO_gg_LP.res + res_higgs_NLO_gg_delta.res )/(res_higgs_NLO_gg_full.res) << endl;
			output << "---------------------------------------------------------------------------" << endl;
	
		}
	}
			
	output.close();
	return 0;
}
