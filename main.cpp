#include <cmath>
#include <iostream>
#include <fstream>
#include "cuba.h"
#include "deriv_pdf.h"
#include "monte_carlo.h"
#include "k_factors_dy.h"
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
	//////////////////////////////////////////////
	
	pdf_sum_qg_charge_weighted(0.5,0.25);
	double z(tau+0.1);
	double eps(0.00001);
	ofstream output;
	string q_str = "output_Q" + to_string(Q) +"_as"+to_string(alphas_Q)+".txt";
	output.open(q_str.c_str()); //.c_str() needed to do a constant string conversion
	output << "z dphi dphi*LP dphi*NLP dphi*NNLP dphi*delta dphi*exact " << endl;
	lumni_params params = {z, pdfs};
	params.z = 0.5;

	/// LO
	cout << "computing LO" << endl;
	results LO_qqbar_full = call_vegas(init_vegas("LO"), params);
	
	/// NLO
	cout << "computing NLO qqbar" << endl;
	results NLO_qqbar_hard = call_vegas(init_vegas("NLO"), params);
	results NLO_qqbar_LP_part1 = call_vegas(init_vegas("NLO","LP"), params);
	results NLO_qqbar_LP_cor = call_vegas(init_vegas("NLO","LP_corr"), params);
	results NLO_qqbar_NLP = call_vegas(init_vegas("NLO","NLP"), params);
	results NLO_qqbar_NNLP = call_vegas(init_vegas("NLO","NNLP"), params);
	results NLO_qqbar_delta = call_vegas(init_vegas("NLO","delta"), params);
	results NLO_qqbar_LP;
	NLO_qqbar_LP.res = NLO_qqbar_LP_part1.res + NLO_qqbar_LP_cor.res;
	NLO_qqbar_LP.err = NLO_qqbar_LP_part1.err + NLO_qqbar_LP_cor.err;
	results NLO_qqbar_full;
	NLO_qqbar_full.res = NLO_qqbar_LP.res + NLO_qqbar_delta.res + NLO_qqbar_hard.res;
	NLO_qqbar_full.err = NLO_qqbar_LP.err + NLO_qqbar_delta.err + NLO_qqbar_hard.err;
	
	
	/// NLO
	cout << "computing NLO qg" << endl;
	results NLO_qg_full = call_vegas(init_vegas("NLO","full","qg"), params);
	
	/// NNLO
	cout << "computing NNLO qqbar" << endl;
	results NNLO_qqbar_hard = call_vegas(init_vegas("NNLO"), params);
	results NNLO_qqbar_LP_part1 = call_vegas(init_vegas("NNLO","LP"), params);
	results NNLO_qqbar_LP_cor = call_vegas(init_vegas("NNLO","LP_corr"), params);
	results NNLO_qqbar_NLP = call_vegas(init_vegas("NNLO","NLP"), params);
	results NNLO_qqbar_NNLP = call_vegas(init_vegas("NNLO","NNLP"), params);
	results NNLO_qqbar_NNNLP = call_vegas(init_vegas("NNLO","NNNLP"), params);
	results NNLO_qqbar_delta = call_vegas(init_vegas("NNLO","delta"), params);
	results NNLO_qqbar_LP;
	NNLO_qqbar_LP.res = NNLO_qqbar_LP_part1.res + NNLO_qqbar_LP_cor.res;
	NNLO_qqbar_LP.err = NNLO_qqbar_LP_part1.err + NNLO_qqbar_LP_cor.err;
	results NNLO_qqbar_full;
	NNLO_qqbar_full.res = NNLO_qqbar_LP.res + NNLO_qqbar_delta.res + NNLO_qqbar_hard.res;
	NNLO_qqbar_full.err = NNLO_qqbar_LP.err + NNLO_qqbar_delta.err + NNLO_qqbar_hard.err;
	
	
	cout << "===========================================================================" << endl;
    cout << endl << "dsigma/dQ^2 for DY" << endl <<endl;
    cout << "===========================================================================" << endl;
    
    cout << endl << "INPUT PARAMETERS" << endl ;
    cout << "..........................................................................." << endl;
    cout << "sqrt[S] = " << S << " GeV" << endl<< "Q = " << Q << " GeV" << endl << "tau = " << tau << endl << "tau*sigma_0 = "<< LO_factor() << endl;
    cout << "alphas_Q " << alphas_Q << endl;
	cout << "alphas_muF " << alphas_muF << endl;
	cout << "alphas_muR " << alphas_muR << endl;
	
    cout << "---------------------------------------------------------------------------" << endl;
    
    cout << endl << "RESULTS" << endl;
    
    cout << "..........................................................................." << endl;
    cout << "Total direct integration (LO): " << LO_qqbar_full.res << " pb/GeV^2 (" <<  LO_qqbar_full.err << ")" << endl;
    cout << "Total direct integration (NLO, qqbar): " << NLO_qqbar_full.res << " pb/GeV^2 (" << NLO_qqbar_full.err << ")" << endl;
    cout << "Total direct integration (NLO, qg): " << NLO_qg_full.res << " pb/GeV^2 (" << NLO_qg_full.err << ")" << endl;
    cout << "Total direct integration (NNLO): " << NNLO_qqbar_full.res << " pb/GeV^2 (" << NNLO_qqbar_full.err << ")" << endl;
    cout << "---------------------------------------------------------------------------" << endl;

    
    cout << endl << "POWER EXPANSION (NLO)" << endl;
    cout << "..........................................................................." << endl;
    cout << "LP         : " << NLO_qqbar_LP.res  << " pb/GeV^2 (" << NLO_qqbar_LP.err  << ")"  << endl;
    cout << "LP (+delta): " << NLO_qqbar_LP.res + NLO_qqbar_delta.res << " pb/GeV^2 (" << NLO_qqbar_LP.err + NLO_qqbar_delta.err << ")"  << "  ; fractional: "<<(NLO_qqbar_LP.res + NLO_qqbar_delta.res)/(NLO_qqbar_full.res) << endl;
    cout << "NLP        : " << NLO_qqbar_NLP.res << " pb/GeV^2 (" << NLO_qqbar_NLP.err << ")"  << "  ; fractional: "<<(NLO_qqbar_NLP.res)/(NLO_qqbar_full.res) << endl;
    cout << "NNLP       : " << NLO_qqbar_NNLP.res << " pb/GeV^2 (" << NLO_qqbar_NNLP.err << ")"  << " ; fractional: "<<(NLO_qqbar_NNLP.res)/(NLO_qqbar_full.res) << endl;
   
	cout << endl << "CUMULATIVE" << endl;
    cout << "..........................................................................." << endl;
    cout << "LP : " << NLO_qqbar_LP.res  << " pb/GeV^2 " << " ; fractional: "<<(NLO_qqbar_LP.res)/(NLO_qqbar_full.res)<< endl;
    cout << "LP + delta : " << NLO_qqbar_LP.res + NLO_qqbar_delta.res  << " pb/GeV^2 " << " ; fractional: "<<(NLO_qqbar_LP.res + NLO_qqbar_delta.res)/(NLO_qqbar_full.res)<< endl;
    cout << "LP + NLP + delta : " << NLO_qqbar_NLP.res+NLO_qqbar_LP.res + NLO_qqbar_delta.res << " pb/GeV^2 " << " ; fractional: "<<(NLO_qqbar_NLP.res+NLO_qqbar_LP.res + NLO_qqbar_delta.res )/(NLO_qqbar_full.res) << endl;
    cout << "LP + NLP + NNLP + delta: " << NLO_qqbar_NNLP.res+NLO_qqbar_NLP.res+NLO_qqbar_LP.res + NLO_qqbar_delta.res << " pb/GeV^2 " << " ; fractional: "<<(NLO_qqbar_NNLP.res+NLO_qqbar_NLP.res+NLO_qqbar_LP.res + NLO_qqbar_delta.res )/(NLO_qqbar_full.res) << endl;
    cout << "---------------------------------------------------------------------------" << endl;
    
	cout << endl << "POWER EXPANSION (NNLO)" << endl;
    cout << "..........................................................................." << endl;
    cout << "LP         : " << NNLO_qqbar_LP.res  << " pb/GeV^2 (" << NNLO_qqbar_LP.err  << ")"  << endl;
    cout << "LP (+delta): " << NNLO_qqbar_LP.res + NNLO_qqbar_delta.res << " pb/GeV^2 (" << NNLO_qqbar_LP.err + NNLO_qqbar_delta.err << ")"  << "  ; fractional: "<<(NNLO_qqbar_LP.res + NNLO_qqbar_delta.res)/(NNLO_qqbar_full.res) << endl;
    cout << "NLP        : " << NNLO_qqbar_NLP.res << " pb/GeV^2 (" << NNLO_qqbar_NLP.err << ")"  << "  ; fractional: "<<(NNLO_qqbar_NLP.res)/(NNLO_qqbar_full.res) << endl;
    cout << "NNLP       : " << NNLO_qqbar_NNLP.res << " pb/GeV^2 (" << NNLO_qqbar_NNLP.err << ")"  << " ; fractional: "<<(NNLO_qqbar_NNLP.res)/(NNLO_qqbar_full.res) << endl;
    cout << "NNLP       : " << NNLO_qqbar_NNNLP.res << " pb/GeV^2 (" << NNLO_qqbar_NNNLP.err << ")"  << " ; fractional: "<<(NNLO_qqbar_NNNLP.res)/(NNLO_qqbar_full.res) << endl;
   
	cout << endl << "CUMULATIVE" << endl;
    cout << "..........................................................................." << endl;
    cout << "LP : " << NNLO_qqbar_LP.res  << " pb/GeV^2 " << " ; fractional: "<<(NNLO_qqbar_LP.res)/(NNLO_qqbar_full.res)<< endl;
    cout << "LP + delta : " << NNLO_qqbar_LP.res + NNLO_qqbar_delta.res  << " pb/GeV^2 " << " ; fractional: "<<(NNLO_qqbar_LP.res + NNLO_qqbar_delta.res)/(NNLO_qqbar_full.res)<< endl;
    cout << "LP + NLP + delta : " << NNLO_qqbar_NLP.res+NNLO_qqbar_LP.res + NNLO_qqbar_delta.res << " pb/GeV^2 " << " ; fractional: "<<(NNLO_qqbar_NLP.res+NNLO_qqbar_LP.res + NNLO_qqbar_delta.res )/(NNLO_qqbar_full.res) << endl;
    cout << "LP + NLP + NNLP + delta: " << NNLO_qqbar_NNLP.res+NNLO_qqbar_NLP.res+NNLO_qqbar_LP.res + NNLO_qqbar_delta.res << " pb/GeV^2 " << " ; fractional: "<<(NNLO_qqbar_NNLP.res+NNLO_qqbar_NLP.res+NNLO_qqbar_LP.res + NNLO_qqbar_delta.res )/(NNLO_qqbar_full.res) << endl;
    cout << "LP + NLP + NNLP + NNNLP + delta: " << NNLO_qqbar_NNNLP.res+NNLO_qqbar_NNLP.res+NNLO_qqbar_NLP.res+NNLO_qqbar_LP.res + NNLO_qqbar_delta.res << " pb/GeV^2 " << " ; fractional: "<<(NNLO_qqbar_NNNLP.res+NNLO_qqbar_NNLP.res+NNLO_qqbar_NLP.res+NNLO_qqbar_LP.res + NNLO_qqbar_delta.res )/(NNLO_qqbar_full.res) << endl;
    cout << "---------------------------------------------------------------------------" << endl;
    
    

   
	
	
	/*
	for(int i = 0; i < 100; i++){
		cout << i << endl;
		z = tau + (i+1.)*increment;
		double zpluseps(z+eps), zmineps(z-eps);
		params.z = zpluseps;
		xl[0] = tau/params.z;
		results resultplus = call_vegas("lumni", xl, xu, params);
		//cout << resultplus.res << "+/-" <<  resultplus.err  << endl;
		params.z = zmineps;
		xl[0] = tau/params.z;
		results resultmin = call_vegas("lumni", xl, xu, params);
		//cout << resultmin.res << "+/-" <<  resultmin.err  << endl;
		double dphi = (1./zpluseps*resultplus.res-1./zmineps*resultmin.res)/(2.*eps);
		double delta = dphi*k_zint_NLO_gg_dy_delta();
		double LP = dphi*k_zint_NLO_gg_dy_LP(z);
		double NLP = dphi*k_zint_NLO_gg_dy_NLP(z);
		double NNLP = dphi*k_zint_NLO_gg_dy_NNLP(z);
		double exact = dphi*k_zint_NLO_gg_dy_exact(z);
		cout << z << " " << 1./zpluseps*resultplus.res << " " << 1./zmineps*resultmin.res <<  " "<< dphi << " " << k_zint_NLO_gg_dy_LP(z) << " " << k_zint_NLO_gg_dy_exact(z) << endl;
		output << z << " " << resultplus.res << " " << resultmin.res << " " << LP+delta << " " << LP+NLP+delta << " " << LP+NLP+NNLP+delta << " " << exact+delta << " " << (LP+delta)/(exact+delta) << " " <<  (LP+NLP+delta)/(exact+delta) << " " <<  (LP+NLP+NNLP+delta)/(exact+delta) << " "<< endl;
	}
	output.close();
  */
	return 0;
}
