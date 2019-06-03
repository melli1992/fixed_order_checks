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

	/*double xl[1] = {0.};
	double xu[1] = {1.};
	results testsig = call_vegas("test1", xl, xu, params);
	
	cout << "Results test1 " << testsig.res << " " << testsig.err << endl;
	params.z = 0.5+0.0001;
	results testsigint1 = call_vegas("test2", xl, xu, params);
	params.z = 0.5-0.0001;
	results testsigint2 = call_vegas("test2", xl, xu, params);
	
	cout << "Results test2 " << (testsigint1.res-testsigint2.res)/(2*0.0001) << " " << testsigint1.err << endl;
	double result = pow(tau-params.z,2.)*(tau+params.z)/(2.*pow(params.z,2.));
	cout <<  result << endl;
	double increment = (1.-tau)/100.;
	*/
	
	
	double xl1[2] = {tau,0.};
	double xu1[2] = {1.,1.};
	results testsig1 = call_vegas("DY_NNLO_LP", xl1, xu1, params);
	xl1[0] = 0;
	xu1[0] = tau;
	results testsig2 = call_vegas("DY_NNLO_LP_corr", xl1, xu1, params);
	double DtermsresNNLO = testsig1.res + testsig2.res;
	double DtermserNNLO = testsig1.err + testsig2.err;
	xl1[0] = tau;
	xu1[0] = 1;
	
	
	
	double xl[1] = {tau};
	double xu[1] = {1.};
	testsig1 = call_vegas("DY_LO", xl, xu, params);
	double resLO = testsig1.res;
	double erLO = testsig1.err;
	cout << resLO << endl;
	cout << erLO << endl;
	testsig1 = call_vegas("sum_PDF", xl, xu, params);
	double resPDF = testsig1.res;
	double erPDF = testsig1.err;
	cout << resPDF << " " << erPDF << endl;
	testsig1 = call_vegas("DY_NNLO_NLP", xl1, xu1, params);
	double NNLO_NLP = testsig1.res;
	double NNLO_NLPer = testsig1.err;
	testsig1 = call_vegas("DY_NNLO_NNLP", xl1, xu1, params);
	double NNLO_NNLP = testsig1.res;
	double NNLO_NNLPer = testsig1.err;
	testsig1 = call_vegas("DY_NNLO_NNNLP", xl1, xu1, params);
	double NNLO_NNNLP = testsig1.res;
	double NNLO_NNNLPer = testsig1.err;
	testsig1 = call_vegas("DY_NNLO_full", xl1, xu1, params);
	double NNLO_Ltermsres = testsig1.res;
	double NNLO_Ltermser = testsig1.err;
	testsig1 = call_vegas("DY_NNLO_delta", xl, xu, params);
	double NNLO_delta = testsig1.res;
	double NNLO_deltaer = testsig1.err;
	
	
	testsig1 = call_vegas("DY_NLO_LP1", xl1, xu1, params);
	xl1[0] = 0;
	xu1[0] = tau;
	testsig2 = call_vegas("DY_NLO_LP_corr", xl1, xu1, params);
	double Dtermsres = testsig1.res + testsig2.res;
	double Dtermser = testsig1.err + testsig2.err;
	xl1[0] = tau;
	xu1[0] = 1;
	testsig1 = call_vegas("DY_NLO_NLP", xl1, xu1, params);
	double NLP = testsig1.res;
	double NLPer = testsig1.err;
	testsig1 = call_vegas("DY_NLO_NNLP", xl1, xu1, params);
	double NNLP = testsig1.res;
	double NNLPer = testsig1.err;
	testsig1 = call_vegas("DY_NLO_full", xl1, xu1, params);
	double Ltermsres = testsig1.res;
	double Ltermser = testsig1.err;
	testsig1 = call_vegas("DY_NLO_qg_full", xl1, xu1, params);
	double qgfullres = testsig1.res;
	double qgfuller = testsig1.err;
	
	xl1[0] = tau;
	xu1[0] = 1;
	testsig1 = call_vegas("DY_NLO_LP_int", xl1, xu1, params);
	double LPint = testsig1.res;
	double LPerint = testsig1.err;
	testsig1 = call_vegas("DY_NNLO_LP_int", xl1, xu1, params);
	double NNLO_LPint = testsig1.res;
	double NNLO_LPerint = testsig1.err;
	testsig1 = call_vegas("DY_NLO_NLP_int", xl1, xu1, params);
	double NLPint = testsig1.res;
	double NLPerint = testsig1.err;
	testsig1 = call_vegas("DY_NLO_NNLP_int", xl1, xu1, params);
	double NNLPint = testsig1.res;
	double NNLPerint = testsig1.err;
	testsig1 = call_vegas("DY_NLO_full_int", xl1, xu1, params);
	double fullresint = testsig1.res;
	double fullerint = testsig1.err;
	testsig1 = call_vegas("DY_NLO_qg_full_int", xl1, xu1, params);
	double qgfullresint = testsig1.res;
	double qgfullerint = testsig1.err;
	testsig1 = call_vegas("DY_NLO_delta", xl, xu, params);
	double zdelta = testsig1.res;
	double deltaer = testsig1.err;
	//testsig1 = call_vegas("LO", xl, xu, params);
	//double LO = testsig1.res;
	//double LOer = testsig1.err;
	
	
	
	cout << "===========================================================================" << endl;
    cout << endl << "dsigma/dQ^2 for DY" << endl <<endl;
    cout << "===========================================================================" << endl;
    
    cout << endl << "INPUT PARAMETERS" << endl ;
    cout << "..........................................................................." << endl;
    cout << "S = " << S2 << " GeV^2" << endl<< "Q = " << Q << " GeV" << endl << "tau = " << tau << endl << "tau*sigma_0 = "<< LO_factor() << endl;
    
    cout << "---------------------------------------------------------------------------" << endl;
    
    cout << endl << "RESULTS" << endl;
    
    cout << "..........................................................................." << endl;
    cout << endl << "Total direct integration (LO): " << resLO << " pb/GeV^2 (" << erLO << ")" << endl;
    cout << endl << "Total direct integration (NLO): " << Dtermsres+Ltermsres+zdelta << " pb/GeV^2 (" << Dtermser+Ltermser+deltaer << ")" << endl;
    cout << endl << "Total direct integration (NNLO): " << DtermsresNNLO+NNLO_Ltermsres+NNLO_delta << " pb/GeV^2 (" << DtermserNNLO+NNLO_Ltermser+NNLO_deltaer << ")" << endl;
    cout << "---------------------------------------------------------------------------" << endl;

    
    cout << endl << "POWER EXPANSION (NLO)" << endl;
    cout << "..........................................................................." << endl;
    cout << "LP         : " << Dtermsres  << " pb/GeV^2 (" << Dtermser  << ")"  << endl;
    cout << "LP (+delta): " << Dtermsres + zdelta << " pb/GeV^2 (" << Dtermser + deltaer << ")"  << "  ; fractional: "<<(Dtermsres + zdelta)/(Dtermsres+Ltermsres+zdelta) << endl;
    
    cout << "NLP        : " << NLP << " pb/GeV^2 (" << NLPer << ")"  << "  ; fractional: "<<(NLP)/(Dtermsres+Ltermsres+zdelta) << endl;
    cout << "NNLP       : " << NNLP << " pb/GeV^2 (" << NNLPer << ")"  << " ; fractional: "<<(NNLP)/(Dtermsres+Ltermsres+zdelta) << endl;
    cout << "qg         : " << qgfullres << " pb/GeV^2 (" << qgfuller << ")"  << endl;
     cout << endl << "Integrated (NLO)" << endl;
    cout << "LP         : " << LPint  << " pb/GeV^2 (" << LPerint  << ")"  << endl;
    cout << "LP (+delta): " << LPint + zdelta << " pb/GeV^2 (" << LPerint + deltaer << ")"  << "  ; fractional: "<<(LPint + zdelta)/(fullresint+zdelta) << endl;
    cout << "NLP        : " << NLPint << " pb/GeV^2 (" << NLPerint << ")"  << "  ; fractional: "<<(NLPint)/(fullresint+zdelta) << endl;
    cout << "NNLP       : " << NNLPint << " pb/GeV^2 (" << NNLPerint << ")"  << " ; fractional: "<<(NNLPint)/(fullresint+zdelta) << endl;
    cout << "qg         : " << qgfullresint << " pb/GeV^2 (" << qgfullerint << ")"  << endl;
    
	cout << endl << "CUMULATIVE" << endl;
    cout << "..........................................................................." << endl;
    cout << "LP  : " << Dtermsres + zdelta  << " pb/GeV^2 " << " ; fractional: "<<(Dtermsres + zdelta)/(Dtermsres+Ltermsres+zdelta) << endl;
    cout << "NLP : " << NLP + Dtermsres + zdelta << " pb/GeV^2 " << " ; fractional: "<<(NLP + Dtermsres + zdelta)/(Dtermsres+Ltermsres+zdelta) << endl;
    cout << "NNLP: " << NNLP + NLP + Dtermsres + zdelta << " pb/GeV^2 " << " ; fractional: "<<(NNLP + NLP + Dtermsres + zdelta)/(Dtermsres+Ltermsres+zdelta) << endl;
    cout << "---------------------------------------------------------------------------" << endl;
    
	cout << endl << "POWER EXPANSION (NNLO)" << endl;
	cout << "..........................................................................." << endl;
    cout << "LP         : " << DtermsresNNLO  << " pb/GeV^2 (" << DtermserNNLO  << ")"  << endl;
    cout << "LP (+delta): " << DtermsresNNLO + NNLO_delta << " pb/GeV^2 (" << DtermserNNLO + NNLO_deltaer << ")"  << "  ; fractional: "<<(DtermsresNNLO + NNLO_delta)/(DtermsresNNLO+NNLO_Ltermsres+NNLO_delta) << endl;
    
    cout << "NLP        : " << NNLO_NLP << " pb/GeV^2 (" << NNLO_NLPer << ")"  << "  ; fractional: "<<(NNLO_NLP)/(DtermsresNNLO+NNLO_Ltermsres+NNLO_delta) << endl;
    cout << "NNLP       : " << NNLO_NNLP << " pb/GeV^2 (" << NNLO_NNLPer << ")"  << " ; fractional: "<<(NNLO_NNLP)/(DtermsresNNLO+NNLO_Ltermsres+NNLO_delta) << endl;
    cout << "NNNLP       : " << NNLO_NNNLP << " pb/GeV^2 (" << NNLO_NNNLPer << ")"  << " ; fractional: "<<(NNLO_NNNLP)/(DtermsresNNLO+NNLO_Ltermsres+NNLO_delta) << endl;
    cout << endl << "Integrated (NNLO)" << endl;
    cout << "LP         : " << NNLO_LPint  << " pb/GeV^2 (" << NNLO_LPerint  << ")"  << endl;
    cout << "LP (+delta): " << NNLO_LPint + NNLO_delta << " pb/GeV^2 (" << NNLO_LPerint + NNLO_deltaer << ")"  << "  ; fractional: "<<(NNLO_LPint + NNLO_delta)/(DtermsresNNLO+NNLO_Ltermsres+NNLO_delta) << endl;
   
    cout << endl << "CUMULATIVE (NNLO)" << endl;
    cout << "..........................................................................." << endl;
    cout << "LP  : " << DtermsresNNLO + NNLO_delta  << " pb/GeV^2 " << " ; fractional: "<<(DtermsresNNLO + NNLO_delta)/(DtermsresNNLO+NNLO_Ltermsres+NNLO_delta) << endl;
    cout << "NLP : " << NNLO_NLP + DtermsresNNLO + NNLO_delta << " pb/GeV^2 " << " ; fractional: "<<(NNLO_NLP + DtermsresNNLO + NNLO_delta)/(DtermsresNNLO+NNLO_Ltermsres+NNLO_delta) << endl;
    cout << "NNLP: " << NNLO_NNLP + NNLO_NLP + DtermsresNNLO + NNLO_delta << " pb/GeV^2 " << " ; fractional: "<<(NNLO_NNLP + NNLO_NLP + DtermsresNNLO + NNLO_delta)/(DtermsresNNLO+NNLO_Ltermsres+NNLO_delta) << endl;
    cout << "NNNLP: " << NNLO_NNNLP + NNLO_NNLP + NNLO_NLP + DtermsresNNLO + NNLO_delta << " pb/GeV^2 " << " ; fractional: "<<(NNLO_NNNLP +NNLO_NNLP + NNLO_NLP + DtermsresNNLO + NNLO_delta)/(DtermsresNNLO+NNLO_Ltermsres+NNLO_delta) << endl;
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
