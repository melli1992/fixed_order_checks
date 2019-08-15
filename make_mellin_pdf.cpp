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
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <string>
#include <sstream>
#include <unordered_map>
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
	LHAPDF::PDFSet setk(setname);
	update_defaults();
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
	double z;
	double eta = 1.5;
	lumni_params params = {z, Q, 2*Q/S, exp(eta), exp(-eta), 0,0,0};
	params.z = 0.5;
	///////////////
	ofstream output, output2;
	ostringstream x_convert;
	x_convert << Q;
	string Qstring  = x_convert.str();
	x_convert << alphas_Q;
	string asstring  = x_convert.str();
	//string q2_str = "fit_pdfs/"+setname+"/muF" + Qstring +"_"+setname;
	string q2_str = "fit_pdfs/"+setname+"/muF" + to_string((float) Q) +"_"+setname;
	
	q2_str = q2_str + "_pdfoutput.txt";
	
	output2.open(q2_str.c_str()); //.c_str() needed to do a constant string conversion
	output2 << "x xf(x) f(x)" << endl;
	//////////////////////////////////////////////
    cout << q2_str << endl;
    for(int flav = -5; flav < 6; flav++)
	{
		cout << flav << endl;
		cout << endl;
		cout << "ErrorType: " << setk.errorType() << endl;
		cout << "ErrorConfLevel: " << setk.errorConfLevel() << endl;


		double MAX = 1;
		output2 << "muF=" << muF <<" GeV, flavor=" << flav << endl;
		for (double i = 1E-10, increment = 1E-11, counter = 1; i <= MAX; i += increment) {
			vector<double> xfx;
			for (size_t imem = 0; imem <= nmem; imem++) {
			xfx.push_back(pdfs[imem]->xfxQ(flav,i,muF));
			//cout << pdfs[imem]->type() << endl;
			}
			const LHAPDF::PDFUncertainty xErr = setk.uncertainty(xfx, -1);
			output2 << i << " " << (xErr.central) << " " << (xErr.errsymm) << " " << xErr.errplus << " " << xErr.errminus << endl;

			if (counter == 100) {
				increment *= 10;
				counter = 1;
			}
			++counter;
		}
	}
	output2.close();
	return 0;
}
