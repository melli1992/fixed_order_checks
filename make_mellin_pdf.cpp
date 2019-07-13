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
#include <string.h>
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

	std::unordered_map<double, vector<vector<double>>> u = {
		{100,{{1},{2},{3}}},
			{500,{{2},{4},{5}}},
			{1000,{{3},{3},{4}}}
			/*  {100,{1,2,3,4}},
				{500,{1,2,3,4}},
				{1000,{1,2,3,4}} */
		};

		cout << u[100][0][0] << endl;
		cout << u[100][1][0] << endl;
    // Iterate and print keys and values of unordered_map
  //  for( const auto& n : u ) {
  //      std::cout << "Key:[" << n.first << "] Value:[" << n.second << "]\n";
  //  }
		//////////////////////////////////////////////
	/// predefinition of everything, setting it up
	//////////////////////////////////////////////
	/*std::unordered_map<double,vector<double> > fitcoefficients = {
        {100,{1,2,3,4}},
        {500,{1,2,3,4}},
        {1000,{10,20,30,40}}
    };
	for( const auto& n : fitcoefficients) {
	        std::cout << "Key:[" << n.first << "] Value:[" << n.second << "]\n";
	    }
	cout << fitcoefficients[100] << endl;*/


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
	double eta = 1.5;
	lumni_params params = {z, Q, 2*Q/S, exp(eta), exp(-eta), 0,0,0};
	params.z = 0.5;
	///////////////
	/*ofstream output, output2;
	string q_str = "fit_pdfs/2muF" + to_string(Q) +"_as"+to_string(alphas_Q)+"_"+setname;
	string q2_str = "fit_pdfs/muF" + to_string(Q) +"_"+setname;
	q2_str = q2_str + "_pdfoutput.txt";
	q_str = q_str + "_coefficients.txt";
	output.open(q_str.c_str());
	output2.open(q2_str.c_str()); //.c_str() needed to do a constant string conversion
	output2 << "x xf(x) f(x)" << endl;
	//////////////////////////////////////////////
    cout << q_str << endl;
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
		output << "muF=" << muF <<" GeV, flavor=" << flav << endl;
		/*for(int i = 0; i < 50; i++)
		{
			cout << "computing a_"<<i << endl;
			results out1;
			params.coefficient = i;
			if(flav!=0){params.flavor =flav;}
			else if(flav==0){params.flavor =21;}
			out1 = call_vegas(init_vegas_coefficients(params.coefficient),params,false,true);
			cout << "a_"<< params.coefficient << " " << out1.res << " " << out1.err << endl;
		}
	}
	output2.close();
	output.close();*/
	return 0;
}
