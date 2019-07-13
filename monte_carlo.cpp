#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include "parameters.h"
#include "monte_carlo.h"
#include "deriv_pdf.h"
#include "k_factors_dy.h"
#include "k_factors_higgs.h"
#include "k_factors_nnlo_dy.h"
#include "k_factors_prompt_photon.h"
#include "mellin_functions.h"
#include "mellin_pdf.h"

using namespace std;

/////////////////////////////////////////////////////////////////////
/// this contains a call to vegas and a functionint to write the results
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////
/// display results of the MC integration
/////////////////////////////////////////
void display_results(string title, double &result, double &error){
	cout << title << endl;
	cout << "result = " << result << endl;
	cout << "sigma = " << error  << endl;
}



/////////////////////////////////////////////////////////
/// initializes the function for the mellin coefficients
/////////////////////////////////////////////////////////
functionint init_vegas_coefficients(int n){
	functionint integrand;
	if(n==0){
					integrand.G.f =&vegas_coefficients_0;
					integrand.G.dim = 1;
					integrand.xl = {-1.};
					integrand.xu = {1.};
				}
	else{
								integrand.G.f =&vegas_coefficients_n;
								integrand.G.dim = 1;
								integrand.xl = {-1.};
								integrand.xu = {1.};
							}
	return integrand;
}

/////////////////////////////////////////////////////////
/// initializes the function to integrate in the photon case
/////////////////////////////////////////////////////////
functionint init_vegas_mellin(std::string process, double cmp_set, double phi_set){
	functionint integrand;
			if (process == "full"){
					integrand.G.f =&vegas_fofx2_full;
					integrand.G.dim = 2;
					integrand.xl = {0.,0.};
					integrand.xu = {1.,1.};
				}
		  else if (process == "deriv"){
						integrand.G.f =&vegas_fofx2_deriv;
						integrand.G.dim = 2;
						integrand.xl = {0.,0.};
						integrand.xu = {1.,1.};
					}
			else if (process == "defor"){
							integrand.G.f =&vegas_fofx2_defor;
							integrand.G.dim = 2;
							integrand.xl = {0.,-1.};
							integrand.xu = {1.,0.};
						}
			else if (process == "nspace"){
							integrand.G.f =&vegas_fofx2_Nspace;
							integrand.G.dim = 1;
							integrand.xl = {0.};
							integrand.xu = {1.};
						}
			else if (process == "LOnul"){
							integrand.G.f =&vegas_sigma0_nomel;
							integrand.G.dim = 1;
							integrand.xl = {tau};
							integrand.xu = {1.};
							//integrand.xl = {-1.,-1.};
							//integrand.xu = {0,0};
						}
			else if (process == "LOderiv"){
							integrand.G.f =&vegas_sigma0_deriv;
							integrand.G.dim = 3;
							integrand.xl = {0.,0.,0.};
							integrand.xu = {1.,1.,1.};
						}
				else if (process == "LOdefor"){
								integrand.G.f =&vegas_sigma0_defor;
								integrand.G.dim = 3;
								integrand.xl = {0.,-1.,-1.};
								integrand.xu = {1.,0.,0.};
							}
					else if (process == "resumdefor"){
									integrand.G.f =&vegas_resum_defor;
									integrand.G.dim = 3;
									integrand.xl = {0.,-1.,-1.};
									integrand.xu = {1.,0.,0.};
								}
	return integrand;
}

/////////////////////////////////////////////////////////
/// initializes the function to integrate in the photon case
/////////////////////////////////////////////////////////
functionint init_vegas_pf(std::string order, std::string power, std::string process, int power_number, bool integrated){
	functionint integrand;

	if (order == "NLO"){
		if (process == "qqbar"){
			if (power_number != 0){
				integrand.G.f =&vegas_FP_qqbartogg_power;
				integrand.G.dim = 2;
				integrand.xl = {0.,0.};
				integrand.xu = {1.,1.};
			}
			else if (power =="LP"){
					integrand.G.f =&vegas_FP_qqbartogg_LP;
					integrand.G.dim = 2;
					integrand.xl = {0.,0.};
					integrand.xu = {1.,1.};
				}
				else if (power =="LP_corr"){
						integrand.G.f =&vegas_FP_qqbartogg_LP_corr;
						integrand.G.dim = 2;
						integrand.xl = {0.,0.};
						integrand.xu = {1.,1.};
					}
			else if (power =="full"){
						integrand.G.f =&vegas_FP_qqbartogg_full;
						integrand.G.dim = 2;
						integrand.xl = {0.,0.};
						integrand.xu = {1.,1.};
					}
			else if (power =="delta"){
							integrand.G.f =&vegas_FP_qqbartogg_delta;
							integrand.G.dim = 1;
							integrand.xl = {0.};
							integrand.xu = {1.};
						}
			else{
				cout << process << " " << integrated << " " << order << " " << power << endl;
				cout << "Wrong process specified" << endl;
				exit(0);
			}
		}
	}
	else{
		cout << process << " " << integrated << " " << order << " " << power << endl;
		cout << "Wrong order specified" << endl;
		exit(0);
	}
	return integrand;
}

/////////////////////////////////////////////////////////
/// initializes the function to integrate in the higgs case
/////////////////////////////////////////////////////////
functionint init_vegas_higgs(std::string order, std::string power, std::string process, int power_number, bool integrated){
	functionint integrand;

	if(order == "LO"){
	    integrand.G.f =&vegas_higgs_LO;
	    integrand.G.dim=1;
		integrand.xl = {tau};
		integrand.xu = {1.};
	}
	else if (order == "NLO"){
		if (process == "gg"){
			if (power_number != 0){
				integrand.G.f =&vegas_higgs_NLO_gg_power;
				integrand.G.dim = 2;
				integrand.xl = {tau,0.};
				integrand.xu = {1.,1.};
			}
			else if (power == "LP"){
				integrand.G.f =&vegas_higgs_NLO_gg_LP;
				integrand.G.dim = 2;
				integrand.xl = {tau,0.};
				integrand.xu = {1.,1.};
			}
			else if (power == "LP_corr"){
				integrand.G.f =&vegas_higgs_NLO_gg_LP_corr;
				integrand.G.dim = 2;
				integrand.xl = {0.,0.};
				integrand.xu = {tau,1.};
			}
			else if (power == "full"){
				integrand.G.f =&vegas_higgs_NLO_gg_full;
				integrand.G.dim = 2;
				integrand.xl = {tau,0.};
				integrand.xu = {1.,1.};
			}
			else if (power == "delta"){
				integrand.G.f =&vegas_higgs_NLO_gg_delta;
				integrand.G.dim = 1;
				integrand.xl = {tau};
				integrand.xu = {1.};
				}
			else{
				cout << process << " " << integrated << " " << order << " " << power << endl;
				cout << "Wrong power (non-int) specified" << endl;
				exit(0);
			}
		}
		else if (process == "qg"){
			if (power_number != 0){
				integrand.G.f =&vegas_higgs_NLO_qg_power;
				integrand.G.dim = 2;
				integrand.xl = {tau,0.};
				integrand.xu = {1.,1.};
			}
			else if (power == "full"){
				integrand.G.f =&vegas_higgs_NLO_qg_full;
				integrand.G.dim = 2;
				integrand.xl = {tau,0.};
				integrand.xu = {1.,1.};
			}
			else{
				cout << process << " " << integrated << " " << order << " " << power << endl;
				cout << "Wrong power (non-int) specified" << endl;
				exit(0);
			}
		}
		else if (process == "qqbar"){
			if (power_number != 0){
				integrand.G.f =&vegas_higgs_NLO_qqbar_power;
				integrand.G.dim = 2;
				integrand.xl = {tau,0.};
				integrand.xu = {1.,1.};
			}
			else if (power == "full"){
				integrand.G.f =&vegas_higgs_NLO_qqbar_full;
				integrand.G.dim = 2;
				integrand.xl = {tau,0.};
				integrand.xu = {1.,1.};
			}
			else{
				cout << process << " " << integrated << " " << order << " " << power << endl;
				cout << "Wrong power (non-int) specified" << endl;
				exit(0);
			}
		}
		else{
			cout << process << " " << integrated << " " << order << " " << power << endl;
			cout << "Wrong process specified" << endl;
			exit(0);
		}
	}
	else{
		cout << process << " " << integrated << " " << order << " " << power << endl;
		cout << "Wrong order specified" << endl;
		exit(0);
	}
	return integrand;
}

/////////////////////////////////////////////////////////
/// initializes the function to integrate in the DY case
/////////////////////////////////////////////////////////
functionint init_vegas_dy(std::string order, std::string power, std::string process, bool integrated){
	functionint integrand;

	if(order == "LO"){
	    integrand.G.f =&vegas_LO;
	    integrand.G.dim=1;
		integrand.xl = {tau};
		integrand.xu = {1.};
	}
	else if (order == "NLO"){
		if (process == "qqbar"){
			if (integrated == false){
				if (power == "LP"){
					integrand.G.f =&vegas_sig_LP_1;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "LP_corr"){
					integrand.G.f =&vegas_sig_LP_correction;
					integrand.G.dim = 2;
					integrand.xl = {0.,0.};
					integrand.xu = {tau,1.};
				}
				else if (power == "NLP"){
					integrand.G.f =&vegas_sig_NLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNLP"){
					integrand.G.f =&vegas_sig_NNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				//else if (power == "NNNLP"){integrand.G.f =&vegas_NNLO_qqbar_NNNLP; integrand.G.dim = 2;}
				else if (power == "full"){
					integrand.G.f =&vegas_sig_full;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "delta"){
					integrand.G.f =&vegas_sig_delta;
					integrand.G.dim = 1;
					integrand.xl = {tau};
					integrand.xu = {1.};
					}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power (non-int) specified" << endl;
					exit(0);
				}
			}
			else{
				if (power == "LP"){
					integrand.G.f =&vegas_sig_LP_int;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NLP"){
					integrand.G.f =&vegas_sig_NLP_int;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNLP"){
					integrand.G.f =&vegas_sig_NNLP_int;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "full"){
					integrand.G.f =&vegas_sig_full_int;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power (int) specified" << endl;
					exit(0);
				}
			}
		}
		else if (process == "qg"){
			if (integrated == false){
				if (power == "full1"){
					integrand.G.f =&vegas_qg_full;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "full"){
					integrand.G.f =&vegas_NLO_qg_full;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NLP"){
					integrand.G.f =&vegas_NLO_qg_NLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNLP"){
					integrand.G.f =&vegas_NLO_qg_NNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNNLP"){
					integrand.G.f =&vegas_NLO_qg_NNNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power (non int) specified" << endl;
					exit(0);
				}
			}
			else{
				if (power == "full"){
					integrand.G.f =&vegas_qg_full_int;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power (int) specified" << endl;
					exit(0);
				}
			}
		}

		else{
			cout << process << " " << integrated << " " << order << " " << power << endl;
			cout << "Wrong process specified" << endl;
			exit(0);
		}
	}
	else if (order == "NNLO"){
		if (integrated == false){
			if (process == "qqbar"){
				if (power == "LP"){
					integrand.G.f =&vegas_NNLO_qqbar_LP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
					}
				else if (power == "LP_corr"){
					integrand.G.f =&vegas_NNLO_qqbar_LP_correction;
					integrand.G.dim = 2;
					integrand.xl = {0,0.};
					integrand.xu = {tau,1.};
					}
				else if (power == "NLP"){
					integrand.G.f =&vegas_NNLO_qqbar_NLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
					}
				else if (power == "NNLP"){
					integrand.G.f =&vegas_NNLO_qqbar_NNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
					}
				else if (power == "NNNLP"){
					integrand.G.f =&vegas_NNLO_qqbar_NNNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
					}
				else if (power == "full"){
					integrand.G.f =&vegas_NNLO_qqbar_full;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
					}
				else if (power == "delta"){
					integrand.G.f =&vegas_NNLO_qqbar_delta;
					integrand.G.dim = 1;
					integrand.xl = {tau};
					integrand.xu = {1.};
					}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power specified" << endl;
					exit(0);
				}
			}
			else if (process == "qg"){
				if (power == "full"){
					integrand.G.f =&vegas_NNLO_qg_full;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NLP"){
					integrand.G.f =&vegas_NNLO_qg_NLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNLP"){
					integrand.G.f =&vegas_NNLO_qg_NNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNNLP"){
					integrand.G.f =&vegas_NNLO_qg_NNNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power specified" << endl;
					exit(0);
				}
			}
			else if (process == "gg"){
				if (power == "full"){
					integrand.G.f =&vegas_NNLO_gg_full;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NLP"){
					integrand.G.f =&vegas_NNLO_gg_NLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNLP"){
					integrand.G.f =&vegas_NNLO_gg_NNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNNLP"){
					integrand.G.f =&vegas_NNLO_gg_NNNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power specified" << endl;
					exit(0);
				}
			}
			else if (process == "qq"){
				if (power == "full"){
					integrand.G.f =&vegas_NNLO_qq_full;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NLP"){
					integrand.G.f =&vegas_NNLO_qq_NLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNLP"){
					integrand.G.f =&vegas_NNLO_qq_NNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNNLP"){
					integrand.G.f =&vegas_NNLO_qq_NNNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power specified" << endl;
					exit(0);
				}
			}
			else if (process == "qqNI"){
				if (power == "full"){
					integrand.G.f =&vegas_NNLO_qqNI_full;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NLP"){
					integrand.G.f =&vegas_NNLO_qqNI_NLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNLP"){
					integrand.G.f =&vegas_NNLO_qqNI_NNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNNLP"){
					integrand.G.f =&vegas_NNLO_qqNI_NNNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power specified" << endl;
					exit(0);
				}
			}
			else if (process == "qqbarNI"){
				if (power == "full"){
					integrand.G.f =&vegas_NNLO_qqbarNI_full;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NLP"){
					integrand.G.f =&vegas_NNLO_qqbarNI_NLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNLP"){
					integrand.G.f =&vegas_NNLO_qqbarNI_NNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNNLP"){
					integrand.G.f =&vegas_NNLO_qqbarNI_NNNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power specified" << endl;
					exit(0);
				}
			}
			else if (process == "qbarqbarNI"){
				if (power == "full"){
					integrand.G.f =&vegas_NNLO_qbarqbarNI_full;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NLP"){
					integrand.G.f =&vegas_NNLO_qbarqbarNI_NLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNLP"){
					integrand.G.f =&vegas_NNLO_qbarqbarNI_NNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNNLP"){
					integrand.G.f =&vegas_NNLO_qbarqbarNI_NNNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power specified" << endl;
					exit(0);
				}
			}
			else if (process == "qbarqbar"){
				if (power == "full"){
					integrand.G.f =&vegas_NNLO_qbarqbar_full;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NLP"){
					integrand.G.f =&vegas_NNLO_qbarqbar_NLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNLP"){
					integrand.G.f =&vegas_NNLO_qbarqbar_NNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else if (power == "NNNLP"){
					integrand.G.f =&vegas_NNLO_qbarqbar_NNNLP;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power specified" << endl;
					exit(0);
				}
			}
			else{
				cout << process << " " << integrated << " " << order << " " << power << endl;
				cout << "Wrong process (non-int) specified" << endl;
				exit(0);
			}
		}
		else{
			if (process == "qqbar"){
				if (power == "LP"){
					integrand.G.f =&vegas_NNLO_LP_int;
					integrand.G.dim = 2;
					integrand.xl = {tau,0.};
					integrand.xu = {1.,1.};
				}
				else{
					cout << process << " " << integrated << " " << order << " " << power << endl;
					cout << "Wrong power (int) specified" << endl;
					exit(0);
				}
			}
			else{
				cout << process << " " << integrated << " " << order << " " << power << endl;
				cout << "Wrong process (int) specified" << endl;
				exit(0);
			}
		}
	}
	else{
		cout << process << " " << integrated << " " << order << " " << power << endl;
		cout << "Wrong order specified" << endl;
		exit(0);
	}
	return integrand;
}


//////////////////////////////////////////////////////////
/// call vegas, store the results in a result struc
/// handles all the declaration
/// input is the functionint to integrate
/// the parameter structure (may change)
/// the number of dimensions
/// the upper and lower bounds of the integrand
//////////////////////////////////////////////////////////
results call_vegas(functionint integrand, lumni_params params, bool verbal, bool high_prec){
	  double res(0),err(0);
	  int MAX_ITER = 10;
	  const gsl_rng_type *T;
	  gsl_rng *r; //the random number generator



	  integrand.G.params = &params;
	  size_t calls = 50000;
	  gsl_rng_env_setup ();
	  T = gsl_rng_default;
	  r = gsl_rng_alloc (T);


	  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (integrand.G.dim);
	  int n_iter = 0;

	  if(high_prec==true)
	  {gsl_monte_vegas_params *params_run = (gsl_monte_vegas_params*)malloc( sizeof(gsl_monte_vegas_params) );
	  gsl_monte_vegas_params_get(s, params_run);
	  params_run->iterations = MAX_ITER;
	  params_run->alpha = 2.0;
	  params_run->stage = 0;
	  gsl_monte_vegas_params_set(s, params_run);

	  gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, calls/50, r, s, &res, &err);
	  params_run->stage = 1;
	  gsl_monte_vegas_params_set(s, params_run);
	  gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, calls, r, s, &res, &err);
	  params_run->stage = 1;
	  gsl_monte_vegas_params_set(s, params_run);
	  if(verbal){display_results ("vegas warm-up", res, err);}
	  do{
			n_iter+=1;
			params_run->stage = 2;
	        params_run->alpha = 1.5;
	        gsl_monte_vegas_params_set(s, params_run);
			gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, calls, r, s, &res, &err);
			if(verbal){ cout << "round " << n_iter << ", result = " << res << " error = " << err << " chisq/dof = " << gsl_monte_vegas_chisq (s) << " err/res = " << err/res << endl;}
			if(n_iter >= MAX_ITER){break;}
		}
	  while ((isnan(res))||(fabs (err/res) > 1E-4));
	  params_run->stage = 1;
	  gsl_monte_vegas_params_set(s, params_run);

	  gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, calls*10, r, s, &res, &err);

	  params_run->stage = 3;
	  gsl_monte_vegas_params_set(s, params_run);
	  n_iter = 0;
	  do{
			n_iter+=1;
			gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, calls/5, r, s, &res, &err);
			if(verbal){ cout << "round " << n_iter << ", result = " << res << " error = " << err << " chisq/dof = " << gsl_monte_vegas_chisq (s) << " err/res = " << err/res << endl;}
			if(n_iter > 100){break;}
		}
	  while ((isnan(res))||(fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.05)||(fabs (err/res) > 1E-4));

	  if(verbal){display_results ("vegas final", res, err);}
	  }
	 if(high_prec==false)
	 {
		 do{
			n_iter+=1;
			gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, calls/5, r, s, &res, &err);
			if(verbal){ cout << "round " << n_iter << ", result = " << res << " error = " << err << " chisq/dof = " << gsl_monte_vegas_chisq (s) << " err/res = " << err/res << endl;}
			if(n_iter > 100){break;}
	  }
	  while ((isnan(res))||(fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.05)||(fabs (err/res) > 1E-4));
      }
	  gsl_monte_vegas_free (s);

	  gsl_rng_free (r);
	  results result = {res,err};
	  return result;

}
