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
#include "k_factors_nnlo_dy.h"

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

functionint init_vegas(std::string order, std::string power, std::string process, bool integrated){
	functionint integrand;
	
	//if(integrand=="lumni"){ G.f = &vegas_sum_pdf_weigthed; ndim=1;}
	//  if(integrand=="valence"){ G.f = &vegas_pdf_up_minus_upbar; ndim=1;}
	//  if(integrand=="momcons"){ G.f = &vegas_pdf_mom_consv; ndim=1;}  
	//  if(integrand=="sum_PDF"){ G.f =&vegas_sum_pdf; ndim=1;}
	
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
				//else if (power == "NNNLP"){integrand.G.f =&vegas_NNLO_qqtogg_NNNLP; integrand.G.dim = 2;}
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
					cout << "Wrong order specified" << endl;
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
					cout << "Wrong order specified" << endl;
					exit(0);
				}
			}
		}
		else if (process == "qg"){
			if (integrated == false){
				if (power == "full"){
					integrand.G.f =&vegas_qg_full; 
					integrand.G.dim = 2;
					integrand.xl = {tau,0.}; 
					integrand.xu = {1.,1.};
				}
				else{		
					cout << process << " " << integrated << " " << order << " " << power << endl;	
					cout << "Wrong order specified" << endl;
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
					cout << "Wrong order specified" << endl;
					exit(0);
				}
			}
		}
		
		else{		
			cout << process << " " << integrated << " " << order << " " << power << endl;	
			cout << "Wrong order specified" << endl;
			exit(0);
		}
	}
	else if (order == "NNLO"){
		if (process == "qqbar"){
			if (power == "LP"){
				integrand.G.f =&vegas_NNLO_qqtogg_LP; 
				integrand.G.dim = 2;
				integrand.xl = {tau,0.}; 
				integrand.xu = {1.,1.};
				}
			else if (power == "LP_corr"){
				integrand.G.f =&vegas_NNLO_qqtogg_LP_correction; 
				integrand.G.dim = 2;
				integrand.xl = {0,0.}; 
				integrand.xu = {tau,1.};
				}
			else if (power == "NLP"){
				integrand.G.f =&vegas_NNLO_qqtogg_NLP; 
				integrand.G.dim = 2;
				integrand.xl = {tau,0.}; 
				integrand.xu = {1.,1.};
				}
			else if (power == "NNLP"){
				integrand.G.f =&vegas_NNLO_qqtogg_NNLP; 
				integrand.G.dim = 2;
				integrand.xl = {tau,0.}; 
				integrand.xu = {1.,1.};
				}
			else if (power == "NNNLP"){
				integrand.G.f =&vegas_NNLO_qqtogg_NNNLP; 
				integrand.G.dim = 2;
				integrand.xl = {tau,0.}; 
				integrand.xu = {1.,1.};
				}
			else if (power == "full"){
				integrand.G.f =&vegas_NNLO_qqtogg_full; 
				integrand.G.dim = 2;
				integrand.xl = {tau,0.}; 
				integrand.xu = {1.,1.};
				}
			else if (power == "delta"){
				integrand.G.f =&vegas_NNLO_qqtogg_delta; 
				integrand.G.dim = 1;
				integrand.xl = {tau}; 
				integrand.xu = {1.};
				}
			else{		
				cout << process << " " << integrated << " " << order << " " << power << endl;	
				cout << "Wrong order specified" << endl;
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
results call_vegas(functionint integrand, /*gsl_monte_functionint G, double *xl, double *xu,*/ lumni_params params, bool verbal, bool high_prec){
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
	  while (fabs (err/res) > 1E-4);
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
	  while ((fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.05));
	 
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
	  while ((fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.05));
      }
	  gsl_monte_vegas_free (s);
  
	  gsl_rng_free (r);
	  results result = {res,err};
	  return result;
	
}
