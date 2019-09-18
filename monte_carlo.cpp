#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include "parameters.h"
#include "monte_carlo.h"
#include "deriv_pdf.h"
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


//////////////////////////////////////////////////////////////////////
/// initializes the function to integrate for the PDFs
//////////////////////////////////////////////////////////////////////
functionint init_vegas_pdfs(std::string lumi){
	functionint integrand;
	std::string process= lumi;
	if(process == "gg"){
			integrand.G.f =&vegas_lumi_gg;
			integrand.G.dim = 1;
			integrand.xl = {tau};
			integrand.xu = {1.};
		}
	else if(process == "qq"){
			integrand.G.f =&vegas_lumi_qq;
			integrand.G.dim = 1;
			integrand.xl = {tau};
			integrand.xu = {1.};
		}
	else if(process == "qqbar"){
			integrand.G.f =&vegas_lumi_qqbar;
			integrand.G.dim = 1;
			integrand.xl = {tau};
			integrand.xu = {1.};
		}
	else if(process == "qqNI"){
			integrand.G.f =&vegas_lumi_qqNI;
			integrand.G.dim = 1;
			integrand.xl = {tau};
			integrand.xu = {1.};
		}
	else if(process == "qg"){
			integrand.G.f =&vegas_lumi_qg;
			integrand.G.dim = 1;
			integrand.xl = {tau};
			integrand.xu = {1.};
		}
	else{cout << "this lumi does not exist" << endl;
exit(0);}
		return integrand;
}

//////////////////////////////////////////////////////////////////////
/// initializes the function to integrate in the resummed cases
//////////////////////////////////////////////////////////////////////
functionint init_vegas_mellin(std::string process, std::string order){
	functionint integrand;
	if(process == "test"){ //test cases for full, derivative method, deformation method and the direct N space formula
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
	}
	else{cout << "Process " << process <<" is not implemented" << endl; exit(0);}


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
