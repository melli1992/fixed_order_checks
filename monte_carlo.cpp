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
/// this contains a call to vegas and a function to write the results
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////
/// display results of the MC integration
/////////////////////////////////////////
void display_results(string title, double &result, double &error){
	cout << title << endl;
	cout << "result = " << result << endl;
	cout << "sigma = " << error  << endl;
}


//////////////////////////////////////////////////////////
/// call vegas, store the results in a result struc
/// handles all the declaration
/// input is the function to integrate
/// the parameter structure (may change)
/// the number of dimensions
/// the upper and lower bounds of the integrand
//////////////////////////////////////////////////////////
results call_vegas(std::string integrand, double *xl, double *xu, lumni_params params, bool verbal, bool high_prec){
	  double res(0),err(0);
	  int MAX_ITER = 10;
	  const gsl_rng_type *T; 
	  int ndim(0);
	  gsl_rng *r; //the random number generator

	  gsl_monte_function G;


	  if(integrand=="lumni"){ G.f = &vegas_sum_pdf_weigthed; ndim=1;}
	  if(integrand=="DYNLO"){ G.f = &vegas_k_NLO_dy_gg; ndim=1;}
	  if(integrand=="valence"){ G.f = &vegas_pdf_up_minus_upbar; ndim=1;}
	  if(integrand=="momcons"){ G.f = &vegas_pdf_mom_consv; ndim=1;}
	  if(integrand=="LO"){ G.f =&vegas_LO; ndim=1;}
	  if(integrand=="DY_NLO_LP1"){ G.f =&vegas_sig_LP_1; ndim=2;}
	  if(integrand=="DY_NLO_LP_corr"){ G.f =&vegas_sig_LP_correction; ndim=2;}
	  if(integrand=="DY_NNLO_LP"){ G.f =&vegas_NNLO_qqtogg_LP; ndim=2;}
	  if(integrand=="DY_NNLO_LP_corr"){ G.f =&vegas_NNLO_qqtogg_LP_correction; ndim=2;}
	  if(integrand=="DY_NLO_NLP"){ G.f =&vegas_sig_NLP; ndim=2;}
	  if(integrand=="DY_NLO_NNLP"){ G.f =&vegas_sig_NNLP; ndim=2;}
	  if(integrand=="DY_NNLO_NLP"){ G.f =&vegas_NNLO_qqtogg_NLP; ndim=2;}
	  if(integrand=="DY_NNLO_NNLP"){ G.f =&vegas_NNLO_qqtogg_NNLP; ndim=2;}
	  if(integrand=="DY_NNLO_NNNLP"){ G.f =&vegas_NNLO_qqtogg_NNNLP; ndim=2;}
	  if(integrand=="DY_LO"){ G.f =&vegas_LO; ndim=1;}
	  if(integrand=="DY_NLO_delta"){ G.f =&vegas_sig_delta; ndim=1;}
	  if(integrand=="sum_PDF"){ G.f =&vegas_sum_pdf; ndim=1;}
	  if(integrand=="DY_NNLO_delta"){ G.f =&vegas_NNLO_qqtogg_delta; ndim=1;}
	  if(integrand=="DY_NLO_full"){ G.f =&vegas_sig_full; ndim=2;}
	  if(integrand=="DY_NNLO_full"){ G.f =&vegas_NNLO_qqtogg_full; ndim=2;}
	  if(integrand=="DY_NLO_qg_full"){ G.f =&vegas_qg_full; ndim=2;}
	  if(integrand=="DY_NLO_LP_int"){ G.f =&vegas_sig_LP_int; ndim=2;}
	  if(integrand=="DY_NNLO_LP_int"){ G.f =&vegas_NNLO_LP_int; ndim=2;}
	  if(integrand=="DY_NLO_NLP_int"){ G.f =&vegas_sig_NLP_int; ndim=2;}
	  if(integrand=="DY_NLO_NNLP_int"){ G.f =&vegas_sig_NNLP_int; ndim=2;}
	  if(integrand=="DY_NLO_full_int"){ G.f =&vegas_sig_full_int; ndim=2;}
	  if(integrand=="DY_NLO_qg_full_int"){ G.f =&vegas_qg_full_int; ndim=2;}
	  if(integrand=="test1"){ G.f =&vegas_sig_testfunction1; ndim=1;}
	  if(integrand=="test2"){ G.f =&vegas_sig_testfunction2; ndim=1;}
	  
	  G.dim = ndim;
	  G.params = &params;
	  size_t calls = 50000;
	  gsl_rng_env_setup ();
	  T = gsl_rng_default;
	  r = gsl_rng_alloc (T);
	  
		  
	  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (ndim);
	  int n_iter = 0;
	 
	  if(high_prec==true)
	  {gsl_monte_vegas_params *params_run = (gsl_monte_vegas_params*)malloc( sizeof(gsl_monte_vegas_params) ); 
	  gsl_monte_vegas_params_get(s, params_run);
	  params_run->iterations = MAX_ITER;
	  params_run->alpha = 2.0; 
	  params_run->stage = 0;
	  gsl_monte_vegas_params_set(s, params_run);
	  gsl_monte_vegas_integrate (&G, xl, xu, ndim, calls/50, r, s, &res, &err);
	  params_run->stage = 1;
	  gsl_monte_vegas_params_set(s, params_run);
	  gsl_monte_vegas_integrate (&G, xl, xu, ndim, calls, r, s, &res, &err);		 
	  params_run->stage = 1;
	  gsl_monte_vegas_params_set(s, params_run);	
	  if(verbal){display_results ("vegas warm-up", res, err);}
	  do{ 
			n_iter+=1;
			params_run->stage = 2;
	        params_run->alpha = 1.5; 
	        gsl_monte_vegas_params_set(s, params_run);	
			gsl_monte_vegas_integrate (&G, xl, xu, ndim, calls, r, s, &res, &err);
			if(verbal){ cout << "round " << n_iter << ", result = " << res << " error = " << err << " chisq/dof = " << gsl_monte_vegas_chisq (s) << " err/res = " << err/res << endl;}
			if(n_iter >= MAX_ITER){break;}
		}
	  while (fabs (err/res) > 1E-4);
	  params_run->stage = 1;
	  gsl_monte_vegas_params_set(s, params_run);	
	  
	  gsl_monte_vegas_integrate (&G, xl, xu, ndim, calls*10, r, s, &res, &err);
	  
	  params_run->stage = 3;
	  gsl_monte_vegas_params_set(s, params_run);	
	  n_iter = 0;
	  do{ 
			n_iter+=1;
			gsl_monte_vegas_integrate (&G, xl, xu, ndim, calls/5, r, s, &res, &err);
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
			gsl_monte_vegas_integrate (&G, xl, xu, ndim, calls/5, r, s, &res, &err);
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
