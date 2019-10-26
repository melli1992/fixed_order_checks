
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "cuba.h"
#include "parameters.h"
#include "mellin_pdf.h"
#include "cuba_integration.h"
#include "resum_functions.h"
#include "k_factors_diboson.h"
#include "k_factors_dihiggs.h"
#include "k_factors_higgs.h"
#include "k_factors_nnlo_higgs.h"
#include "k_factors_dy.h"
#include "k_factors_nnlo_dy.h"
#include "deriv_pdf.h"

using namespace std;

int NDIM =3; //number of dimensions in integral
int NCOMP = 2; //number of components in integral
// integrand should be a function int integrand(ndim, x, ncomp, f, userdata, nvec, core)
// x(ndim, nvec), f(ncomp, nvec)
// nvec = number of samples x that is received, fill array f with integrand values
// returned value of integrand is irrelevant, unless it is -999, then should the integration be aborted
// nvec, userdata and core are optional in integrand function, not in the cuhre and vegas functions
void* USERDATA = NULL;
int NVEC = 1;
double EPSREL = 1.e-5;
double EPSABS = 1.e-32;
int VERBOSE = 0; //bit 0-1 is verbosity (so 0,1,2,3) , bit 2 = samples, bit 3 can improve smoothing, bit 4 = state file, bit 5 = vegas specific, bit 8 and higher is random number
int LAST = 4;
int SEED = 0; //other choices of seed
int MINEVAL = 1000; //min number of evaluations integrand
int MAXEVAL = 10000000; //max number of evaluations integrand
char* STATEFILE = NULL; //filename for storing internal state
int* SPIN = NULL;

//Vegas specific
int NSTART = 1000; //number evaluations to start with
int NINCREASE = 500; //increase in number of integrands
int NBATCH = 1000; //batch size for sampling
int GRIDNO = 0; //keep grids during one integratoin for the next one (if integrands are similar) if gridno > 0


//Cuhre
int KEY = 0; //cubature rule of degree?
int FUN = 3;
std::string order = "LO";
std::string process = "gg";
std::string power = "0";
bool fitted_pdfs = true;
int MAX_POW = 10;


//for dy
static int dy_integrand(const int *ndim, const cubareal xx[],const int *ncomp, cubareal ff[], void *userdata) {

  if(order == "LO"){
      if(!fitted_pdfs) ff[0] = DY_LO_factor()*(pdf_sum_qqbar_charge_weighted(tau+(1.-tau)*xx[0],tau))*(1.-tau);
      if(fitted_pdfs) ff[0] = DY_LO_factor()*(real(fit_sum_qqbar_charge_weighted(tau+(1.-tau)*xx[0],tau)))*(1.-tau);
  }
  else if(order == "LOfullZ"){ // Z full xsec
	  double q2 = xx[0]*xx[1]*S2;
	  double sum_pdf(0);
	  double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	  for(int i = 1; i <=5; i++){
		sum_pdf+= eq[i-1]*eq[i-1]*(pdfs[0]->xfxQ(i,xx[0],muF)*pdfs[0]->xfxQ(-i,xx[1],muF)+pdfs[0]->xfxQ(i,xx[1],muF)*pdfs[0]->xfxQ(-i,xx[0],muF));
		}
		double gammaZll = alphaEM*sqrt(mZ2)*(1.+pow(1.-4.*pow(sinw,2),2))/(48.*pow(sinw,2)*pow(cosw,2));
		double result = S2*pbunits*M_PI*alphaEM*gammaZll/(4.*sqrt(mZ2)*pow(sinw,2)*pow(cosw,2)*(pow(q2-mZ2,2)+mZ2*pow(2.4952,2)))*1./CA*(sum_pdf);
		if(isnan(result)){ff[0] = 0;}
		if(q2 < mZ2){ff[0] = 0;}
        else{ff[0] = result;}
  }

  else if(order == "test"){
	  phiMP = M_PI/2;
	  CMP = 0.;
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
      double result =  2.*imag(1./(2*M_PI)*Njac*pow(1.5,-Nint)*-1.*pow(5.,2.+Nint)/(2.+Nint));
      if(isnan(result)){ff[0] = 0;}
      if(xx[0]==0){ff[0]=0;}
      else{ff[0] = result;}
      }
  else if (order == "NLO"){
    if (process == "qqbar"){

      if (power == "reg"){
        double z =  tau+(1.-tau)*xx[0];
      	double jacobian = (1.-tau)*(1.-tau/z);
      	double x = tau/z+xx[1]*(1.-tau/z);
        double result = 0;
        if(!fitted_pdfs) result = DY_LO_factor()*(DY_NLO_qqbar_reg(z)*jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z+DY_NLO_qqbar_plus(z)*(jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z-pow(1.-tau,2)*real(pdf_sum_qqbar_charge_weighted(tau+xx[1]*(1.-tau),tau))));
        if(fitted_pdfs) result = DY_LO_factor()*(DY_NLO_qqbar_reg(z)*jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z+DY_NLO_qqbar_plus(z)*(jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z-pow(1.-tau,2)*real(fit_sum_qqbar_charge_weighted(tau+xx[1]*(1.-tau),tau))));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
      }
      else if (power == "delta"){
        double jacobian = (1.-tau);
        double x = tau+xx[0]*(1.-tau);
        double result = 0;
        if(!fitted_pdfs) result = DY_LO_factor()*DY_NLO_qqbar_delta()*real(pdf_sum_qqbar_charge_weighted(x,tau));
        if(fitted_pdfs) result = DY_LO_factor()*DY_NLO_qqbar_delta()*real(fit_sum_qqbar_charge_weighted(x,tau));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "LP"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
      	double result =  0;
        if(!fitted_pdfs) result = DY_LO_factor()*(DY_NLO_qqbar_plus(z))*(jacobian*pdf_sum_qqbar_charge_weighted(x,tau/z)/z - pow(1.-tau,2)*pdf_sum_qqbar_charge_weighted(tau+xx[1]*(1.-tau),tau));
        if(fitted_pdfs) result = DY_LO_factor()*(DY_NLO_qqbar_plus(z))*(jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z - pow(1.-tau,2)*real(fit_sum_qqbar_charge_weighted(tau+xx[1]*(1.-tau),tau)));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "LP_corr"){
        double jacobian = tau*(1.-tau);
        double z = tau*xx[0];
        double x = tau+xx[1]*(1.-tau);
      	double result =  0;
        if(!fitted_pdfs) result = DY_LO_factor()*(DY_NLO_qqbar_plus(z))*jacobian*(-pdf_sum_qqbar_charge_weighted(x,tau));
        if(fitted_pdfs) result = DY_LO_factor()*(DY_NLO_qqbar_plus(z))*jacobian*(-real(fit_sum_qqbar_charge_weighted(x,tau)));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "exp"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
      	double result =  0;
        for(int i = 0; i < NCOMP; i++){
          if(!fitted_pdfs) result = DY_LO_factor()*(DY_NLO_qqbar_expansion(z, i+1))*jacobian*(pdf_sum_qqbar_charge_weighted(x,tau/z)/z);
          if(fitted_pdfs) result = DY_LO_factor()*(DY_NLO_qqbar_expansion(z, i+1))*jacobian*(real(fit_sum_qqbar_charge_weighted(x,tau/z)/z));
          if(isnan(result)){ff[i] = 0;}
          else{ff[i] = result;}
          }
        }
      else{
        cout << process  << " " << order << " " << power  << endl;
        cout << "Wrong power specified" << endl;
        exit(0);
        }
      }
    else if (process == "qg"){
      if (power == "full"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result =  0;
        if(!fitted_pdfs) result = DY_LO_factor()*(DY_NLO_qg_full(z))*jacobian*(pdf_sum_qg_charge_weighted(x,tau/z)/z);
        if(fitted_pdfs) result = DY_LO_factor()*(DY_NLO_qg_full(z))*jacobian*(real(fit_sum_qg_charge_weighted(x,tau/z)/z));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "exp"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result =  0;
        for(int i = 0; i < NCOMP; i++){
          if(!fitted_pdfs) result = DY_LO_factor()*(DY_NLO_qg_expansion(z, i+1))*jacobian*(pdf_sum_qg_charge_weighted(x,tau/z)/z);
          if(fitted_pdfs) result = DY_LO_factor()*(DY_NLO_qg_expansion(z, i+1))*jacobian*(real(fit_sum_qg_charge_weighted(x,tau/z)/z));
          if(isnan(result)){ff[i] = 0;}
          else{ff[i] = result;}
          }
        }
      else{
        cout << process << " " << order << " " << power << endl;
        cout << "Wrong power specified" << endl;
        exit(0);
        }
      }
    else{
      cout << process << " " << order << " " << power << endl;
      cout << "Wrong process specified" << endl;
      exit(0);
      }
  }
  else if (order == "NNLO"){
    if (process == "qqbar"){
      if (power == "reg"){
        double z =  tau+(1.-tau)*xx[0];
        double jacobian = (1.-tau)*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result = 0;
        if(!fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_qqbar_NS(z)*jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z
      										+ DY_NNLO_BB_full(z)*jacobian*real(pdf_sum_qqbar_charge_unweighted(x,tau/z))/z
      										+2.*DY_NNLO_BC_full(z)*jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z
      										+DY_NNLO_qqbar_plus(z)*(jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z-pow(1.-tau,2)*real(pdf_sum_qqbar_charge_weighted(tau+xx[1]*(1.-tau),tau)))
      									);
        if(fitted_pdfs) result = result = DY_LO_factor()*(DY_NNLO_qqbar_NS(z)*jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z
      										+ DY_NNLO_BB_full(z)*jacobian*real(fit_sum_qqbar_charge_unweighted(x,tau/z))/z
      										+2.*DY_NNLO_BC_full(z)*jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z
      										+DY_NNLO_qqbar_plus(z)*(jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z-pow(1.-tau,2)*real(fit_sum_qqbar_charge_weighted(tau+xx[1]*(1.-tau),tau)))
      									);
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
      }
      else if (power == "delta"){
        double jacobian = (1.-tau);
        double x = tau+xx[0]*(1.-tau);
        double result = 0;
        if(!fitted_pdfs) result = jacobian*DY_LO_factor()*DY_NNLO_qqbar_delta()*real(pdf_sum_qqbar_charge_weighted(x,tau));
        if(fitted_pdfs) result = jacobian*DY_LO_factor()*DY_NNLO_qqbar_delta()*real(fit_sum_qqbar_charge_weighted(x,tau));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "LP"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result =  0;
        if(!fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_qqbar_plus(z))*(jacobian*pdf_sum_qqbar_charge_weighted(x,tau/z)/z - pow(1.-tau,2)*pdf_sum_qqbar_charge_weighted(tau+xx[1]*(1.-tau),tau));
        if(fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_qqbar_plus(z))*(jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z - pow(1.-tau,2)*real(fit_sum_qqbar_charge_weighted(tau+xx[1]*(1.-tau),tau)));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "LP_corr"){
        double jacobian = tau*(1.-tau);
        double z = tau*xx[0];
        double x = tau+xx[1]*(1.-tau);
        double result =  0;
        if(!fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_qqbar_plus(z))*jacobian*(-pdf_sum_qqbar_charge_weighted(x,tau));
        if(fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_qqbar_plus(z))*jacobian*(-real(fit_sum_qqbar_charge_weighted(x,tau)));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "exp"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result =  0;
        for(int i = 0; i < NCOMP; i++){
          if(!fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_qqbar_NS_expansion(z,i+1)*jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z
        										+ DY_NNLO_BB_expansion(z,i+1)*jacobian*real(pdf_sum_qqbar_charge_unweighted(x,tau/z))/z
        										+2.*DY_NNLO_BC_expansion(z,i+1)*jacobian*real(pdf_sum_qqbar_charge_weighted(x,tau/z))/z
        									);
          if(fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_qqbar_NS_expansion(z,i+1)*jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z
        										+ DY_NNLO_BB_expansion(z,i+1)*jacobian*real(fit_sum_qqbar_charge_unweighted(x,tau/z))/z
        										+2.*DY_NNLO_BC_expansion(z,i+1)*jacobian*real(fit_sum_qqbar_charge_weighted(x,tau/z))/z
        									);
          if(isnan(result)){ff[i] = 0;}
          else{ff[i] = result;}
          }
        }
      else{
        cout << process  << " " << order << " " << power  << endl;
        cout << "Wrong power specified" << endl;
        exit(0);
        }
      }
    else if (process == "qg"){
      if (power == "full"){
          double jacobian = (1.-tau);
          double z = tau+(1.-tau)*xx[0];
          jacobian = jacobian*(1.-tau/z);
          double x = tau/z+xx[1]*(1.-tau/z);
          double result =  0;
          if(!fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_qg_full(z))*jacobian*(pdf_sum_qg_charge_weighted(x,tau/z)/z);
          if(fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_qg_full(z))*jacobian*(real(fit_sum_qg_charge_weighted(x,tau/z)/z));
          if(isnan(result)){ff[0] = 0;}
          else{ff[0] = result;}
          }
      else if (power == "exp"){
          double jacobian = (1.-tau);
          double z = tau+(1.-tau)*xx[0];
          jacobian = jacobian*(1.-tau/z);
          double x = tau/z+xx[1]*(1.-tau/z);
          double result =  0;
          for(int i = 0; i < NCOMP; i++){
            if(!fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_qg_expansion(z, i+1))*jacobian*(pdf_sum_qg_charge_weighted(x,tau/z)/z);
            if(fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_qg_expansion(z, i+1))*jacobian*(real(fit_sum_qg_charge_weighted(x,tau/z)/z));
            if(isnan(result)){ff[i] = 0;}
            else{ff[i] = result;}
            }
          }
      else{
        cout << process << " " << order << " " << power << endl;
        cout << "Wrong power specified" << endl;
        exit(0);
        }
      }
    else if (process == "gg"){
      if (power == "full"){
          double jacobian = (1.-tau);
          double z = tau+(1.-tau)*xx[0];
          jacobian = jacobian*(1.-tau/z);
          double x = tau/z+xx[1]*(1.-tau/z);
          double result =  0;
          if(!fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_gg_full(z))*jacobian*(pdf_sum_gg_charge_weighted(x,tau/z)/z);
          if(fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_gg_full(z))*jacobian*(real(fit_sum_gg_charge_weighted(x,tau/z)/z));
          if(isnan(result)){ff[0] = 0;}
          else{ff[0] = result;}
          }
      else if (power == "exp"){
          double jacobian = (1.-tau);
          double z = tau+(1.-tau)*xx[0];
          jacobian = jacobian*(1.-tau/z);
          double x = tau/z+xx[1]*(1.-tau/z);
          double result =  0;
          for(int i = 0; i < NCOMP; i++){
            if(!fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_gg_expansion(z, i+1))*jacobian*(pdf_sum_gg_charge_weighted(x,tau/z)/z);
            if(fitted_pdfs) result = DY_LO_factor()*(DY_NNLO_gg_expansion(z, i+1))*jacobian*(real(fit_sum_gg_charge_weighted(x,tau/z)/z));
            if(isnan(result)){ff[i] = 0;}
            else{ff[i] = result;}
            }
          }
      else{
        cout << process << " " << order << " " << power << endl;
        cout << "Wrong power specified" << endl;
        exit(0);
        }
      }
    else if (process == "qq"){
      if (power == "full"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result =  0;

        if(!fitted_pdfs) result = DY_LO_factor()*jacobian/z*(DY_NNLO_CC_full(z)*real(pdf_sum_qq_charge_weighted_double(x,tau/z))
      												+ DY_NNLO_CD_full(z)*real(pdf_sum_qq_charge_weighted_double_vivj(x,tau/z))
      												+ DY_NNLO_CE_full(z)*real(pdf_sum_qq_charge_weighted_single(x,tau/z))
      												+ DY_NNLO_CF_full(z)*real(pdf_sum_qq_charge_weighted_single_vivi(x,tau/z))
                            ); //factor of 1/2 is taken care of in DeltaCF!
        if(fitted_pdfs) result = DY_LO_factor()*jacobian/z*(DY_NNLO_CC_full(z)*real(fit_sum_qq_charge_weighted_double(x,tau/z))
      												+ DY_NNLO_CD_full(z)*real(fit_sum_qq_charge_weighted_double_vivj(x,tau/z))
      												+ DY_NNLO_CE_full(z)*real(fit_sum_qq_charge_weighted_single(x,tau/z))
      												+ DY_NNLO_CF_full(z)*real(fit_sum_qq_charge_weighted_single_vivi(x,tau/z))
      												);
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "exp"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result =  0;
        for(int i = 0; i < NCOMP; i++){
            if(!fitted_pdfs) result = DY_LO_factor()*jacobian/z*(DY_NNLO_CC_expansion(z, i+1)*real(pdf_sum_qq_charge_weighted_double(x,tau/z))
          												+ DY_NNLO_CD_expansion(z, i+1)*real(pdf_sum_qq_charge_weighted_double_vivj(x,tau/z))
          												+ DY_NNLO_CE_expansion(z, i+1)*real(pdf_sum_qq_charge_weighted_single(x,tau/z))
          												+ DY_NNLO_CF_expansion(z, i+1)*real(pdf_sum_qq_charge_weighted_single_vivi(x,tau/z))
          												);
            if(fitted_pdfs) result = DY_LO_factor()*jacobian/z*(DY_NNLO_CC_expansion(z, i+1)*real(fit_sum_qq_charge_weighted_double(x,tau/z))
          												+ DY_NNLO_CD_expansion(z, i+1)*real(fit_sum_qq_charge_weighted_double_vivj(x,tau/z))
          												+ DY_NNLO_CE_expansion(z, i+1)*real(fit_sum_qq_charge_weighted_single(x,tau/z))
          												+ DY_NNLO_CF_expansion(z, i+1)*real(fit_sum_qq_charge_weighted_single_vivi(x,tau/z))
          												);
            if(isnan(result)){ff[i] = 0;}
            else{ff[i] = result;}
            }
        }
        else{
          cout << process << " " << order << " " << power << endl;
          cout << "Wrong power specified" << endl;
          exit(0);
          }
        }
    else{
        cout << process << " " << order << " " << power << endl;
        cout << "Wrong process specified" << endl;
        exit(0);
      }
  }
  else if (order == "resum"){
    if(power == "test"){
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
      double result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*DY_LO_factor()*fit_mellin_pdf_sum_qqbar_charge_weighted(Nint-1.));
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
      }
    else if(power == "full"){
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(phiMP*I);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(phiMP*I);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
    	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*DY_LO_factor()*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda)))*fit_mellin_pdf_sum_qqbar_charge_weighted(Nint-1.));
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
      }
    else if(power == "NLO"){
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(phiMP*I);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(phiMP*I);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
      double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*DY_LO_factor()*NLOmatch_DY(Nint)*fit_mellin_pdf_sum_qqbar_charge_weighted(Nint-1.));
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
      }
    else if(power == "NNLO"){
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(phiMP*I);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(phiMP*I);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
      double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*DY_LO_factor()*NNLOmatch_DY(Nint)*fit_mellin_pdf_sum_qqbar_charge_weighted(Nint-1.));
  	  if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
      }
    }
  else{
    cout << process << " " << order << " " << power << endl;
    cout << "Wrong order specified" << endl;
    exit(0);
  }
  return 0;
}

// for higgs
static int higgs_integrand(const int *ndim, const cubareal xx[],const int *ncomp, cubareal ff[], void *userdata) {

  if(order == "LO"){
      if(!fitted_pdfs) ff[0] = higgs_LO_factor()*(pdf_sum_gg(tau+(1.-tau)*xx[0],tau))*(1.-tau);
      if(fitted_pdfs) ff[0] = higgs_LO_factor()*(real(fit_sum_gg(tau+(1.-tau)*xx[0],tau)))*(1.-tau);
  }
  else if (order == "NLO"){
    if (process == "gg"){

      if (power == "reg"){
        double z =  tau+(1.-tau)*xx[0];
      	double jacobian = (1.-tau)*(1.-tau/z);
      	double x = tau/z+xx[1]*(1.-tau/z);
        double result = 0;
        if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_gg_reg(z)*jacobian*real(pdf_sum_gg(x,tau/z))/z+higgs_NLO_gg_plus(z)*(jacobian*real(pdf_sum_gg(x,tau/z))/z-pow(1.-tau,2)*real(pdf_sum_gg(tau+xx[1]*(1.-tau),tau))));
        if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_gg_reg(z)*jacobian*real(fit_sum_gg(x,tau/z))/z+higgs_NLO_gg_plus(z)*(jacobian*real(fit_sum_gg(x,tau/z))/z-pow(1.-tau,2)*real(fit_sum_gg(tau+xx[1]*(1.-tau),tau))));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
      }
      else if (power == "delta"){
        double jacobian = (1.-tau);
        double x = tau+xx[0]*(1.-tau);
        double result = 0;
        if(!fitted_pdfs) result = higgs_LO_factor()*higgs_NLO_gg_delta()*real(pdf_sum_gg(x,tau));
        if(fitted_pdfs) result = higgs_LO_factor()*higgs_NLO_gg_delta()*real(fit_sum_gg(x,tau));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "LP"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
      	double result =  0;
        if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_gg_plus(z))*(jacobian*pdf_sum_gg(x,tau/z)/z - pow(1.-tau,2)*pdf_sum_gg(tau+xx[1]*(1.-tau),tau));
        if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_gg_plus(z))*(jacobian*real(fit_sum_gg(x,tau/z))/z - pow(1.-tau,2)*real(fit_sum_gg(tau+xx[1]*(1.-tau),tau)));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "LP_corr"){
        double jacobian = tau*(1.-tau);
        double z = tau*xx[0];
        double x = tau+xx[1]*(1.-tau);
      	double result =  0;
        if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_gg_plus(z))*jacobian*(-pdf_sum_gg(x,tau));
        if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_gg_plus(z))*jacobian*(-real(fit_sum_gg(x,tau)));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "exp"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
      	double result =  0;
        for(int i = 0; i < NCOMP; i++){
          if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_gg_expansion(z, i+1))*jacobian*(pdf_sum_gg(x,tau/z)/z);
          if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_gg_expansion(z, i+1))*jacobian*(real(fit_sum_gg(x,tau/z)/z));
          if(isnan(result)){ff[i] = 0;}
          else{ff[i] = result;}
          }
        }
      else{
        cout << process  << " " << order << " " << power  << endl;
        cout << "Wrong power specified" << endl;
        exit(0);
        }
      }
    else if (process == "qg"){
      if (power == "full"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result =  0;
        if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_qg_full(z))*jacobian*(pdf_sum_qg(x,tau/z)/z);
        if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_qg_full(z))*jacobian*(real(fit_sum_qg(x,tau/z)/z));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "exp"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result =  0;
        for(int i = 0; i < NCOMP; i++){
          if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_qg_expansion(z, i+1))*jacobian*(pdf_sum_qg(x,tau/z)/z);
          if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_qg_expansion(z, i+1))*jacobian*(real(fit_sum_qg(x,tau/z)/z));
          if(isnan(result)){ff[i] = 0;}
          else{ff[i] = result;}
          }
        }
      else{
        cout << process << " " << order << " " << power << endl;
        cout << "Wrong power specified" << endl;
        exit(0);
        }
      }
    else if (process == "qqbar"){
      if (power == "full"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result =  0;
        if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_qqbar_full(z))*jacobian*(pdf_sum_qqbar(x,tau/z)/z);
        if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_qqbar_full(z))*jacobian*(real(fit_sum_qqbar(x,tau/z)/z));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "exp"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result =  0;
        for(int i = 0; i < NCOMP; i++){
          if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_qqbar_expansion(z, i+1))*jacobian*(pdf_sum_qqbar(x,tau/z)/z);
          if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NLO_qqbar_expansion(z, i+1))*jacobian*(real(fit_sum_qqbar(x,tau/z)/z));
          if(isnan(result)){ff[i] = 0;}
          else{ff[i] = result;}
          }
        }
      else{
        cout << process << " " << order << " " << power << endl;
        cout << "Wrong power specified" << endl;
        exit(0);
        }
      }
    else{
      cout << process << " " << order << " " << power << endl;
      cout << "Wrong process specified" << endl;
      exit(0);
      }
  }
  else if (order == "NNLO"){
    if (process == "gg"){
      if (power == "reg"){
        double z =  tau+(1.-tau)*xx[0];
        double jacobian = (1.-tau)*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result = 0;
      	if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_gg_reg(z)*jacobian*real(pdf_sum_gg(x,tau/z))/z+higgs_NNLO_gg_plus(z)*(jacobian*real(pdf_sum_gg(x,tau/z))/z-pow(1.-tau,2)*real(pdf_sum_gg(tau+xx[1]*(1.-tau),tau))));
        if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_gg_reg(z)*jacobian*real(fit_sum_gg(x,tau/z))/z+higgs_NNLO_gg_plus(z)*(jacobian*real(fit_sum_gg(x,tau/z))/z-pow(1.-tau,2)*real(fit_sum_gg(tau+xx[1]*(1.-tau),tau))));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
      }
      else if (power == "delta"){
        double jacobian = (1.-tau);
        double x = tau+xx[0]*(1.-tau);
        double result = 0;
        if(!fitted_pdfs) result = higgs_LO_factor()*higgs_NNLO_gg_delta()*real(pdf_sum_gg(x,tau));
        if(fitted_pdfs) result = higgs_LO_factor()*higgs_NNLO_gg_delta()*real(fit_sum_gg(x,tau));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "LP"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result =  0;
        if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_gg_plus(z))*(jacobian*pdf_sum_gg(x,tau/z)/z - pow(1.-tau,2)*pdf_sum_gg(tau+xx[1]*(1.-tau),tau));
        if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_gg_plus(z))*(jacobian*real(fit_sum_gg(x,tau/z))/z - pow(1.-tau,2)*real(fit_sum_gg(tau+xx[1]*(1.-tau),tau)));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "LP_corr"){
        double jacobian = tau*(1.-tau);
        double z = tau*xx[0];
        double x = tau+xx[1]*(1.-tau);
        double result =  0;
        if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_gg_plus(z))*jacobian*(-pdf_sum_gg(x,tau));
        if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_gg_plus(z))*jacobian*(-real(fit_sum_gg(x,tau)));
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
      else if (power == "exp"){
        double jacobian = (1.-tau);
        double z = tau+(1.-tau)*xx[0];
        jacobian = jacobian*(1.-tau/z);
        double x = tau/z+xx[1]*(1.-tau/z);
        double result =  0;
        for(int i = 0; i < NCOMP; i++){
          if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_gg_expansion(z, i+1))*jacobian*(pdf_sum_gg(x,tau/z)/z);
          if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_gg_expansion(z, i+1))*jacobian*(real(fit_sum_gg(x,tau/z)/z));
          if(isnan(result)){ff[i] = 0;}
          else{ff[i] = result;}
          }
        }
      else{
        cout << process  << " " << order << " " << power  << endl;
        cout << "Wrong power specified" << endl;
        exit(0);
        }
      }
    else if (process == "qg"){
        if (power == "full"){
          double jacobian = (1.-tau);
          double z = tau+(1.-tau)*xx[0];
          jacobian = jacobian*(1.-tau/z);
          double x = tau/z+xx[1]*(1.-tau/z);
          double result =  0;
          if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qg_reg(z))*jacobian*(pdf_sum_qg(x,tau/z)/z);
          if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qg_reg(z))*jacobian*(real(fit_sum_qg(x,tau/z)/z));
          if(isnan(result)){ff[0] = 0;}
          else{ff[0] = result;}
          }
        else if (power == "exp"){
          double jacobian = (1.-tau);
          double z = tau+(1.-tau)*xx[0];
          jacobian = jacobian*(1.-tau/z);
          double x = tau/z+xx[1]*(1.-tau/z);
          double result =  0;
          for(int i = 0; i < NCOMP; i++){
            if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qg_expansion(z, i+1))*jacobian*(pdf_sum_qg(x,tau/z)/z);
            if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qg_expansion(z, i+1))*jacobian*(real(fit_sum_qg(x,tau/z)/z));
            if(isnan(result)){ff[i] = 0;}
            else{ff[i] = result;}
            }
          }
        else{
          cout << process << " " << order << " " << power << endl;
          cout << "Wrong power specified" << endl;
          exit(0);
          }
        }
    else if (process == "qq"){
        if (power == "full"){
          double jacobian = (1.-tau);
          double z = tau+(1.-tau)*xx[0];
          jacobian = jacobian*(1.-tau/z);
          double x = tau/z+xx[1]*(1.-tau/z);
          double result =  0;
          if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qq_reg(z))*jacobian*(pdf_sum_qq(x,tau/z)/z);
          if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qq_reg(z))*jacobian*(real(fit_sum_qq(x,tau/z)/z));
          if(isnan(result)){ff[0] = 0;}
          else{ff[0] = result;}
          }
        else if (power == "exp"){
          double jacobian = (1.-tau);
          double z = tau+(1.-tau)*xx[0];
          jacobian = jacobian*(1.-tau/z);
          double x = tau/z+xx[1]*(1.-tau/z);
          double result =  0;
          for(int i = 0; i < NCOMP; i++){
            if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qq_expansion(z, i+1))*jacobian*(pdf_sum_qq(x,tau/z)/z);
            if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qq_expansion(z, i+1))*jacobian*(real(fit_sum_qq(x,tau/z)/z));
            if(isnan(result)){ff[i] = 0;}
            else{ff[i] = result;}
            }
          }
        else{
          cout << process << " " << order << " " << power << endl;
          cout << "Wrong power specified" << endl;
          exit(0);
          }
        }
    else if (process == "qqp"){
        if (power == "full"){
          double jacobian = (1.-tau);
          double z = tau+(1.-tau)*xx[0];
          jacobian = jacobian*(1.-tau/z);
          double x = tau/z+xx[1]*(1.-tau/z);
          double result =  0;
          if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qqp_reg(z))*jacobian*(pdf_sum_qqNI(x,tau/z)/z);
          if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qqp_reg(z))*jacobian*(real(fit_sum_qqNI(x,tau/z)/z));
          if(isnan(result)){ff[0] = 0;}
          else{ff[0] = result;}
          }
        else if (power == "exp"){
          double jacobian = (1.-tau);
          double z = tau+(1.-tau)*xx[0];
          jacobian = jacobian*(1.-tau/z);
          double x = tau/z+xx[1]*(1.-tau/z);
          double result =  0;
          for(int i = 0; i < NCOMP; i++){
            if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qqp_expansion(z, i+1))*jacobian*(pdf_sum_qqNI(x,tau/z)/z);
            if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qqp_expansion(z, i+1))*jacobian*(real(fit_sum_qqNI(x,tau/z)/z));
            if(isnan(result)){ff[i] = 0;}
            else{ff[i] = result;}
            }
          }
        else{
          cout << process << " " << order << " " << power << endl;
          cout << "Wrong power specified" << endl;
          exit(0);
          }
        }
    else if (process == "qqbar"){
        if (power == "full"){
          double jacobian = (1.-tau);
          double z = tau+(1.-tau)*xx[0];
          jacobian = jacobian*(1.-tau/z);
          double x = tau/z+xx[1]*(1.-tau/z);
          double result =  0;
          if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qqbar_reg(z))*jacobian*(pdf_sum_qqbar(x,tau/z)/z);
          if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qqbar_reg(z))*jacobian*(real(fit_sum_qqbar(x,tau/z)/z));
          if(isnan(result)){ff[0] = 0;}
          else{ff[0] = result;}
          }
        else if (power == "exp"){
          double jacobian = (1.-tau);
          double z = tau+(1.-tau)*xx[0];
          jacobian = jacobian*(1.-tau/z);
          double x = tau/z+xx[1]*(1.-tau/z);
          double result =  0;
          for(int i = 0; i < NCOMP; i++){
            if(!fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qqbar_expansion(z, i+1))*jacobian*(pdf_sum_qqbar(x,tau/z)/z);
            if(fitted_pdfs) result = higgs_LO_factor()*(higgs_NNLO_qqbar_expansion(z, i+1))*jacobian*(real(fit_sum_qqbar(x,tau/z)/z));
            if(isnan(result)){ff[i] = 0;}
            else{ff[i] = result;}
            }
          }
        else{
          cout << process << " " << order << " " << power << endl;
          cout << "Wrong power specified" << endl;
          exit(0);
          }
        }
    else{
        cout << process << " " << order << " " << power << endl;
        cout << "Wrong process specified" << endl;
        exit(0);
      }
  }
  else if (order == "resum"){
    if(power == "test"){
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
      double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*higgs_LO_factor()/tau*fit_mellin_pdf_sum_gg(Nint));
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
      }
    else if(power == "full"){
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(phiMP*I);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(phiMP*I);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
      double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*higgs_LO_factor()*exp(ISNNLL*alphas_muR*wideangle(D2higgs,lambda)+2.*(1./alphas_muR*ISLL*h0g(lambda)+ISNLL*h1g(lambda)+ISNNLL*alphas_muR*h2g(lambda)+ISNLP*h0gNLP(Nint,lambda)))*fit_mellin_pdf_sum_gg(Nint-1.));
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
      }
    else if(power == "NLO"){
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(phiMP*I);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(phiMP*I);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
      double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*higgs_LO_factor()*NLOmatch_higgs(Nint)*fit_mellin_pdf_sum_gg(Nint-1.));
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
      }
    else if(power == "NNLO"){
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(phiMP*I);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(phiMP*I);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
      double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*higgs_LO_factor()*NNLOmatch_higgs(Nint)*fit_mellin_pdf_sum_gg(Nint-1.));
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
      }
    }
  else{
    cout << process << " " << order << " " << power << endl;
    cout << "Wrong order specified" << endl;
    exit(0);
  }
  return 0;
}

//for di-higgs
static int dihiggs_integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {
  if(process == "diff"){
    if(order == "SM"){
    double jacobian = (1.-tau)*2;
    double x = tau+(1.-tau)*xx[0];
    double ctheta = -1.+2.*xx[1];
    double result = 0;
    if(!fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_LO_factor_SM(tau*S2,ctheta)*(pdf_sum_gg(x,tau)/S2);
    if(fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_LO_factor_SM(tau*S2,ctheta)*real(fit_sum_gg(x,tau)/S2);
    if(tau*S2 < 4.*mH2){ff[0] = 0;}
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }

  else if(order == "SUSY_hh"){
    double jacobian = (1.-tau)*2;
    double x = tau+(1.-tau)*xx[0];
    double ctheta = -1.+2.*xx[1];
    double result = 0;
    if(!fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_hh(tau*S2, ctheta)*(pdf_sum_gg(x,tau)/S2);
    if(fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_hh(tau*S2, ctheta)*real(fit_sum_gg(x,tau)/S2);
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }
  else if(order == "SUSY_hH"){
    double jacobian = (1.-tau)*2;
    double x = tau+(1.-tau)*xx[0];
    double ctheta = -1.+2.*xx[1];
    double result = 0;
    if(!fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_hH(tau*S2, ctheta)*(pdf_sum_gg(x,tau)/S2);
    if(fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_hH(tau*S2, ctheta)*real(fit_sum_gg(x,tau)/S2);
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }else if(order == "SUSY_HH"){
    double jacobian = (1.-tau)*2;
    double x = tau+(1.-tau)*xx[0];
    double ctheta = -1.+2.*xx[1];
    double result = 0;
    if(!fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_HH(tau*S2, ctheta)*(pdf_sum_gg(x,tau)/S2);
    if(fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_HH(tau*S2, ctheta)*real(fit_sum_gg(x,tau)/S2);
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }else if(order == "SUSY_Ah"){
    double jacobian = (1.-tau)*2;
    double x = tau+(1.-tau)*xx[0];
    double ctheta = -1.+2.*xx[1];
    double result = 0;
    if(!fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_Ah(tau*S2, ctheta)*(pdf_sum_gg(x,tau)/S2);
    if(fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_Ah(tau*S2, ctheta)*real(fit_sum_gg(x,tau)/S2);
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }else if(order == "SUSY_AH"){
    double jacobian = (1.-tau)*2;
    double x = tau+(1.-tau)*xx[0];
    double ctheta = -1.+2.*xx[1];
    double result = 0;
    if(!fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_AH(tau*S2, ctheta)*(pdf_sum_gg(x,tau)/S2);
    if(fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_AH(tau*S2, ctheta)*real(fit_sum_gg(x,tau)/S2);
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }else if(order == "SUSY_AA"){
    double jacobian = (1.-tau)*2;
    double x = tau+(1.-tau)*xx[0];
    double ctheta = -1.+2.*xx[1];
    double result = 0;
    if(!fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_AA(tau*S2, ctheta)*(pdf_sum_gg(x,tau)/S2);
    if(fitted_pdfs) result = 2.*jacobian*sqrt(tau*S2)*dihiggs_AA(tau*S2, ctheta)*real(fit_sum_gg(x,tau)/S2);
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }
}
  else if(process== "full"){
  if(order == "SM"){
    double z = tau+(1.-tau)*xx[0];
  	double x = z+(1.-z)*xx[1];
    double jacobian = (1.-z)*(1.-tau)*2;
    double ctheta = -1.+2.*xx[2];
    double result = 0;
    if(!fitted_pdfs) result = jacobian*dihiggs_LO_factor_SM(z*S2,ctheta)*(pdf_sum_gg(x,z));
    if(fitted_pdfs) result = jacobian*dihiggs_LO_factor_SM(z*S2,ctheta)*real(fit_sum_gg(x,z));
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }
  else if(order == "SUSY_hh"){
    double z = tau+(1.-tau)*xx[0];
  	double x = z+(1.-z)*xx[1];
    double jacobian = (1.-z)*(1.-tau)*2;
    double ctheta = -1.+2.*xx[2];
    double result = 0;
    if(!fitted_pdfs) result = jacobian*dihiggs_hh(z*S2,ctheta)*(pdf_sum_gg(x,z));
    if(fitted_pdfs) result = jacobian*dihiggs_hh(z*S2,ctheta)*real(fit_sum_gg(x,z));
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }
  else if(order == "SUSY_hH"){
    double z = tau+(1.-tau)*xx[0];
  	double x = z+(1.-z)*xx[1];
    double jacobian = (1.-z)*(1.-tau)*2;
    double ctheta = -1.+2.*xx[2];
    double result = 0;
    if(!fitted_pdfs) result = jacobian*dihiggs_hH(z*S2,ctheta)*(pdf_sum_gg(x,z));
    if(fitted_pdfs) result = jacobian*dihiggs_hH(z*S2,ctheta)*real(fit_sum_gg(x,z));
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }
  else if(order == "SUSY_HH"){
    double z = tau+(1.-tau)*xx[0];
  	double x = z+(1.-z)*xx[1];
    double jacobian = (1.-z)*(1.-tau)*2;
    double ctheta = -1.+2.*xx[2];
    double result = 0;
    if(!fitted_pdfs) result = jacobian*dihiggs_HH(z*S2,ctheta)*(pdf_sum_gg(x,z));
    if(fitted_pdfs) result = jacobian*dihiggs_HH(z*S2,ctheta)*real(fit_sum_gg(x,z));
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }
  else if(order == "SUSY_Ah"){
    double z = tau+(1.-tau)*xx[0];
  	double x = z+(1.-z)*xx[1];
    double jacobian = (1.-z)*(1.-tau)*2;
    double ctheta = -1.+2.*xx[2];
    double result = 0;
    if(!fitted_pdfs) result = jacobian*dihiggs_Ah(z*S2,ctheta)*(pdf_sum_gg(x,z));
    if(fitted_pdfs) result = jacobian*dihiggs_Ah(z*S2,ctheta)*real(fit_sum_gg(x,z));
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }
  else if(order == "SUSY_AH"){
    double z = tau+(1.-tau)*xx[0];
  	double x = z+(1.-z)*xx[1];
    double jacobian = (1.-z)*(1.-tau)*2;
    double ctheta = -1.+2.*xx[2];
    double result = 0;
    if(!fitted_pdfs) result = jacobian*dihiggs_AH(z*S2,ctheta)*(pdf_sum_gg(x,z));
    if(fitted_pdfs) result = jacobian*dihiggs_AH(z*S2,ctheta)*real(fit_sum_gg(x,z));
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }
  else if(order == "SUSY_AA"){
    double z = tau+(1.-tau)*xx[0];
  	double x = z+(1.-z)*xx[1];
    double jacobian = (1.-z)*(1.-tau)*2;
    double ctheta = -1.+2.*xx[2];
    double result = 0;
    if(!fitted_pdfs) result = jacobian*dihiggs_AA(z*S2,ctheta)*(pdf_sum_gg(x,z));
    if(fitted_pdfs) result = jacobian*dihiggs_AA(z*S2,ctheta)*real(fit_sum_gg(x,z));
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }
  else if(order == "LOapprox"){
    double z = tau+(1.-tau)*xx[0];
  	double x = z+(1.-z)*xx[1];
    double jacobian = (1.-z)*(1.-tau)*2;
    double ctheta = -1.+2.*xx[2];
    double result = 0;
    if(!fitted_pdfs) result = jacobian*dihiggs_LO_factor_approx(z*S2,ctheta)*(pdf_sum_gg(x,z));
    if(fitted_pdfs) result = jacobian*dihiggs_LO_factor_approx(z*S2,ctheta)*real(fit_sum_gg(x,z));
    if(isnan(result)){ff[0] = 0;}
    else{ff[0] = result;}
  }
}
  else if (process == "resumdiff"){
    if(order == "SMinverse"){
      double jacobian = 2.;
      double ctheta = -1.+2.*xx[1];
      double result = 0;
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
      result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*jacobian*sqrt(tau*S2)*dihiggs_LO_factor_SM(tau*S2,ctheta)/tau*fit_mellin_pdf_sum_gg(Nint)/S2);
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
      }
    else if(order == "LOinversefull"){ //does not work
        CMP = 2.4;
        phiMP = 1./2.*M_PI;
        complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
        complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
          // z = 4mH2/s => s = 4mH2/z;
        double scale2 = 4.*mH2/xx[1];
        // extra factpr of 2 is jacobian

        double ctheta = -1.+2.*xx[2];
        double result = 2.*2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*pow(xx[1],Nint-1.)*dihiggs_LO_factor_approx(scale2,ctheta)*fit_mellin_pdf_sum_gg(Nint));
        if(scale2 >= S2){ff[0]=0;}
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
      } //does not work
    else if(order == "fullapprox"){
        CMP = 2.1;
        phiMP = 1./2.*M_PI;
        complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
        complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
        double scale2 = 4.*mH2/xx[1];
        complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
        double ctheta = -1.+2.*xx[2];
        double result = 2.*2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*pow(xx[1],Nint-1.)*dihiggs_LO_factor_approx(scale2,ctheta)*exp(ISNNLL*alphas_muR*wideangle(D2higgs,lambda)+2.*(1./alphas_muR*ISLL*h0g(lambda)+ISNLL*h1g(lambda)+ISNNLL*alphas_muR*h2g(lambda)+ISNLP*h0gNLP(Nint,lambda)))*fit_mellin_pdf_sum_gg(Nint));
        if(scale2 >= S2){ff[0]=0;}
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        } //does not work
    else if(order == "full"){
        CMP = 2.1;
        phiMP = 1./2.*M_PI;
        complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
        complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
              // z = 4mH2/s => s = 4mH2/z;
        double scale2 = 4.*mH2/xx[1];
            // extra factpr of 2 is jacobian
        complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
        double ctheta = -1.+2.*xx[2];
        double result = 2.*2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*pow(xx[1],Nint-1.)*dihiggs_LO_factor_SM(scale2,ctheta)*exp(ISNNLL*alphas_muR*wideangle(D2higgs,lambda)+2.*(1./alphas_muR*ISLL*h0g(lambda)+ISNLL*h1g(lambda)+ISNNLL*alphas_muR*h2g(lambda)+ISNLP*h0gNLP(Nint,lambda)))*fit_mellin_pdf_sum_gg(Nint));
        if(scale2 >= S2){ff[0]=0;}
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
    else if(order == "SM"){
          double jacobian = 2.;
          double ctheta = -1.+2.*xx[1];
          double result = 0;
          complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
          complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
          complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
          result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*jacobian*sqrt(tau*S2)*dihiggs_LO_factor_SM(tau*S2,ctheta)/tau*exp(ISNNLL*alphas_muR*wideangle(D2higgs,lambda)+2.*(1./alphas_muR*ISLL*h0g(lambda)+ISNLL*h1g(lambda)+ISNNLL*alphas_muR*h2g(lambda)+ISNLP*h0gNLP(Nint,lambda)))*fit_mellin_pdf_sum_gg(Nint)/S2);
          if(tau*S2 < 4.*mH2){ff[0] = 0;}
          if(isnan(result)){ff[0] = 0;}
          else{ff[0] = result;}
          }
    else if(order == "SUSY_hh"){
      double jacobian = 2.;
      double ctheta = -1.+2.*xx[1];
      double result = 0;
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
      result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*jacobian*sqrt(tau*S2)*dihiggs_hh(tau*S2,ctheta)/tau*exp(ISNNLL*alphas_muR*wideangle(D2higgs,lambda)+2.*(1./alphas_muR*ISLL*h0g(lambda)+ISNLL*h1g(lambda)+ISNNLL*alphas_muR*h2g(lambda)+ISNLP*h0gNLP(Nint,lambda)))*fit_mellin_pdf_sum_gg(Nint)/S2);
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
    }
    else if(order == "SUSY_hH"){
      double jacobian = 2.;
      double ctheta = -1.+2.*xx[1];
      double result = 0;
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
      result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*jacobian*sqrt(tau*S2)*dihiggs_hH(tau*S2,ctheta)/tau*exp(ISNNLL*alphas_muR*wideangle(D2higgs,lambda)+2.*(1./alphas_muR*ISLL*h0g(lambda)+ISNLL*h1g(lambda)+ISNNLL*alphas_muR*h2g(lambda)+ISNLP*h0gNLP(Nint,lambda)))*fit_mellin_pdf_sum_gg(Nint)/S2);
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
    }
    else if(order == "SUSY_HH"){
      double jacobian = 2.;
      double ctheta = -1.+2.*xx[1];
      double result = 0;
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
      result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*jacobian*sqrt(tau*S2)*dihiggs_HH(tau*S2,ctheta)/tau*exp(ISNNLL*alphas_muR*wideangle(D2higgs,lambda)+2.*(1./alphas_muR*ISLL*h0g(lambda)+ISNLL*h1g(lambda)+ISNNLL*alphas_muR*h2g(lambda)+ISNLP*h0gNLP(Nint,lambda)))*fit_mellin_pdf_sum_gg(Nint)/S2);
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
    }
    else if(order == "SUSY_Ah"){
      double jacobian = 2.;
      double ctheta = -1.+2.*xx[1];
      double result = 0;
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
      result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*jacobian*sqrt(tau*S2)*dihiggs_Ah(tau*S2,ctheta)/tau*exp(ISNNLL*alphas_muR*wideangle(D2higgs,lambda)+2.*(1./alphas_muR*ISLL*h0g(lambda)+ISNLL*h1g(lambda)+ISNNLL*alphas_muR*h2g(lambda)+ISNLP*h0gNLP(Nint,lambda)))*fit_mellin_pdf_sum_gg(Nint)/S2);
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
    }
    else if(order == "SUSY_AH"){
      double jacobian = 2.;
      double ctheta = -1.+2.*xx[1];
      double result = 0;
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
      result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*jacobian*sqrt(tau*S2)*dihiggs_AH(tau*S2,ctheta)/tau*exp(ISNNLL*alphas_muR*wideangle(D2higgs,lambda)+2.*(1./alphas_muR*ISLL*h0g(lambda)+ISNLL*h1g(lambda)+ISNNLL*alphas_muR*h2g(lambda)+ISNLP*h0gNLP(Nint,lambda)))*fit_mellin_pdf_sum_gg(Nint)/S2);
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
    }
    else if(order == "SUSY_AA"){
      double jacobian = 2.;
      double ctheta = -1.+2.*xx[1];
      double result = 0;
      complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
      complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
      result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*jacobian*sqrt(tau*S2)*dihiggs_AA(tau*S2,ctheta)/tau*exp(ISNNLL*alphas_muR*wideangle(D2higgs,lambda)+2.*(1./alphas_muR*ISLL*h0g(lambda)+ISNLL*h1g(lambda)+ISNNLL*alphas_muR*h2g(lambda)+ISNLP*h0gNLP(Nint,lambda)))*fit_mellin_pdf_sum_gg(Nint)/S2);
      if(isnan(result)){ff[0] = 0;}
      else{ff[0] = result;}
    }
    }
  else{
    cout << process << " " << order << " " << power << endl;
    cout << "Wrong order specified" << endl;
    exit(0);
  }
  return 0;
}

//for diboson
static int diboson_integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {
  if(process == "WW"){
	  if(order == "LOdiff"){
			double jacobian = (1.-tau);
			double x = tau+(1.-tau)*xx[0];
			double result = 0;
			if(!fitted_pdfs){
					result = 2.*jacobian*sqrt(tau*S2)/S2*(partonic_down_wpwm(tau*S2)*(pdf_sum_qqbarDOWN(x,tau)));
					result += 2.*jacobian*sqrt(tau*S2)/S2*(partonic_up_wpwm(tau*S2)*(pdf_sum_qqbarUP(x,tau)));
					}
			if(fitted_pdfs){
					result = 2.*jacobian*sqrt(tau*S2)/S2*(partonic_down_wpwm(tau*S2)*real(fit_sum_qqbarDOWN(x,tau)));
					result += 2.*jacobian*sqrt(tau*S2)/S2*(partonic_up_wpwm(tau*S2)*real(fit_sum_qqbarUP(x,tau)));
					}
      if(tau*S2 < 4.*mW2){ff[0] =0;}
			if(isnan(result)){ff[0] = 0;}
			else{ff[0] = result;}
		}
	  else if(order == "LOfull"){

		double z = tau+(1.-tau)*xx[0];
		double x = z+(1.-z)*xx[1];
		double jacobian = (1.-z)*(1.-tau);
		double result = 0;
		if(!fitted_pdfs){  result = jacobian*(partonic_down_wpwm(z*S2)*(pdf_sum_qqbarDOWN(x,z)));
						   result += jacobian*(partonic_up_wpwm(z*S2)*(pdf_sum_qqbarUP(x,z)));
						}
		if(fitted_pdfs) {  result = jacobian*(partonic_down_wpwm(z*S2)*real(fit_sum_qqbarDOWN(x,z)));
						   result += jacobian*(partonic_up_wpwm(z*S2)*real(fit_sum_qqbarUP(x,z)));
						}
		if(isnan(result)){ff[0] = 0;}
		else{ff[0] = result;}
	  }
	  else if(order == "LOinversefullQ"){
        double result = 0;
		complex<double> Nint = 1.+CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
		complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
		double z = tau+(1.-tau)*xx[1];
		double jacz = (1.-tau);
		result = 2.*imag(1./(2*M_PI)*jacz*Njac*pow(z,-Nint)*partonic_down_wpwm(z*S2)*fit_mellin_pdf_sum_qqbarDOWN(Nint-1.));
		result += 2.*imag(1./(2*M_PI)*jacz*Njac*pow(z,-Nint)*partonic_up_wpwm(z*S2)*fit_mellin_pdf_sum_qqbarUP(Nint-1.));
		if(isnan(result)){ff[0] = 0;}
		else{ff[0] = result;}

        }
	  else if(order == "LOinversefullM"){
        complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
        complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
        double xvar = -1.+xx[1];
	    double wr = xvar/(1.+xvar);
	    double wjac = 1./(pow(1.+xvar,2));
	    complex<double> z = exp(wr/Nint);
        complex<double> scale2 = 4.*mW2/z;
	    double result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac*exp(wr)/Nint*cpartonic_down_wpwm(scale2)*fit_mellin_pdf_sum_qqbarDOWN(Nint));
        result += 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac*exp(wr)/Nint*cpartonic_up_wpwm(scale2)*fit_mellin_pdf_sum_qqbarUP(Nint));
        if(abs(scale2) >= S2){ff[0]=0;}
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
    else if(order == "LOinversediff"){ // works
		  double result = 0;
		  complex<double> Nint = 1.+CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
		  complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
		  result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*sqrt(tau*S2)*partonic_down_wpwm(tau*S2)*fit_mellin_pdf_sum_qqbarDOWN(Nint-1.)/S2);
		  result += 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*sqrt(tau*S2)*partonic_up_wpwm(tau*S2)*fit_mellin_pdf_sum_qqbarUP(Nint-1.)/S2);
      if(tau*S2 < 4.*mW2){ff[0] =0;}
      if(isnan(result)){ff[0] = 0;}
		  else{ff[0] = result;}
        }
    else if(order == "resumfullQ"){ // check if the wide angle is correct
		  double result = 0;
		  complex<double> Nint = 1.+CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
		  complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
          complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
		  double z = tau+(1.-tau)*xx[1];
		  double jacz = (1.-tau);
		  result = 2.*imag(1./(2*M_PI)*jacz*Njac*pow(z,-Nint)*partonic_down_wpwm(z*S2)*fit_mellin_pdf_sum_qqbarDOWN(Nint-1.)*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda))));
		  result += 2.*imag(1./(2*M_PI)*jacz*Njac*pow(z,-Nint)*partonic_up_wpwm(z*S2)*fit_mellin_pdf_sum_qqbarUP(Nint-1.)*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda))));
		  if(isnan(result)){ff[0] = 0;}
		  else{ff[0] = result;}
        }
    else if(order == "resumfullM"){ // check if the wide angle is correct
			complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
			complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
			  // z = 4mH2/s => s = 4mH2/z;
			double xvar = -1.+xx[1];
			double wr = xvar/(1.+xvar);
			double wjac = 1./(pow(1.+xvar,2));
			complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
		    complex<double> z = exp(wr/Nint);
			complex<double> scale2 = 4.*mW2/z;
			double result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac*exp(wr)/Nint*cpartonic_down_wpwm(scale2)*fit_mellin_pdf_sum_qqbarDOWN(Nint)*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda))));
			result += 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac*exp(wr)/Nint*cpartonic_up_wpwm(scale2)*fit_mellin_pdf_sum_qqbarUP(Nint)*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda))));
			if(abs(scale2) >= S2){ff[0]=0;}
			if(isnan(result)){ff[0] = 0;}
			else{ff[0] = result;}

        }
    else if(order == "resumdiff"){
		  double result = 0;
		  complex<double> Nint = 1.+CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
		  complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
      complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
		  result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*sqrt(tau*S2)*partonic_down_wpwm(tau*S2)*fit_mellin_pdf_sum_qqbarDOWN(Nint-1.)/S2*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda))));
		  result += 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*sqrt(tau*S2)*partonic_up_wpwm(tau*S2)*fit_mellin_pdf_sum_qqbarUP(Nint-1.)/S2*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda))));
      if(tau*S2 < 4.*mW2){ff[0] =0;}
      if(isnan(result)){ff[0] = 0;}
		  else{ff[0] = result;}
        }
	}
  else if(process == "ZZ"){
	  if(order == "LOdiff"){
			double jacobian = (1.-tau);
			double x = tau+(1.-tau)*xx[0];
			double result = 0;
			if(!fitted_pdfs){
					result = 2.*jacobian*sqrt(tau*S2)/S2*(partonic_down_zz(tau*S2)*(pdf_sum_qqbarDOWN(x,tau)));
					result += 2.*jacobian*sqrt(tau*S2)/S2*(partonic_up_zz(tau*S2)*(pdf_sum_qqbarUP(x,tau)));
					}
			if(fitted_pdfs){
					result = 2.*jacobian*sqrt(tau*S2)/S2*(partonic_down_zz(tau*S2)*real(fit_sum_qqbarDOWN(x,tau)));
					result += 2.*jacobian*sqrt(tau*S2)/S2*(partonic_up_zz(tau*S2)*real(fit_sum_qqbarUP(x,tau)));
					}
			if(isnan(result)){ff[0] = 0;}
			else{ff[0] = result;}
		}
	  else if(order == "LOfull"){

		double z = tau+(1.-tau)*xx[0];
		double x = z+(1.-z)*xx[1];
		double jacobian = (1.-z)*(1.-tau);
		double result = 0;
		if(!fitted_pdfs){  result = jacobian*(partonic_down_zz(z*S2)*(pdf_sum_qqbarDOWN(x,z)));
						   result += jacobian*(partonic_up_zz(z*S2)*(pdf_sum_qqbarUP(x,z)));
						}
		if(fitted_pdfs) {  result = jacobian*(partonic_down_zz(z*S2)*real(fit_sum_qqbarDOWN(x,z)));
						   result += jacobian*(partonic_up_zz(z*S2)*real(fit_sum_qqbarUP(x,z)));
						}
		if(isnan(result)){ff[0] = 0;}
		else{ff[0] = result;}
	  }
	  else if(order == "LOinversefullQ"){
        double result = 0;
		complex<double> Nint = 1.+CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
		complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
		double z = tau+(1.-tau)*xx[1];
		double jacz = (1.-tau);
		result = 2.*imag(1./(2*M_PI)*jacz*Njac*pow(z,-Nint)*partonic_down_zz(z*S2)*fit_mellin_pdf_sum_qqbarDOWN(Nint-1.));
		result += 2.*imag(1./(2*M_PI)*jacz*Njac*pow(z,-Nint)*partonic_up_zz(z*S2)*fit_mellin_pdf_sum_qqbarUP(Nint-1.));
		if(isnan(result)){ff[0] = 0;}
		else{ff[0] = result;}

        }
	  else if(order == "LOinversefullM"){
        complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
        complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
          // z = 4mH2/s => s = 4mH2/z;
        double xvar = -1.+xx[1];
	    double wr = xvar/(1.+xvar);
	    double wjac = 1./(pow(1.+xvar,2));
	    complex<double> z = exp(wr/Nint);
        complex<double> scale2 = 4.*mZ2/z;
	    double result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac*exp(wr)/Nint*cpartonic_down_zz(scale2)*fit_mellin_pdf_sum_qqbarDOWN(Nint));
        result += 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac*exp(wr)/Nint*cpartonic_up_zz(scale2)*fit_mellin_pdf_sum_qqbarUP(Nint));
        if(abs(scale2) >= S2){ff[0]=0;}
        if(isnan(result)){ff[0] = 0;}
        else{ff[0] = result;}
        }
    else if(order == "LOinversediff"){ // works
		  double result = 0;
		  complex<double> Nint = 1.+CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
		  complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
		  result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*sqrt(tau*S2)*partonic_down_zz(tau*S2)*fit_mellin_pdf_sum_qqbarDOWN(Nint-1.)/S2);
		  result += 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*sqrt(tau*S2)*partonic_up_zz(tau*S2)*fit_mellin_pdf_sum_qqbarUP(Nint-1.)/S2);
      if(tau*S2 < 4.*mZ2){ff[0] =0;}
      if(isnan(result)){ff[0] = 0;}
		  else{ff[0] = result;}
        }
    else if(order == "resumfullQ"){ // check if the wide angle is correct
		  double result = 0;
		  complex<double> Nint = 1.+CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
		  complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
          complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
		  double z = tau+(1.-tau)*xx[1];
		  double jacz = (1.-tau);
		  result = 2.*imag(1./(2*M_PI)*jacz*Njac*pow(z,-Nint)*partonic_down_zz(z*S2)*fit_mellin_pdf_sum_qqbarDOWN(Nint-1.)*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda))));
		  result += 2.*imag(1./(2*M_PI)*jacz*Njac*pow(z,-Nint)*partonic_up_zz(z*S2)*fit_mellin_pdf_sum_qqbarUP(Nint-1.)*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda))));
		  if(isnan(result)){ff[0] = 0;}
		  else{ff[0] = result;}
        }
    else if(order == "resumfullM"){ // check if the wide angle is correct
			complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
			complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
			  // z = 4mH2/s => s = 4mH2/z;
			double xvar = -1.+xx[1];
			double wr = xvar/(1.+xvar);
			double wjac = 1./(pow(1.+xvar,2));
			complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
		    complex<double> z = exp(wr/Nint);
			complex<double> scale2 = 4.*mZ2/z;
			double result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac*exp(wr)/Nint*cpartonic_down_zz(scale2)*fit_mellin_pdf_sum_qqbarDOWN(Nint)*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda))));
			result += 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac*exp(wr)/Nint*cpartonic_up_zz(scale2)*fit_mellin_pdf_sum_qqbarUP(Nint)*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda))));
			if(abs(scale2) >= S2){ff[0]=0;}
			if(isnan(result)){ff[0] = 0;}
			else{ff[0] = result;}

        }
    else if(order == "resumdiff"){
		  double result = 0;
		  complex<double> Nint = 1.+CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
		  complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
          complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
		  result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*sqrt(tau*S2)*partonic_down_zz(tau*S2)*fit_mellin_pdf_sum_qqbarDOWN(Nint-1.)/S2*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda))));
		  result += 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*2.*sqrt(tau*S2)*partonic_up_zz(tau*S2)*fit_mellin_pdf_sum_qqbarUP(Nint-1.)/S2*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda))));
      if(tau*S2 < 4.*mZ2){ff[0] =0;}
      if(isnan(result)){ff[0] = 0;}
		  else{ff[0] = result;}
        }
  }
  else{
    cout << process << " " << order << " " << power << endl;
    cout << "Wrong order specified" << endl;
    exit(0);
  }
  return 0;
}


//for test
static int test_integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {
  if(process == "defor"){// note that xx[1] needs to go from -1 to 0!
	  complex<double> Nint = CMP+xx[0]/(1.-xx[0])*I;
	  complex<double> Njac = 1./pow(1.-xx[0],2)*I;
	  double xvar = -1.+xx[1];
	  double wr = xvar/(1.+xvar);
	  double wjac = 1./(pow(1.+xvar,2));
	  complex<double> x = exp(wr/Nint);
	  double result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac*exp(wr)/Nint*x*(1.-x));
	  if(isnan(result)){ff[0] = 0;}
	  else{ff[0] = result;}
  }
  else if(process == "deriv"){
	  complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
	  complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
	  double result = -2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*pow(xx[1],Nint)/Nint*(1.-2.*xx[1]));
	  if(isnan(result)){ff[0] = 0;}
	  else{ff[0] = result;}
  }
  else if(process == "normal"){
	  complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
	  complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);
	  double result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*pow(xx[1],Nint)*(1.-xx[1]));
	  if(isnan(result)){ff[0] = 0;}
	  else{ff[0] = result;}
	}
  else if(process == "nspace"){
	  complex<double> Nint = CMP+xx[0]/(1.-xx[0])*exp(I*phiMP);
	  complex<double> Njac = 1./pow(1.-xx[0],2)*exp(I*phiMP);

	  double result= 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*1./((Nint*Nint+3.*Nint+2.)));
	  if(isnan(result)){ff[0] = 0;}
	  else{ff[0] = result;}
	}
  else{
    cout << process << " " << order << " " << power << endl;
    cout << "Wrong order specified" << endl;
    exit(0);
  }
  return 0;
}


vector<results_c> call_cuhre_dy(std::string orde, std::string chan, std::string pow, bool fitted, int maxpower, int verbose){
	int comp, nregions, neval, fail;
	NDIM = 2;
  order = orde;
  process = chan;
  power = pow;
  NCOMP = maxpower;
  fitted_pdfs=fitted;
  cout << "orde=" << orde << " chan=" <<chan << " pow=" << pow << " NCOMP " << NCOMP << endl;
	cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
	Cuhre(NDIM, NCOMP, dy_integrand, USERDATA, NVEC,
    EPSREL, EPSABS, verbose | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
	vector<results_c> result;
  for(comp = 0; comp < NCOMP; comp++){result.push_back({integral[comp], error[comp], prob[comp]});};
	return result;
}



vector<results_c> call_cuhre_higgs(std::string orde, std::string chan, std::string pow, bool fitted, int maxpower, int verbose){
	int comp, nregions, neval, fail;
	NDIM = 2;
  order = orde;
  process = chan;
  power = pow;
  NCOMP = maxpower;
  fitted_pdfs=fitted;
	cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
	Cuhre(NDIM, NCOMP, higgs_integrand, USERDATA, NVEC,
    EPSREL, EPSABS, verbose | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
	vector<results_c> result;
  for(comp = 0; comp < NCOMP+1; comp++){result.push_back({integral[comp], error[comp], prob[comp]});};
	return result;
}



vector<results_c> call_cuhre_dihiggs(std::string orde, std::string chan, bool fitted, int verbose){
	int comp, nregions, neval, fail;
	NDIM = 3;
	  order = orde; //indicate SM, SUSY_hh etc
	  process = chan; //indicate diff, LO or resum
	  NCOMP = 1;
	  fitted_pdfs=fitted;
	cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
	Cuhre(NDIM, NCOMP, dihiggs_integrand, USERDATA, NVEC,
    EPSREL, EPSABS, verbose | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
	vector<results_c> result;
  for(comp = 0; comp < NCOMP+1; comp++){result.push_back({integral[comp], error[comp], prob[comp]});};
	return result;
}


vector<results_c> call_cuhre_diboson(std::string orde, std::string chan, bool fitted, int verbose){
	int comp, nregions, neval, fail;
	NDIM = 3;
  order = orde;
  process = chan;
  NCOMP = 1;
  fitted_pdfs=fitted;
	cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
	Cuhre(NDIM, NCOMP, diboson_integrand, USERDATA, NVEC,
    EPSREL, EPSABS, verbose | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
	vector<results_c> result;
  for(comp = 0; comp < NCOMP+1; comp++){result.push_back({integral[comp], error[comp], prob[comp]});};
	return result;
}



vector<results_c> call_cuhre_test(std::string orde, std::string chan, bool fitted, int verbose){
	int comp, nregions, neval, fail;
	NDIM = 3;
    order = orde;
    process = chan;
    NCOMP = 1;
    fitted_pdfs=fitted;
	cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
	Cuhre(NDIM, NCOMP, test_integrand, USERDATA, NVEC,
    EPSREL, EPSABS, verbose | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
	vector<results_c> result;
    for(comp = 0; comp < NCOMP+1; comp++){result.push_back({integral[comp], error[comp], prob[comp]});};
	return result;
}
