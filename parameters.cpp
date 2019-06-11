#include <iostream> //for stirng
#include <fstream>
#include <sstream>
#include <unistd.h> //for getting options
#include <gsl/gsl_math.h>
#include <stdlib.h> // for exit()
#include <getopt.h> //for getting options
#include <vector>
#include "LHAPDF/LHAPDF.h"

using namespace std;

//////////////////////////////////////////////////////////////////
///
/// parameters.cpp contains the functions to read in the parameters
/// that are used in the code and unchanged
/// contains stuff like invariant mass, com energy, tau, color 
/// factors and the electromagnetic coupling constant
///
//////////////////////////////////////////////////////////////////

double Q(500.);
double Q2(pow(Q,2.));
double S(13000.);
double S2(pow(S,2.));
double tau(Q2/S2);
double muF(Q);
double muR(Q);
double muF2(muF*muF);
double muR2(muR*muR);
double CF(4./3.);
double CA(3.);
double TF(1./2.);
double alphaEM(1./132.507);
double GF(1.166379e-5);
double alphas_muF(0);
double alphas_muR(0);
double alphas_Q(0);
double nF(5);
double zeta2(1.6449340668482264);
double zeta3(1.2020569031595942);
double pbunits(0.38937966*pow(10.,9.));
string setname("PDF4LHC15_nnlo_100");
std::vector<LHAPDF::PDF*> pdfs;
double xmin_pdfs(0.), xmax_pdfs(0.); //min x, max x and alphas


///////////////////////////////////////
/// usage of the read argument function
///////////////////////////////////////
void print_usage() {
  cout << endl
       << "Usage:" << endl
       << " -p   <set>     specify PDFs being used (default \"PDF4LHC15_nnlo_100\")" << endl
       << "                installed: NNPDF31_nnlo_as_0118, CT14nnlo, PDF4LHC15_nnlo_100" << endl
       << "                           NNPDF31_nlo_as_0118, CT14nlo, PDF4LHC15_nlo_100 " << endl
	   << "							  MSTW2008nnlo90cl, MMHT2014nlo68cl, GRVPI0" << endl
       << " -Q   <Q>      set transferred momentum Q [GeV]" << endl
       << " -r   <muR>         set renormalization scale [GeV]" << endl
       << " -f   <muF>         set factorization scale [GeV]" << endl
       << " -S   <com> 	    set C.O.M. energy [TeV]" << endl
       << endl;
  exit(0);
  return;
}

///////////////////////////////////////
/// defaults of the programm
///////////////////////////////////////
void print_defaults(){
cout << "The values are set as follows: " << endl
	<< "PDFset: 				" << setname << endl
	<< "Momentum [GeV]: 			" << Q << endl
	<< "Renormalization scale [GeV]: 		" << muR << endl
	<< "Factorization scale [GeV]:		" << muF << endl
	<< "Center of Mass energy [TeV]: 		" << S/1000. << endl
	<< "First set Q if muR!=muF!=Q:"  << endl;
}


///////////////////////////////////////
/// read arguments of command line
///////////////////////////////////////

void read_arguments(int argc, char* argv[]) {
  const char* const short_options = "hp:f:r:Q:S:"; //: is needed as there is a variable after the -p, -S etc
  int next_option;
  ostringstream sset;
  while ((next_option = getopt(argc, argv, short_options)) != -1)
  { switch (next_option) {
    case 'h':
      print_usage();
    case 'p':
      sset << optarg;
      setname = sset.str();
      break;
    case 'f':
      muF = strtod(optarg, NULL);
      muF2 = muF*muF;
      break;
    case 'r':
      muR = strtod(optarg, NULL);
      muR2 = muR*muR;
      break;
    case 'Q':
      Q = strtod(optarg, NULL);
      Q2 = pow(Q,2);
      muF = Q;
      muR = Q; 
      muF2 = muF*muF;
      muR2 = muR*muR;
      break;
    case 'S':
       S = strtod(optarg, NULL)*1000.;
       S2 = pow(S,2);
      break;
    case -1: break;
    case ':': abort();
    }
  }
  tau = Q2/S2; //reset tau
  return;
}

