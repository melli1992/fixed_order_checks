#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <vector>
#include "polygamma.h"
#include <iostream>

//////////////////////////////////////////////////////
///
/// Implementation Li3 - M. van Beekveld
/// and S_{1,2} (Nielsen generalized polylog)
///
//////////////////////////////////////////////////////
/// contains implementation of the Li3(x) function
/// uses expansion in Chebychev polynomials 
/// since there is an anoying log(x) around x = 0 
/// for the Li3(1.-x) taylor expansion, we need to split
/// the transformation up. First we transform to 
/// the Chebychev domain, which is from -1 to 1 so 
/// we need y = 2x-1. But the problem is with x = 0
/// so in this domain lets add the troublemaker which
/// is 1/2*log(x)*log(1-x)^2 and then do the Chebychev
/// fit. Then of course set an upper bound, which we pick
/// to be 1/2. So y = 4x-1. For the other one we take
/// from 1/2 to 1 (x), so y = 4x-3. to guarantee the
/// Chebychev domain. Then we make a fit to the functions
/// in each domain and evaluate it, this is in Li3.
/// the precision is around 3E-17 in the domain x[0,1]
/// in mathematica code polylog_numerical_cpp_implementation.nb are the fits. 
//////////////////////////////////////////////////////

using namespace std;
vector<double> coeff_with_log{0.7866376838445367,-0.41620952725593785,-0.00034228344812004474,0.0005198347859432284,0.00008186319779763126,0.000011446588102007482,1.5892062886999848e-6,2.2338423252224516e-7,3.189188783892597e-8,4.621454681117568e-9,6.786958679408646e-10,1.0084751549203312e-10,1.5139491991787295e-11,2.293313233042111e-12,3.501500411878991e-13,5.38379303454814e-14,8.32972771014032e-15,1.295974110879938e-15};
vector<double> coeff_without_log{0.26351741692148545,-0.2683653003435577,0.005072463662441,-0.00023987125391778635,0.00001657773235414073,-1.4104850004471674e-6,1.3681696988221392e-7,-1.4522249304080988e-8,1.6460488323922186e-9,-1.9613017167800043e-10,2.4305026460713415e-11,-3.1088309014127186e-12,4.081483566621797e-13,-5.476708301362134e-14,7.48645184976289e-15,-1.0398244951753513e-15};

vector<double> coeff_s12_with_log{0.9385394862381089,0.2683653003435577,-0.005072463662441,0.00023987125391778635,-0.00001657773235414073,1.4104850004471674e-6,-1.3681696988221392e-7,1.4522249304080988e-8,-1.6460488323922186e-9,1.9613017167800043e-10,-2.4305026460713415e-11,3.1088309014127186e-12,-4.081483566621797e-13,5.476708301362134e-14,-7.48645184976289e-15,1.0398244951753513e-15,-1.4644180329607417e-16};
vector<double> coeff_s12_without_log{0.03290390726967671,0.045835603887253114,0.014296807139225137,0.001519409672618663,0.00017302133861581628,0.000021118341141779037,2.715324482155722e-6,3.63009111307922e-7,5.00065688981209e-8,7.053481051513469e-9,1.0140718221176674e-9,1.4810107288084014e-10,2.191595620940357e-11,3.2795413058193226e-12,4.95489600432486e-13,7.548780560145603e-14,1.1584919786276463e-14,1.7894276899206954e-15,2.7799171469275844e-16,4.340961820471936e-17};


double Li2(double x){
	return gsl_sf_dilog(x);
}

double Li3(double x){
	x = 1.-x;
	if(x <= 0.5){return clenshaw(coeff_with_log, 4.*x-1.)-1./2.*log(x)*pow(log(1.-x),2);}
	if(x > 0.5){return clenshaw(coeff_without_log,4.*x-3.);}
}


double S12(double x){
	if(x <= 0.5){return clenshaw(coeff_s12_without_log, 4.*x-1.);}
	if(x > 0.5){return clenshaw(coeff_s12_with_log, 4.*x-3.) + 1./2.*log(x)*pow(log(1.-x),2) + gsl_sf_dilog(1.-x)*log(1.-x);}
}
////////////////////////////////////////////////////////
/// implementation of clenshaw algorithm 
/// sum of the chebychevs is quickly done this way
///////////////////////////////////////////////////////
double clenshaw(vector<double> coeff, double x){
	double alpha = 2.*x, beta = -1;
	double bkpp = 0,bkp = 0,bk = 0;
	for(int i =coeff.size(); i > 0 ; i--){
		bkpp = bkp;
		bkp = bk;
		bk = coeff[i]+2.*x*bkp-bkpp;
	}
	return coeff[0]+x*bk-bkp;
	
}
