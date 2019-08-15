#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_gamma.h>
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
vector<double> coeff_without_log_min{-0.68938870688796914238692004929684,0.21441761459841800098124390853181,0.0023166683910024487036808719881751,0.000054733182590597519438754734161676,1.7745288538231597645987933696917e-6,6.8437045907106193566996452577016e-8,2.947060744092889090738342640066e-9,1.3699793323770418612719644976808e-10,6.7373903781187947913604791558209e-12,3.459730584201987021974909266821e-13,1.8385357724131998544633785938459e-14,1.0045988359647687312715120064408e-15,5.6174207252965252608445164258392e-17,3.2027951145461002617937971690034e-18,1.8566924200310092367418945019891e-19,1.0919217233917323797174597235267e-20,6.5026881927623498871751017577074e-22,3.9155737492442172164934104388177e-23,2.3810030155046108776745485269144e-24,1.4606085438163192870420809235499e-25,9.0309435773785635793682538973758e-27};
vector<double> coeff_without_log{0.016016180449195828736706919847563,0.50364244007530129181209541016961,0.016150992430500258887454469519295,0.0012440242104244936126561052441311,0.00013757218124461673921971996409271,0.000018563818526037733164867951831292,2.841735345154404159345057900388e-6,4.7459967905789937221638951390084e-7,8.4480385437812006760918194736415e-8,1.5787671240434005432468704745206e-8,3.0657620713990379812888900401962e-9,6.1407921728125845808062188948123e-10,1.2618830243156719690872483510663e-10,2.6493141798196099571267833100346e-11,5.6647085463642592615881193424544e-12,1.2304111577958111751746743425556e-12,2.7093457836786768143959997910272e-13,6.0380264633837012791971169067897e-14,1.3600089939957496823523631172345e-14,3.0924474063185687585526105660708e-15,7.0917249609207158220307479772166e-16,1.6388083639226002471134216524917e-16,3.8134643501689946130066172536341e-17,8.9301073961181165576802169067821e-18,2.1033134159935941622136823325517e-18,4.9802988416537865766736420915489e-19};
vector<double> coeff_with_log{0.78663768384453670033888320677278,-0.41620952725593784667953842242806,-0.00034228344812004476202652239655774,0.00051983478594322839878400723976638,0.000081863197797631253790804177909117,0.000011446588102007482416312605876182,1.5892062886999847941775628695067e-6,2.2338423252224515068231007415833e-7,3.1891887838925967071719982063812e-8,4.6214546811175686710007173549982e-9,6.7869586794086462849346346283628e-10,1.0084751549203312137600111811472e-10,1.5139491991787295463756565311691e-11,2.2933132330421112748609261090515e-12,3.5015004118789907486820109950181e-13,5.3837930345481399119494586411224e-14,8.3297277101403199084766307019966e-15,1.295974110879937883479179368404e-15,2.02647336035796002185391428377e-16,3.1831312832495903503796521239796e-17,5.0205993519207247124651461342787e-18,7.9485113709665495616862163483289e-19,1.2627194849864269987683678605277e-19,2.0123203632470790939829582792006e-20,3.216243943999222314415843339999e-21,5.1542680037470152358191187964105e-22};

vector<double> coeff_s12_with_log{0.9385394862381089,0.2683653003435577,-0.005072463662441,0.00023987125391778635,-0.00001657773235414073,1.4104850004471674e-6,-1.3681696988221392e-7,1.4522249304080988e-8,-1.6460488323922186e-9,1.9613017167800043e-10,-2.4305026460713415e-11,3.1088309014127186e-12,-4.081483566621797e-13,5.476708301362134e-14,-7.48645184976289e-15,1.0398244951753513e-15,-1.4644180329607417e-16};
vector<double> coeff_s12_without_log1{0.034296726151657783361374008695853,0.017504968761145416714930074438557,0.035378271506331131823948097105364,0.0061764940588794180312724154968855,0.0011236340536774410116394874204879,0.00021709883386103752822378794837679,0.000043993484179921723358142112342373,9.2464559757059895682022284137595e-6,1.9994029939890986498758167338896e-6,4.4221989210810433667110595297177e-7,9.9619406400580314632018274759815e-8,2.2784411669105567190827818360172e-8,5.2779501098238905073537446151731e-9,1.2359581244924201622798706024144e-9,2.9214610462529939471612096933494e-10,6.9618951124252817929018579857696e-11,1.6709174462086701331047031786919e-11,4.0357650125532259749441847387231e-12,9.8026035243514964973993495418354e-13,2.3930358725244763009084009143736e-13,5.8685912034624695689183297920517e-14,1.4451412028329030387130003813529e-14,3.5720556278617025235785976016295e-15,8.8596785384604692094778377947014e-16,2.2043830440802058007210049995229e-16,5.5006945399993942077919884420482e-17};
vector<double> coeff_s12_without_log2{0.09624881942560648691872123300995,-0.05179476600063473873028904405946,0.0023408173375369381855706922015049,0.00013404708785011296890145212630355,7.1253925300384430673629607489391e-6,3.8888347751213459456336084666111e-7,2.1907454262409345624139479564837e-8,1.2683162437034123863591145754445e-9,7.5099312688131694231522842728167e-11,4.5298972640520130753860525636915e-12,2.7747662262692536993345278636863e-13,1.72182143550158569572003823856e-14,1.0802800301114477364713842648252e-15,6.8423177791080219817491879745696e-17,4.369674956948432011556695628731e-18,2.8108171009124870437535744690189e-19};

double Li2(double x){
	return gsl_sf_dilog(x);
}

double Li3(double x){
	if(x < -1 || x > 1){cout << "OOB Li3" << endl; exit(0);}
	if(x > -0.5 && x < 0.5){return clenshaw(coeff_without_log, 2.*x);}
	if(x >= 0.5){return clenshaw(coeff_with_log, 4.*(1.-x)-1.)-1./2.*log(1.-x)*pow(log(x),2);}
	if(x <= -0.5){return clenshaw(coeff_without_log_min, 4.*x+3.);}
}


double S12(double x){
	if(x < -1  || x > 1){cout << "OOB S12" << endl; exit(0);}
	if(x <= -0.5){return clenshaw(coeff_s12_without_log2, 4.*x+3.);}
	if(x>-0.5 && x <= 0.5){return clenshaw(coeff_s12_without_log1, 2.*x);}
	if(x > 0.5){return clenshaw(coeff_s12_with_log, 4.*x-3.) + 1./2.*log(x)*pow(log(1.-x),2) + gsl_sf_dilog(1.-x)*log(1.-x);}
}
////////////////////////////////////////////////////////
/// implementation of clenshaw algorithm 
/// sum of the chebychevs is quickly done this way
///////////////////////////////////////////////////////
double clenshaw(vector<double> coeff, double x){
	double bkpp = 0.,bkp = 0.,bk = 0.;
	for(int i = coeff.size()-1; i > 0 ; i--){
		bkpp = bkp;
		bkp = bk;
		bk = coeff[i]+2.*x*bkp-bkpp;
	}
	
	return coeff[0]+x*bk-bkp;
	
}

complex<double> Gamma(complex<double> x){
	gsl_sf_result lnr;
	gsl_sf_result arg;
	complex<double> I(0,1);
	gsl_sf_lngamma_complex_e(real(x),imag(x),&lnr, &arg);
	return exp(lnr.val+I*arg.val);
}
