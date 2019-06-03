#include <iostream>
#include <cmath>
#include <string>
#include "DY-num-v01.h"
#include "LHAPDF/LHAPDF.h"

using namespace std;
using namespace LHAPDF;

#define PI 3.141592653589793
#define E 2.718281828459045
#define Nc 3.

/////////////////////////////////////////////////////////////////////////////////////

double S;           // center of mass energy
double Q,Qfix,tau;  // invariant mass 

double mur; // renormalization scale
double muf; // factorization scale

int member;

double Minf,Msup; // integration limits

std::string setname;
size_t nmem;
std::vector<LHAPDF::PDF*> pdfs;
std::vector<int> pids;

/////////////////////////////////////////////////////////////////////////////////////

int main(){

   //centre of mass energy
   S = pow(14000.,2);

   // fix one value of invariant mass
   Qfix = 500.;
   tau = Qfix*Qfix/S;
  
   //which PDFmemeber to use
   member = 0;

   // define some tring that may be needed
   std::string nlo = "NLO";
   std::string nnlo = "NNLO";

   std::string Vegas = "vegas";
   std::string Suave = "suave";
   std::string Cuhre = "cuhre";
   std::string setname = "MSTW2008nnlo90cl";
   LHAPDF::PDFSet setk(setname);
	nmem = setk.size()-1;
	pdfs = setk.mkPDFs();

	pids = pdfs[0]->flavors();
   // load PDF
	   
   cout << "HI" << endl;
   initPDFSet(setname, LHGRID);

   // set factorisation and renormalisation scales
   muf = Qfix;
   mur = muf;

   //evaluate alpha_s at different scales (for checks)
   double asmur = pdfs[0]->alphasQ(mur);
   double asmuf = pdfs[0]->alphasQ(muf);
   double asmuZ = pdfs[0]->alphasQ(91.1876);
   //evaluate invariant mass prefactor (for checks)
   double pre = Prefactor(S,Qfix);

   // evaluate luminosity at specific point (for checks)

   double tfix = 0.09;
   double zfix = 0.699999999;
   double lum1= lum(tfix,tau/zfix,muf,pdfs[0]);

   tfix = 0.20;
   zfix = 0.899999999;
   double lum2= lum(tfix,tau/zfix,muf,pdfs[0]);

   double aux = 1.;   

   double sigmaDlo = fsigmaD(tfix,aux,Qfix,S,muf,mur,pdfs[0],0);
   double sigmaDnlo = fsigmaD(tfix,aux,Qfix,S,muf,mur,pdfs[0],1);
   double sigmaZnlo = fsigmaZ(tfix,zfix,Qfix,S,muf,mur,pdfs[0],1);

   double* TabsigmaDlo = new double[100];
   
   for(int i = 0; i < 100; i++){
		  
      tfix = 0.0 + 0.01*i;
      TabsigmaDlo[i] = fsigmaD(tfix,aux,Qfix,S,muf,mur,pdfs[0],0);  
   }
      
   double* MassD0 = suaveMassD(Qfix,S,muf,mur,pdfs[0],0,1);
   double* MassZ0 = suaveMassZ(Qfix,S,muf,mur,pdfs[0],0,1);
  
   double* MassD1 = suaveMassD(Qfix,S,muf,mur,pdfs[0],1,1);
   double* MassZ1 = suaveMassZ(Qfix,S,muf,mur,pdfs[0],1,1);

 /* 
   double* Mass0v = fMass(Qfix,S,muf,mur,member,0,"vegas");
   double* Mass0s = fMass(Qfix,S,muf,mur,member,0,"suave");
   double* Mass0c = fMass(Qfix,S,muf,mur,member,0,"cuhre");

   double* Mass1v = fMass(Qfix,S,muf,mur,member,1,"vegas");
   double* Mass1s = fMass(Qfix,S,muf,mur,member,1,"suave");
   double* Mass1c = fMass(Qfix,S,muf,mur,member,1,"cuhre");


///////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// print result on the terminal
///
///////////////////////////////////////////////////////////////////////////////////////////////////////

   // set format for terminal ourput 
   cout.precision(8);

   cout << "********************************" << endl  << endl; 

   cout << "mur = " << mur  << endl; 
   cout << "muf = " << muf  << endl  << endl; 

   cout << "********************************" << endl  << endl; 

   cout << "asmur = " << asmur  << endl; 
   cout << "asmuf = " << asmuf  << endl; 
   cout << "asmuZ = " << asmuZ  << endl  << endl; 

   cout << "********************************" << endl  << endl; 

   cout << "Prefactor = " << pre  << endl  << endl; 

   cout << "********************************" << endl  << endl; 

   cout << "Luminosity test 1 = " << lum1  << endl; 
   cout << "Luminosity test 2 = " << lum2 << endl  << endl; 

   cout << "********************************" << endl  << endl; 

   cout << "Invariant mass integrand LO = " << sigmaDlo  << endl; 
   cout << "Invariant mass integrand NLO delta term = " << sigmaDnlo  << endl; 
   cout << "Invariant mass integrand NLO z term = " << sigmaZnlo << endl  << endl;

   for(int i = 0; i < 100; i++){
	   
      tfix = 0.0 + 0.01*i;

      cout << "Invariant mass integrand LO for t = " << tfix  << " : " << TabsigmaDlo[i] << endl;     
   }

   cout << "********************************" << endl  << endl; 

   cout << "Invariant mass LO D = " << MassD0[0] << " += " << MassD0[1] << endl; 
   cout << "Invariant mass LO Z = " << MassZ0[0] << " += " << MassZ0[1] << endl << endl; 

   cout << "********************************" << endl  << endl; 
   
   cout << "Invariant mass NLO D = " << MassD1[0] << " += " << MassD1[1] << endl; 
   cout << "Invariant mass NLO Z = " << MassZ1[0] << " += " << MassZ1[1] << endl << endl; 
      
   cout << "********************************" << endl  << endl; 

   cout << "Invariant mass LO vegas = " << Mass0v[0] << " += " << Mass0v[1] << endl; 
   cout << "Invariant mass LO suave = " << Mass0s[0] << " += " << Mass0s[1] << endl; 
   cout << "Invariant mass LO cuhre = " << Mass0c[0] << " += " << Mass0c[1] << endl << endl; 

   cout << "********************************" << endl  << endl; 
   
   cout << "Invariant mass NLO vegas = " << Mass1v[0] << " += " << Mass1v[1] << endl; 
   cout << "Invariant mass NLO suave = " << Mass1s[0] << " += " << Mass1s[1] << endl; 
   cout << "Invariant mass NLO cuhre = " << Mass1c[0] << " += " << Mass1c[1] << endl << endl; 
      
   cout << "********************************" << endl  << endl; 
   */
   
   return 0;
}
