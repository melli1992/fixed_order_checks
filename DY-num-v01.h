
#include "LHAPDF/LHAPDF.h"

double Prefactor(double &S, double &Q);


double lum(double t, double y, double &mu, LHAPDF::PDF* pdf);

double* lumi(double t, double y, double &mu, LHAPDF::PDF* pdf);


double fsigmaD(double t, double aux, double &Q, double &S, double &muf, double &mur, LHAPDF::PDF* pdf, int ord);

double fsigmaZ(double t, double z, double &Q, double &S, double &muf, double &mur, LHAPDF::PDF* pdf, int ord);


double fmassIntD(double t, double aux, void *params);

double fmassIntZ(double t, double z, void *params);


double* vegasMassD( double &Q, double &S, double &muf, double &mur, int &member, int ord, int verbose);

double* suaveMassD( double &Q, double &S, double &muf, double &mur, LHAPDF::PDF* pdf, int ord, int verbose);

double* cuhreMassD( double &Q, double &S, double &muf, double &mur, int &member, int ord, int verbose);


double* vegasMassZ( double &Q, double &S, double &muf, double &mur, int &member, int ord, int verbose);

double* suaveMassZ( double &Q, double &S, double &muf, double &mur, LHAPDF::PDF* pdf, int ord, int verbose);

double* cuhreMassZ( double &Q, double &S, double &muf, double &mur, int &member, int ord, int verbose);


double* fMass( double &Q, double &S, double &muf, double &mur, int &member, int ord, std::string cuba_method);
