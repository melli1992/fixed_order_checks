#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORNNLO_H
#define KFACTORNNLO_H

double vegas_NNLO_qqtogg_LP(double *k, size_t dim, void *params);
double vegas_NNLO_qqtogg_LP_correction(double *k, size_t dim, void *params);
double vegas_NNLO_qqtogg_NLP(double *k, size_t dim, void *params);
double vegas_NNLO_qqtogg_NNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qqtogg_NNNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qqtogg_delta(double *k, size_t dim, void *params);
double vegas_NNLO_qqtogg_full(double *k, size_t dim, void *params);


double vegas_NNLO_LP_int(double *k, size_t dim, void *params);

#endif
