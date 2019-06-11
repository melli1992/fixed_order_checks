#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTOR_H
#define KFACTOR_H

double LO_factor();
double k_zint_NLO_gg_dy_LP(double z);
double k_zint_NLO_gg_dy_NLP(double z);
double k_zint_NLO_gg_dy_NNLP(double z);
double k_zint_NLO_gg_dy_exact(double z);
double k_zint_NLO_gg_dy_delta();
double k_NLO_dy_gg_delta();
double k_NLO_dy_gg_nonconst(double z);
double vegas_sum_pdf(double *k, size_t dim, void *params);
double vegas_k_NLO_dy_gg(double *k, size_t dim, void *params);
double k_zint_NLO_qg_dy_NLP(double z);
double k_zint_NLO_qg_dy_NNLP(double z);
double k_zint_NLO_qg_dy_exact(double z);

double vegas_LO(double *k, size_t dim, void *params);

double vegas_sig_LP_1(double *k, size_t dim, void *params);
double vegas_sig_LP_correction(double *k, size_t dim, void *params);
double vegas_sig_NLP(double *k, size_t dim, void *params);
double vegas_sig_NNLP(double *k, size_t dim, void *params);
double vegas_sig_delta(double *k, size_t dim, void *params);
double vegas_sig_full(double *k, size_t dim, void *params);
double vegas_qg_full(double *k, size_t dim, void *params);

double vegas_sig_LP_int(double *k, size_t dim, void *params);
double vegas_sig_NLP_int(double *k, size_t dim, void *params);
double vegas_sig_NNLP_int(double *k, size_t dim, void *params);
double vegas_sig_full_int(double *k, size_t dim, void *params);
double vegas_qg_full_int(double *k, size_t dim, void *params);

double vegas_sig_testfunction1(double *k, size_t dim, void *params);
double vegas_sig_testfunction2(double *k, size_t dim, void *params);



//needed for qg
double vegas_NLO_qg_NLP(double *k, size_t dim, void *params);
double vegas_NLO_qg_NNLP(double *k, size_t dim, void *params);
double vegas_NLO_qg_NNNLP(double *k, size_t dim, void *params);
double vegas_NLO_qg_full(double *k, size_t dim, void *params);
double NLO_qg_NLP(double x);
double NLO_qg_NNLP(double x);
double NLO_qg_NNNLP(double x);
double NLO_qg_full(double x);

#endif
