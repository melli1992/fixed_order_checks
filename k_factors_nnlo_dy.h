#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORNNLO_H
#define KFACTORNNLO_H
//needed for qqbar channel (identical quarks)
double vegas_NNLO_qqbar_LP(double *k, size_t dim, void *params);
double vegas_NNLO_qqbar_LP_correction(double *k, size_t dim, void *params);
double vegas_NNLO_qqbar_NLP(double *k, size_t dim, void *params);
double vegas_NNLO_qqbar_NNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qqbar_NNNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qqbar_delta(double *k, size_t dim, void *params);
double vegas_NNLO_qqbar_full(double *k, size_t dim, void *params);
double NNLO_qqbar_A2(double x);
double NNLO_qqbar_B2(double x);
double NNLO_qqbar_BC(double x);
double NNLO_qqbar_AC(double x);
double NNLO_qqbar_HCA(double x);
double NNLO_qqbar_HCF(double x);
double NNLO_qqbar_b0_nonconst(double x);
double NNLO_qqbar_b0_const();
double NNLO_qqbar_delta();
double NNLO_qqbar_LP(double z);
double NNLO_qqbar_b0_LP(double x);
double NNLO_qqbar_NLP(double z);
double NNLO_qqbar_NNLP(double x);
double NNLO_qqbar_NNNLP(double x);
// also an integration routine to check the result
double vegas_NNLO_LP_int(double *k, size_t dim, void *params);

//needed for qg
double vegas_NNLO_qg_NLP(double *k, size_t dim, void *params);
double vegas_NNLO_qg_NNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qg_NNNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qg_full(double *k, size_t dim, void *params);
double NNLO_qg_NLP(double x);
double NNLO_qg_NNLP(double x);
double NNLO_qg_NNNLP(double x);
double NNLO_qg_full(double x);

//needed for gg
double vegas_NNLO_gg_NLP(double *k, size_t dim, void *params);
double vegas_NNLO_gg_NNLP(double *k, size_t dim, void *params);
double vegas_NNLO_gg_NNNLP(double *k, size_t dim, void *params);
double vegas_NNLO_gg_full(double *k, size_t dim, void *params);
double NNLO_gg_NLP(double x);
double NNLO_gg_NNLP(double x);
double NNLO_gg_NNNLP(double x);
double NNLO_gg_full(double x);

//needed for qq, qqbar, qbarqbar (non-identical)
double vegas_NNLO_qqNI_NLP(double *k, size_t dim, void *params);
double vegas_NNLO_qqNI_NNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qqNI_NNNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qqNI_full(double *k, size_t dim, void *params);
double vegas_NNLO_qqbarNI_NLP(double *k, size_t dim, void *params);
double vegas_NNLO_qqbarNI_NNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qqbarNI_NNNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qqbarNI_full(double *k, size_t dim, void *params);
double vegas_NNLO_qbarqbarNI_NLP(double *k, size_t dim, void *params);
double vegas_NNLO_qbarqbarNI_NNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qbarqbarNI_NNNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qbarqbarNI_full(double *k, size_t dim, void *params);
double NNLO_qqNI_NLP(double x);
double NNLO_qqNI_NNLP(double x);
double NNLO_qqNI_NNNLP(double x);
double NNLO_qqNI_full(double x);
double NNLO_qqbarNI_NLP(double x);
double NNLO_qqbarNI_NNLP(double x);
double NNLO_qqbarNI_NNNLP(double x);
double NNLO_qqbarNI_full(double x);

//needed for qq, qbarqbar (identical)
double vegas_NNLO_qq_NLP(double *k, size_t dim, void *params);
double vegas_NNLO_qq_NNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qq_NNNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qq_full(double *k, size_t dim, void *params);
double vegas_NNLO_qbarqbar_NLP(double *k, size_t dim, void *params);
double vegas_NNLO_qbarqbar_NNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qbarqbar_NNNLP(double *k, size_t dim, void *params);
double vegas_NNLO_qbarqbar_full(double *k, size_t dim, void *params);
double NNLO_qq_NLP(double x);
double NNLO_qq_NNLP(double x);
double NNLO_qq_NNNLP(double x);
double NNLO_qq_full(double x);
#endif
