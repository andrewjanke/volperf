/* gamma_fit.h */

#ifndef GAMMA_FIT_H
#define GAMMA_FIT_H

#include "perf_util.h"

/* calculate the area under a time curve */
void     TimeCurveArea(gsl_vector * conc, double tr, double *area, double *gparm);
void     TimeCurveAreaAif(gsl_vector * aif, gsl_vector * conc, double tr,
                          double te, double *area, double *gparm, int *start, int *rest);

#endif
