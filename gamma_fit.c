#include <math.h>
#include "gamma_fit.h"

/* gamma variate fit structure */
typedef struct {
   gsl_vector *signal;
   gsl_vector *conc;
   gsl_vector *fit;
   gsl_vector *gparm;
   double   te;
   double   tr;
   } Gamma_Struct;

/* function prototypes */
double   gamma_f(const gsl_vector * v, void *params);
void     gamma_df(const gsl_vector * v, void *params, gsl_vector * g);
void     gamma_fdf(const gsl_vector * x, void *params, double *f, gsl_vector * df);

double   gamma_f_aif(const gsl_vector * v, void *params);
void     gamma_fdf_aif(const gsl_vector * x, void *params, double *f, gsl_vector * df);

/************************************************************/

void TimeCurveArea(gsl_vector * conc, double tr, double *area, double *gparm)
{
   int      i, status, peak;

   const gsl_multimin_fdfminimizer_type *T;
   gsl_multimin_fdfminimizer *s;
   gsl_vector *gamma_vary;
   gsl_multimin_function_fdf func;
   Gamma_Struct gamma_const;

   /* set up "func" data and parameters */
   func.f = &gamma_f;
   func.df = &gamma_df;
   func.fdf = &gamma_fdf;
   func.n = 4;

   gamma_const.conc = conc;
   gamma_const.tr = tr;
   func.params = &gamma_const;

   if(gsl_vector_max(conc) == 0){
      *area = 0.0;
      for(i = 0; i < GAM_NUM; i++){
         gparm[i] = 0.0;
         }
      return;
      }

/* 
 *	predict the location of the global minima
 */
   peak = gsl_vector_max_index(conc);
   if(peak == 0){
      *area = 0.0;
      for(i = 0; i < GAM_NUM; i++){
         gparm[i] = 0.0;
         }
      return;
      }

   gamma_vary = gsl_vector_alloc(4);
   gsl_vector_set(gamma_vary, 1, gsl_vector_get(conc, peak));
   gsl_vector_set(gamma_vary, 2, peak * tr);

   i = peak;
   while((gsl_vector_get(conc, i)
          > 0.5 * (gsl_vector_get(gamma_vary, 1))) && (i > 0)){
      i--;
      }
   gsl_vector_set(gamma_vary, 0, (2.0 * i - peak) * tr);
   gsl_vector_set(gamma_vary, 3, gsl_vector_get(gamma_vary, 1) / 10.0);
   gamma_const.fit = gsl_vector_alloc(conc->size);

/* 
 *	do funky wierd GSL stuff 
 */
   T = gsl_multimin_fdfminimizer_vector_bfgs;
   s = gsl_multimin_fdfminimizer_alloc(T, 4);
   gsl_multimin_fdfminimizer_set(s, &func, gamma_vary, 1.0e-1, 1.0e-2);

   i = 0;
   status = GSL_CONTINUE;
   while(status == GSL_CONTINUE && i++ < 100){
      if(gsl_multimin_fdfminimizer_iterate(s)){
         break;
         }
      status = gsl_multimin_test_gradient(s->gradient, 0.5);
      }

/* 
 *	area under curve 
 */
   ConvertParm(s->x, gparm);
   *area = gparm[GAM_FAC] * pow(gparm[GAM_BETA], gparm[GAM_ALPHA] + 1)
      * gsl_sf_gamma(gparm[GAM_ALPHA] + 1);

   gsl_vector_free(gamma_vary);
   gsl_vector_free(gamma_const.fit);
   gsl_multimin_fdfminimizer_free(s);
   }

/************************************************************/

double gamma_f(const gsl_vector * v, void *params)
{
   Gamma_Struct *p = (Gamma_Struct *) params;
   double   gparm[GAM_NUM];

   ConvertParm(v, gparm);
   eval_gamma(p->fit, gparm, p->tr);

   return (vector_distance(p->fit, p->conc));
   }

/************************************************************/

void gamma_fdf(const gsl_vector * x, void *params, double *f, gsl_vector * df)
{
   *f = gamma_f(x, params);
   gamma_df(x, params, df);
   }

/************************************************************/

void TimeCurveAreaAif(gsl_vector * signal, gsl_vector * conc,
                      double tr, double te, double *area, double *gparm, int *start,
                      int *rest)
{
   const gsl_multimin_fdfminimizer_type *T;
   gsl_multimin_fdfminimizer *s;
   gsl_vector *gamma_vary;
   gsl_multimin_function_fdf func;
   Gamma_Struct gamma_const;
   int      i, status, peak, N;
   double   baseline;

   /* set up "func" data and parameters */
   func.f = &gamma_f_aif;
   func.df = &gamma_df;
   func.fdf = &gamma_fdf_aif;
   func.n = 4;

   N = signal->size;
   gamma_const.signal = signal;
   gamma_const.conc = conc;
   gamma_const.tr = tr;
   gamma_const.te = te;
   func.params = &gamma_const;

   if(gsl_vector_min(signal) == 0){
      fprintf(stderr, "Error in Aif\n");
      exit(EXIT_FAILURE);
      }

/* 
 *	convert the signal to a concentration
 */
   baseline = gsl_vector_max(signal);
   SignalToConc(signal, conc, baseline, te);
   peak = gsl_vector_max_index(conc);
   if(peak == 0){
      gsl_vector_free(conc);
      fprintf(stderr, "Error in Aif\n");
      exit(EXIT_FAILURE);
      }

/*
 *	predict the location of the global minima
 */
   gamma_vary = gsl_vector_alloc(4);
   gsl_vector_set(gamma_vary, 1, gsl_vector_get(conc, peak));
   gsl_vector_set(gamma_vary, 2, peak * tr);

   i = peak;
   while((gsl_vector_get(conc, i)
          > 0.5 * (gsl_vector_get(gamma_vary, 1))) && (i >= 0)){
      i--;
      }
   gsl_vector_set(gamma_vary, 0, (2.0 * i - peak) * tr);
   gsl_vector_set(gamma_vary, 3, gsl_vector_get(gamma_vary, 1) / 10.0);
   gamma_const.fit = gsl_vector_alloc(N);
   gamma_const.gparm = gsl_vector_alloc(GAM_NUM);

/* 
 *	do funky wierd GSL stuff 
 */
   T = gsl_multimin_fdfminimizer_vector_bfgs;
   s = gsl_multimin_fdfminimizer_alloc(T, 4);
   gsl_multimin_fdfminimizer_set(s, &func, gamma_vary, 0.01, 1.0e-4);

   i = 0;
   status = GSL_CONTINUE;
   while(status == GSL_CONTINUE && i++ < 100){
      if(gsl_multimin_fdfminimizer_iterate(s)){
         break;
         }
      status = gsl_multimin_test_gradient(s->gradient, 1.0e-6);
      }

/* 
 *	area under curve 
 */
   ConvertParm(s->x, gparm);
   *area = gparm[GAM_FAC] * pow(gparm[GAM_BETA], gparm[GAM_ALPHA] + 1)
      * gsl_sf_gamma(gparm[GAM_ALPHA] + 1);
   *start = rint(gparm[GAM_BAT] / tr);
   *rest = signal->size - *start;

   gsl_vector_free(gamma_vary);
   gsl_vector_free(gamma_const.fit);
   gsl_vector_free(gamma_const.gparm);
   gsl_multimin_fdfminimizer_free(s);
   }

/************************************************************/

double gamma_f_aif(const gsl_vector * v, void *params)
{
   Gamma_Struct *p = (Gamma_Struct *) params;
   double   gparm[GAM_NUM], baseline;

   ConvertParm(v, gparm);
   baseline = gsl_stats_mean(p->signal->data, 1, gparm[GAM_BAT] / (p->tr));
   SignalToConc(p->signal, p->conc, baseline, p->te);
   eval_gamma(p->fit, gparm, p->tr);
   return (vector_distance(p->fit, p->conc));
   }

/************************************************************/

void gamma_fdf_aif(const gsl_vector * x, void *params, double *f, gsl_vector * df)
{
   *f = gamma_f_aif(x, params);
   gamma_df(x, params, df);
   }

/************************************************************/

void gamma_df(const gsl_vector * v, void *params, gsl_vector * df)
{
   Gamma_Struct *p = (Gamma_Struct *) params;
   size_t   i;
   double   fit_i, t;
   double   gparm[GAM_NUM];
   double   bat, fac, alpha, beta;
   double   dfac, dbat, dalpha, dbeta;
   double   dfac2, dalpha2, dbeta2;

   gsl_vector_set_zero(df);
   ConvertParm(v, gparm);
   bat = gparm[GAM_BAT];
   fac = gparm[GAM_FAC];
   alpha = gparm[GAM_ALPHA];
   beta = gparm[GAM_BETA];

   dfac = 0.0;
   dbat = 0.0;
   dalpha = 0.0;
   dbeta = 0.0;

   for(i = 0; i < p->conc->size; i++){
      /* retrieve the value of fit_i (set previously in gamma_f) */
      fit_i = gsl_vector_get(p->fit, i);
      t = i * (p->tr) - bat;

      if(t > 0){
         dfac += 2 * (fit_i - gsl_vector_get(p->conc, i)) * fit_i / fac;
         if(t != 0.0){
            dbat += 2 * (fit_i - gsl_vector_get(p->conc, i)) * fit_i
               * ((1 / beta) - (alpha / t));
            dalpha += 2 * (fit_i - gsl_vector_get(p->conc, i)) * fit_i * log(t);
            dbeta += 2 * (fit_i - gsl_vector_get(p->conc, i)) * fit_i * t / (beta * beta);
            }
         }
      }
   gsl_vector_set(df, 0, dbat + (dfac * fac / beta) - (dbeta / alpha));

   dalpha2 = 1.0 / (gsl_vector_get(v, 1) * (1.0 - log(2.0)));
   dbeta2 = -dalpha2 * beta / alpha;
   dfac2 = pow(alpha * beta, -alpha) * exp(alpha);
   dfac2 = dfac2 * (1.0 + gsl_vector_get(v, 1) *
                    (-log(alpha * beta) * dalpha2 - (alpha / beta) * dbeta2));
   gsl_vector_set(df, 1, dfac * dfac2 + dalpha * dalpha2 + dbeta * dbeta2);

   dfac2 = pow(alpha * beta, -alpha) * exp(alpha);
   dfac2 = -dfac2 * gsl_vector_get(v, 1) / beta;
   gsl_vector_set(df, 2, dfac * dfac2 + dbeta / alpha);

   dalpha2 = 1.0 / (gsl_vector_get(v, 3) * (log(2.0) - 1.0));
   dbeta2 = -dalpha2 * beta / alpha;
   dfac2 = pow(alpha * beta, -alpha) * exp(alpha);
   dfac2 = dfac2 * gsl_vector_get(v, 1) *
      (-log(alpha * beta) * dalpha2 - (alpha / beta) * dbeta2);
   gsl_vector_set(df, 3, dfac * dfac2 + dalpha * dalpha2 + dbeta * dbeta2);
   }
