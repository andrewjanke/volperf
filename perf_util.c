#include <math.h>
#include "perf_util.h"
#include "minc_vector_io.h"

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

void convolve_vector(gsl_vector * data, gsl_vector * kernel)
{
   gsl_vector *tmp;
   size_t   half_ks;
   size_t   i;
   int      j;

   if((kernel->size % 2) == 0){
      fprintf(stderr, "convolve_vector passed an even kernel! [%d]\n", kernel->size);
      exit(EXIT_FAILURE);
      }

   tmp = gsl_vector_calloc(data->size);

/* 
 *	do the convolution 
 */
   half_ks = (kernel->size - 1) / 2;
   for(i = half_ks; i < (data->size - half_ks); i++){
      for(j = -half_ks; j <= (int)half_ks; j++){
         gsl_vector_set(tmp, i, gsl_vector_get(tmp, i)
                        + (gsl_vector_get(kernel, half_ks + j)
                           * gsl_vector_get(data, i - j)));
         }
      }

   gsl_vector_memcpy(data, tmp);
   gsl_vector_free(tmp);
   }

/************************************************************/

/* convolves a vector using zero padding */
void convolve_vector_zero(gsl_vector * data, gsl_vector * kernel)
{
   fprintf(stderr, "NOT IMPLEMENTED YET, dying...\n");
   exit(EXIT_FAILURE);
   convolve_vector(data, kernel);
   }

/************************************************************/

/* convolves a vector using averaging of end points */
void convolve_vector_avg(gsl_vector * data, gsl_vector * kernel)
{
   size_t   half_ks;
   double   sum0, sum1;
   size_t   i;

/* 
 *	calc avg value for ends 
 */
   half_ks = (kernel->size - 1) / 2;

   sum0 = sum1 = 0.0;
   for(i = 0; i < half_ks; i++){
      sum0 += gsl_vector_get(data, i);
      sum1 += gsl_vector_get(data, data->size - i - 1);
      }

/* 
 *	now convolve the vector 
 */
   convolve_vector(data, kernel);

/* 
 *	set the ends 
 */
   sum0 = (sum0 == 0.0) ? 0.0 : sum0 / half_ks;
   sum1 = (sum1 == 0.0) ? 0.0 : sum1 / half_ks;
   for(i = 0; i < half_ks; i++){
      gsl_vector_set(data, i, sum0);
      gsl_vector_set(data, data->size - i - 1, sum1);
      }
   }

/************************************************************/

/*
 *  create_filter_kernel -
 *
 *	create a low-pass filter
 */
gsl_vector *create_filter_kernel(void)
{
   gsl_vector *v;

   v = gsl_vector_alloc(11);
   gsl_vector_set(v, 0, -0.00390625);
   gsl_vector_set(v, 1, -0.0195312);
   gsl_vector_set(v, 2, -0.0195312);
   gsl_vector_set(v, 3, 0.078125);
   gsl_vector_set(v, 4, 0.273438);
   gsl_vector_set(v, 5, 0.382812);
   gsl_vector_set(v, 6, 0.273438);
   gsl_vector_set(v, 7, 0.078125);
   gsl_vector_set(v, 8, -0.0195312);
   gsl_vector_set(v, 9, -0.0195312);
   gsl_vector_set(v, 10, -0.00390625);

   return (v);
   }

/************************************************************/

/* 
 *	reads in a Art_IF from a file
 */
Status input_art_if(const char *aif_file, Art_IF * aif)
{
   MINC_Vector *tmp;
   Status   status;
   int      i;

   tmp = new_MINC_Vector(0);

   if(input_MINC_Vector(aif_file, tmp) == OK){
      aif->signal = gsl_vector_alloc(tmp->size);
      aif->conc = gsl_vector_alloc(tmp->size);
      for(i = 0; i < tmp->size; i++){
         gsl_vector_set(aif->signal, i, tmp->V[i]);
         }
      status = OK;
      }
   else {
      status = ERROR;
      }
   free_MINC_Vector(tmp);

   return (status);
   }

/************************************************************/

void print_vector(const char *name, gsl_vector * v)
{
   size_t   i;

   fprintf(stdout, "%s[%d]: ", name, v->size);
   for(i = 0; i < v->size; i++){
      fprintf(stdout, "%.10g ", gsl_vector_get(v, i));
      }
   fprintf(stdout, "\n");
   }

/************************************************************/

void print_matrix(const char *name, gsl_matrix * m)
{
   size_t   i, j;

   fprintf(stdout, "%s[%d,%d]:\n", name, m->size1, m->size2);
   for(i = 0; i < m->size1; i++){
      for(j = 0; j < m->size2; j++){
         fprintf(stdout, "%.10g ", gsl_matrix_get(m, i, j));
         }
      fprintf(stdout, "\n");
      }
   fprintf(stdout, "\n");
   }

/************************************************************/

void SignalToConc(gsl_vector * s, gsl_vector * c, double baseline, double te)
{
   int      i;

/*
 *      convert from an intensity to a concentration
 */
   for(i = 0; i < c->size; i++){
      gsl_vector_set(c, i, -log(gsl_vector_get(s, i) / baseline) / te);
      }

/*
 *      Make sure nothing is negative
 */
   for(i = 0; i < c->size; i++){
      if(gsl_vector_get(c, i) < 0.0){
         gsl_vector_set(c, i, 0.0);
         }
      }
   }

/************************************************************/

void ConvertParm(const gsl_vector * v, double *gparm)
{
   double   bat, fac, alpha, beta;

   bat = gsl_vector_get(v, 0);
   if((gsl_vector_get(v, 3) > 0.0) && (gsl_vector_get(v, 1) > 0.0)){
      alpha = (log(gsl_vector_get(v, 3)) - log(gsl_vector_get(v, 1)))
         / (log(2.0) - 1.0);
      beta = (gsl_vector_get(v, 2) - bat) / alpha;
      if((alpha > 0.0) && (beta > 0.0)){
         fac = gsl_vector_get(v, 1) / (pow(alpha * beta, alpha)
                                       * exp(-alpha));
         }
      else {
         fac = 1.0;
         }
      }
   else {
      alpha = 1.0;
      beta = 1.0;
      fac = 1.0;
      }

   gparm[GAM_BAT] = bat;
   gparm[GAM_FAC] = fac;
   gparm[GAM_ALPHA] = alpha;
   gparm[GAM_BETA] = beta;
   }

/************************************************************/

/*
 *  eval_gamma -
 *
 *	create a gamma variate function based upon the four
 *	parameters
 */
void eval_gamma(gsl_vector * v, double *gparm, double step)
{
   double   t;
   int      i;

   for(i = 0; i < v->size; i++){
      t = i * step - gparm[GAM_BAT];
      if(t <= 0.0){
         gsl_vector_set(v, i, 0.0);
         }
      else {
         gsl_vector_set(v, i, gparm[GAM_FAC] * powf(t, gparm[GAM_ALPHA])
                        * exp(-t / gparm[GAM_BETA]));
         }
      }
   }

/************************************************************/

/*
 *  fit_gamma -
 *
 *      given a concentration profile and an estimated bat
 *      determine the fac, alpha, and beta parameters
 */
/*
void fit_gamma(gsl_vector *log_conc, gsl_vector *t, double bat,
 double *gamma, gsl_matrix *mwork, gsl_matrix *covwork,
 gsl_multifit_linear_workspace *fwork){
  double t, chisq;
  int i, j;
*/

/*
 *      setup the data
 */
/*
  for (i = 0; i < log_conc->size; i++){
    t = gsl_vector_get(t, i) - bat;
    gsl_matrix_set(mwork, i, GAM_FAC, 1.0);
    gsl_matrix_set(mwork, i, GAM_ALPHA, -t);
    gsl_matrix_set(mwork, i, GAM_BETA, log(t));
                          }
*/

/*
 *      linear regression
 */
//  gsl_multifit_linear(mwork, log_conc, gamma, covwork, &chisq, fwork);

/*
 *      convert log(fac) to fac
 */
/*
  gamma[GAM_FAC] = exp(gamma[GAM_FAC]);
                        }
*/

/************************************************************/

double vector_distance(gsl_vector * a, gsl_vector * b)
{
   double   dist = 0.0;
   int      i;

   if(a->size != b->size){
      fprintf(stderr,
              "ERROR - program trying to calculate the distance between vectors of different sizes\n");
      exit(EXIT_FAILURE);
      }

   for(i = 0; i < a->size; i++){
      dist += powf(gsl_vector_get(a, i) - gsl_vector_get(b, i), 2.0);
      }

   return (sqrt(dist));
   }

/************************************************************/

double vector_sum(gsl_vector * v, int a, int b)
{
   double   sum = 0.0;
   int      i;

   if(a > b){
      fprintf(stderr, "ERROR - can not sum a vector from %d to %d\n", a, b);
      exit(EXIT_FAILURE);
      }
   if(b >= v->size){
      fprintf(stderr, "ERROR - can not sum a vector past it's length\n");
      exit(EXIT_FAILURE);
      }

   for(i = a; i < b; i++){
      sum += gsl_vector_get(v, i);
      }

   return (sum);
   }

/************************************************************/

double bolus_arrival_time(gsl_vector * conc, double tr, Bat_info * b)
{
   int      peakx, batx, i;
   double   peaky, baty, tempy, m;
   double   arrival;

   int      mina, minb, tmp;
   double   vala = DBL_MAX, valb = DBL_MAX;
   double   euc;
   double   m1, m2, c1, c2;
   double   width, scale;

   switch (b->type){
   case BAT_GAMMA:
      arrival = b->gparm[GAM_BAT];
      break;

   case BAT_SLOPE:
      peakx = gsl_vector_max_index(conc);
      peaky = gsl_vector_get(conc, peakx);
      m = peaky / peakx;
      batx = 0;
      baty = 0.0;
      for(i = 1; i < peakx; i++){
         tempy = m * i - gsl_vector_get(conc, i);
         if(tempy > baty){
            batx = i;
            baty = tempy;
            }
         }
      arrival = tr * batx;
      break;

   case BAT_CUTOFF:
      i = gsl_vector_max_index(conc);
      peaky = gsl_vector_get(conc, i);
      while((i > 0) && (gsl_vector_get(conc, i) > b->cutoff * peaky)){
         i--;
         }
      if(i > 0){
         arrival = tr * (i + ((b->cutoff * peaky) - gsl_vector_get(conc, i))
                         / (gsl_vector_get(conc, i + 1) - gsl_vector_get(conc, i)));
         }
      else {
         arrival = 0.0;
         }
      break;
   
   /* see the attached jpg to understand what the heck is going on here (min_dist.jpg) */
   case BAT_MIN_DIST:

      /* get the x and y co-ordinate of the peak */
      peakx = gsl_vector_max_index(conc);
      
      /* special case with only one possible point */
      if(peakx < 2){
         arrival = mina;
         }
      
      else{
         peaky = gsl_vector_get(conc, peakx);
      
      
         /* get the half-width of the peak */
         i = peakx;
         while(i > 0 && gsl_vector_get(conc, i) > peaky /2){
            i--;
            }
         m1 = gsl_vector_get(conc, i+1) - gsl_vector_get(conc, i);
         c1 = gsl_vector_get(conc, i) - (m1 * i);
         width = peakx - (((peaky/2) - c1)/m1);
      
         /* get the scale */
         scale = peaky / width * 2;
      
         /* find the two points with the minimum distance */
         mina = minb = 0;
         vala = valb = sqrt(SQR(peakx * tr * scale) + SQR(gsl_vector_get(conc, 0) * 2));
      
         for(i = 0; i < peakx; i++){

            /* calc euc distance */
            euc = sqrt(SQR((peakx - i) * tr * scale) + SQR(gsl_vector_get(conc, i) * 2));
            
            /* don't even ask what's happening here */
            if(euc < vala){
               valb = vala;
               minb = mina;

               vala = euc;
               mina = i;
               }
            else if(euc < valb){
               valb = euc;
               minb = i;
               }
            }
      
         /* make sure the points are ordered */
         if(mina > minb){
            tmp = minb;
            minb = mina;
            mina = tmp;
            }
         
         /* find the closest point on the line mina:minb to our point */
         m1 = (gsl_vector_get(conc, mina) - gsl_vector_get(conc, minb))
            / (mina - minb) * scale;
         c1 = (gsl_vector_get(conc, mina) / scale) - (m1 * mina);
         m2 = -1 / m1;
         c2 = -peaky - (m1 * peakx);
         arrival = (c2 - c1) / (m1 - m2);

         /* check that we aren't out of bounds */
         if(arrival < mina){
            arrival = mina;
            }
         else if(arrival > minb){
            arrival = minb;
            }

         /* scale arrival index to the time-series */
         arrival *= tr;
         }
      
      break;

   default:
      fprintf(stderr, "ERROR - invalid choice of method for bolus arrival time\n");
      exit(EXIT_FAILURE);
      }

   return (arrival);
   }

/************************************************************/

/*
 *  CalcAifMatrix -
 *
 *    Set up the AIF matrix for a standard deconvolution (Ostergaard 1996)
 */
void CalcAifMatrix(gsl_vector * v, gsl_matrix * m, double fac, int Nstart, int Nrest)
{
   double   temp;
   size_t   i, j;

/* 
 *	diagonal elements 
 */
   for(j = 0; j < Nrest; j++){
      gsl_matrix_set(m, j, j, fac * gsl_vector_get(v, Nstart));
      }

/* 
 *	off diagonal elements 
 */
   for(i = 1; i < Nrest - 1; i++){
      temp = fac * (gsl_vector_get(v, Nstart + (i - 1))
                    + (4.0 * gsl_vector_get(v, Nstart + i))
                    + (gsl_vector_get(v, Nstart + (i + 1)))) / 6.0;
      for(j = 0; j < Nrest - i; j++){
         gsl_matrix_set(m, j + i, j, temp);
         }
      }

/* 
 *	do the last element in the corner 
 */
   gsl_matrix_set(m, Nrest - 1, 0, fac * gsl_vector_get(v, Nstart + Nrest - 1));
   }

/************************************************************/

/*
 *  CalcCircAifMatrix -
 *
 *    Set up the AIF matrix for circular deconvolution (Wu 2003)
 */
void CalcCircAifMatrix(gsl_vector * v, gsl_matrix * m, double fac, int Nstart, int Nrest)
{
   double   temp;
   size_t   i, j;

/* 
 *	Ca(t0)
 */
   for(j = 0; j < 2 * Nrest; j++){
      gsl_matrix_set(m, j, j, fac * gsl_vector_get(v, Nstart));
      }

/* 
 *	Ca(t1) ... Ca(tn-2) 
 */
   for(i = 1; i < Nrest - 1; i++){
      temp = fac * (gsl_vector_get(v, Nstart + (i - 1))
                    + (4.0 * gsl_vector_get(v, Nstart + i))
                    + (gsl_vector_get(v, Nstart + (i + 1)))) / 6.0;
      for(j = 0; j < 2 * Nrest; j++){
         gsl_matrix_set(m, j, (j - i + 2 * Nrest) % (2 * Nrest), temp);
         }
      }

/* 
 *	Ca(tn-1)
 */
   for(j = 0; j < 2 * Nrest; j++){
      gsl_matrix_set(m, j, (j + Nrest + 1) % (2 * Nrest),
                     fac * gsl_vector_get(v, Nstart));
      }
   }

/************************************************************/

/*
 * svd_oscillation_index -
 *
 *     calculate svd oscillation index (as a double) 
 */
double svd_oscillation_index(gsl_vector * v)
{
   int      i;
   double   oi = 0.0;

   for(i = 1; i < v->size - 1; i++){
      oi +=
         fabs(gsl_vector_get(v, i - 1) - 2 * gsl_vector_get(v, i) +
              gsl_vector_get(v, i + 1));
      }
   return (oi / (v->size * gsl_vector_max(v)));
   }


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
