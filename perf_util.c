/* perf_util.c */

#include "perf_util.h"
#include "minc_vector_io.h"

extern int verbose;

/* gamma variate fit structure */
typedef struct {
   gsl_vector *tcurve;
   int      tcurve_size;
   int      start_tpoint;
   double   lower[4];
   double   upper[4];
   gsl_vector *fit_i;
} Gamma_Struct;

/* function prototypes */
void     convolve_vector(gsl_vector * data, gsl_vector * kernel);
double   gamma_f(const gsl_vector * v, void *params);
void     gamma_df(const gsl_vector * v, void *params, gsl_vector * g);
void     gamma_fdf(const gsl_vector * x, void *params, double *f, gsl_vector * df);


void do_math(void *caller_data, long num_voxels, int input_num_buffers, int input_vector_length,
             double *input_data[], int output_num_buffers, int output_vector_length,
             double *output_data[], Loop_Info * loop_info)
{
   /* Get pointer to window info */
   Math_Data *md = (Math_Data *) caller_data;

   size_t   i, j;
   long     ivox;
   double   raw_baseline, raw_avg;
   double   tmp, sum;

   gsl_vector *conc_ts;
   gsl_vector *residue;

   /* calc_perf SVD vector and matricies */
   conc_ts = gsl_vector_alloc(md->rem_tpoints);
   residue = gsl_vector_alloc(md->rem_tpoints);

   for(ivox = 0; ivox < num_voxels * input_vector_length; ivox++){

      /* calc baseline for the raw time-series */
      raw_baseline = 0.0;
      for(i = 0; i < md->start_tpoint; i++){
         raw_baseline += input_data[i][ivox];
         }
      raw_baseline /= md->start_tpoint;

      /* calc average value for the entire series */
      raw_avg = 0.0;
      for(i = 0; i < (size_t) input_num_buffers; i++){
         raw_avg += input_data[i][ivox];
         }
      raw_avg /= input_num_buffers;

      /* check if we have an abberant timeseries */
      if((raw_baseline <= 0) || (raw_avg / md->aif->baseline < md->cutoff)){
         for(i = 0; i < (size_t) output_num_buffers; i++){
            output_data[i][ivox] = 0.0;
            }
         }

      /* else do some real work */
      else{

         /* convert signal intensity to concentration */
         for(i = md->start_tpoint; i < (size_t) input_num_buffers; i++){
            if(input_data[i][ivox] == 0){
               tmp = 0.0;
               }
            else{
               tmp = -log(input_data[i][ivox] / raw_baseline) / md->te;
               if(tmp < 0.0){
                  tmp = 0.0;
                  }
               }

            gsl_vector_set(conc_ts, i - md->start_tpoint, tmp);
            }


         /* filter the time-series */
         if(md->filter){
            convolve_vector_avg(conc_ts, md->kernel);
            }

         /* usd SVD to evaluate the residue function   */
         /* from the AIF and concentration time-curve  */
         gsl_linalg_SV_solve(md->svd_u, md->svd_v, md->svd_s, conc_ts, residue);

         /* chi-squared fit */
         output_data[0][ivox] = 0.0;
         for(i = 0; i < md->rem_tpoints; i++){

            sum = 0.0;
            for(j = 0; j < md->rem_tpoints; j++){
               sum += gsl_vector_get(residue, j) * gsl_matrix_get(md->svd_u, i, j);
               }

            tmp = gsl_vector_get(conc_ts, i) - sum;
            output_data[0][ivox] += tmp * tmp;
            }

         /* CBF */
         output_data[1][ivox] = gsl_vector_max(residue);

         /* CBV */
         tmp = 0.5 * (gsl_vector_get(conc_ts, 0) + gsl_vector_get(conc_ts, md->rem_tpoints - 1));
         for(i = 1; i < md->rem_tpoints - 1; i++){
            tmp += gsl_vector_get(conc_ts, i);
            }
         output_data[2][ivox] = tmp * md->tr / md->aif->area;

         /* MTT - should range between 0.0 and 1.0 */
         output_data[3][ivox] = (output_data[1][ivox] <= 0.0) ?
            0.0 : output_data[2][ivox] / output_data[1][ivox];

         /* BAT */
//         print_vector("conc_ts", conc_ts);
//         time_curve_area(conc_ts, 0, md->tr, &tmp, &output_data[4][ivox], 0);
//         fprintf(stderr, "BAT: %g\n", output_data[4][ivox]);
         }
      }

   /* free math_data memory */
   gsl_vector_free(conc_ts);
   gsl_vector_free(residue);
   return;
   }


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

   /* do the convolution */
   half_ks = (kernel->size - 1) / 2;
   for(i = half_ks; i < (data->size - half_ks); i++){
      for(j = -half_ks; j <= (int)half_ks; j++){
         gsl_vector_set(tmp, i,
                        gsl_vector_get(tmp, i) +
                        (gsl_vector_get(kernel, half_ks + j) * gsl_vector_get(data, i - j)));
         }
      }

   gsl_vector_memcpy(data, tmp);
   gsl_vector_free(tmp);
   }


/* convolves a vector using zero padding */
void convolve_vector_zero(gsl_vector * data, gsl_vector * kernel)
{

   fprintf(stderr, "NOT IMPLEMENTED YET, dying...\n");
   exit(EXIT_FAILURE);

   convolve_vector(data, kernel);
   }


/* convolves a vector using averaging of end points */
void convolve_vector_avg(gsl_vector * data, gsl_vector * kernel)
{
   size_t   half_ks;
   double   sum0, sum1;
   size_t   i;

   /* calc avg value for ends */
   half_ks = (kernel->size - 1) / 2;

   sum0 = sum1 = 0.0;
   for(i = 0; i < half_ks; i++){
      sum0 += gsl_vector_get(data, i);
      sum1 += gsl_vector_get(data, data->size - i - 1);
      }

   /* now convolve the vector */
   convolve_vector(data, kernel);

   /* set the ends */
   sum0 = (sum0 == 0.0) ? 0.0 : sum0 / half_ks;
   sum1 = (sum1 == 0.0) ? 0.0 : sum1 / half_ks;
   for(i = 0; i < half_ks; i++){
      gsl_vector_set(data, i, sum0);
      gsl_vector_set(data, data->size - i - 1, sum1);
      }
   }


void time_curve_area(gsl_vector * tcurve, size_t start_tpoint, double tr,
                     double *area, double *bat, int flg)
{
   int      iter, status;

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

   gamma_const.tcurve = tcurve;
   gamma_const.start_tpoint = start_tpoint;
   gamma_const.lower[0] = 0.0;
   gamma_const.lower[1] = -1.0;
   gamma_const.lower[2] = 0.01;
   gamma_const.lower[3] = 0.01;

   gamma_const.upper[0] = 100.0;
   gamma_const.upper[1] = 3.0;
   gamma_const.upper[2] = 30.0;
   gamma_const.upper[3] = 30.0;

   gamma_const.fit_i = gsl_vector_alloc(tcurve->size);
   func.params = &gamma_const;

   if(gsl_vector_max(tcurve) == 0){
      *area = 0.0;
      *bat = 0.0;
      return;
      }

   /* set up the guess parameters -- see Eq 19 of Ostergaard 1996 */
   gamma_vary = gsl_vector_alloc(4);
   gsl_vector_set(gamma_vary, 1, 0.0); /* arrival time */
   gsl_vector_set(gamma_vary, 2, gsl_vector_max_index(tcurve) / 2.0); /* alpha        */
   gsl_vector_set(gamma_vary, 3, 2.0); /* beta         */

   if(gsl_vector_get(gamma_vary, 2) == 0){
      gsl_vector_set(gamma_vary, 0, 0.0);
      }
   else{
      gsl_vector_set(gamma_vary, 0, gsl_vector_max(tcurve) /
                     pow(2.0 * gsl_vector_get(gamma_vary, 2) / 2.71,
                         gsl_vector_get(gamma_vary, 2)));
      }

   /* do funky wierd GSL stuff */
   T = gsl_multimin_fdfminimizer_conjugate_fr;
   s = gsl_multimin_fdfminimizer_alloc(T, 4);
   gsl_multimin_fdfminimizer_set(s, &func, gamma_vary, 0.01, 1.0e-4);

   iter = 0;
   status = GSL_CONTINUE;
   while (status == GSL_CONTINUE && iter++ < 500){

      if(gsl_multimin_fdfminimizer_iterate(s)){
         break;
         }

//      fprintf(stdout, "[%03d]  %g %g %g %g  %d\n", iter, 
//                       gsl_vector_get(s->x, 0),
//                       gsl_vector_get(s->x, 1),
//                       gsl_vector_get(s->x, 2),
//                       gsl_vector_get(s->x, 3),
//                       status);

      status = gsl_multimin_test_gradient(s->gradient, 1.0e-6);
      }

   if(flg){
      /* area under curve */
      *area = gsl_vector_get(s->x, 0) *
         pow(gsl_vector_get(s->x, 3), gsl_vector_get(s->x, 2) + 1) *
         gsl_sf_gamma(gsl_vector_get(s->x, 2) + 1) * tr;
      }

   else{
      *area = 0.0;
      }

   /* BAT */
   *bat = gsl_vector_get(s->x, 1);

//   print_vector("X:", s->x);
//   fprintf(stdout, "BAT: %g\n", *bat);

   gsl_multimin_fdfminimizer_free(s);
   gsl_vector_free(gamma_vary);
   gsl_vector_free(gamma_const.fit_i);
   }


double gamma_f(const gsl_vector * v, void *params)
{
   size_t   i;
   Gamma_Struct *p = (Gamma_Struct *) params;

   double   fit = 0.0;
   double   fit_i;
   gsl_vector *v_cp;

   v_cp = gsl_vector_alloc(v->size);
   gsl_vector_memcpy(v_cp, v);

   for(i = 0; i < 4; i++){
      if(gsl_vector_get(v, i) <= p->lower[i]){
         fit += 1000000 * pow(gsl_vector_get(v, i) - p->lower[i], 2.0);
         gsl_vector_set(v_cp, i, p->lower[i]);
         }
      else if(gsl_vector_get(v, i) >= p->upper[i]){
         fit += 1000000 * pow(gsl_vector_get(v, i) - p->upper[i], 2.0);
         gsl_vector_set(v_cp, i, p->upper[i]);
         }
      }

   for(i = p->start_tpoint; i < p->tcurve->size; i++){
      if(i > gsl_vector_get(v_cp, 1)){
         if(gsl_vector_get(v_cp, 3) == 0){
            fit_i = 1000000;
            }
         else{
            fit_i = gsl_vector_get(v_cp, 0) * pow((i - gsl_vector_get(v_cp, 1)),
                                                  gsl_vector_get(v_cp,
                                                                 2)) * exp(-(i -
                                                                             gsl_vector_get(v_cp,
                                                                                            1)) /
                                                                           gsl_vector_get(v_cp, 3));
            }
         }
      else{
         fit_i = 0.0;
         }

      /* keep the value of fit_i for gamma_df */
      gsl_vector_set(p->fit_i, i, fit_i);

      fit += pow(fit_i - gsl_vector_get(p->tcurve, i), 2.0);
      }

   gsl_vector_free(v_cp);
   return fit;
   }


void gamma_df(const gsl_vector * v, void *params, gsl_vector * df)
{
   Gamma_Struct *p = (Gamma_Struct *) params;
   size_t   i;
   double   fit_i;
   gsl_vector *v_cp;
   int      flg;

   gsl_vector_set_zero(df);

   flg = 0;
   for(i = 0; i < 4; i++){
      if(gsl_vector_get(v, i) <= p->lower[i]){
         gsl_vector_set(df, i, gsl_vector_get(df, i) +
                        2000000 * (gsl_vector_get(v, i) - p->lower[i]));
         flg = 1;
         }
      if(gsl_vector_get(v, i) >= p->upper[i]){
         gsl_vector_set(df, i, gsl_vector_get(df, i) +
                        2000000 * (gsl_vector_get(v, i) - p->upper[i]));
         flg = 1;
         }
      }

   if(!flg){
      v_cp = gsl_vector_alloc(v->size);
      gsl_vector_memcpy(v_cp, v);

      for(i = p->start_tpoint; i < p->tcurve->size; i++){

         /* retrieve the value of fit_i (set previously in gamma_f) */
         fit_i = gsl_vector_get(p->fit_i, i);

         if(i > gsl_vector_get(v, 1)){

            gsl_vector_set(df, 0, gsl_vector_get(df, 0) +
                           2 * (fit_i - gsl_vector_get(p->tcurve, i)) *
                           fit_i / gsl_vector_get(v_cp, 0));

            gsl_vector_set(df, 1, gsl_vector_get(df, 1) +
                           2 *
                           (fit_i - gsl_vector_get(p->tcurve, i)) *
                           fit_i *
                           ((1 / gsl_vector_get(v_cp, 3)) -
                            (gsl_vector_get(v_cp, 2) / (i - gsl_vector_get(v_cp, 1)))
                           )
               );

            gsl_vector_set(df, 2, gsl_vector_get(df, 2) +
                           2 * (fit_i - gsl_vector_get(p->tcurve, i)) *
                           fit_i * log(i - gsl_vector_get(v_cp, 1)));

            gsl_vector_set(df, 3, gsl_vector_get(df, 3) -
                           2 * (fit_i - gsl_vector_get(p->tcurve, i)) *
                           fit_i * (i - gsl_vector_get(v_cp, 1)) /
                           (gsl_vector_get(v_cp, 3) * gsl_vector_get(v_cp, 3)));
            }
         }

      gsl_vector_free(v_cp);
      }
   }


void gamma_fdf(const gsl_vector * x, void *params, double *f, gsl_vector * df)
{

   *f = gamma_f(x, params);
   gamma_df(x, params, df);
   }


void calc_aif_matrix(Art_IF * aif, size_t start_tpoint, size_t rem_tpoints,
                     double tr, double te, gsl_matrix * aif_m)
{

   size_t   i, j;
   double   sum, tmp;

   /* calculate the AIF baseline */
   sum = 0.0;
   for(i = 0; i < start_tpoint; i++){
      sum += gsl_vector_get(aif->AIF, i);
      }
   aif->baseline = sum / start_tpoint;

   /* convert from an intensity to a concentration */
   for(i = 0; i < aif->AIF->size; i++){
      tmp = -log(gsl_vector_get(aif->AIF, i) / aif->baseline) / te;

      gsl_vector_set(aif->AIF, i, (gsl_finite(tmp)) ? tmp : 0.0);
      }

   /* Make sure nothing is negative */
   for(i = 0; i < aif->AIF->size; i++){
      if(gsl_vector_get(aif->AIF, i) < 0.0){
         gsl_vector_set(aif->AIF, i, 0.0);
         }
      }

   /* diagonal elements */
   for(j = 0; j < rem_tpoints; j++){
      gsl_matrix_set(aif_m, j, j, gsl_vector_get(aif->AIF, start_tpoint));
      }

   /* off diagonal elements */
   for(i = 1; i < rem_tpoints - 1; i++){
      for(j = 0; j < rem_tpoints - i; j++){
         gsl_matrix_set(aif_m, j + i, j, ((gsl_vector_get(aif->AIF, start_tpoint + (i - 1)))
                                          + (4 * gsl_vector_get(aif->AIF, start_tpoint + i))
                                          + (gsl_vector_get(aif->AIF, start_tpoint + (i + 1))))
                        / 6);
         }
      }

   /* do the last element in the corner */
   gsl_matrix_set(aif_m, rem_tpoints - 1, 0, gsl_vector_get(aif->AIF, aif->AIF->size - 1));

   /* multiply by the TR */
   for(i = 1; i < rem_tpoints - 1; i++){
      for(j = 0; j < rem_tpoints - i; j++){
         gsl_matrix_set(aif_m, i, j, gsl_matrix_get(aif_m, i, j) * tr);
         }
      }
   }


/* reads in a Art_IF from a file                        */
Status input_art_if(const char *aif_file, Art_IF * aif)
{
   MINC_Vector *tmp;
   Status   status;
   int      i;

   tmp = new_MINC_Vector(0);

   if(input_MINC_Vector(aif_file, tmp) == OK){
      aif->AIF = gsl_vector_alloc(tmp->size);
      for(i = 0; i < tmp->size; i++){
         gsl_vector_set(aif->AIF, i, tmp->V[i]);
         }

      status = OK;
      }
   else{
      status = ERROR;
   };

   free_MINC_Vector(tmp);

   return (status);
   }


void print_vector(const char *name, gsl_vector * v)
{
   size_t   i;

   fprintf(stdout, "%s[%d]: ", name, v->size);
   for(i = 0; i < v->size; i++){
      fprintf(stdout, "%g ", gsl_vector_get(v, i));
      }
   fprintf(stdout, "\n");
   }

void setup_math_data(Math_Data ** md)
{

   int      i;

   /* setup the kernel */
   ((*md))->kernel = gsl_vector_alloc(11);
   gsl_vector_set((*md)->kernel, 0, -0.00390625);
   gsl_vector_set((*md)->kernel, 1, -0.0195312);
   gsl_vector_set((*md)->kernel, 2, -0.0195312);
   gsl_vector_set((*md)->kernel, 3, 0.078125);
   gsl_vector_set((*md)->kernel, 4, 0.273438);
   gsl_vector_set((*md)->kernel, 5, 0.382812);
   gsl_vector_set((*md)->kernel, 6, 0.273438);
   gsl_vector_set((*md)->kernel, 7, 0.078125);
   gsl_vector_set((*md)->kernel, 8, -0.0195312);
   gsl_vector_set((*md)->kernel, 9, -0.0195312);
   gsl_vector_set((*md)->kernel, 10, -0.00390625);

   /* filter the AIF */
   if(((*md))->filter){
      convolve_vector_avg((*md)->aif->AIF, (*md)->kernel);
      }

   /* Calculate the AIF matrix from the raw AIF */
   (*md)->svd_u = gsl_matrix_calloc((*md)->rem_tpoints, (*md)->rem_tpoints);
   calc_aif_matrix((*md)->aif, (*md)->start_tpoint, (*md)->rem_tpoints,
                   (*md)->tr, (*md)->te, (*md)->svd_u);

   /* setup SVD fit matrices */
   (*md)->svd_v = gsl_matrix_alloc((*md)->rem_tpoints, (*md)->rem_tpoints);
   (*md)->svd_s = gsl_vector_alloc((*md)->rem_tpoints);
   (*md)->svd_work = gsl_vector_alloc((*md)->rem_tpoints);

   gsl_linalg_SV_decomp((*md)->svd_u, (*md)->svd_v, (*md)->svd_s, (*md)->svd_work);
   for(i = 0; i < (*md)->rem_tpoints; i++){
      if(gsl_vector_get((*md)->svd_s, i) < (*md)->svd_tol){
         gsl_vector_set((*md)->svd_s, i, 0.0);
         }
      }

   /* time curve area and arrival time of AIF */
   time_curve_area((*md)->aif->AIF, (*md)->start_tpoint, (*md)->tr,
                   &(*md)->aif->area, &(*md)->aif->arrival_time, 1);
   }
