/* perf_util.h */

#ifndef PERF_UTIL_H
#define PERF_UTIL_H

#include <float.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics_double.h>

#include <volume_io.h>

#define SQR(x) ((x) * (x))

/************************************************************/

typedef enum { SHIFT_NONE = 0, SHIFT_AIF, SHIFT_AIF_EXACT, SHIFT_CONC, SHIFT_CIRC
   } SHIFT_enum;
typedef enum { BAT_NONE = 0, BAT_GAMMA, BAT_SLOPE, BAT_CUTOFF, BAT_MIN_DIST } BAT_enum;
typedef enum { GAM_BAT, GAM_FAC, GAM_ALPHA, GAM_BETA, GAM_NUM } GAMMA_enum;
typedef enum { SVD_TOL, SVD_OI } SVD_enum;

typedef enum { OUT_CBF, OUT_CBV, OUT_MTT, NUM_OUT_BASIC } OUT_BASIC_ENUM;
typedef enum { OUT_CHI, NUM_OUT_CHI } OUT_CHI_ENUM;
typedef enum { OUT_BAT, OUT_FAC, OUT_ALPHA, OUT_BETA, NUM_OUT_GAMMA } OUT_GAMMA_ENUM;
typedef enum { OUT_ARRIVAL, NUM_OUT_ARRIVAL } OUT_ARRIVAL_ENUM;
typedef enum { OUT_DELAY, NUM_OUT_DELAY } OUT_DELAY_ENUM;

/************************************************************/

/* Structure for BAT information */
typedef struct {
   BAT_enum type;
   double   cutoff;
   double   min_dist;

   double   gparm[GAM_NUM];
   } Bat_info;

/* Structure for Art_IF information */
typedef struct {
   gsl_vector *signal;
   gsl_vector *conc;
   double   baseline;
   double   area;
   double   arrival_time;

   Bat_info bat;
   } Art_IF;

/* Structure for math information */
typedef struct {

   /* arterial input function */
   Art_IF  *aif;
   int      aif_fit;
   size_t   Nstart;
   size_t   Nrest;
   size_t   length;

   /* number of input data files */
   size_t   n_datafiles;

   /* mask information */
   int      mask;
   int      mask_idx;

   /* kernel to filter with */
   gsl_vector *kernel;

   /* SVD fit bits */
   SVD_enum svd_type;
   gsl_matrix *svd_u;
   gsl_matrix *svd_v;
   gsl_vector *svd_s;
   gsl_vector *svd_s_cp;
   gsl_vector *svd_work;
   double   svd_tol;
   double   svd_oi;

   /* voxel concentration */
   gsl_vector *signal;
   gsl_vector *conc;
   double   gparm[GAM_NUM];
   gsl_vector *residue;
   int      conc_fit;

   /* tracer delay */
   SHIFT_enum shift_type;
   Bat_info vox_bat;

   double   tr;
   double   te;
   int      filter;
   double   cutoff;
   double   area;
   double   normalization;
   double   dosage;

   int      output_chi;
   int      output_gamma;
   int      output_arrival;
   int      output_delay;
   int      output_residue;
   } Math_Data;

/************************************************************/

void     convolve_vector(gsl_vector * data, gsl_vector * kernel);

/* convolves a vector using zero padding */
void     convolve_vector_zero(gsl_vector * data, gsl_vector * kernel);

/* convolves a vector using averaging of end points */
void     convolve_vector_avg(gsl_vector * data, gsl_vector * kernel);

gsl_vector *create_filter_kernel(void);

/* reads an allocated Art_IF from a file (pointer) */
Status   input_art_if(const char *aif_file, Art_IF * aif);

/* debugging function to print a vector */
void     print_vector(const char *name, gsl_vector * v);
void     print_matrix(const char *name, gsl_matrix * m);

void     SignalToConc(gsl_vector * s, gsl_vector * c, double baseline, double te);
void     ConvertParm(const gsl_vector * v, double *gparm);
void     eval_gamma(gsl_vector * v, double *gparm, double step);

/*
void fit_gamma(gsl_vector *log_conc, gsl_vector *t, double bat,
 double *gamma, gsl_matrix *mwork, gsl_matrix *covwork,
 gsl_multifit_linear_workspace *fwork);
*/
double   vector_distance(gsl_vector * a, gsl_vector * b);
double   vector_sum(gsl_vector * v, int a, int b);
double   bolus_arrival_time(gsl_vector * conc, double tr, Bat_info * b);

/* calculate the AIF matrix */
void     CalcAifMatrix(gsl_vector * v, gsl_matrix * m, double fac, int Nstart, int Nrest);
void     CalcCircAifMatrix(gsl_vector * v, gsl_matrix * m, double fac, int Nstart,
                           int Nrest);

/* calculate svd oscillation index (as a double) */
double   svd_oscillation_index(gsl_vector * v);

#endif
