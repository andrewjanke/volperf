/* perf_util.h */

#include <volume_io.h>
#include <voxel_loop.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
/* #include <gsl/gsl_linalg.h> */
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>

/* Structure for Art_IF information */
typedef struct {
   gsl_vector *AIF;
   double   baseline;
   double   area;
   double   arrival_time;
} Art_IF;

/* Structure for math information */
typedef struct {
   Art_IF  *aif;
   size_t   start_tpoint;
   size_t   rem_tpoints;

   /* kernel to filter with */
   gsl_vector *kernel;

   /* SVD fit bits */
   gsl_matrix *svd_u;
   gsl_matrix *svd_v;
   gsl_vector *svd_s;
   gsl_vector *svd_work;

   double   tr;
   double   te;
   int      filter;
   double   svd_tol;
   double   cutoff;
   int      calc_bat;
} Math_Data;

/* voxel_loop function */
void     do_math(void *caller_data, long num_voxels, int input_num_buffers, int input_vector_length,
                 double *input_data[], int output_num_buffers, int output_vector_length,
                 double *output_data[], Loop_Info * loop_info);

/* reads an allocated Art_IF from a file (pointer) */
Status   input_art_if(const char *aif_file, Art_IF * aif);

/* convolves a vector using zero padding */
void     convolve_vector_zero(gsl_vector * data, gsl_vector * kernel);

/* convolves a vector using averaging of end points */
void     convolve_vector_avg(gsl_vector * data, gsl_vector * kernel);

/* calculate the AIF matrix */
void     calc_aif_matrix(Art_IF * aif, size_t start_tpoint, size_t rem_tpoints,
                         double tr, double te, gsl_matrix * aif_m);

/* calculate the area under a time curve */
void time_curve_area(gsl_vector * tcurve, size_t start_tpoint, double tr,
                     double *area, double *bat, int flg);

/* debugging function to print a vector */
void     print_vector(const char *name, gsl_vector * v);

/* setup the math data structure */
void     setup_math_data(Math_Data ** md);
