/* volperf.c                                                                 */
/*                                                                           */
/* The MINC perfusion calculator                                             */
/*                                                                           */
/* Mark Griffin - mark.griffin@cmr.uq.edu.au                                 */
/* Andrew Janke - rotor@cmr.uq.edu.au                                        */
/* Center for Magnetic Resonance                                             */
/* University of Queensland                                                  */
/*                                                                           */
/* Copyright Andrew Janke & Mark Griffin The University of Queensland.       */
/* Permission to use, copy, modify, and distribute this software and its     */
/* documentation for any purpose and without fee is hereby granted,          */
/* provided that the above copyright notice appear in all copies.  The       */
/* author and the University of Queensland make no representations about the */
/* suitability of this software for any purpose.  It is provided "as is"     */
/* without express or implied warranty.                                      */
/*                                                                           */
/* Thu Dec  6 22:33:12 EST 2001 - initial version (that didn't work)         */
/* Wed Jul  3 15:59:46 EST 2002 - major revisions (to make it work)          */
/* Mon Jul 28 18:15:24 EST 2003 - added bolus delay correction               */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <ParseArgv.h>
#include <time_stamp.h>
#include <voxel_loop.h>

#include "perf_util.h"
#include "gamma_fit.h"

#define DEFAULT_BOOL -1

/* function prototypes */
void     do_math(void *caller_data, long num_voxels, int input_num_buffers,
                 int input_vector_length, double *input_data[], int output_num_buffers,
                 int output_vector_length, double *output_data[], Loop_Info * loop_info);

/************************************************************/
char    *SHIFT_names[] = { "shift_none", "shift_aif", "shift_aif_exact", "shift_conc" };
char    *BAT_names[] = { "bat_none", "bat_gamma", "bat_slope", "bat_cutoff" };
char    *outfile_basic[] = { "CBF", "CBV", "MTT", NULL };
char    *outfile_chi[] = { "chi", NULL };
char    *outfile_gamma[] = { "bat", "fac", "alpha", "beta", NULL };
char    *outfile_arrival[] = { "arrival", NULL };
char    *outfile_delay[] = { "delay", NULL };
char     outfile_residue[] = "residue_";

/************************************************************/

/* Argument variables and table */
int      verbose = FALSE;
int      quiet = FALSE;
int      clobber = FALSE;
int      max_buffer = 4 * 1024;
int      copy_header = DEFAULT_BOOL;
int      is_signed = FALSE;
nc_type  dtype = NC_UNSPECIFIED;

char    *aif_file = NULL;
int      aif_fit = FALSE;
int      conc_fit = FALSE;

SHIFT_enum shift_type = SHIFT_NONE;

int      aif_bat_gamma = FALSE;
int      aif_bat_slope = FALSE;
double   aif_bat_cutoff = -1;
int      vox_bat_gamma = FALSE;
int      vox_bat_slope = FALSE;
double   vox_bat_cutoff = -1;
int      bat_gamma = FALSE;
int      bat_slope = FALSE;
double   bat_cutoff = -1;
double   tr = -1;
double   te = -1;
int      filter = FALSE;
double   svd_tol = 0.2;
double   svd_oi = -1;
double   cutoff = 0.4;
double   normalization = 1.0;
double   dosage = 0.2;
double   tstart = -1.0;
int      output_chi = FALSE;
int      output_gamma = FALSE;
int      output_arrival = FALSE;
int      output_delay = FALSE;
int      output_residue = FALSE;

//Math_Data md = { 
//   NULL, 0, 0, 0,
//   NULL,
//   NULL, NULL, NULL, NULL, 0.2,
//   NULL, NULL, NULL, NULL, FALSE,
//   SHIFT_NONE, BAT_NONE, 0.0,
//   -1, -1, FALSE, 0.4, 0.0, 1.0, 0.2,
//   FALSE,  FALSE, FALSE, FALSE, FALSE
//   }

/************************************************************/

ArgvInfo argTable[] = {
   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "General options:"},
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
    "Print out extra information."},
   {"-quiet", ARGV_CONSTANT, (char *)TRUE, (char *)&quiet,
    "Be very quiet."},
   {"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&clobber,
    "Clobber existing files."},
   {"-max_buffer", ARGV_INT, (char *)1, (char *)&max_buffer,
    "maximum size of buffers (in kbytes)"},

   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, "\nOutfile Options:"},
   {"-copy_header", ARGV_CONSTANT, (char *)TRUE, (char *)&copy_header,
    "Copy header from 1st file (default)"},
   {"-nocopy_header", ARGV_CONSTANT, (char *)FALSE, (char *)&copy_header,
    "Don't copy header from 1st file."},
   {"-filetype", ARGV_CONSTANT, (char *)NC_UNSPECIFIED, (char *)&dtype,
    "Use data type of 1st file (default)"},
   {"-byte", ARGV_CONSTANT, (char *)NC_BYTE, (char *)&dtype,
    "Write out byte data."},
   {"-short", ARGV_CONSTANT, (char *)NC_SHORT, (char *)&dtype,
    "Write out short integer data."},
   {"-long", ARGV_CONSTANT, (char *)NC_LONG, (char *)&dtype,
    "Write out long integer data."},
   {"-float", ARGV_CONSTANT, (char *)NC_FLOAT, (char *)&dtype,
    "Write out single-precision data."},
   {"-double", ARGV_CONSTANT, (char *)NC_DOUBLE, (char *)&dtype,
    "Write out double-precision data."},
   {"-signed", ARGV_CONSTANT, (char *)TRUE, (char *)&is_signed,
    "Write signed integer data."},
   {"-unsigned", ARGV_CONSTANT, (char *)FALSE, (char *)&is_signed,
    "Write unsigned integer data."},
   {"-chi", ARGV_CONSTANT, (char *)TRUE, (char *)&output_chi,
    "Ouput the Chi-squared fit"},
   {"-gamma", ARGV_CONSTANT, (char *)TRUE, (char *)&output_gamma,
    "Ouput the Gamma function variables for each voxel"},
   {"-arrival", ARGV_CONSTANT, (char *)TRUE, (char *)&output_arrival,
    "Output the arrival time (according to the method chosen)"},
   {"-delay", ARGV_CONSTANT, (char *)TRUE, (char *)&output_delay,
    "Output the delay time (similar to arrival time)"},
   {"-residue", ARGV_CONSTANT, (char *)TRUE, (char *)&output_residue,
    "Ouput the Residue Function"},

   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, "\nPerfusion Options:"},
   {"-aif", ARGV_STRING, (char *)1, (char *)&aif_file,
    "<art_if.aif> File containing the Arterial input Function"},
   {"-aif_fit", ARGV_CONSTANT, (char *)TRUE, (char *)&aif_fit,
    "Fit a gamma variate function to the aif"},
   {"-conc_fit", ARGV_CONSTANT, (char *)TRUE, (char *)&conc_fit,
    "Fit a gamma variate function to the voxel"},
   {"-start", ARGV_FLOAT, (char *)1, (char *)&tstart,
    "Time at the start of the aif bolus"},

/*
 *	shifting concentration profiles to remove delay
 */
   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, "\nShift types:"},
   {"-shift_none", ARGV_CONSTANT, (char *)SHIFT_NONE, (char *)&shift_type,
    "Don't do any shifting (Default)"},
   {"-shift_aif", ARGV_CONSTANT, (char *)SHIFT_AIF, (char *)&shift_type,
    "Shift the aif to the voxel concentration (resolution of tr)"},
   {"-shift_aif_exact", ARGV_CONSTANT, (char *)SHIFT_AIF_EXACT,
    (char *)&shift_type,
    "Shift the aif to the voxel concentration (exact resolution)"},
   {"-shift_conc", ARGV_CONSTANT, (char *)SHIFT_CONC, (char *)&shift_type,
    "Shift the concentration to the aif"},
   {"-shift_circ", ARGV_CONSTANT, (char *)SHIFT_CIRC, (char *)&shift_type,
    "Correct for delay (shift) using Circular deconvolution"},

/*
 *	methods for calculating arrival time for both aif and voxel
 */
   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, "\nBAT/delay types:"},
   {"-aif_bat_gamma", ARGV_CONSTANT, (char *)TRUE, (char *)&aif_bat_gamma,
    "The aif arrival is the bat of the gamma variate function (Default)"},
   {"-aif_bat_slope", ARGV_CONSTANT, (char *)TRUE, (char *)&aif_bat_slope,
    "The aif arrival is the slope of the concentration profile"},
   {"-aif_bat_cutoff", ARGV_FLOAT, (char *)1, (char *)&aif_bat_cutoff,
    "The aif arrival is a percentage of the peak concentration"},
   {"-vox_bat_gamma", ARGV_CONSTANT, (char *)TRUE, (char *)&vox_bat_gamma,
    "The voxel arrival is the bat of the gamma variate function (Default)"},
   {"-vox_bat_slope", ARGV_CONSTANT, (char *)TRUE, (char *)&vox_bat_slope,
    "The voxel arrival is the slope of the concentration profile"},
   {"-vox_bat_cutoff", ARGV_FLOAT, (char *)1, (char *)&vox_bat_cutoff,
    "The voxel arrival is a percentage of the peak concentration"},

   {"-bat_gamma", ARGV_CONSTANT, (char *)TRUE, (char *)&bat_gamma,
    "Synonym for '-aif_bat_gamma -vox_bat_gamma' (Default)"},
   {"-bat_slope", ARGV_CONSTANT, (char *)TRUE, (char *)&bat_slope,
    "Synonym for '-aif_bat_slope -vox_bat_slope'"},
   {"-bat_cutoff", ARGV_FLOAT, (char *)1, (char *)&bat_cutoff,
    "Synonym for '-aif_bat_cutoff -vox_bat_cutoff'"},

   {"-tr", ARGV_FLOAT, (char *)1, (char *)&tr,
    "TR of the perfusion experiment (default: read this from the input files)"},
   {"-te", ARGV_FLOAT, (char *)1, (char *)&te,
    "TE of the perfusion experiment (default: read this from the input files)"},
   {"-filter", ARGV_CONSTANT, (char *)TRUE, (char *)&filter,
    "Filter the data before fitting (low-pass gaussian)"},
   {"-svd_tolerance", ARGV_FLOAT, (char *)1, (char *)&svd_tol,
    "[0..1] Tolerance to use for the SVD fitting"},
   {"-svd_oi", ARGV_FLOAT, (char *)1, (char *)&svd_oi,
    "[0..4] Oscillation index for SVD fitting (Suggest 0.09 -- Ona Wu MRM 2003)"},

   {"-cutoff", ARGV_FLOAT, (char *)1, (char *)&cutoff,
    "[0..1] Minimum signal change required at a voxel for a perfusion calculation"},
   {"-normalization", ARGV_FLOAT, (char *)1, (char *)&normalization,
    "Normalization factor, default 0.87 for SE-EPI in humans"},
   {"-dosage", ARGV_FLOAT, (char *)1, (char *)&dosage,
    "dosage, default 0.2 mmol/kg"},

   {NULL, ARGV_HELP, NULL, NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
   };

/************************************************************/

int main(int argc, char *argv[])
{
   char   **infiles;
   char   **outfiles;
   int      n_infiles, n_outfiles;
   int      nbasic, nchi, ngamma, narrival, ndelay, nresidue;
   char    *arg_string;
   Loop_Options *loop_opts;
   double   fac;
   size_t   i, j;
   BAT_enum aif_bat_type = BAT_GAMMA;
   int      n_bat, n_aif_bat, n_vox_bat;
   Math_Data *md;

/* 
 *	Save time stamp and args 
 */
   arg_string = time_stamp(argc, argv);

/* 
 *	Get arguments 
 */
   if(ParseArgv(&argc, argv, argTable, 0) || (argc < 2)){
      fprintf(stdout,
              "\nUsage: %s [options] <in1.mnc> <in2.mnc> ... <out_base>\n", argv[0]);
      fprintf(stdout, "       %s -help\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }
   md = (Math_Data *) malloc(sizeof(Math_Data));

/* 
 *	check for an input TR and TE 
 */
   if(tr == -1 || te == -1){
      fprintf(stdout, "%s: You need to input a tr and a te\n", argv[0]);
      fprintf(stdout, "%s -help  for more info\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }
   md->tr = tr;
   md->te = te;

   /* check input parameters */

   /* SVD options */
   md->svd_oi = svd_oi;
   md->svd_tol = svd_tol;
   if(md->svd_oi != -1){
      if(md->svd_oi > 4 || md->svd_oi < 0){
         fprintf(stderr, "%s: -svd_oi must be in the range 0..4\n\n", argv[0]);
         exit(EXIT_FAILURE);
         }

      md->svd_type = SVD_OI;
      }
   else {
      if(md->svd_tol > 1 || md->svd_tol < 0){
         fprintf(stderr, "%s: -svd_tol must be in the range 0..1\n\n", argv[0]);
         exit(EXIT_FAILURE);
         }
      md->svd_type = SVD_TOL;
      }

/* 
 *	Initialise and read in the AIF 
 */
   if(access(aif_file, F_OK) != 0){
      fprintf(stderr, "%s: Couldn't find Arterial input file: %s\n", argv[0], aif_file);
      exit(EXIT_FAILURE);
      }
   md->aif = (Art_IF *) malloc(sizeof(Art_IF));
   if(input_art_if(aif_file, md->aif) != OK){
      fprintf(stdout, "%s: Died whilst reading in AIF file: %s\n", argv[0], aif_file);
      exit(EXIT_FAILURE);
      }
   md->aif_fit = aif_fit || shift_type == SHIFT_AIF_EXACT;
   md->conc_fit = conc_fit;

/*
 *	shift type
 */
   md->shift_type = shift_type;

/*
 *	check the arrival time measures
 */
   n_bat = (bat_gamma + bat_slope + (bat_cutoff > 0));
   if(n_bat > 1){
      fprintf(stdout,
              "%s: Only one arrival method required for aif and voxel\n", argv[0]);
      exit(EXIT_FAILURE);
      }

   n_aif_bat = (aif_bat_gamma + aif_bat_slope + (aif_bat_cutoff > 0));
   if(n_aif_bat > 1){
      fprintf(stdout, "%s: Only one arrival method required for aif\n", argv[0]);
      exit(EXIT_FAILURE);
      }

   n_vox_bat = (vox_bat_gamma + vox_bat_slope + (vox_bat_cutoff > 0));
   if(n_vox_bat > 1){
      fprintf(stdout, "%s: Only one arrival method required for voxel\n", argv[0]);
      exit(EXIT_FAILURE);
      }

   if(n_aif_bat != n_vox_bat){
      fprintf(stdout, "%s: Must either specify arrival methods for both\n", argv[0]);
      fprintf(stdout, "  the aif and the voxel or for neither of them\n");
      exit(EXIT_FAILURE);
      }

   if(n_bat + n_aif_bat > 1){
      fprintf(stdout, "%s: Must specify arrival methods for the\n", argv[0]);
      fprintf(stdout, "  aif and the voxel either together or\n");
      fprintf(stdout, "  separately\n");
      exit(EXIT_FAILURE);
      }

   if((md->shift_type != SHIFT_NONE && md->shift_type != SHIFT_CIRC)
      && (n_bat + n_aif_bat == 0)){
      fprintf(stdout, "%s: Must specify an arrival marker\n", argv[0]);
      fprintf(stdout, "  order to account for tracer delay\n");
      exit(EXIT_FAILURE);
      }

   md->vox_bat_cutoff = vox_bat_cutoff;
   if(bat_gamma){
      aif_bat_type = BAT_GAMMA;
      md->vox_bat_type = BAT_GAMMA;
      }
   else if(bat_slope){
      aif_bat_type = BAT_SLOPE;
      md->vox_bat_type = BAT_SLOPE;
      }
   else if(bat_cutoff > 0.0){
      aif_bat_type = BAT_CUTOFF;
      md->vox_bat_type = BAT_CUTOFF;
      aif_bat_cutoff = bat_cutoff;
      md->vox_bat_cutoff = bat_cutoff;
      }
   else {
      md->vox_bat_type = BAT_GAMMA;
      }

   if(aif_bat_gamma){
      aif_bat_type = BAT_GAMMA;
      }
   else if(aif_bat_slope){
      aif_bat_type = BAT_SLOPE;
      }
   else if(aif_bat_cutoff > 0.0){
      aif_bat_type = BAT_CUTOFF;
      }

   if(vox_bat_gamma){
      md->vox_bat_type = BAT_GAMMA;
      }
   else if(vox_bat_slope){
      md->vox_bat_type = BAT_SLOPE;
      }
   else if(vox_bat_cutoff > 0.0){
      md->vox_bat_type = BAT_CUTOFF;
      md->vox_bat_cutoff = vox_bat_cutoff;
      }

   if((aif_bat_cutoff > 1.0) || (md->vox_bat_cutoff > 1.0)){
      fprintf(stdout, "%s: Arrival cutoffs must be between 0 and 1\n", argv[0]);
      exit(EXIT_FAILURE);
      }

/*
 *	check the arrival outputs
 */
   md->output_gamma = output_gamma;
   md->output_arrival = output_arrival;
   md->output_delay = output_delay;

/* 
 *	Set up math data structure 
 */
   md->filter = filter;
   md->cutoff = cutoff;
   md->normalization = normalization;
   md->dosage = dosage;
   md->output_chi = output_chi;
   md->output_residue = output_residue;

/*
 *	gamma variate fit to the arterial input function
 */
   md->kernel = create_filter_kernel();
   if(md->filter){
      convolve_vector_avg(md->aif->signal, md->kernel);
      }
   if(tstart < -0.5){
      TimeCurveAreaAif(md->aif->signal, md->aif->conc, md->tr, md->te,
                       &(md->aif->area), md->aif->gparm, &(md->Nstart), &(md->Nrest));
      md->aif->baseline = 0.0;
      }
   else {
      md->Nstart = tstart / md->tr;
      md->Nrest = md->aif->signal->size - md->Nstart;
      }
   for(i = 0; i < md->Nstart; i++){
      md->aif->baseline += gsl_vector_get(md->aif->signal, i);
      }
   md->aif->baseline /= md->Nstart;
   SignalToConc(md->aif->signal, md->aif->conc, md->aif->baseline, md->te);
   md->aif->arrival_time = bolus_arrival_time(md->aif->conc, md->tr,
                                              md->aif->gparm, aif_bat_type,
                                              aif_bat_cutoff);
   if(tstart > -0.5){
      TimeCurveArea(md->aif->conc, md->tr, &(md->aif->area), md->aif->gparm);
      }
   md->aif->arrival_time -= md->Nstart * md->tr;
   md->aif->gparm[GAM_BAT] -= md->Nstart * md->tr;

/* 
 *	set up in and outfile counts 
 */
   n_infiles = argc - 2;
   infiles = &argv[1];
   nbasic = NUM_OUT_BASIC;
   nchi = output_chi * NUM_OUT_CHI;
   ngamma = output_gamma * NUM_OUT_GAMMA;
   narrival = output_arrival * NUM_OUT_ARRIVAL;
   ndelay = output_delay * NUM_OUT_DELAY;
   nresidue = output_residue * md->Nrest;
   n_outfiles = nbasic + nchi + ngamma + narrival + ndelay + nresidue;

/* 
 *	check the infiles match and exist 
 */
   if(n_infiles != (int)md->aif->signal->size){
      fprintf(stdout, "%s: %d infiles != %d AIF time-points\n", argv[0],
              n_infiles, md->aif->signal->size);
      exit(EXIT_FAILURE);
      }

   if(verbose){
      fprintf(stdout, "\n==== Infiles / AIF ====\n");
      }
   for(i = 0; i < (size_t) n_infiles; i++){
      if(verbose){
         fprintf(stdout, "[%02d]: %s  |  %g\n", i, infiles[i],
                 gsl_vector_get(md->aif->signal, i));
         }
      if(access(infiles[i], F_OK) != 0){
         fprintf(stdout, "%s: Couldn't find %s\n", argv[0], infiles[i]);
         exit(EXIT_FAILURE);
         }
      }

/* 
 *	create and check for the outfiles 
 */
   outfiles = (char **)malloc(sizeof(char *) * n_outfiles);
   if(verbose){
      fprintf(stdout, "\n==== Outfiles ====\n");
      }
   j = 0;
   for(i = 0; i < (size_t) (nbasic); i++){
      outfiles[j] = (char *)malloc(sizeof(char)
                                   * (strlen(argv[n_infiles + 1]) +
                                      strlen(outfile_basic[i]) + 5));
      sprintf(outfiles[j], "%s.%s.mnc", argv[n_infiles + 1], outfile_basic[i]);
      j++;
      }
   for(i = 0; i < (size_t) (nchi); i++){
      outfiles[j] = (char *)malloc(sizeof(char)
                                   * (strlen(argv[n_infiles + 1]) +
                                      strlen(outfile_chi[i]) + 5));
      sprintf(outfiles[j], "%s.%s.mnc", argv[n_infiles + 1], outfile_chi[i]);
      j++;
      }
   for(i = 0; i < (size_t) (ngamma); i++){
      outfiles[j] = (char *)malloc(sizeof(char)
                                   * (strlen(argv[n_infiles + 1]) +
                                      strlen(outfile_gamma[i]) + 5));
      sprintf(outfiles[j], "%s.%s.mnc", argv[n_infiles + 1], outfile_gamma[i]);
      j++;
      }
   for(i = 0; i < (size_t) (narrival); i++){
      outfiles[j] = (char *)malloc(sizeof(char)
                                   * (strlen(argv[n_infiles + 1]) +
                                      strlen(outfile_arrival[i]) + 5));
      sprintf(outfiles[j], "%s.%s.mnc", argv[n_infiles + 1], outfile_arrival[i]);
      j++;
      }
   for(i = 0; i < (size_t) (ndelay); i++){
      outfiles[j] = (char *)malloc(sizeof(char)
                                   * (strlen(argv[n_infiles + 1]) +
                                      strlen(outfile_delay[i]) + 5));
      sprintf(outfiles[j], "%s.%s.mnc", argv[n_infiles + 1], outfile_delay[i]);
      j++;
      }
   if(output_residue){
      for(i = 0; (i < (size_t) (10)) && (i < (size_t) (nresidue)); i++){
         outfiles[j] = (char *)malloc(sizeof(char)
                                      * (strlen(argv[n_infiles + 1]) +
                                         strlen(outfile_residue) + 8));
         sprintf(outfiles[j], "%s.%s0%d.mnc", argv[n_infiles + 1], outfile_residue, i);
         j++;
         }
      for(i = 10; i < (size_t) (nresidue); i++){
         outfiles[j] = (char *)malloc(sizeof(char)
                                      * (strlen(argv[n_infiles + 1]) +
                                         strlen(outfile_residue) + 8));
         sprintf(outfiles[j], "%s.%s%d.mnc", argv[n_infiles + 1], outfile_residue, i);
         j++;
         }
      }

   for(i = 0; i < n_outfiles; i++){
      if(verbose){
         fprintf(stdout, "[%02d]: %s\n", i, outfiles[i]);
         }
      if(access(outfiles[i], F_OK) == 0 && !clobber){
         fprintf(stdout, "%s: %s exists, use -clobber to overwrite\n",
                 argv[0], outfiles[i]);
         exit(EXIT_FAILURE);
         }
      }

   if(verbose){
      fprintf(stdout, "====Parameters======\n");
      fprintf(stdout, "| aif->baseline:     %g\n", md->aif->baseline);
      fprintf(stdout, "| aif->area:         %g\n", md->aif->area);
      fprintf(stdout, "| aif->arrival_time: %g\n", md->aif->arrival_time);

      fprintf(stdout, "| aif->signal           ");
      for(i = 0; i < md->aif->signal->size; i++){
         if((i % 10) == 0 && i != 0){
            fprintf(stdout, "\n|                    ");
            }
         fprintf(stdout, "%6.3f ", gsl_vector_get(md->aif->signal, i));
         }
      fprintf(stdout, "\n");

      fprintf(stdout, "| Nstart:            %d\n", md->Nstart);
      fprintf(stdout, "| Nrest:             %d\n", md->Nrest);
      fprintf(stdout, "| svd_tol:           %g\n", md->svd_tol);
      fprintf(stdout, "| svd_oi:            %g\n", md->svd_oi);
      fprintf(stdout, "| svd_type:          %d\n", md->svd_type);

      fprintf(stdout, "| shift_type:        %s\n", SHIFT_names[md->shift_type]);

      fprintf(stdout, "| TR:                %g\n", md->tr);
      fprintf(stdout, "| TE:                %g\n", md->te);
      fprintf(stdout, "| filter:            %s\n", (md->filter) ? "TRUE" : "FALSE");
      fprintf(stdout, "| cutoff:            %g\n", md->cutoff);
      }

/*
 *	Create the aif matrix
 */
   if(md->aif_fit){
      eval_gamma(md->aif->conc, md->aif->gparm, md->tr);
      }
   else {
      for(i = 0; i < md->Nrest; i++){
         gsl_vector_set(md->aif->conc, i, gsl_vector_get(md->aif->conc, i + md->Nstart));
         }
      }

   switch (md->shift_type){
   case SHIFT_NONE:

   case SHIFT_CONC:
      md->length = md->Nrest;
      md->svd_u = gsl_matrix_alloc(md->length, md->length);
      md->svd_v = gsl_matrix_alloc(md->length, md->length);
      md->svd_s = gsl_vector_alloc(md->length);
      md->svd_s_cp = md->svd_s;
      md->svd_work = gsl_vector_alloc(md->length);
      CalcAifMatrix(md->aif->conc, md->svd_u, md->tr, 0, md->length);
      gsl_linalg_SV_decomp(md->svd_u, md->svd_v, md->svd_s, md->svd_work);
      fac = md->svd_tol * gsl_vector_max(md->svd_s);
      for(i = 0; i < md->length; i++){
         if(gsl_vector_get(md->svd_s, i) < fac){
            gsl_vector_set(md->svd_s, i, 0.0);
            }
         }
      md->residue = gsl_vector_alloc(md->length);
      md->signal = gsl_vector_alloc(md->aif->signal->size);
      md->conc = gsl_vector_alloc(md->length);
      break;

   case SHIFT_AIF:
   case SHIFT_AIF_EXACT:
      md->length = md->Nrest;
      md->signal = gsl_vector_alloc(md->aif->signal->size);
      md->conc = gsl_vector_alloc(md->length);
      break;

   case SHIFT_CIRC:
      md->length = 2 * md->Nrest;
      md->svd_u = gsl_matrix_alloc(md->length, md->length);
      md->svd_v = gsl_matrix_alloc(md->length, md->length);
      md->svd_s = gsl_vector_alloc(md->length);
      md->svd_s_cp = gsl_vector_alloc(md->length);
      md->svd_work = gsl_vector_alloc(md->length);
      CalcCircAifMatrix(md->aif->conc, md->svd_u, md->tr, 0, md->Nrest);
      gsl_linalg_SV_decomp(md->svd_u, md->svd_v, md->svd_s, md->svd_work);
      fac = md->svd_tol * gsl_vector_max(md->svd_s);
      for(i = 0; i < md->length; i++){
         if(gsl_vector_get(md->svd_s, i) < fac){
            gsl_vector_set(md->svd_s, i, 0.0);
            }
         }
      md->residue = gsl_vector_alloc(md->length);
      md->signal = gsl_vector_alloc(md->aif->signal->size);
      md->conc = gsl_vector_alloc(md->length);
      break;

   default:
      fprintf(stderr, "ERROR - shift method not implemented\n");
      exit(EXIT_FAILURE);
      }

/* 
 *	Set up and run the Voxel Loop 
 */
   loop_opts = create_loop_options();
   set_loop_copy_all_header(loop_opts, copy_header);
   set_loop_verbose(loop_opts, !quiet);
   set_loop_clobber(loop_opts, clobber);
   set_loop_datatype(loop_opts, dtype, is_signed, 0.0, 0.0);
   set_loop_buffer_size(loop_opts, (long)1024 * max_buffer);

   voxel_loop(n_infiles, infiles, n_outfiles, outfiles, arg_string,
              loop_opts, do_math, (void *)md);

/* 
 *	be tidy 
 */
   if((md->shift_type == SHIFT_CONC) || (md->shift_type == SHIFT_NONE)
      || md->shift_type == SHIFT_CIRC){
      gsl_matrix_free(md->svd_u);
      gsl_matrix_free(md->svd_v);
      gsl_vector_free(md->svd_s);
      if(md->shift_type == SHIFT_CIRC){
         gsl_vector_free(md->svd_s_cp);
         }
      gsl_vector_free(md->svd_work);
      gsl_vector_free(md->conc);
      gsl_vector_free(md->residue);
      }
   for(i = 0; i < n_outfiles; i++){
      free(outfiles[i]);
      }
   free(outfiles);
   free_loop_options(loop_opts);
   gsl_vector_free(md->conc);
   gsl_vector_free(md->signal);
   gsl_vector_free(md->aif->signal);
   gsl_vector_free(md->aif->conc);
   free(md->aif);
   free(md);

   return EXIT_SUCCESS;
   }

/************************************************************/
void do_math(void *caller_data, long num_voxels, int input_num_buffers,
             int input_vector_length, double *input_data[], int output_num_buffers,
             int output_vector_length, double *output_data[], Loop_Info * loop_info)
{

/* 
 *	Get pointer to window info 
 */
   Math_Data *md = (Math_Data *) caller_data;
   gsl_vector_view view;
   gsl_vector *conc, *aif;
   size_t   i, n_outputs;
   long     ivox, Nivox;
   double   raw_baseline, raw_avg;
   double   arrival, delay, temp;
   int      length, delayi, residue_shift;
   int      svd_s_len;

   /* shut the compiler up */
   (void)loop_info;
   (void)output_vector_length;

/* 
 *	calc_perf SVD vector and matrices 
 */
   Nivox = num_voxels * input_vector_length;
   for(ivox = 0; ivox < Nivox; ivox++){

/* 
 *	calc baseline for the raw time-series 
 */
      raw_baseline = 0.0;
      for(i = 0; i < md->Nstart; i++){
         raw_baseline += input_data[i][ivox];
         }
      raw_baseline /= md->Nstart;

/* 
 *	calc average value for the entire series 
 */
      raw_avg = 0.0;
      for(i = 0; i < input_num_buffers; i++){
         raw_avg += input_data[i][ivox];
         }
      raw_avg /= input_num_buffers;

/* 
 *	check if we have an abberant timeseries 
 */
      if((raw_baseline <= 0) || (raw_avg / md->aif->baseline < md->cutoff)){
         for(i = 0; i < (size_t) output_num_buffers; i++){
            output_data[i][ivox] = 0.0;
            }
         }

/* 
 *	else do some real work 
 */
      else {

/* 
 *	convert signal intensity to concentration 
 */
         for(i = 0; i < (size_t) input_num_buffers; i++){
            gsl_vector_set(md->signal, i, input_data[i][ivox]);
            }
         SignalToConc(md->signal, md->signal, raw_baseline, md->te);

/* 
 *	filter the time-series 
 */
         if(md->filter){
            convolve_vector_avg(md->signal, md->kernel);
            }
         for(i = md->Nstart; i < (size_t) input_num_buffers; i++){
            gsl_vector_set(md->conc, i - md->Nstart, gsl_vector_get(md->signal, i));
            }

/*
 *	fit a gamma variate function to the voxel if required
 */
         if((md->conc_fit) || (md->output_gamma)
            || (((md->output_arrival) || (md->shift_type != SHIFT_NONE))
                && (md->vox_bat_type == BAT_GAMMA))){
            TimeCurveArea(md->conc, md->tr, &(md->area), md->gparm);
            }
         if(md->conc_fit){
            eval_gamma(md->conc, md->gparm, md->tr);
            }

/*
 *	shift the concentration and aif as required  
 *    -- yes, there is some duplication here, but it's easier to read
 */
         switch (md->shift_type){
         case SHIFT_NONE:
            conc = md->conc;
            break;

         case SHIFT_AIF:
            arrival = bolus_arrival_time(md->conc, md->tr, md->gparm,
                                         md->vox_bat_type, md->vox_bat_cutoff);
            delay = arrival - md->aif->arrival_time;
            delayi = rint(delay / md->tr);
            length = md->conc->size - delayi;

            view = gsl_vector_subvector(md->conc, delayi, length);
            conc = &(view.vector);
            md->residue = gsl_vector_alloc(conc->size);

            view = gsl_vector_subvector(md->aif->conc, md->aif->arrival_time, length);
            aif = &(view.vector);
            break;

         case SHIFT_AIF_EXACT:
            arrival = bolus_arrival_time(md->conc, md->tr, md->gparm,
                                         md->vox_bat_type, md->vox_bat_cutoff);
            delay = arrival - md->aif->arrival_time;
            delayi = rint(delay / md->tr);
            length = md->conc->size - delayi;

            view = gsl_vector_subvector(md->conc, delayi, length);
            conc = &(view.vector);
            md->residue = gsl_vector_alloc(conc->size);

            aif = gsl_vector_alloc(length);
            md->aif->gparm[GAM_BAT] += delay;
            eval_gamma(aif, md->aif->gparm, md->tr);
            md->aif->gparm[GAM_BAT] -= delay;
            break;

         case SHIFT_CONC:
            arrival = bolus_arrival_time(md->conc, md->tr, md->gparm,
                                         md->vox_bat_type, md->vox_bat_cutoff);
            delay = arrival - md->aif->arrival_time;
            delayi = rint(delay / md->tr);
            length = md->conc->size - delayi;

            if((delayi > 0) && (delayi < length)){
               for(i = 0; i < length; i++){
                  gsl_vector_set(md->conc, i, gsl_vector_get(md->conc, i + delayi));
                  }
               for(i = length; i < md->conc->size; i++){
                  gsl_vector_set(md->conc, i, 0.0);
                  }
               }
            conc = md->conc;

         case SHIFT_CIRC:
            conc = md->conc;
            break;

         default:
            fprintf(stderr, "ERROR - shift method not implemented\n");
            exit(EXIT_FAILURE);
            }

/*
 *	calculate a new aif matrix if required
 */
         if((md->shift_type == SHIFT_AIF || md->shift_type == SHIFT_AIF_EXACT)){
            md->svd_u = gsl_matrix_alloc(length, length);
            md->svd_v = gsl_matrix_alloc(length, length);
            md->svd_s = gsl_vector_alloc(length);
            md->svd_work = gsl_vector_alloc(length);
            CalcAifMatrix(aif, md->svd_u, md->tr, 0, length);
            gsl_linalg_SV_decomp(md->svd_u, md->svd_v, md->svd_s, md->svd_work);
            for(i = 0; i < md->Nrest; i++){
               if(gsl_vector_get(md->svd_s, i) < temp){
                  gsl_vector_set(md->svd_s, i, 0.0);
                  }
               }
            }

/* 
 *	use SVD to evaluate the residue function from the AIF and 
 *	concentration time-curve  
 */
         switch (md->svd_type){
         case SVD_TOL:
            gsl_linalg_SV_solve(md->svd_u, md->svd_v, md->svd_s, conc, md->residue);
            break;

         case SVD_OI:
            gsl_vector_set_zero(md->svd_s_cp);
            svd_s_len = -1;
            do {
               svd_s_len++;
               gsl_vector_set(md->svd_s_cp, svd_s_len,
                              gsl_vector_get(md->svd_s, svd_s_len));
               gsl_linalg_SV_solve(md->svd_u, md->svd_v, md->svd_s_cp, conc, md->residue);
               } while((svd_s_len < md->length - 1)
                    && (svd_oscillation_index(md->residue) < md->svd_oi));
            break;

         default:
            fprintf(stderr, "ERROR - SVD type is undefined (%d)\n", md->svd_type);
            exit(EXIT_FAILURE);
            }

/* 
 *	output CBF, CBV and MTT
 */
         temp = 0.5 * (gsl_vector_get(conc, 0)
                       + gsl_vector_get(conc, md->Nrest - 1));
         temp += vector_sum(conc, 1, md->Nrest - 1);

         output_data[OUT_CBV][ivox] = temp * md->tr / md->aif->area;
         if(output_data[OUT_CBV][ivox] < 0.0){
            output_data[OUT_CBV][ivox] = 0.0;
            }
         output_data[OUT_CBF][ivox] = gsl_vector_max(md->residue);
         if((output_data[OUT_CBF][ivox] > 0.0)
            && (output_data[OUT_CBV][ivox] > 0.0)){
            output_data[OUT_MTT][ivox] = output_data[OUT_CBV][ivox]
               / (output_data[OUT_CBF][ivox]);
            }
         else {
            output_data[OUT_MTT][ivox] = 0.0;
            }
         if(output_data[OUT_MTT][ivox] > md->tr * md->Nrest){
            output_data[OUT_MTT][ivox] = md->tr * md->Nrest;
            }
         n_outputs = NUM_OUT_BASIC;

/*
 *	correct for normalization and dosage from Ostergaard
 */
         output_data[OUT_CBV][ivox] *= md->normalization / md->dosage;
         output_data[OUT_CBF][ivox] *= 60.0 * md->normalization / md->dosage;

/* 
 *	output chi-squared fit 
 */
         if(md->output_chi){
            gsl_vector_memcpy(md->svd_work, conc);
            gsl_blas_dgemv(CblasNoTrans, 1.0, md->svd_u, md->residue, -1.0, md->svd_work);
            temp = gsl_blas_dnrm2(md->svd_work);
            output_data[n_outputs + OUT_CHI][ivox] = temp * temp;
            n_outputs += NUM_OUT_CHI;
            }

/*
 *	output gamma
 */
         if(md->output_gamma){
            output_data[n_outputs + OUT_BAT][ivox] = md->gparm[GAM_BAT];
            output_data[n_outputs + OUT_FAC][ivox] = md->gparm[GAM_FAC];
            output_data[n_outputs + OUT_ALPHA][ivox] = md->gparm[GAM_ALPHA];
            output_data[n_outputs + OUT_BETA][ivox] = md->gparm[GAM_BETA];
            n_outputs += NUM_OUT_GAMMA;
            }

         if(md->shift_type == SHIFT_CIRC){
            residue_shift = gsl_vector_max_index(md->residue);
            delay = residue_shift * md->tr;
            if(delay > md->length / 2.0){
               delay -= md->length;
               }
            arrival = delay + md->aif->arrival_time;
            }
         else {
            residue_shift = 0;
            }

/*
 *	output arrival
 */
         if(md->output_arrival){
            output_data[n_outputs + OUT_ARRIVAL][ivox] = arrival;
            n_outputs += NUM_OUT_ARRIVAL;
            }

/*
 *	output delay
 */
         if(md->output_delay){
            output_data[n_outputs + OUT_DELAY][ivox] = delay;
            n_outputs += NUM_OUT_DELAY;
            }

/*
 *	output residue
 */
         if(md->output_residue){
            for(i = 0; i < md->Nrest; i++){
               output_data[n_outputs + i][ivox]
                  = gsl_vector_get(md->residue, (i + residue_shift) % md->length);
               }
            }

/*
 *	free memory
 */
         if(md->shift_type == SHIFT_AIF_EXACT){
            gsl_vector_free(aif);
            }
         if((md->shift_type == SHIFT_AIF || md->shift_type == SHIFT_AIF_EXACT)){
            gsl_matrix_free(md->svd_u);
            gsl_matrix_free(md->svd_v);
            gsl_vector_free(md->svd_s);
            gsl_vector_free(md->svd_work);
            gsl_vector_free(md->residue);
            }
         }
      }

   return;
   }
