/* volperf.c                                                                 */
/*                                                                           */
/* The MINC perfusion calculator                                             */
/*                                                                           */
/* Andrew Janke - rotor@cmr.uq.edu.au                                        */
/* Mark Griffin - mark.griffin@cmr.uq.edu.au                                 */
/* Center for Magnetic Resonance                                             */
/* University of Queensland                                                  */
/*                                                                           */
/* Copyright Andrew Janke & Mark Griffin, The University of Queensland.      */
/* Permission to use, copy, modify, and distribute this software and its     */
/* documentation for any purpose and without fee is hereby granted,          */
/* provided that the above copyright notice appear in all copies.  The       */
/* author and the University of Queensland make no representations about the */
/* suitability of this software for any purpose.  It is provided "as is"     */
/* without express or implied warranty.                                      */
/*                                                                           */
/* Thu Dec  6 22:33:12 EST 2001 - initial version (that didn't work)         */
/* Wed Jul  3 15:59:46 EST 2002 - major revisions (to make it work)          */
/* Wed Aug 21 14:03:47 EST 2002 - more revisions (and now it does work)      */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <float.h>
#include <math.h>

#include <ParseArgv.h>
#include <time_stamp.h>
#include "perf_util.h"

#define DEFAULT_BOOL -1

/* Output files                      */
/* chi  - Chi Squared fit            */
/* CBF  - Cerebral Blood Flow        */
/* CBV  - Cerebral Blood Volume      */
/* MTT  - Mean Transit Time          */
/* BAT  - Bolus Arrival Time         */
char     outfile_names[][11] = {
   "chi",
   "CBF",
   "CBV",
   "MTT",
   "BAT",
   };

/* Argument variables and table */
int      verbose = FALSE;
int      quiet = FALSE;
int      clobber = FALSE;
int      max_buffer = 4 * 1024;
int      copy_header = DEFAULT_BOOL;
int      is_signed = FALSE;
nc_type  dtype = NC_UNSPECIFIED;

char    *aif_file = NULL;
int      start_tpoint = -1;
double   tr = -1;
double   te = -1;
int      filter = FALSE;
double   svd_tol = 0.2;
double   cutoff = 0.4;
int      calc_bat = FALSE;

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

   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, "\nPerfusion Options:"},
   {"-aif", ARGV_STRING, (char *)1, (char *)&aif_file,
    "<art_if.aif> File containing the Arterial input Function"},
   {"-start", ARGV_INT, (char *)1, (char *)&start_tpoint,
    "Time point from which to start fitting (indexed from 0)"},
   {"-bat", ARGV_CONSTANT, (char *)TRUE, (char *)&calc_bat,
    "Calculate and correct for the Bolus Arrival Time. (-start not required)"},
   {"-tr", ARGV_FLOAT, (char *)1, (char *)&tr,
    "TR of the perfusion experiment (units)"},
   {"-te", ARGV_FLOAT, (char *)1, (char *)&te,
    "TE of the perfusion experiment (units)"},
   {"-filter", ARGV_CONSTANT, (char *)TRUE, (char *)&filter,
    "Filter the data before fitting (low-pass gaussian)"},
   {"-tolerance", ARGV_FLOAT, (char *)1, (char *)&svd_tol,
    "[0..1] Tolerance to use for the SVD fitting"},
   {"-cutoff", ARGV_FLOAT, (char *)1, (char *)&cutoff,
    "[0..1] Minimum signal change required at a voxel for a perfusion calculation"},

   {NULL, ARGV_HELP, NULL, NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
   };

int main(int argc, char *argv[])
{
   char   **infiles;
   char   **outfiles;
   int      n_infiles;
   int      n_outfiles;
   char    *arg_string;
   size_t   i, j;

   Loop_Options *loop_opts;
   Math_Data *md;

   /* Save time stamp and args */
   arg_string = time_stamp(argc, argv);

   /* Get arguments */
   if(ParseArgv(&argc, argv, argTable, 0) || (argc < 2)){
      fprintf(stdout, "\nUsage: %s [options] <in1.mnc> <in2.mnc> ... <out_base>\n",
              argv[0]);
      fprintf(stdout, "       %s -help\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }

   /* check for an input TR and TE */
   if(tr == -1 || te == -1){
      fprintf(stdout, "%s: You need to input a tr and a te\n", argv[0]);
      fprintf(stdout, "%s -help  for more info\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }

   /* check for an input start time-point (if we are not using the BAT) */
   if(!calc_bat && start_tpoint == -1){
      fprintf(stdout, "%s: A -start argument is required if not using -bat\n", argv[0]);
      fprintf(stdout, "%s -help  for more info\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }

   /* Initialise and read in the AIF */
   if(access(aif_file, F_OK) != 0){
      fprintf(stdout, "%s: Couldn't find Arterial input file: %s\n", argv[0], aif_file);
      exit(EXIT_FAILURE);
      }

   md = (Math_Data *) malloc(sizeof(Math_Data));
   md->aif = (Art_IF *) malloc(sizeof(Art_IF));
   if(input_art_if(aif_file, md->aif) != OK){
      fprintf(stdout, "%s: Died whilst reading in AIF file: %s\n", argv[0], aif_file);
      exit(EXIT_FAILURE);
      }

   /* set up in and outfile counts */
   n_infiles = argc - 2;
   infiles = &argv[1];
   n_outfiles = sizeof(outfile_names) / sizeof(outfile_names[0]);

   /* check the infiles match and exist */
   if(n_infiles != (int)md->aif->AIF->size){
      fprintf(stdout, "%s: %d infiles != %d AIF time-points\n",
              argv[0], n_infiles, md->aif->AIF->size);
      exit(EXIT_FAILURE);
      }

   if(verbose){
      fprintf(stdout, "\n==== Infiles / AIF ====\n");
      }
   for(i = 0; i < (size_t) n_infiles; i++){
      if(verbose){
         fprintf(stdout, "[%02d]: %s  |  %g\n", i, infiles[i],
                 gsl_vector_get(md->aif->AIF, i));
         }
      if(access(infiles[i], F_OK) != 0){
         fprintf(stdout, "%s: Couldn't find %s\n", argv[0], infiles[i]);
         exit(EXIT_FAILURE);
         }
      }

   /* create and check for the outfiles */
   outfiles = (char **)malloc(sizeof(char *) * n_outfiles);
   if(verbose){
      fprintf(stdout, "\n==== Outfiles ====\n");
      }
   for(i = 0; i < (size_t) n_outfiles; i++){
      outfiles[i] = (char *)malloc(sizeof(char) * (strlen(argv[n_infiles + 1]) +
                                                   strlen(outfile_names[i]) + 5));
      sprintf(outfiles[i], "%s.%s.mnc", argv[n_infiles + 1], outfile_names[i]);
      if(verbose){
         fprintf(stdout, "[%02d]: %s\n", i, outfiles[i]);
         }
      if(access(outfiles[i], F_OK) == 0 && !clobber){
         fprintf(stdout, "%s: %s exists, use -clobber to overwrite\n", argv[0],
                 outfiles[i]);
         exit(EXIT_FAILURE);
         }
      }

   /* Set up math data structure */
   md->tr = tr;
   md->te = te;
   md->filter = filter;
   md->svd_tol = svd_tol;
   md->cutoff = cutoff;
   md->calc_bat = calc_bat;

   md->start_tpoint = (size_t) start_tpoint;
   md->rem_tpoints = md->aif->AIF->size - md->start_tpoint;

   /* do calculations required to set up math_data */
   setup_math_data(&md);

   if(verbose){

      fprintf(stdout, "====Parameters======\n");
      fprintf(stdout, "| aif->baseline:     %g\n", md->aif->baseline);
      fprintf(stdout, "| aif->area:         %g\n", md->aif->area);
      fprintf(stdout, "| aif->arrival_time: %g\n", md->aif->arrival_time);

      fprintf(stdout, "| aif->AIF           ");
      for(i = 0; i < md->aif->AIF->size; i++){
         if((i % 10) == 0 && i != 0){
            fprintf(stdout, "\n|                    ");
            }
         fprintf(stdout, "%6.3f ", gsl_vector_get(md->aif->AIF, i));
         }
      fprintf(stdout, "\n");

      fprintf(stdout, "| start_tpoint:      %d\n", md->start_tpoint);
      fprintf(stdout, "| rem_tpoints:       %d\n", md->rem_tpoints);
      fprintf(stdout, "| TR:                %g\n", md->tr);
      fprintf(stdout, "| TE:                %g\n", md->te);
      fprintf(stdout, "| filter:            %s\n", (md->filter) ? "TRUE" : "FALSE");
      fprintf(stdout, "| svd_tol:           %g\n", md->svd_tol);
      fprintf(stdout, "| cutoff:            %g\n", md->cutoff);

      fprintf(stdout, "====svd_u======\n");
      for(i = 0; i < md->rem_tpoints; i++){
         fprintf(stdout, "| %02d: ", i);
         for(j = 0; j < md->rem_tpoints; j++){
            fprintf(stdout, "%6.3f ", gsl_matrix_get(md->svd_u, i, j));
            }
         fprintf(stdout, "\n");
         }
      }

   /* Set up and run the Voxel Loop */
   loop_opts = create_loop_options();
   set_loop_copy_all_header(loop_opts, copy_header);
   set_loop_verbose(loop_opts, !quiet);
   set_loop_clobber(loop_opts, clobber);
   set_loop_datatype(loop_opts, dtype, is_signed, 0.0, 0.0);
   set_loop_buffer_size(loop_opts, (long)1024 * max_buffer);

   voxel_loop(n_infiles, infiles,
              n_outfiles, outfiles, arg_string, loop_opts, do_math, (void *)md);

   /* be tidy */
   free_loop_options(loop_opts);

   gsl_vector_free(md->aif->AIF);
   free(md->aif);

   gsl_matrix_free(md->svd_u);
   gsl_matrix_free(md->svd_v);
   gsl_vector_free(md->svd_s);
   gsl_vector_free(md->svd_work);
   free(md);

   return (EXIT_SUCCESS);
   }
