/* find_art.c                                                           */
/*                                                                      */
/* The MINC perfusion aif finder                                        */
/*                                                                      */
/* Andrew Janke - rotor@cmr.uq.edu.au                                   */
/* Center for Magnetic Resonance                                        */
/* University of Queensland                                             */
/*                                                                      */
/* Wed Jul 17 16:17:59 EST 2002 - initial version                       */


#include <stdlib.h>
#include <unistd.h>
#include <ParseArgv.h>
#include <time_stamp.h>
#include "find_art.h"
#include "interface.h"


Main_Info *mi;

int main(int argc, char *argv[])
{
   int      c, i;

   volume_input_struct input_info;
   Volume   volume;
   char    *axis_order[3] = { MIxspace, MIyspace, MIzspace };

   Real     tmp_steps[MAX_VAR_DIMS];
   int      tmp_sizes[MAX_VAR_DIMS];

   /* alloc space for main info struct */
   mi = (Main_Info *) malloc(sizeof(Main_Info));
   
   /* Argument variables and table */
   mi->verbose = FALSE;
   mi->clobber = FALSE;
   mi->aif_file = NULL;
   mi->c_slice = 0;
   mi->c_frame = 0;
   mi->p_slice = -1;
   mi->p_frame = -1;

   ArgvInfo argTable[] = {
      {NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
       "General options:"},
      {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&mi->verbose,
       "Print out extra information."},
      {"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&mi->clobber,
       "Clobber existing files."},
      {"-slice", ARGV_FLOAT, (char *)1, (char *)&mi->c_slice,
       "Slice to start at"},
      {"-frame", ARGV_FLOAT, (char *)1, (char *)&mi->c_frame,
       "Frame to start at"},
      {"-aif", ARGV_STRING, (char *)1, (char *)&mi->aif_file,
       "<art_if.aif> Output file for the Arterial Input Function"},
      {NULL, ARGV_HELP, NULL, NULL, ""},
      {NULL, ARGV_END, NULL, NULL, NULL}
   };
   
   /* Save time stamp and args */
   mi->arg_string = time_stamp(argc, argv);

   /* Get arguments */
   if(ParseArgv(&argc, argv, argTable, 0) || (argc < 2)){
      fprintf(stdout, "\nUsage: %s [options] <in1.mnc> <in2.mnc> ...\n", argv[0]);
      fprintf(stdout, "       %s -help\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }

   /* set up infiles array and count */
   mi->n_infiles = argc - 1;
   mi->infiles = &argv[1];

   /* create and check for the aif_file */
   if(mi->aif_file == NULL){
      g_warning("%s: No output aif file specified!", argv[0]);
      }
   
   if(mi->verbose){
      g_message("outfile (AIF): %s", mi->aif_file);
      }
   if(access(mi->aif_file, F_OK) == 0 && !mi->clobber){
      g_error("%s: %s exists, use -clobber to overwrite", argv[0], mi->aif_file);
      }

   /* read in and check infiles */
   if(mi->verbose){
      g_message("==== infiles ====");
      }
   for(c = 0; c < mi->n_infiles; c++){
      if(mi->verbose){
         g_message("[%02d]: %s", c, mi->infiles[c]);
         }
      if(access(mi->infiles[c], F_OK) != 0){
         g_error("%s: Couldn't find %s", argv[0], mi->infiles[c]);
         }

      /* get start and step info */
      start_volume_input(mi->infiles[c], 3, axis_order, NC_UNSPECIFIED,
                         TRUE, 0.0, 0.0, TRUE, &volume, (minc_input_options *) NULL, &input_info);
      get_volume_real_range(volume, &mi->data_min[c], &mi->data_max[c]);
      get_volume_sizes(volume, tmp_sizes);
      get_volume_separations(volume, tmp_steps);
      cancel_volume_input(volume, &input_info);

      /* check against first volume */
      for(i = 0; i < 3; i++){
         if(c == 0){
            mi->steps[i] = tmp_steps[i];
            mi->sizes[i] = tmp_sizes[i];
            }
         else{
            if(mi->steps[i] != tmp_steps[i]){
               g_error("%s: Step of %s doesn't equal that of %s  [%g != %g]", argv[0],
                       mi->infiles[c], mi->infiles[0], mi->steps[i], tmp_steps[i]);
               }
            if(mi->sizes[i] != tmp_sizes[i]){
               g_error("%s: Size of %s doesn't equal that of %s  [%d != %d]", argv[0],
                       mi->infiles[c], mi->infiles[0], mi->sizes[i], tmp_sizes[i]);
               }
            }

         }

      /* get mincid */
      mi->minc_id[c] = miopen(mi->infiles[c], NC_NOWRITE);

      /* create the icv's and cdfid's */
      mi->icv[c] = miicv_create();
      miicv_setint(mi->icv[c], MI_ICV_TYPE, NC_DOUBLE);
      miicv_setint(mi->icv[c], MI_ICV_DO_NORM, TRUE);

      /* Attach image variable */
      mi->img_id[c] = ncvarid(mi->minc_id[c], MIimage);
      miicv_attach(mi->icv[c], mi->minc_id[c], mi->img_id[c]);

      /* calloc space for data */
      mi->data[c] = calloc(mi->sizes[0] * mi->sizes[1], sizeof(double));
      }

   /* calloc space for image data */
   mi->image_data = calloc(mi->sizes[0] * mi->sizes[1], sizeof(unsigned char));

   /* check a few args */
   if(mi->c_frame > (double)mi->n_infiles - 1){
      mi->c_frame = (double)mi->n_infiles - 1;
      }
   if(mi->c_frame < 0){
      mi->c_frame = 0;
      }
   if(mi->c_slice > (double)mi->sizes[2] - 1){
      mi->c_slice = (double)mi->sizes[2] - 1;
      }
   if(mi->c_slice < 0){
      mi->c_slice = 0;
      }

   mi->c_point[0] = (double)mi->sizes[0] / 2.0;
   mi->c_point[1] = (double)mi->sizes[1] / 2.0;

   mi->roi[0] = 0.0;
   mi->roi[1] = 0.0;
   mi->roi[2] = 0.0;
   mi->roi[3] = 0.0;
   
   /* set up the various vectors */
   mi->t_vector = new_MINC_Vector(mi->n_infiles);
   mi->x_vector = new_MINC_Vector(mi->sizes[0]);
   mi->y_vector = new_MINC_Vector(mi->sizes[1]);
   
   mi->graph_max = -DBL_MAX;
   mi->window_size = 512;
   mi->scale_fac = mi->window_size / mi->sizes[0];

   /* read in data */
   if(!load_minc_data((int)mi->c_slice, mi)){
      g_error("%s: Failed loading data for slice %d", argv[0], (int)mi->c_slice);
      }

   /* begin GTK code */
   gtk_set_locale();
   gtk_init(&argc, &argv);

   /* make the main window */
   mi->w.main_window = GTK_WIDGET(create_main_window(mi));
   gtk_widget_show(mi->w.main_window);
   gtk_main();

   
   fprintf(stdout, "GOT OUT!!!\n");
   /* clean up */
   clean_up(mi);

   return 0;
   }

void clean_up(Main_Info * mi)
{
   int      c;

   if(mi->verbose){
      g_message("Cleaning up...");
      }

   /* free tmp data store and icv's */
   for(c = 0; c < mi->n_infiles; c++){

      free(mi->data[c]);
      ncclose(mi->minc_id[c]);
      miicv_free(mi->icv[c]);
      }

   /* free image data store */
   free(mi->image_data);
   
      /* free the vectors */
      free_MINC_Vector(mi->t_vector);
      free_MINC_Vector(mi->x_vector);
      free_MINC_Vector(mi->y_vector);

   free(mi);
   }

int load_minc_data(int slice, Main_Info * mi)
{
   long     coord[3] = { 0.0, 0.0, 0.0 };
   long     count[3];
   int      c;

   count[0] = 1.0;
   count[1] = mi->sizes[1];
   count[2] = mi->sizes[0];

   /* read in data */
   for(c = 0; c < mi->n_infiles; c++){

      coord[0] = slice;
      if(miicv_get(mi->icv[c], coord, count, mi->data[c]) == MI_ERROR){
         return FALSE;
         }

      if(mi->verbose){
         g_message("[%d] Got slice %g [%dx%d]", c, mi->c_slice, mi->sizes[0], mi->sizes[1]);
         }
      }

   return TRUE;
   }


Main_Info *get_main_info_ptr(void){
   return(mi);
   }
