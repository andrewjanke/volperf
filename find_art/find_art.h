/* find_art.h */

#ifndef FIND_ART_H
#define FIND_ART_H

#define VERSION "find_art 1.0"

#include <stdio.h>
#include <gtk/gtk.h>
#include <gtkgl/gdkgl.h>
#include <gtkgl/gtkglarea.h>
#include <volume_io.h>
#include "../minc_vector_io.h"

typedef struct {
   GtkWidget *main_window;
   GtkWidget *window_vbox;

   GtkWidget *menubar_hbox;
   GtkWidget *menubar;
   GtkWidget *toolbar;

   GtkWidget *statusbar;
   
   GtkWidget *progressbar;

   GtkWidget *main_hbox;

   /* controls */
   GtkWidget *controls_table;
   GtkObject *slice_adjustment;
   GtkObject *frame_adjustment;
   GtkWidget *slice_label;
   GtkWidget *frame_label;
   GtkWidget *slice_vscale;
   GtkWidget *frame_vscale;

   /* Open GL windows */
   GtkWidget *slice_gtk_glarea;

   GtkWidget *graph_table;
   GtkWidget *graph_gtk_glarea;
   GtkWidget *CBV_gtk_glarea;
   GtkWidget *CBF_gtk_glarea;
   GtkWidget *MTT_gtk_glarea;
   GtkWidget *CBF_label;
   GtkWidget *CBV_label;
   GtkWidget *MTT_label;
} Widget_Info;


typedef struct {
   Widget_Info w;

   int      verbose;
   int      clobber;
   char    *arg_string;
   int      n_infiles;
   char   **infiles;
   char    *aif_file;

   double   c_slice;
   double   c_frame;
   double   p_slice;
   double   p_frame;

   /* current ROI and point */
   double   roi[4];
   double   c_point[2];

   /* minc info */
   int      minc_id[100];
   int      img_id[100];
   int      icv[100];
   Real     steps[MAX_VAR_DIMS];
   int      sizes[MAX_VAR_DIMS];
   Real     data_min[100];
   Real     data_max[100];
   Real     calc_min[100];
   Real     calc_max[100];
   double  *data[100];
   
   /* data stores */
   MINC_Vector *t_vector;
   MINC_Vector *y_vector;
   MINC_Vector *x_vector;
   unsigned char *image_data;
   double    graph_max;
   
   double  *maps[5];
   unsigned char *CBF_image;
   unsigned char *CBV_image;
   unsigned char *MTT_image;

   /* OpenGL */
   double   begin_x;
   double   begin_y;
   double   scale_fac;
   int      window_size;

} Main_Info;

/* public functions */
void     output_vector_file(Main_Info * mi, char *filename, MINC_Vector * mv);
int      load_minc_data(int slice, Main_Info * mi);
void     clean_up(Main_Info * mi);

#endif
