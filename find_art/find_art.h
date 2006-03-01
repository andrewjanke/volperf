/* find_art.h */

#ifndef FIND_ART_H
#define FIND_ART_H

#define VERSION "find_art 1.0"

#define MAX_NUM_FILES 1000

#include <stdio.h>
#include <gtk/gtk.h>
#include <gtkgl/gdkgl.h>
#include <gtkgl/gtkglarea.h>
#include <volume_io.h>

#undef X
#undef Y
#undef Z

#include "../minc_vector_io.h"

typedef struct {
   GtkWidget *main_window;
   GtkWidget *window_vbox;

   GtkWidget *menubar_hbox;
   GtkWidget *menubar;

   GtkWidget *statusbar;

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
   int      minc_id[MAX_NUM_FILES];
   int      img_id[MAX_NUM_FILES];
   int      icv[MAX_NUM_FILES];
   Real     steps[MAX_VAR_DIMS];
   int      sizes[MAX_VAR_DIMS];
   Real     data_min[MAX_NUM_FILES];
   Real     data_max[MAX_NUM_FILES];
   Real     calc_min[MAX_NUM_FILES];
   Real     calc_max[MAX_NUM_FILES];
   double  *data[MAX_NUM_FILES];
   
   /* data stores */
   MINC_Vector *t_vector;
   MINC_Vector *y_vector;
   MINC_Vector *x_vector;
   unsigned char *image_data;
   double    graph_max;

   /* OpenGL */
   double   begin_x;
   double   begin_y;
   double   scale_fac;
   int      window_size;

} Main_Info;

/* public functions */
int      load_minc_data(int slice, Main_Info * mi);
void     clean_up(Main_Info * mi);
Main_Info *get_main_info_ptr(void);

#endif
