/* gtk_gl.h */

#ifndef GTK_GL_H
#define GTK_GL_H

#include "find_art.h"

/* structure for map gtk_glarea widgets */
typedef struct {
   unsigned char **data;
   int     *xsize, *ysize;
   Real    *xstep, *ystep;
} Map_Data;

/* public functions */
GtkWidget *create_slice_glarea(Main_Info * mi);
GtkWidget *create_graph_glarea(Main_Info * mi);
GtkWidget *create_map_glarea(Map_Data * md);

int load_maps(Main_Info *mi);

#endif
