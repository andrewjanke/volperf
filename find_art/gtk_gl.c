/* gtkglarea_calls.c */

#include <math.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include "gtk_gl.h"
#include "find_art.h"
#include "callbacks.h"
#include "../perf_util.h"

/* internal function prototypes */
int      load_image(Main_Info * mi);
void     load_t_vector(Main_Info * mi, double *max);

gint     glarea_button_press(GtkWidget *, GdkEventButton *, gpointer func_data);
gint     glarea_motion_notify(GtkWidget *, GdkEventMotion *, gpointer func_data);
gint     glarea_button_release(GtkWidget *, GdkEventButton *, gpointer func_data);

gint     glarea_slice_draw(GtkWidget *, GdkEventExpose *, gpointer func_data);
gint     glarea_graph_draw(GtkWidget *, GdkEventExpose *, gpointer func_data);
gint     glarea_reshape(GtkWidget *, GdkEventConfigure *, gpointer func_data);

GtkWidget *create_slice_glarea(Main_Info * mi)
{
   int      attrlist[13];
   GtkWidget *glarea;

   /* From the glXChooseVisual manpage:                        */
   /* glXChooseVisual returns a pointer to an XVisualInfo      */
   /* structure describing the visual that best meets a        */
   /* minimum specification.                                   */
   attrlist[0] = GDK_GL_RGBA;
   attrlist[1] = GDK_GL_DOUBLEBUFFER;
   attrlist[2] = GDK_GL_RED_SIZE;
   attrlist[3] = 1;
   attrlist[4] = GDK_GL_GREEN_SIZE;
   attrlist[5] = 1;
   attrlist[6] = GDK_GL_BLUE_SIZE;
   attrlist[7] = 1;
   attrlist[8] = GDK_GL_ALPHA_SIZE;
   attrlist[9] = 1;
   attrlist[10] = GDK_GL_DEPTH_SIZE;
   attrlist[11] = 1;
   attrlist[12] = GDK_GL_NONE;

   if(gdk_gl_query() == FALSE){
      g_error("OpenGL not supported!\n");
      return NULL;
      }

   if((glarea = gtk_gl_area_new(attrlist)) == NULL){
      g_error("Error creating GtkGLArea!\n");
      return NULL;
      }

   /* Check out gdk/gdktypes.h in your include directory for a */
   /* complete list of event masks that you can use.           */
   gtk_widget_set_events(GTK_WIDGET(glarea), GDK_EXPOSURE_MASK
                         | GDK_BUTTON_PRESS_MASK
                         | GDK_BUTTON_RELEASE_MASK
                         | GDK_POINTER_MOTION_MASK | GDK_POINTER_MOTION_HINT_MASK);

   gtk_signal_connect(GTK_OBJECT(glarea), "button_press_event",
                      GTK_SIGNAL_FUNC(glarea_button_press), mi);
   gtk_signal_connect(GTK_OBJECT(glarea), "motion_notify_event",
                      GTK_SIGNAL_FUNC(glarea_motion_notify), mi);
   gtk_signal_connect(GTK_OBJECT(glarea), "button_release_event",
                      GTK_SIGNAL_FUNC(glarea_button_release), mi);
   gtk_signal_connect(GTK_OBJECT(glarea), "expose_event", GTK_SIGNAL_FUNC(glarea_slice_draw), mi);
   gtk_signal_connect(GTK_OBJECT(glarea), "configure_event", GTK_SIGNAL_FUNC(glarea_reshape), mi);

   return (glarea);
   }

GtkWidget *create_graph_glarea(Main_Info * mi)
{
   int      attrlist[] = {
      GDK_GL_RGBA,
      GDK_GL_DOUBLEBUFFER,
      GDK_GL_RED_SIZE, 1,
      GDK_GL_GREEN_SIZE, 1,
      GDK_GL_BLUE_SIZE, 1,
      GDK_GL_ALPHA_SIZE, 1,
      GDK_GL_NONE
   };
   GtkWidget *glarea;

   if(gdk_gl_query() == FALSE){
      g_error("OpenGL not supported!\n");
      return NULL;
      }

   if((glarea = gtk_gl_area_new(attrlist)) == NULL){
      g_error("Error creating GtkGLArea!\n");
      return NULL;
      }

   /* Check out gdk/gdktypes.h in your include directory for a */
   /* complete list of event masks that you can use.           */
   gtk_widget_set_events(GTK_WIDGET(glarea), GDK_EXPOSURE_MASK);
   gtk_signal_connect(GTK_OBJECT(glarea), "expose_event", GTK_SIGNAL_FUNC(glarea_graph_draw), mi);
   gtk_signal_connect(GTK_OBJECT(glarea), "configure_event", GTK_SIGNAL_FUNC(glarea_reshape), mi);

   return (glarea);
   }

GtkWidget *create_map_glarea(Map_Data * md)
{
   int      attrlist[] = {
      GDK_GL_RGBA,
      GDK_GL_DOUBLEBUFFER,
      GDK_GL_RED_SIZE, 1,
      GDK_GL_GREEN_SIZE, 1,
      GDK_GL_BLUE_SIZE, 1,
      GDK_GL_ALPHA_SIZE, 1,
      GDK_GL_NONE
   };
   GtkWidget *glarea;

   if(gdk_gl_query() == FALSE){
      g_error("OpenGL not supported!\n");
      return NULL;
      }

   if((glarea = gtk_gl_area_new(attrlist)) == NULL){
      g_error("Error creating GtkGLArea!\n");
      return NULL;
      }

   /* Check out gdk/gdktypes.h in your include directory for a */
   /* complete list of event masks that you can use.           */
   gtk_widget_set_events(GTK_WIDGET(glarea), GDK_EXPOSURE_MASK);
   gtk_signal_connect(GTK_OBJECT(glarea), "configure_event", GTK_SIGNAL_FUNC(glarea_reshape), NULL);

   return (glarea);
   }

gint glarea_button_press(GtkWidget *widget, GdkEventButton * event, gpointer func_data)
{
   Main_Info *mi = (Main_Info *) func_data;
   gchar    buf[128];

   if(gdk_pointer_grab(GTK_WIDGET(widget)->window, FALSE,
                        GDK_POINTER_MOTION_HINT_MASK |
                        GDK_BUTTON1_MOTION_MASK |
                        GDK_BUTTON2_MOTION_MASK |
                        GDK_BUTTON3_MOTION_MASK |
                        GDK_BUTTON_RELEASE_MASK, NULL, NULL, event->time) == 0){
      }

   mi->begin_x = event->x / mi->scale_fac;
   mi->begin_y = (widget->allocation.height - event->y) / mi->scale_fac;

   switch (event->button){
   default:
   case 1:
      if(event->state & GDK_SHIFT_MASK){
         }
      else{
         mi->c_point[0] = mi->begin_x;
         mi->c_point[1] = mi->begin_y;

         g_snprintf(buf, 128, "Current Point: [%d:%d]", (int)mi->begin_x, (int)mi->begin_y);
         push_statusbar(mi, buf);

         gtk_widget_queue_draw(GTK_WIDGET(mi->w.graph_gtk_glarea));
         gtk_widget_queue_draw(GTK_WIDGET(mi->w.slice_gtk_glarea));
         }
      break;

   case 2:
      if(event->state & GDK_SHIFT_MASK){
         ;
         }
      else{
         }
      break;

   case 3:
      if(event->state & GDK_SHIFT_MASK){
         }
      else{
         mi->roi[0] = mi->roi[2] = mi->begin_x;
         mi->roi[1] = mi->roi[3] = mi->begin_y;

         gtk_widget_queue_draw(GTK_WIDGET(mi->w.slice_gtk_glarea));
         }
      break;
      }

   return TRUE;
   }

/* This function handles motion events for the GtkGLArea */
gint glarea_motion_notify(GtkWidget *widget, GdkEventMotion * event, gpointer func_data)
{
   Main_Info *mi = (Main_Info *) func_data;

   GdkModifierType state;
   double   delta_x, delta_y;
   int      x, y;
   double   tmp;
   gchar    buf[128];

   if(event->is_hint){
      gdk_window_get_pointer(event->window, &x, &y, &state);
      }
   else{
      x = event->x;
      y = event->y;
      state = event->state;
      }

   x = event->x / mi->scale_fac;
   y = (widget->allocation.height - event->y) / mi->scale_fac;

   delta_x = x - mi->begin_x;
   delta_y = y - mi->begin_y;

   if(state & GDK_BUTTON1_MASK){
      if(state & GDK_SHIFT_MASK){
         }
      else{
         mi->c_point[0] = x;
         mi->c_point[1] = y;

         g_snprintf(buf, 128, "Current Point: [%d:%d]", x, y);
         push_statusbar(mi, buf);

         gtk_widget_queue_draw(GTK_WIDGET(mi->w.graph_gtk_glarea));
         gtk_widget_queue_draw(GTK_WIDGET(mi->w.slice_gtk_glarea));
         }
      }
   else if(state & GDK_BUTTON2_MASK){
      if(state & GDK_SHIFT_MASK){

         tmp = mi->c_slice - delta_y / 3;
         if(tmp < 0.0){
            tmp = 0.0;
            }
         if(tmp > mi->sizes[2]){
            tmp = mi->sizes[2];
            }
         gtk_adjustment_set_value(GTK_ADJUSTMENT(mi->w.slice_adjustment), tmp);
         }
      else{

         tmp = mi->c_frame - delta_y / 3;
         if(tmp < 0.0){
            tmp = 0.0;
            }
         if(tmp > (double)(mi->n_infiles - 1)){
            tmp = (double)(mi->n_infiles - 1);
            }
         gtk_adjustment_set_value(GTK_ADJUSTMENT(mi->w.frame_adjustment), tmp);
         }
      }
   else if(state & GDK_BUTTON3_MASK){
      if(state & GDK_SHIFT_MASK){
         ;
         }
      else{
         mi->roi[2] = x;
         mi->roi[3] = y;

         gtk_widget_queue_draw(GTK_WIDGET(mi->w.slice_gtk_glarea));
         }
      }

   mi->begin_x += delta_x;
   mi->begin_y += delta_y;
   return TRUE;
   }

/* This function handles button-release events for the GtkGLArea */
gint glarea_button_release(GtkWidget *widget, GdkEventButton * event, gpointer func_data)
{
   Main_Info *mi = (Main_Info *) func_data;
   int      x, y;

   x = event->x / mi->scale_fac;
   y = (widget->allocation.height - event->y) / mi->scale_fac;

   if(event->button == 3){

      mi->roi[2] = x;
      mi->roi[3] = y;

      gtk_widget_queue_draw(GTK_WIDGET(mi->w.slice_gtk_glarea));
      gtk_widget_queue_draw(GTK_WIDGET(mi->w.graph_gtk_glarea));
      }

   gdk_pointer_ungrab(event->time);
   return TRUE;
   }

int load_image(Main_Info * mi)
{

   int      c;
   double   min, max;

   int      c_frame = (int)mi->c_frame;

   min = 1000000000;
   max = -1000000000;
   for(c = 0; c < mi->sizes[0] * mi->sizes[1]; c++){
      if(mi->data[c_frame][c] < min){
         min = mi->data[c_frame][c];
         }
      if(mi->data[c_frame][c] > max){
         max = mi->data[c_frame][c];
         }
      }

   /* rescale the data */
   for(c = 0; c < mi->sizes[0] * mi->sizes[1]; c++){
      mi->image_data[c] =
         (unsigned char)((mi->data[c_frame][c] - min) / (max - min) * 255);
      }

   return 1;
   }

/* slice window draw function */
gint glarea_slice_draw(GtkWidget *widget, GdkEventExpose * event, gpointer func_data)
{
   Main_Info *mi = (Main_Info *) func_data;

   /* Draw only on the last expose event. */
   if(event->count > 0){
      return TRUE;
      }

   if(gtk_gl_area_make_current(GTK_GL_AREA(widget))){

      glClearColor(0.0, 0.0, 0.0, 1.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      /* load image if we have to */
      if(mi->p_frame != mi->c_frame || mi->p_slice != mi->c_slice){
         if(!load_image(mi)){
            g_message("Failed loading image F:%g  S:%g\n", mi->c_frame, mi->c_slice);
            }
         }

      glRasterPos2f(0.0, 0.0);
      glPixelZoom(mi->scale_fac, mi->scale_fac);
      glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
      glDrawPixels(mi->sizes[0], mi->sizes[1], GL_LUMINANCE, GL_UNSIGNED_BYTE, mi->image_data);


      /* draw the "decorations" */
      glPushMatrix();
      glScaled(mi->scale_fac, mi->scale_fac, mi->scale_fac);

      /* ROI */
      glColor3d(1.0, 0.0, 0.0);
      glPointSize(2.0);
      glBegin(GL_LINE_LOOP);
      glVertex2d(mi->roi[0], mi->roi[1]);
      glVertex2d(mi->roi[2], mi->roi[1]);
      glVertex2d(mi->roi[2], mi->roi[3]);
      glVertex2d(mi->roi[0], mi->roi[3]);
      glEnd();

      /* cross-hair */
      glTranslatef(mi->c_point[0], mi->c_point[1], 0.0);
      glColor3f(1.0, 1.0, 0.0);
      glBegin(GL_LINES);
      glVertex2d(-2.5, 2.5);
      glVertex2d(2.5, -2.5);

      glVertex2d(-2.5, -2.5);
      glVertex2d(2.5, 2.5);
      glEnd();
      glPopMatrix();

      /* update previous counters */
      mi->p_frame = mi->c_frame;
      mi->p_slice = mi->c_slice;
      gtk_gl_area_swap_buffers(GTK_GL_AREA(widget));
      }

   return TRUE;
   }

/* load the t_vector */
void load_t_vector(Main_Info * mi, double *max)
{

   int      c;
   int      offset;

   /* calc offset */
   offset = mi->sizes[1] * mi->c_point[1] + mi->c_point[0];

   /* set up the vector */
   *max = mi->data[0][offset];
   for(c = 0; c < mi->n_infiles; c++){
      if(mi->data[c][offset] > *max){
         *max = mi->data[c][offset];
         }
      mi->t_vector->V[c] = mi->data[c][offset];
      }
   }

/* graph window draw function */
gint glarea_graph_draw(GtkWidget *widget, GdkEventExpose * event, gpointer func_data)
{
   Main_Info *mi = (Main_Info *) func_data;
   int      i;
   float    x_div, y_div;
   double   max;

   double   pad = 0.05;

   /* Draw only on the last expose event. */
   if(event->count > 0){
      return TRUE;
      }

   if(gtk_gl_area_make_current(GTK_GL_AREA(widget))){

      load_t_vector(mi, &max);
      if(max > mi->graph_max){
         mi->graph_max = max;
         }

      x_div = (float)widget->allocation.width / mi->n_infiles;
      y_div = (float)widget->allocation.height / mi->graph_max;

      glClearColor(0.0, 0.0, 0.0, 1.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      glPushMatrix();

      /* a bit of tom-foolery to center the plot */
      glTranslated((double)widget->allocation.width * (pad / 2.0),
                   (double)widget->allocation.height * (pad / 2.0), 0);
      glScaled(1.0 - pad, 1.0 - pad, 1.0);

      /* draw the axes */
      glPointSize(2.0);
      glColor4d(1.0, 1.0, 0.0, 1.0);
      glBegin(GL_LINE_STRIP);
      glVertex2i(0, widget->allocation.height);
      glVertex2i(0, 0);
      glVertex2i(widget->allocation.width, 0);
      glEnd();

      /* draw the line */
      glLineWidth(1.0);
      glColor4d(1.0, 1.0, 1.0, 1.0);
      glBegin(GL_LINE_STRIP);
      for(i = 0; i < mi->n_infiles; i++){
         glVertex2f(i * x_div, mi->t_vector->V[i] * y_div);
         }
      glEnd();

      /* draw the points */
      glPointSize(2.0);
      glColor3d(0.0, 1.0, 0.0);
      glBegin(GL_POINTS);
      for(i = 0; i <= mi->n_infiles; i++){
         glVertex2f(i * x_div, mi->t_vector->V[i] * y_div);
         }
      glEnd();

      /* draw the vertical tick marks */
      glLineWidth(1.0);
      glColor4d(1.0, 1.0, 0.0, 0.8);
      for(i = 0; i < mi->graph_max; i++){
         glBegin(GL_LINES);
         if(i % 100 == 0){
            glVertex2d(-5.0, i * y_div);
            glVertex2d(0.0, i * y_div);
            }
         else if(i % 20 == 0){
            glVertex2d(-2.5, i * y_div);
            glVertex2d(0.0, i * y_div);
            }
         glEnd();
         }

      glPopMatrix();

      gtk_gl_area_swap_buffers(GTK_GL_AREA(widget));
      }

   return (TRUE);
   }


/* This should be called whenever the size of the area changes   */
gint glarea_reshape(GtkWidget *widget, GdkEventConfigure * event, gpointer func_data)
{
   GLint    pix_w, pix_h;

   if(gtk_gl_area_make_current(GTK_GL_AREA(widget))){
      pix_w = (GLint) widget->allocation.width;
      pix_h = (GLint) widget->allocation.height;

      glLoadIdentity();
      glMatrixMode(GL_PROJECTION);
      glViewport(0, 0, pix_w, pix_h);
      gluOrtho2D(0, pix_w, 0, pix_h);
      glMatrixMode(GL_MODELVIEW);
      }

   return TRUE;
   }
