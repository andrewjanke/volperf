/* callbacks.h */

#include "callbacks.h"
#include "interface.h"
#include "gtk_gl.h"

void slice_adjustment_value_changed(GtkAdjustment *adjustment, gpointer user_data){
   Main_Info *mi = (Main_Info*)user_data;
   gchar    buf[128];

   
   mi->c_slice = (double)adjustment->value;
   
   /* reset the current slice max counter */
   mi->graph_max = -DBL_MAX;
   
   if(!load_minc_data(mi->c_slice, mi)){
      g_snprintf(buf, 128, "Error loading data for slice %g[%d]", 
                           mi->c_slice, (int)mi->c_slice);
      push_statusbar(mi, buf);
      return;
      }
   
   gtk_widget_queue_draw(GTK_WIDGET(mi->w.slice_gtk_glarea));
   gtk_widget_queue_draw(GTK_WIDGET(mi->w.graph_gtk_glarea));
   }

void frame_adjustment_value_changed(GtkAdjustment *adjustment, gpointer user_data){
   Main_Info *mi = (Main_Info*)user_data;
   
   mi->c_frame = (double)adjustment->value;

   gtk_widget_queue_draw(GTK_WIDGET(mi->w.slice_gtk_glarea));
   }

/* updates the statusbar */
void push_statusbar(Main_Info *mi, char *buf)
{
   if(buf == NULL){
      g_message("push_statusbar passed a NULL pointer!");
      return;
      }
   else{
      gtk_statusbar_pop(GTK_STATUSBAR(mi->w.statusbar), 1);
      gtk_statusbar_push(GTK_STATUSBAR(mi->w.statusbar), 1, buf);
      }

   }

/* toolbar buttons */
void     on_calc_button_clicked(GtkButton * button, gpointer user_data){
   Main_Info *mi = (Main_Info*)user_data;
   
   gchar    buf[128];
   
   /* calc perfusion maps */
 
   g_snprintf(buf, 128, "Calculating maps for AIF[%g:%g]", mi->c_point[0], mi->c_point[1]);
   push_statusbar(mi, buf);
 
   if(load_maps(mi)){
      gtk_widget_queue_draw(GTK_WIDGET(mi->w.CBV_gtk_glarea));
      gtk_widget_queue_draw(GTK_WIDGET(mi->w.CBF_gtk_glarea));
      gtk_widget_queue_draw(GTK_WIDGET(mi->w.MTT_gtk_glarea));
      }
   
   }

void     on_save_button_clicked(GtkButton * button, gpointer user_data){
   Main_Info *mi = (Main_Info*)user_data;
   
   output_vector_file(mi, mi->aif_file, mi->t_vector);
   }

void     on_quit_button_clicked(GtkButton * button, gpointer user_data){
   gtk_main_quit();
   }
