/* callbacks.h */

#include <gtk/gtk.h>
#include "find_art.h"

/* convenience functions */
void     push_statusbar(Main_Info * mi, char *buf);

void     slice_adjustment_value_changed(GtkAdjustment * adjustment, gpointer user_data);
void     frame_adjustment_value_changed(GtkAdjustment * adjustment, gpointer user_data);

void     on_calc_button_clicked(GtkButton * button, gpointer user_data);
void     on_save_button_clicked(GtkButton * button, gpointer user_data);
void     on_quit_button_clicked(GtkButton * button, gpointer user_data);
