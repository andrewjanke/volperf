/* interface.c */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

#include <gdk/gdkkeysyms.h>
#include <gtk/gtk.h>

#include "find_art.h"
#include "callbacks.h"
#include "interface.h"
#include "gtk_gl.h"

/* function prototypes */
GtkWidget *create_menus(GtkWidget *window, Main_Info * mi);


/* This is the GtkItemFactoryEntry structure used to generate new menus.
   Item 1: The menu path. The letter after the underscore indicates an
           accelerator key once the menu is open.
   Item 2: The accelerator key for the entry
   Item 3: The callback function.
   Item 4: The callback action.  This changes the parameters with
           which the function is called.  The default is 0.
   Item 5: The item type, used to define what kind of an item it is.
           Here are the possible values:

           NULL               -> "<Item>"
           ""                 -> "<Item>"
           "<Title>"          -> create a title item
           "<Item>"           -> create a simple item
           "<CheckItem>"      -> create a check item
           "<ToggleItem>"     -> create a toggle item
           "<RadioItem>"      -> create a radio item
           <path>             -> path of a radio item to link against
           "<Separator>"      -> create a separator
           "<Branch>"         -> create an item to hold sub items (optional)
           "<LastBranch>"     -> create a right justified branch 
*/
static GtkItemFactoryEntry menu_items[] = {
   {"/_File", NULL, NULL, 0, "<Branch>"},
   {"/File/_Save AIF", "<control>S", on_save_button_clicked, 1, "<Item>"},
   {"/File/Quit", "<control>Q", on_quit_button_clicked, 0, "<Item>"},
};


GtkWidget *create_menus(GtkWidget *window, Main_Info * mi)
{
   GtkItemFactory *item_factory;
   GtkAccelGroup *accel_group;

   accel_group = gtk_accel_group_new();

   /* This function initializes the item factory.
      Param 1: The type of menu - can be GTK_TYPE_MENU_BAR, GTK_TYPE_MENU,
      or GTK_TYPE_OPTION_MENU.
      Param 2: The path of the menu.
      Param 3: A pointer to a gtk_accel_group.  The item factory sets up
      the accelerator table while generating menus.          */
   item_factory = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<main>", accel_group);

   /* This function generates the menu items. Pass the item factory,
      the number of items in the array, the array itself, and any
      callback data for the the menu items. */
   gtk_item_factory_create_items(item_factory,
                                 sizeof(menu_items) / sizeof(menu_items[0]), 
                                 menu_items, 
                                 mi);

   /* Attach the new accelerator group to the window. */
   gtk_window_add_accel_group(GTK_WINDOW(window), accel_group);

   return (gtk_item_factory_get_widget(item_factory, "<main>"));
   }

GtkWidget *create_main_window(Main_Info * mi)
{

   /* create main window & vbox */
   mi->w.main_window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
   gtk_window_set_title(GTK_WINDOW(mi->w.main_window), VERSION);

   mi->w.window_vbox = gtk_vbox_new(FALSE, 0);
   gtk_widget_show(mi->w.window_vbox);
   gtk_container_add(GTK_CONTAINER(mi->w.main_window), mi->w.window_vbox);

   mi->w.menubar_hbox = gtk_hbox_new(FALSE, 0);
   gtk_widget_show(mi->w.menubar_hbox);
   gtk_box_pack_start(GTK_BOX(mi->w.window_vbox), mi->w.menubar_hbox, TRUE, TRUE, 0);

   /* menubar */
   mi->w.menubar = create_menus(mi->w.main_window, mi);
   gtk_widget_show(mi->w.menubar);
   gtk_box_pack_start(GTK_BOX(mi->w.menubar_hbox), mi->w.menubar, TRUE, TRUE, 0);

   /* main part of window */
   mi->w.main_hbox = gtk_hbox_new(FALSE, 0);
   gtk_widget_show(mi->w.main_hbox);
   gtk_box_pack_start(GTK_BOX(mi->w.window_vbox), mi->w.main_hbox, TRUE, TRUE, 0);


   /* controls table */
   mi->w.controls_table = gtk_table_new(2, 2, FALSE);
   gtk_widget_show(mi->w.controls_table);
   gtk_box_pack_start(GTK_BOX(mi->w.main_hbox), mi->w.controls_table, FALSE, FALSE, 0);

   /* slice adjustment, label and vscale */
   mi->w.slice_label = gtk_label_new("slice");
   gtk_widget_show(mi->w.slice_label);
   gtk_table_attach(GTK_TABLE(mi->w.controls_table), mi->w.slice_label, 0, 1, 0, 1,
                    (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

   mi->w.slice_adjustment = gtk_adjustment_new((gfloat) mi->c_slice, 0.0,
                                               (gfloat) mi->sizes[2] - 1, 1.0, 10.0, 1.0);
   gtk_signal_connect(GTK_OBJECT(mi->w.slice_adjustment), "value-changed",
                      GTK_SIGNAL_FUNC(slice_adjustment_value_changed), mi);

   mi->w.slice_vscale = gtk_vscale_new(GTK_ADJUSTMENT(mi->w.slice_adjustment));
   gtk_widget_show(mi->w.slice_vscale);
   gtk_table_attach(GTK_TABLE(mi->w.controls_table), mi->w.slice_vscale, 0, 1, 1, 2,
                    (GtkAttachOptions) (GTK_FILL),
                    (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), 0, 0);
   gtk_scale_set_value_pos(GTK_SCALE(mi->w.slice_vscale), GTK_POS_BOTTOM);
   gtk_scale_set_digits(GTK_SCALE(mi->w.slice_vscale), 0);
   gtk_range_set_update_policy(GTK_RANGE(mi->w.slice_vscale), GTK_UPDATE_DELAYED);

   mi->w.frame_label = gtk_label_new("frame");
   gtk_widget_show(mi->w.frame_label);
   gtk_table_attach(GTK_TABLE(mi->w.controls_table), mi->w.frame_label, 1, 2, 0, 1,
                    (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

   mi->w.frame_adjustment = gtk_adjustment_new((gfloat) mi->c_frame, 0.0,
                                               (gfloat) mi->n_infiles, 1.0, 10.0, 1.0);
   gtk_signal_connect(GTK_OBJECT(mi->w.frame_adjustment), "value-changed",
                      GTK_SIGNAL_FUNC(frame_adjustment_value_changed), mi);

   mi->w.frame_vscale = gtk_vscale_new(GTK_ADJUSTMENT(mi->w.frame_adjustment));
   gtk_widget_show(mi->w.frame_vscale);
   gtk_table_attach(GTK_TABLE(mi->w.controls_table), mi->w.frame_vscale, 1, 2, 1, 2,
                    (GtkAttachOptions) (GTK_FILL),
                    (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), 0, 0);
   gtk_scale_set_value_pos(GTK_SCALE(mi->w.frame_vscale), GTK_POS_BOTTOM);
   gtk_scale_set_digits(GTK_SCALE(mi->w.frame_vscale), 0);

   /* create our slice gtkglarea window */
   mi->w.slice_gtk_glarea = create_slice_glarea(mi);
   gtk_widget_show(GTK_WIDGET(mi->w.slice_gtk_glarea));
   gtk_widget_set_usize(GTK_WIDGET(mi->w.slice_gtk_glarea), mi->window_size, mi->window_size);
   gtk_box_pack_start(GTK_BOX(mi->w.main_hbox), GTK_WIDGET(mi->w.slice_gtk_glarea), TRUE, TRUE, 0);

   /* main graph window */
   mi->w.graph_gtk_glarea = create_graph_glarea(mi);
   gtk_widget_show(GTK_WIDGET(mi->w.graph_gtk_glarea));
   gtk_widget_set_usize(GTK_WIDGET(mi->w.graph_gtk_glarea), mi->window_size, mi->window_size);
   gtk_box_pack_start(GTK_BOX(mi->w.main_hbox), GTK_WIDGET(mi->w.graph_gtk_glarea), TRUE, TRUE, 0);

   /* statusbar */
   mi->w.statusbar = gtk_statusbar_new();
   gtk_widget_show(mi->w.statusbar);
   gtk_box_pack_start(GTK_BOX(mi->w.window_vbox), mi->w.statusbar, FALSE, FALSE, 0);

   return mi->w.main_window;
   }
