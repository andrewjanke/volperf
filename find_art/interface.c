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
   {"/File/_Save", "<control>S", on_save_button_clicked, 1, "<Item>"},
   {"/File/Save _As", NULL, NULL, 0, "<Item>"},
   {"/File/sep1", NULL, NULL, 0, "<Separator>"},
   {"/File/Quit", "<control>Q", on_quit_button_clicked, 0, "<Item>"},

   {"/_Tools", NULL, NULL, 0, "<Branch>"},
   {"/Tools/_Calculate Maps", "<control>K", on_calc_button_clicked, 0, "<Item>"},

   {"/_Options", NULL, NULL, 0, "<Branch>"},
   {"/Options/Test", NULL, NULL, 0, "<Item>"},

   {"/_Help", NULL, NULL, 0, "<LastBranch>"},
   {"/_Help/About", NULL, NULL, 0, "<Item>"},
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
   /* toolbar */
   GtkWidget *tb_calc_button;
   GtkWidget *tb_save_button;
   GtkWidget *tb_quit_button;

   /* structures for gtk_glareas */
   Map_Data *cbf_data = (Map_Data *) malloc(sizeof(Map_Data));
   Map_Data *cbv_data = (Map_Data *) malloc(sizeof(Map_Data));
   Map_Data *mtt_data = (Map_Data *) malloc(sizeof(Map_Data));

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

   /* toolbar */
   mi->w.toolbar = gtk_toolbar_new(GTK_ORIENTATION_HORIZONTAL, GTK_TOOLBAR_ICONS);
   gtk_widget_show(mi->w.toolbar);
   gtk_box_pack_start(GTK_BOX(mi->w.menubar_hbox), mi->w.toolbar, FALSE, FALSE, 0);
   gtk_toolbar_set_space_size(GTK_TOOLBAR(mi->w.toolbar), 6);

   tb_calc_button = gtk_button_new_with_label("Calc Perf Maps");
   gtk_toolbar_append_widget(GTK_TOOLBAR(mi->w.toolbar), tb_calc_button, NULL, NULL);
   gtk_widget_show(tb_calc_button);
   gtk_signal_connect(GTK_OBJECT(tb_calc_button), "clicked",
                      GTK_SIGNAL_FUNC(on_calc_button_clicked), mi);

   tb_save_button = gtk_button_new_with_label("Save AIF");
   gtk_toolbar_append_widget(GTK_TOOLBAR(mi->w.toolbar), tb_save_button, NULL, NULL);
   gtk_widget_show(tb_save_button);
   gtk_signal_connect(GTK_OBJECT(tb_save_button), "clicked",
                      GTK_SIGNAL_FUNC(on_save_button_clicked), mi);

   tb_quit_button = gtk_button_new_with_label("Quit");
   gtk_toolbar_append_widget(GTK_TOOLBAR(mi->w.toolbar), tb_quit_button, NULL, NULL);
   gtk_widget_show(tb_quit_button);
   gtk_signal_connect(GTK_OBJECT(tb_quit_button), "clicked",
                      GTK_SIGNAL_FUNC(on_quit_button_clicked), mi);
   
   /* progressbar */
   mi->w.progressbar = gtk_progress_bar_new();
   gtk_widget_show (mi->w.progressbar);
   gtk_box_pack_start(GTK_BOX(mi->w.menubar_hbox), mi->w.progressbar, TRUE, TRUE, 0);

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

   /* create our graph table and gtkglarea windows */
   mi->w.graph_table = gtk_table_new(3, 3, FALSE);
   gtk_widget_show(mi->w.graph_table);
   gtk_box_pack_start(GTK_BOX(mi->w.main_hbox), mi->w.graph_table, TRUE, TRUE, 0);
   gtk_container_set_border_width(GTK_CONTAINER(mi->w.graph_table), 1);
   gtk_table_set_row_spacings(GTK_TABLE(mi->w.graph_table), 1);
   gtk_table_set_col_spacings(GTK_TABLE(mi->w.graph_table), 1);

   /* main graph window */
   mi->w.graph_gtk_glarea = create_graph_glarea(mi);
   gtk_widget_show(GTK_WIDGET(mi->w.graph_gtk_glarea));
   gtk_table_attach(GTK_TABLE(mi->w.graph_table), mi->w.graph_gtk_glarea, 0, 3, 0, 1,
                    (GtkAttachOptions) (GTK_EXPAND | GTK_FILL),
                    (GtkAttachOptions) (GTK_FILL), 0, 0);
   gtk_widget_set_usize(GTK_WIDGET(mi->w.graph_gtk_glarea), mi->window_size, mi->window_size * 0.6);

   /* CBF, CBV, MTT gtk_glareas */
   cbf_data->xsize = &mi->sizes[0];
   cbf_data->ysize = &mi->sizes[1];
   cbf_data->data = &mi->CBF_image;
   mi->w.CBF_gtk_glarea = create_map_glarea(cbf_data);
   gtk_widget_show(mi->w.CBF_gtk_glarea);
   gtk_table_attach(GTK_TABLE(mi->w.graph_table), mi->w.CBF_gtk_glarea, 0, 1, 1, 2,
                    (GtkAttachOptions) (GTK_EXPAND | GTK_FILL),
                    (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), 0, 0);

   cbv_data->xsize = &mi->sizes[0];
   cbv_data->ysize = &mi->sizes[1];
   cbv_data->data = &mi->CBV_image;
   mi->w.CBV_gtk_glarea = create_map_glarea(cbv_data);
   gtk_widget_show(mi->w.CBV_gtk_glarea);
   gtk_table_attach(GTK_TABLE(mi->w.graph_table), mi->w.CBV_gtk_glarea, 1, 2, 1, 2,
                    (GtkAttachOptions) (GTK_EXPAND | GTK_FILL),
                    (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), 0, 0);

   mtt_data->xsize = &mi->sizes[0];
   mtt_data->ysize = &mi->sizes[1];
   mtt_data->data = &mi->MTT_image;
   mi->w.MTT_gtk_glarea = create_map_glarea(mtt_data);
   gtk_widget_show(mi->w.MTT_gtk_glarea);
   gtk_table_attach(GTK_TABLE(mi->w.graph_table), mi->w.MTT_gtk_glarea, 2, 3, 1, 2,
                    (GtkAttachOptions) (GTK_FILL),
                    (GtkAttachOptions) (GTK_EXPAND | GTK_FILL), 0, 0);

   /* CBF, CBV, MTT labels */
   mi->w.CBF_label = gtk_label_new("CBF");
   gtk_widget_show(mi->w.CBF_label);
   gtk_table_attach(GTK_TABLE(mi->w.graph_table), mi->w.CBF_label, 0, 1, 2, 3,
                    (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

   mi->w.CBV_label = gtk_label_new("CBV");
   gtk_widget_show(mi->w.CBV_label);
   gtk_table_attach(GTK_TABLE(mi->w.graph_table), mi->w.CBV_label, 1, 2, 2, 3,
                    (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

   mi->w.MTT_label = gtk_label_new("MTT");
   gtk_widget_show(mi->w.MTT_label);
   gtk_table_attach(GTK_TABLE(mi->w.graph_table), mi->w.MTT_label, 2, 3, 2, 3,
                    (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 0, 0);

   /* statusbar */
   mi->w.statusbar = gtk_statusbar_new();
   gtk_widget_show(mi->w.statusbar);
   gtk_box_pack_start(GTK_BOX(mi->w.window_vbox), mi->w.statusbar, FALSE, FALSE, 0);

   return mi->w.main_window;
   }
