/*
 #
 #  File        : gmic4gimp.cpp
 #                ( C++ source file )
 #
 #  Description : G'MIC for GIMP - A plug-in to allow the use
 #                of G'MIC commands in GIMP.
 #                This file is a part of the CImg Library project.
 #                ( http://cimg.sourceforge.net )
 #
 #  Copyright   : David Tschumperle (GREYCstoration API)
 #
 #  License     : CeCILL v2.0
 #                ( http://www.cecill.info/licences/Licence_CeCILL_V2-en.html )
 #
 #  This software is governed by the CeCILL  license under French law and
 #  abiding by the rules of distribution of free software.  You can  use,
 #  modify and/ or redistribute the software under the terms of the CeCILL
 #  license as circulated by CEA, CNRS and INRIA at the following URL
 #  "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and  rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and, more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL license and that you accept its terms.
 #
*/

// Include necessary header files.
//--------------------------------
#define cimg_display_type 0
#include "gmic.h"
#include "gmic4gimp_def.h"
#include <pthread.h>
#include <locale>
#include <gtk/gtk.h>
#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
using namespace cimg_library;

// Define plug-in global variables.
//---------------------------------
CImgList<char> gmic_entries;           // The list of recognized G'MIC menu entries (stored as 'char*' strings).
CImgList<char> gmic_commands;          // The list of corresponding G'MIC commands to process the image.
CImgList<char> gmic_preview_commands;  // The list of corresponding G'MIC commands to preview the image.
CImgList<char> gmic_arguments;         // The list of corresponding needed filter arguments.
char *gmic_macros;                     // The array of customized G'MIC macros.
GtkTreeStore *filter_store;            // A list of available filter entries (used by the GtkTreeView).
bool return_create_dialog;             // Return value of the 'create_gui_dialog()' function (set by events handlers).
void **event_infos;                    // Infos that are passed to the GUI callback functions.
char *path_home;                       // The path where configuration files are looked for.

// Replace '[]' by '()' in a C-string.
//------------------------------------
void strparenthesis(char *const s) {
  for (char *ns = s; *ns; ++ns) if (*ns=='[') *ns = '('; else if (*ns==']') *ns = ')';
}

// Set/get plug-in global variables in GIMP.
//------------------------------------------
void set_current_filter(const unsigned int current_filter) {
  const unsigned int ncurrent_filter = current_filter>gmic_entries.size?0:current_filter;
  gimp_set_data("gmic_current_filter",&ncurrent_filter,sizeof(unsigned int));
}

unsigned int get_current_filter() {
  unsigned int current_filter = 0;
  gimp_get_data("gmic_current_filter",&current_filter);
  if (current_filter>gmic_entries.size) current_filter = 0;
  return current_filter;
}

void set_filter_nbparams(const unsigned int filter, const unsigned int nbparams) {
  char s_tmp[256] = { 0 };
  std::sprintf(s_tmp,"gmic_filter%u_nbparams",filter);
  gimp_set_data(s_tmp,&nbparams,sizeof(unsigned int));
}

unsigned int get_filter_nbparams(const unsigned int filter) {
  char s_tmp[256] = { 0 };
  std::sprintf(s_tmp,"gmic_filter%u_nbparams",filter);
  unsigned int nbparams = 0;
  gimp_get_data(s_tmp,&nbparams);
  return nbparams;
}

void set_filter_parameter(const unsigned int filter, const unsigned int n, const char *const param) {
  char s_tmp[256] = { 0 };
  std::sprintf(s_tmp,"gmic_filter%u_parameter%u",filter,n);
  gimp_set_data(s_tmp,param,cimg::strlen(param)+1);
}

const char *get_filter_parameter(const unsigned int filter, const unsigned int n) {
  char s_tmp[256] = { 0 };
  std::sprintf(s_tmp,"gmic_filter%u_parameter%u",filter,n);
  static char res[4096] = { 0 };
  res[0] = 0;
  gimp_get_data(s_tmp,res);
  return res;
}

unsigned int get_verbosity_level() {
  unsigned int verbosity = 0;
  gimp_get_data("gmic_verbosity",&verbosity);
  return verbosity;
}

void set_verbosity_level(const unsigned int verbosity) {
  gimp_set_data("gmic_verbosity",&verbosity,sizeof(unsigned int));
}

// Return G'MIC command line needed to run the selected filter.
//--------------------------------------------------------------
const char* get_commandline(const bool preview) {
  const unsigned int
    verbosity_level = get_verbosity_level(),
    filter = get_current_filter(),
    nbparams = get_filter_nbparams(filter);
  if (!filter) return 0;

  static CImg<char> res;

  CImgList<char> lres;
  switch (verbosity_level) {
  case 0: lres.insert(CImg<char>("-v- -",5)); break;
  case 1: lres.insert(CImg<char>("-",1)); break;
  default: lres.insert(CImg<char>("-v+ -debug -",12));
  }

  const unsigned int N = filter - 1;
  const CImg<char> &command_item = (preview?gmic_preview_commands[N]:gmic_commands[N]);
  if (command_item) {
    lres.insert(command_item);
    if (nbparams) {
      lres[1].last() = ' ';
      for (unsigned int p = 0; p<nbparams; ++p) {
        const char *const param = get_filter_parameter(filter,p);
        if (param) lres.insert(CImg<char>(param,cimg::strlen(param)+1)).last().last() = ',';
      }
    }
    (res = lres.get_append('x')).last() = 0;
  }
  return res.ptr();
}

// Process image region with G'MIC.
//---------------------------------

// Define structure to store the arguments needed by the processing thread.
struct st_process_thread {
  pthread_t thread;
  CImgList<float> images;
  const char *commandline;
  unsigned int verbosity_level;
  pthread_mutex_t is_running;
};

// Thread that does the image processing part (call the G'MIC library).
void *process_thread(void *arg) {
  st_process_thread &spt = *(st_process_thread*)arg;
  try {
    if (spt.verbosity_level>0)
      std::fprintf(stderr,"\n*** Plug-in 'gmic4gimp' : Running G'MIC to process the image, with command : %s\n",spt.commandline);
    std::setlocale(LC_NUMERIC,"C");
    gmic(spt.commandline,spt.images,gmic_macros,false);
    if (spt.verbosity_level>0)
      std::fprintf(stderr,"\n*** Plug-in 'gmic4gimp' : G'MIC successfully returned !\n");
  } catch (gmic_exception &e) {
    if (spt.verbosity_level>0)
      std::fprintf(stderr,"\n*** Plug-in 'gmic4gimp' : Error encountered when running G'MIC :\n*** %s\n",e.message);
    spt.images.assign();
  }
  pthread_mutex_unlock(&spt.is_running);
  pthread_exit(0);
  return 0;
}

// Routine called to process the current GIMP image.
void process_image(GimpDrawable *drawable, const char *last_commandline) {
  const unsigned int filter = get_current_filter();
  if (!last_commandline && !filter) return;
  const char *commandline = last_commandline?last_commandline:get_commandline(false);
  if (!commandline || !cimg::strcmp(commandline,"-v- -nop")) return;
  gimp_progress_init_printf(" G'MIC Toolbox : %s...",gmic_entries[filter-1].ptr());

  // Read GIMP image region data and make a CImg<float> instance from it.
  GimpPixelRgn src_region;
  gint x1, y1, x2, y2;
  gimp_drawable_mask_bounds(drawable->drawable_id,&x1,&y1,&x2,&y2);  // Get coordinates of the current layer selection.
  const gint width  = x2 - x1, height = y2 - y1, channels = drawable->bpp;
  gimp_pixel_rgn_init(&src_region,drawable,x1,y1,width,height,false,false);
  guchar *const src_row = g_new(guchar,width*channels);
  CImg<float> img(width,height,1,channels);
  cimg_forY(img,y) {
    gimp_pixel_rgn_get_row(&src_region,src_row,x1,y1+y,width);
    const guchar *ptrs = src_row;
    cimg_forX(img,x) cimg_forV(img,k) img(x,y,k) = (float)*(ptrs++);
  }
  g_free(src_row);

  // Call G'MIC interpreter on the CImg<float> image in a new thread.
  st_process_thread spt;
  spt.images.assign(1);
  img.transfer_to(spt.images[0]);
  spt.commandline = commandline;
  spt.verbosity_level = get_verbosity_level();
  pthread_mutex_init(&spt.is_running,0);
  pthread_mutex_lock(&spt.is_running);
  pthread_create(&(spt.thread),0,process_thread,(void*)&spt);

  // Do a small animation with the progress bar, while waiting for
  // the termination of the processing thread.
  while (pthread_mutex_trylock(&spt.is_running)) { gimp_progress_pulse(); cimg::wait(500); }
  gimp_progress_update(1.0);
  pthread_join(spt.thread,0);
  pthread_mutex_unlock(&spt.is_running);
  pthread_mutex_destroy(&spt.is_running);

  // Force the resulting images to have all the same 2D GRAY, GRAYA, RGB or RGBA format.
  if (!spt.images) { gimp_progress_end(); return; }
  unsigned int max_width = 0, max_height = 0, max_channels = 0;
  cimglist_for(spt.images,p) {
    const CImg<float>& img = spt.images[p];
    if (img.width>max_width) max_width = img.width;
    if (img.height>max_height) max_height = img.height;
    if (img.dim>max_channels) max_channels = img.dim;
  }
  if (max_channels>4) max_channels = 4;
  cimglist_apply(spt.images,resize)(-100,-100,1,max_channels);

  // Transfer the result image back into GIMP.
  if (spt.images.size==1 && (int)max_width==width && (int)max_height==height && (int)max_channels==channels) {

    // When the result image has same dimensions than the source :
    // Replace the selected region of the original GIMP image.
    CImg<float> &res = spt.images[0];
    GimpPixelRgn dest_region;
    guchar *const dest_row = g_new(guchar,res.dimx()*res.dimv());
    gimp_pixel_rgn_init(&dest_region,drawable,0,0,drawable->width,drawable->height,true,true);
    cimg_forY(res,y) {
      guchar *ptrd = dest_row;
      cimg_forX(res,x) cimg_forV(res,k) *(ptrd++) = (guchar)res(x,y,k);
      gimp_pixel_rgn_set_row(&dest_region,dest_row,x1,y1+y,width);
    }
    g_free(dest_row);
    spt.images.assign();
    gimp_drawable_flush(drawable);
    gimp_drawable_merge_shadow(drawable->drawable_id,true);
    gimp_drawable_update(drawable->drawable_id,x1,y1,x2-x1,y2-y1);
    gimp_displays_flush();
  } else {

    // When the result image has different dimensions than the source :
    // Returns a new GIMP image.
    gint id_img = gimp_image_new(max_width,max_height,max_channels<=2?GIMP_GRAY:GIMP_RGB);
    gimp_image_undo_group_start(id_img);

    cimglist_for(spt.images,p) {
      CImg<float> &res = spt.images[p];
      gint id_layer = gimp_layer_new(id_img,"image",res.dimx(),res.dimy(),
                                     res.dimv()==1?GIMP_GRAY_IMAGE:
                                     res.dimv()==2?GIMP_GRAYA_IMAGE:
                                     res.dimv()==3?GIMP_RGB_IMAGE:
                                     GIMP_RGBA_IMAGE,
                                     100.0,GIMP_NORMAL_MODE);
      gimp_image_add_layer(id_img,id_layer,0);
      GimpDrawable *ndrawable = gimp_drawable_get(id_layer);

      GimpPixelRgn dest_region;
      guchar *const dest_row = g_new(guchar,res.dimx()*res.dimv());
      gimp_pixel_rgn_init(&dest_region,ndrawable,0,0,ndrawable->width,ndrawable->height,true,true);
      cimg_forY(res,y) {
        guchar *ptrd = dest_row;
        cimg_forX(res,x) cimg_forV(res,k) *(ptrd++) = (guchar)res(x,y,k);
        gimp_pixel_rgn_set_row(&dest_region,dest_row,0,y,res.dimx());
      }
      g_free(dest_row);
      res.assign();
      gimp_drawable_flush(ndrawable);
      gimp_drawable_merge_shadow(ndrawable->drawable_id,true);
      gimp_drawable_update(ndrawable->drawable_id,0,0,ndrawable->width,ndrawable->height);
      gimp_drawable_detach(ndrawable);
    }
    gimp_display_new(id_img);
    gimp_image_undo_group_end(id_img);
    gimp_displays_flush();
  }
  gimp_progress_end();
}

// Process preview with G'MIC.
//-----------------------------
void process_preview(GimpPreview *preview) {
  const unsigned int filter = get_current_filter();
  if (!filter) return;
  const char *const commandline = get_commandline(true);
  if (!commandline || !cimg::strcmp(commandline,"-v- -nop")) return;

  // Read GIMP image preview and make a CImg<float> instance from it.
  gint width, height, channels;
  guchar *const ptr0 = gimp_zoom_preview_get_source(GIMP_ZOOM_PREVIEW(preview),&width,&height,&channels), *ptrs = ptr0;
  CImg<float> img(width,height,1,channels);
  cimg_forXY(img,x,y) cimg_forV(img,k) img(x,y,k) = (float)*(ptrs++);

  // Call G'MIC interpreter on the preview image.
  CImgList<float> gmic_images(1);
  img.transfer_to(gmic_images[0]);
  try {
    if (get_verbosity_level()>0)
      std::fprintf(stderr,"\n*** Plug-in 'gmic4gimp' : Running G'MIC to process the preview, with command : %s\n",commandline);
    std::setlocale(LC_NUMERIC,"C");
    gmic(commandline,gmic_images,gmic_macros,false);
    if (get_verbosity_level()>0)
      std::fprintf(stderr,"\n*** Plug-in 'gmic4gimp' : G'MIC successfully returned !\n");
  } catch (gmic_exception &e) {
    if (get_verbosity_level()>0)
      std::fprintf(stderr,"\n*** Plug-in 'gmic4gimp' : Error encountered when running G'MIC :\n*** %s\n",e.message);
    gmic_images.assign();
  }

  // Get current image preview from the processed data.
  if (gmic_images.size && gmic_images[0]) {
    CImg<float>& res = gmic_images[0];
    if (res.width>res.height) {
      const unsigned int _nheight = res.height*width/res.width, nheight = _nheight?_nheight:1;
      res.resize(width,nheight,1,-100,2);
    } else {
      const unsigned int _nwidth = res.width*height/res.height, nwidth = _nwidth?_nwidth:1;
      res.resize(nwidth,height,1,-100,2);
    }
    if (res.dimx()!=width || res.dimy()!=height) res.resize(width,height,1,-100,0,0,1);
    switch (channels) {
    case 1:
      switch (res.dim) {
      case 1: break;
      case 2: res.channel(0); break;
      case 3: res.channel(0); break;
      case 4: res.channel(0); break;
      default: res.channel(0);
      } break;
    case 2:
      switch (res.dim) {
      case 1: res.resize(-100,-100,1,2,0).get_shared_channel(1).fill(255); break;
      case 2: break;
      case 3: res.channels(0,1).get_shared_channel(1).fill(255); break;
      case 4: res.get_shared_channel(1) = res.get_shared_channel(3); res.channels(0,1); break;
      default: res.channels(0,1).get_shared_channel(1).fill(255);
      } break;
    case 3:
      switch (res.dim) {
      case 1: res.resize(-100,-100,1,3); break;
      case 2: res.channel(0).resize(-100,-100,1,3); break;
      case 3: break;
      case 4: res.channels(0,2); break;
      default: res.channels(0,2);
      } break;
    case 4:
      switch (res.dim) {
      case 1: res.resize(-100,-100,1,4).get_shared_channel(3).fill(255); break;
      case 2:
        res.resize(-100,-100,1,4,0);
        res.get_shared_channel(3) = res.get_shared_channel(1);
        res.get_shared_channel(1) = res.get_shared_channel(0);
        res.get_shared_channel(2) = res.get_shared_channel(0);
        break;
      case 3: res.resize(-100,-100,1,4,0).get_shared_channel(3).fill(255); break;
      case 4: break;
      default: res.resize(-100,-100,1,4,0);
      } break;
    }
    guchar *ptrd = ptr0;
    cimg_forXY(res,x,y) cimg_forV(res,k) *(ptrd++) = (guchar)res(x,y,k);
    gimp_preview_draw_buffer(preview,ptr0,width*channels);
    g_free(ptr0);
  }
}

// Define event functions for GUI.
//--------------------------------

// Handle responses to the parameter widgets.
void on_float_parameter_changed(GtkAdjustment *scale, gpointer user_data) {
  const unsigned int arg = *(unsigned int*)user_data;
  double value = 0;
  gimp_double_adjustment_update(scale,&value);
  char s_value[1024] = { 0 };
  std::sprintf(s_value,"%g",value);
  set_filter_parameter(get_current_filter(),arg,s_value);
  return_create_dialog = true;
}

void on_int_parameter_changed(GtkAdjustment *scale, gpointer user_data) {
  const unsigned int arg = *(unsigned int*)user_data;
  int value = 0;
  gimp_int_adjustment_update(scale,&value);
  char s_value[1024] = { 0 };
  std::sprintf(s_value,"%d",value);
  set_filter_parameter(get_current_filter(),arg,s_value);
  return_create_dialog = true;
}

void on_bool_parameter_changed(GtkCheckButton *checkbutton, gpointer user_data) {
  const unsigned int arg = *(unsigned int*)user_data;
  int value = 0;
  g_object_get(checkbutton,"active",&value,NULL);
  char s_value[1024] = { 0 };
  std::sprintf(s_value,"%d",value?1:0);
  set_filter_parameter(get_current_filter(),arg,s_value);
  return_create_dialog = true;
}

void on_list_parameter_changed(GtkComboBox *combobox, gpointer user_data) {
  const unsigned int arg = *(unsigned int*)user_data;
  int value = 0;
  g_object_get(combobox,"active",&value,NULL);
  char s_value[1024] = { 0 };
  std::sprintf(s_value,"%d",value);
  set_filter_parameter(get_current_filter(),arg,s_value);
  return_create_dialog = true;
}

void on_text_parameter_changed(GtkButton *button, gpointer user_data) {
  button = 0;
  const unsigned int arg = *(unsigned int*)user_data;
  GtkWidget *entry = *((GtkWidget**)user_data+1);
  const char *s_value = gtk_entry_get_text(GTK_ENTRY(entry));
  set_filter_parameter(get_current_filter(),arg,s_value);
  return_create_dialog = true;
}

void on_file_parameter_changed(GtkFileChooserButton *widget, gpointer user_data){
  const unsigned int arg = *(unsigned int*)user_data;
  const char
    *const filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(widget)),
    *s_value = filename?filename:"";
  set_filter_parameter(get_current_filter(),arg,s_value);
  return_create_dialog = true;
}

void on_color_parameter_changed(GtkColorButton *widget, gpointer user_data){
  const unsigned int arg = *(unsigned int*)user_data;
  GdkColor color;
  gtk_color_button_get_color(GTK_COLOR_BUTTON(widget),&color);
  char s_value[1024] = { 0 };
  if (gtk_color_button_get_use_alpha(GTK_COLOR_BUTTON(widget)))
    std::sprintf(s_value,"%d,%d,%d,%d",
                 color.red>>8,color.green>>8,color.blue>>8,gtk_color_button_get_alpha(GTK_COLOR_BUTTON(widget))>>8);
  else std::sprintf(s_value,"%d,%d,%d",
                    color.red>>8,color.green>>8,color.blue>>8);
  set_filter_parameter(get_current_filter(),arg,s_value);
  return_create_dialog = true;
}

// Create parameter GUI for specific chosen filter.
//--------------------------------------------------
void create_parameters_gui(const bool reset) {
  const unsigned int filter = get_current_filter();

  // Remove widget in the current frame if necessary.
  GtkWidget *frame = 0;
  gimp_get_data("gmic_gui_frame",&frame);
  if (frame) {
    GtkWidget *child = GTK_WIDGET(gtk_bin_get_child(GTK_BIN(frame)));
    if (child) gtk_container_remove(GTK_CONTAINER(frame),child);
  }

  GtkWidget *table = 0;
  if (!filter) {  // No filter selected -> Default message.
    table = gtk_table_new(1,1,false);
    gtk_widget_show(table);
    GtkWidget *label = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(label),"<i>Select a filter...</i>");
    gtk_widget_show(label);
    gtk_table_attach(GTK_TABLE(table),label,0,1,0,1,
                     (GtkAttachOptions)(GTK_EXPAND),(GtkAttachOptions)(GTK_EXPAND),0,0);
    gtk_misc_set_alignment (GTK_MISC(label),0,0.5);
    gtk_frame_set_label(GTK_FRAME(frame),NULL);
  } else { // Filter selected -> Build parameter table.
    GtkWidget *preview = 0;
    gimp_get_data("gmic_gui_preview",&preview);
    const unsigned int N = filter - 1;
    char nlabel[4096] = { 0 };
    std::sprintf(nlabel,"<b>  %s : </b>",gmic_entries[N].ptr());
    GtkWidget *frame_title = gtk_label_new(NULL);
    gtk_widget_show(frame_title);
    gtk_label_set_markup(GTK_LABEL(frame_title),nlabel);
    gtk_frame_set_label_widget(GTK_FRAME(frame),frame_title);

    char argname[4096] = { 0 }, argtype[4096] = { 0 }, argarg[4096] = { 0 };
    unsigned int nb_arguments = 0;
    for (const char *argument = gmic_arguments[N].ptr(); *argument; ) {
      if (std::sscanf(argument,"%4095[^=]=%4095[^(](%4095[^)]",argname,argtype,&(argarg[0]=0))>=2) {
        argument += cimg::strlen(argname) + cimg::strlen(argtype) + cimg::strlen(argarg) + 3;
        if (*argument) ++argument;
        ++nb_arguments;
      } else break;
    }

    if (!nb_arguments) { // Selected filter has no parameters -> Default message.
      table = gtk_table_new(1,1,false);
      gtk_widget_show(table);
      GtkWidget *label = gtk_label_new(NULL);
      gtk_label_set_markup(GTK_LABEL(label),"<i>No parameters to set...</i>");
      gtk_widget_show(label);
      gtk_table_attach(GTK_TABLE(table),label,0,1,0,1,
                       (GtkAttachOptions)(GTK_EXPAND),(GtkAttachOptions)(GTK_EXPAND),0,0);
      gtk_misc_set_alignment (GTK_MISC(label),0,0.5);
    } else { // Selected filter has parameters -> Create parameter table.

      // Create new table for putting parameters inside.
      table = gtk_table_new(3,nb_arguments,false);
      gtk_widget_show(table);
      gtk_table_set_row_spacings(GTK_TABLE(table),6);
      gtk_table_set_col_spacings(GTK_TABLE(table),6);
      gtk_container_set_border_width(GTK_CONTAINER(table),8);

      // Parse arguments list and add recognized one to the table.
      event_infos = new void*[2*nb_arguments];
      int current_parameter = 0, current_line = 0;
      for (const char *argument = gmic_arguments[N].ptr(); *argument; ) {
        if (std::sscanf(argument,"%4095[^=]=%4095[^(](%4095[^)]",argname,argtype,&(argarg[0]=0))>=2) {
          argument += cimg::strlen(argname) + cimg::strlen(argtype) + cimg::strlen(argarg) + 3;
          if (*argument) ++argument;
          cimg::strclean(argname);
          cimg::strclean(argtype);
          const char *const s_value = get_filter_parameter(filter,current_parameter);

          // Check for a float-valued parameter -> Create GtkAdjustment.
          bool found_valid_item = false;
          if (!found_valid_item && !cimg::strcasecmp(argtype,"float")) {
            float initial_value = 0, min_value = 0, max_value = 100;
            std::setlocale(LC_NUMERIC,"C");
            std::sscanf(argarg,"%f%*c%f%*c%f",&initial_value,&min_value,&max_value);
            if (!reset && std::sscanf(s_value,"%f",&initial_value)) {}
            GtkObject *scale = gimp_scale_entry_new(GTK_TABLE(table),0,current_line,argname,100,6,
                                                    (gdouble)initial_value,(gdouble)min_value,(gdouble)max_value,
                                                    0.1,0.1,2,true,0,0,0,0);
            event_infos[2*current_parameter] = (void*)current_parameter;
            event_infos[2*current_parameter+1] = (void*)0;
            on_float_parameter_changed(GTK_ADJUSTMENT(scale),(void*)(event_infos+2*current_parameter));
            g_signal_connect(scale,"value_changed",G_CALLBACK(on_float_parameter_changed),
                             (void*)(event_infos+2*current_parameter));
            g_signal_connect_swapped(scale,"value_changed",G_CALLBACK(gimp_preview_invalidate),preview);
            found_valid_item = true;
            ++current_parameter;
          }

          // Check for an int-valued parameter -> Create GtkAdjustment.
          if (!found_valid_item && !cimg::strcasecmp(argtype,"int")) {
            float initial_value = 0, min_value = 0, max_value = 100;
            std::setlocale(LC_NUMERIC,"C");
            std::sscanf(argarg,"%f%*c%f%*c%f",&initial_value,&min_value,&max_value);
            if (!reset && std::sscanf(s_value,"%f",&initial_value)) {}
            GtkObject *scale = gimp_scale_entry_new(GTK_TABLE(table),0,current_line,argname,100,6,
                                                    (gdouble)(int)initial_value,(gdouble)(int)min_value,
                                                    (gdouble)(int)max_value,
                                                    1,1,0,true,0,0,0,0);
            event_infos[2*current_parameter] = (void*)current_parameter;
            event_infos[2*current_parameter+1] = (void*)0;
            on_int_parameter_changed(GTK_ADJUSTMENT(scale),(void*)(event_infos+2*current_parameter));
            g_signal_connect(scale,"value_changed",G_CALLBACK(on_int_parameter_changed),
                             (void*)(event_infos+2*current_parameter));
            g_signal_connect_swapped(scale,"value_changed",G_CALLBACK(gimp_preview_invalidate),preview);
            found_valid_item = true;
            ++current_parameter;
          }

          // Check for a bool-valued parameter -> Create GtkCheckButton.
          if (!found_valid_item && !cimg::strcasecmp(argtype,"bool")) {
            unsigned int initial_value = 0;
            std::sscanf(argarg,"%u",&initial_value);
            if (!reset && std::sscanf(s_value,"%u",&initial_value)) {}
            GtkWidget *checkbutton = gtk_check_button_new_with_label(argname);
            gtk_widget_show(checkbutton);
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton),initial_value?true:false);
            gtk_table_attach(GTK_TABLE(table),checkbutton,0,2,current_line,current_line+1,GTK_FILL,GTK_SHRINK,0,0);
            event_infos[2*current_parameter] = (void*)current_parameter;
            event_infos[2*current_parameter+1] = (void*)0;
            on_bool_parameter_changed(GTK_CHECK_BUTTON(checkbutton),(void*)(event_infos+2*current_parameter));
            g_signal_connect(checkbutton,"toggled",G_CALLBACK(on_bool_parameter_changed),
                             (void*)(event_infos+2*current_parameter));
            g_signal_connect_swapped(checkbutton,"toggled",G_CALLBACK(gimp_preview_invalidate),preview);
            found_valid_item = true;
            ++current_parameter;
          }

          // Check for a list-valued parameter -> Create GtkComboBox.
          if (!found_valid_item && !cimg::strcasecmp(argtype,"choice")) {
            GtkWidget *label = gtk_label_new(argname);
            gtk_widget_show(label);
            gtk_table_attach(GTK_TABLE(table),label,0,1,current_line,current_line+1,GTK_FILL,GTK_SHRINK,0,0);
            gtk_misc_set_alignment(GTK_MISC(label),0,0.5);
            GtkWidget *combobox = gtk_combo_box_new_text();
            gtk_widget_show(combobox);
            char s_entry[4096] = { 0 }, end = 0; int err2 = 0;
            unsigned int initial_value = 0;
            const char *entries = argarg;
            if (std::sscanf(entries,"%u",&initial_value)==1) entries+=std::sprintf(s_entry,"%u",initial_value) + 1;
            while (*entries) {
              if ((err2 = std::sscanf(entries,"%4095[^,]%c",s_entry,&end))>0) {
                entries += cimg::strlen(s_entry) + (err2==2?1:0);
                cimg::strclean(s_entry);
                strparenthesis(s_entry);
                gtk_combo_box_append_text(GTK_COMBO_BOX(combobox),s_entry);
              } else break;
            }
            if (!reset && std::sscanf(s_value,"%u",&initial_value)) {}
            gtk_combo_box_set_active(GTK_COMBO_BOX(combobox),initial_value);
            gtk_table_attach(GTK_TABLE(table),combobox,1,3,current_line,current_line+1,
                             (GtkAttachOptions)(GTK_EXPAND | GTK_FILL),(GtkAttachOptions)(GTK_FILL),0,0);
            event_infos[2*current_parameter] = (void*)current_parameter;
            event_infos[2*current_parameter+1] = (void*)0;
            on_list_parameter_changed(GTK_COMBO_BOX(combobox),(void*)(event_infos+2*current_parameter));
            g_signal_connect(combobox,"changed",G_CALLBACK(on_list_parameter_changed),
                             (void*)(event_infos+2*current_parameter));
            g_signal_connect_swapped(combobox,"changed",G_CALLBACK(gimp_preview_invalidate),preview);
            found_valid_item = true;
            ++current_parameter;
          }

          // Check for a text-valued parameter -> Create GtkEntry.
          if (!found_valid_item && !cimg::strcasecmp(argtype,"text")) {
            GtkWidget *label = gtk_label_new(argname);
            gtk_widget_show(label);
            gtk_table_attach(GTK_TABLE(table),label,0,1,current_line,current_line+1,GTK_FILL,GTK_SHRINK,0,0);
            gtk_misc_set_alignment(GTK_MISC(label),0,0.5);
            GtkWidget *entry = gtk_entry_new_with_max_length(4095);
            gtk_widget_show(entry);
            cimg::strclean(argarg);
            if (!reset && *s_value) gtk_entry_set_text(GTK_ENTRY(entry),s_value);
            else gtk_entry_set_text(GTK_ENTRY(entry),argarg);
            gtk_table_attach(GTK_TABLE(table),entry,1,2,current_line,current_line+1,
                             (GtkAttachOptions)(GTK_EXPAND | GTK_FILL),(GtkAttachOptions)0,0,0);
            GtkWidget *button = gtk_button_new_with_label("Update");
            gtk_widget_show(button);
            gtk_table_attach(GTK_TABLE(table),button,2,3,current_line,current_line+1,GTK_FILL,GTK_SHRINK,0,0);
            event_infos[2*current_parameter] = (void*)current_parameter;
            event_infos[2*current_parameter+1] = (void*)entry;
            on_text_parameter_changed(GTK_BUTTON(button),(void*)(event_infos+2*current_parameter));
            g_signal_connect(button,"clicked",G_CALLBACK(on_text_parameter_changed),
                             (void*)(event_infos+2*current_parameter));
            g_signal_connect_swapped(button,"clicked",G_CALLBACK(gimp_preview_invalidate),preview);
            found_valid_item = true;
            ++current_parameter;
          }

          // Check for a filename parameter -> Create GtkFileChooserButton.
          if (!found_valid_item && !cimg::strcasecmp(argtype,"file")) {
            GtkWidget *label = gtk_label_new(argname);
            gtk_widget_show(label);
            gtk_table_attach(GTK_TABLE(table),label,0,1,current_line,current_line+1,GTK_FILL,GTK_SHRINK,0,0);
            gtk_misc_set_alignment(GTK_MISC(label),0,0.5);
            GtkWidget *filechooser = gtk_file_chooser_button_new(argname,GTK_FILE_CHOOSER_ACTION_OPEN);
            gtk_widget_show(filechooser);
            cimg::strclean(argarg);
            if (!reset && *s_value) gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(filechooser),s_value);
            else gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(filechooser),argarg);
            gtk_table_attach(GTK_TABLE(table),filechooser,1,3,current_line,current_line+1,
                             (GtkAttachOptions)(GTK_EXPAND | GTK_FILL),(GtkAttachOptions)0,0,0);
            event_infos[2*current_parameter] = (void*)current_parameter;
            event_infos[2*current_parameter+1] = (void*)0;
            on_file_parameter_changed(GTK_FILE_CHOOSER_BUTTON(filechooser),(void*)(event_infos+2*current_parameter));
            g_signal_connect(filechooser,"file-set",G_CALLBACK(on_file_parameter_changed),
                             (void*)(event_infos+2*current_parameter));
            g_signal_connect_swapped(filechooser,"file-set",G_CALLBACK(gimp_preview_invalidate),preview);
            found_valid_item = true;
            ++current_parameter;
          }

          // Check for a color -> Create GtkColorButton.
          if (!found_valid_item && !cimg::strcasecmp(argtype,"color")) {
            GtkWidget *hbox = gtk_hbox_new(false,6);
            gtk_widget_show(hbox);
            gtk_table_attach(GTK_TABLE(table),hbox,0,2,current_line,current_line+1,GTK_FILL,GTK_SHRINK,0,0);
            GtkWidget *label = gtk_label_new(argname);
            gtk_widget_show(label);
            gtk_box_pack_start(GTK_BOX(hbox),label,false,false,0);
            GtkWidget *colorchooser = gtk_color_button_new();
            gtk_widget_show(colorchooser);
            gtk_color_button_set_title(GTK_COLOR_BUTTON(colorchooser),argname);
            gtk_box_pack_start(GTK_BOX(hbox),colorchooser,false,false,0);
            event_infos[2*current_parameter] = (void*)current_parameter;
            event_infos[2*current_parameter+1] = (void*)0;
            cimg::strclean(argarg);
            unsigned int red = 0, green = 0, blue = 0, alpha = 255;
            const int err = std::sscanf(argarg,"%u%*c%u%*c%u%*c%u",&red,&green,&blue,&alpha);
            if (!reset && std::sscanf(s_value,"%u%*c%u%*c%u%*c%u",&red,&green,&blue,&alpha)==err) {}
            GdkColor col;
            col.pixel = 0; col.red = red<<8; col.green = green<<8; col.blue = blue<<8;
            gtk_color_button_set_color(GTK_COLOR_BUTTON(colorchooser),&col);
            if (err==4) {
              gtk_color_button_set_use_alpha(GTK_COLOR_BUTTON(colorchooser),true);
              gtk_color_button_set_alpha(GTK_COLOR_BUTTON(colorchooser),alpha<<8);
            } else gtk_color_button_set_use_alpha(GTK_COLOR_BUTTON(colorchooser),false);
            on_color_parameter_changed(GTK_COLOR_BUTTON(colorchooser),(void*)(event_infos+2*current_parameter));
            g_signal_connect(colorchooser,"color-set",G_CALLBACK(on_color_parameter_changed),
                             (void*)(event_infos+2*current_parameter));
            g_signal_connect_swapped(colorchooser,"color-set",G_CALLBACK(gimp_preview_invalidate),preview);
            found_valid_item = true;
            ++current_parameter;
          }

          // Check for a note -> Create GtkLabel.
          if (!found_valid_item && !cimg::strcasecmp(argtype,"note")) {
            cimg::strclean(argarg);
            GtkWidget *label = gtk_label_new(NULL);
            cimg::strescape(argarg);
            strparenthesis(argarg);
            gtk_label_set_markup(GTK_LABEL(label),argarg);
            gtk_label_set_line_wrap(GTK_LABEL(label),true);
            gtk_widget_show(label);
            gtk_table_attach(GTK_TABLE(table),label,0,3,current_line,current_line+1,GTK_FILL,GTK_SHRINK,0,0);
            gtk_misc_set_alignment(GTK_MISC(label),0,0.5);
            found_valid_item = true;
          }

          if (!found_valid_item) {
            if (get_verbosity_level()>0)
              std::fprintf(stderr,"\n*** Plug-in 'gmic4gimp' : Found invalid parameter type '%s' for argument '%s'.\n",argtype,argname);
          } else ++current_line;
        } else break;
      }
      set_filter_nbparams(filter,current_parameter);
    }
  }
  gtk_container_add(GTK_CONTAINER(frame),table);
}

// Called when the selected filter changed (in the combo-box).
void on_filter_changed(GtkTreeSelection *selection, gpointer user_data) {
  user_data = 0;
  GtkTreeIter iter;
  GtkTreeModel *model;
  unsigned int choice = 0;
  if (gtk_tree_selection_get_selected(selection,&model,&iter))
    gtk_tree_model_get(model,&iter,0,&choice,-1);
  set_current_filter(choice);
  create_parameters_gui(false);
  return_create_dialog = true;
}

// Handle responses to the dialog window buttons.
void on_verbosity_level_changed(GtkComboBox *combobox, gpointer user_data) {
  user_data = 0;
  int value = 0;
  g_object_get(combobox,"active",&value,NULL);
  set_verbosity_level(value);
}

void on_dialog_reset_clicked(GtkButton *widget, gpointer data) {
  widget = 0; data = 0;
  create_parameters_gui(true);
  return_create_dialog = true;
}

void on_dialog_cancel_clicked(GtkButton *widget, gpointer data) {
  widget = 0; data = 0;
  return_create_dialog = false;
  gtk_main_quit();
}

void on_dialog_apply_clicked(GtkButton *widget, gpointer data) {
  widget = 0;
  GimpDrawable *drawable = (GimpDrawable*)data;
  process_image(drawable,0);
  return_create_dialog = false;
}

void on_dialog_ok_clicked(GtkButton *widget, gpointer data) {
  widget = 0; data = 0;
  gtk_main_quit();
}

void on_update_button_clicked(GtkButton *widget, gpointer data) {
  widget = 0;
  GtkWidget *dialog = (GtkWidget*)data;
  char update_filename[1024] = { 0 }, update_command[1024] = { 0 }, src_filename[1024] = { 0 }, dest_filename[1024] = { 0 };
  const char
    *const update_url = "http://www.greyc.ensicaen.fr/~dtschump",
    *const path_tmp = cimg::temporary_path();
  std::sprintf(update_filename,"gmic4gimp_def.%d",gmic_version);
  std::sprintf(src_filename,"%s/%s",path_tmp,update_filename);
  std::sprintf(dest_filename,"%s/.%s",path_home,update_filename);
  if (get_verbosity_level()>0) {
    std::sprintf(update_command,"wget %s/%s -O %s",update_url,update_filename,src_filename);
    std::fprintf(stderr,"\n*** Plug-in 'gmic4gimp' : Running update procedure, with command : %s\n",update_command);
  } else std::sprintf(update_command,"wget --quiet %s/%s -O %s",update_url,update_filename,src_filename);
  int status = cimg::system(update_command);
  status = 0;
  std::FILE *file_s = std::fopen(src_filename,"r");
  bool succeed = false;
  if (file_s) {
    unsigned int size_s = 0;
    std::fseek(file_s,0,SEEK_END);
    size_s = (unsigned int)std::ftell(file_s);
    std::rewind(file_s);
    if (size_s) {
      std::FILE *file_d = std::fopen(dest_filename,"w");
      char *buffer = new char[size_s], sep = 0;
      if (file_d &&
          std::fread(buffer,sizeof(char),size_s,file_s)==size_s &&
          std::sscanf(buffer,"#@gim%c",&sep)==1 && sep=='p' &&
          std::fwrite(buffer,sizeof(char),size_s,file_d)==size_s) { succeed = true; std::fclose(file_d); }
      delete[] buffer;
    }
    std::fclose(file_s);
  }
  if (!succeed) {
    GtkWidget *message = gtk_message_dialog_new_with_markup(GTK_WINDOW(dialog),GTK_DIALOG_MODAL,GTK_MESSAGE_ERROR,GTK_BUTTONS_OK,
                                                            "<b>Filters update failed !</b>\n\n"
                                                            "A valid version of the update file :\n\n"
                                                            "<i>%s/%s</i>\n\n"
                                                            "  ...could not be retrieved from the G'MIC server.\n\n"
                                                            "Please check your Internet connection or\n"
                                                            "try a manual update instead.",update_url,update_filename);
    gtk_widget_show(message);
    gtk_dialog_run(GTK_DIALOG(message));
    gtk_widget_destroy(message);
  } else {
    GtkWidget *message = gtk_message_dialog_new_with_markup(GTK_WINDOW(dialog),GTK_DIALOG_MODAL,GTK_MESSAGE_INFO,GTK_BUTTONS_OK,
                                                "<b>Filters update succeed !</b>\n\n"
                                                "The G'MIC Toolbox must be restarted now.");
    gtk_widget_show(message);
    gtk_dialog_run(GTK_DIALOG(message));
    gtk_widget_destroy(message);
    return_create_dialog = false;
    set_current_filter(0);
    gtk_main_quit();
  }
}

// Create main plug-in dialog window and wait for a response.
//-----------------------------------------------------------
bool create_dialog_gui(GimpDrawable *drawable) {

  // Init GUI_specific variables
  gimp_ui_init("gmic",true);
  event_infos = 0;

  // Create main plug-in dialog window.
  GtkWidget
    *dialog = gimp_dialog_new("The G'MIC Toolbox","gmic",0,(GtkDialogFlags)0,gimp_standard_help_func,"gmic",NULL),
    *cancel_button = gimp_dialog_add_button(GIMP_DIALOG(dialog),GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL),
    *reset_button = gimp_dialog_add_button(GIMP_DIALOG(dialog),GIMP_STOCK_RESET,1),
    *apply_button = gimp_dialog_add_button(GIMP_DIALOG(dialog),GTK_STOCK_APPLY,GTK_RESPONSE_APPLY),
    *ok_button = gimp_dialog_add_button(GIMP_DIALOG(dialog),GTK_STOCK_OK,GTK_RESPONSE_OK);
  gimp_window_set_transient(GTK_WINDOW(dialog));
  g_signal_connect(dialog,"close",G_CALLBACK(on_dialog_cancel_clicked),0);
  g_signal_connect(dialog,"delete-event",G_CALLBACK(on_dialog_cancel_clicked),0);
  g_signal_connect(cancel_button,"clicked",G_CALLBACK(on_dialog_cancel_clicked),0);
  g_signal_connect(apply_button,"clicked",G_CALLBACK(on_dialog_apply_clicked),drawable);
  g_signal_connect(ok_button,"clicked",G_CALLBACK(on_dialog_ok_clicked),0);

  GtkWidget *dialog_hbox = gtk_hbox_new(false,0);
  gtk_widget_show(dialog_hbox);
  gtk_container_add(GTK_CONTAINER(GTK_DIALOG(dialog)->vbox),dialog_hbox);

  // Create the left pane, containing the preview, the show commmand and the update buttons and the author name.
  GtkWidget *left_pane = gtk_vbox_new(false,4);
  gtk_widget_show(left_pane);
  gtk_box_pack_start(GTK_BOX(dialog_hbox),left_pane,true,true,0);

  GtkWidget *preview = gimp_zoom_preview_new(drawable);
  gtk_widget_show(preview);
  gtk_box_pack_start(GTK_BOX(left_pane),preview,true,true,0);
  gimp_set_data("gmic_gui_preview",&preview,sizeof(GtkWidget*));
  g_signal_connect(preview,"invalidated",G_CALLBACK(process_preview),0);
  g_signal_connect_swapped(apply_button,"clicked",G_CALLBACK(gimp_preview_invalidate),preview);

  GtkWidget *verbosity_hbuttonbox = gtk_hbutton_box_new();
  gtk_widget_show(verbosity_hbuttonbox);
  gtk_box_pack_start(GTK_BOX(left_pane),verbosity_hbuttonbox,false,false,0);

  GtkWidget *verbosity_combobox = gtk_combo_box_new_text();
  gtk_widget_show(verbosity_combobox);
  gtk_combo_box_append_text(GTK_COMBO_BOX(verbosity_combobox),"Quiet mode");
  gtk_combo_box_append_text(GTK_COMBO_BOX(verbosity_combobox),"Verbose mode");
  gtk_combo_box_append_text(GTK_COMBO_BOX(verbosity_combobox),"Debug mode");
  gtk_combo_box_set_active(GTK_COMBO_BOX(verbosity_combobox),get_verbosity_level());
  gtk_container_add(GTK_CONTAINER(verbosity_hbuttonbox),verbosity_combobox);
  g_signal_connect(verbosity_combobox,"changed",G_CALLBACK(on_verbosity_level_changed),0);

  GtkWidget *update_hbuttonbox = gtk_hbutton_box_new();
  gtk_widget_show(update_hbuttonbox);
  gtk_box_pack_start(GTK_BOX(left_pane),update_hbuttonbox,false,false,0);
  GtkWidget
    *tmp_button = gtk_button_new_from_stock(GTK_STOCK_REFRESH),
    *update_image = gtk_button_get_image(GTK_BUTTON(tmp_button)),
    *update_button = gtk_button_new_with_mnemonic("_Update filters");
  gtk_button_set_image(GTK_BUTTON(update_button),update_image);
  gtk_widget_show(update_button);
  gtk_container_add(GTK_CONTAINER(update_hbuttonbox),update_button);
  g_signal_connect(update_button,"clicked",G_CALLBACK(on_update_button_clicked),(void*)dialog);

  GtkWidget *about_label = gtk_label_new(NULL);
  gtk_label_set_markup(GTK_LABEL(about_label),
                       "\n<span color=\"#666666\"><small>"
                       "<b>G'MIC</b> is proposed to you\n"
                       "   by <i>David Tschumperle</i>"
                       "</small></span>");
  gtk_widget_show(about_label);
  gtk_box_pack_start(GTK_BOX(left_pane),about_label,false,false,0);

  const unsigned int logo_width = 102, logo_height = 22;
  GdkPixbuf *pixbuf = gdk_pixbuf_new_from_data(data_logo,GDK_COLORSPACE_RGB,false,8,
                                               logo_width,logo_height,3*logo_width,0,0);
  GtkWidget *image = gtk_image_new_from_pixbuf(pixbuf);
  gtk_widget_show(image);
  gtk_box_pack_start(GTK_BOX(left_pane),image,false,false,0);

  // Create the middle pane, which contains the filters treeview.
  GtkWidget *middle_pane = gtk_frame_new(NULL);
  gtk_widget_show(middle_pane);
  gtk_container_set_border_width(GTK_CONTAINER(middle_pane),4);
  gtk_widget_set_size_request(middle_pane,250,-1);
  gtk_box_pack_start(GTK_BOX(dialog_hbox),middle_pane,false,false,0);

  GtkWidget *scrolledwindow = gtk_scrolled_window_new(NULL,NULL);
  gtk_widget_show(scrolledwindow);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolledwindow),GTK_POLICY_AUTOMATIC,GTK_POLICY_AUTOMATIC);
  gtk_container_add(GTK_CONTAINER(middle_pane),scrolledwindow);

  GtkWidget *treeview = gtk_tree_view_new_with_model(GTK_TREE_MODEL(filter_store));
  GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
  GtkTreeViewColumn *column = gtk_tree_view_column_new_with_attributes(" Available filters :",renderer,"text",1,NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(treeview),column);

  GtkTreeSelection *select = gtk_tree_view_get_selection(GTK_TREE_VIEW(treeview));
  gtk_tree_selection_set_mode(select,GTK_SELECTION_SINGLE);
  g_signal_connect(G_OBJECT(select),"changed",G_CALLBACK(on_filter_changed),0);
  g_signal_connect_swapped(select,"changed",G_CALLBACK(gimp_preview_invalidate),preview);
  gtk_widget_show(treeview);
  gtk_container_add(GTK_CONTAINER(scrolledwindow),treeview);
  g_signal_connect(reset_button,"clicked",G_CALLBACK(on_dialog_reset_clicked),select);
  g_signal_connect_swapped(reset_button,"clicked",G_CALLBACK(gimp_preview_invalidate),preview);

  // Create the right pane which contains the parameters frame.
  GtkWidget *parameters_frame = gtk_frame_new(NULL);
  gtk_widget_show(parameters_frame);
  gtk_container_set_border_width(GTK_CONTAINER(parameters_frame),4);
  gtk_widget_set_size_request(parameters_frame,450,-1);
  gtk_box_pack_start(GTK_BOX(dialog_hbox),parameters_frame,false,false,0);
  gimp_set_data("gmic_gui_frame",&parameters_frame,sizeof(GtkWidget*));
  create_parameters_gui(false);

  // Show dialog window and wait for user response.
  gtk_widget_show(dialog);
  gtk_main();

  // Destroy dialog box widget and free resources.
  gtk_widget_destroy(dialog);
  gtk_widget_destroy(tmp_button);
  if (event_infos) delete[] event_infos;
  return return_create_dialog;
}

// 'Run' function needed by GIMP plug-in API.
//-------------------------------------------
void gmic_run(const gchar *name, gint nparams, const GimpParam *param, gint *nreturn_vals, GimpParam **return_vals) {

  // Init plug-in variables.
  static GimpParam values[1];
  values[0].type = GIMP_PDB_STATUS;
  *return_vals  = values;
  *nreturn_vals = 1;
  name = 0;
  nparams = 0;
  GimpRunMode run_mode;
  run_mode = (GimpRunMode)param[0].data.d_int32;
  if (run_mode==GIMP_RUN_NONINTERACTIVE) {
    std::fprintf(stderr,"\n*** Plug-in 'gmic4gimp' : ERROR, this plug-in cannot be run in non-interactive mode.\n");
    values[0].data.d_status = GIMP_PDB_CALLING_ERROR;
    return;
  }
  gmic_macros = 0;
  filter_store = 0;
  return_create_dialog = true;
  path_home = getenv(cimg_OS!=2?"HOME":"APPDATA");

  // Check that no instance of the plug-in is already running.
  bool is_existing_instance = 0;
  gimp_get_data("gmic_instance",&is_existing_instance);
  if (is_existing_instance) {
    std::fprintf(stderr,"\n*** Plug-in 'gmic4gimp' : Existing instance of the plug-in is already running.\n");
    return;
  }
  is_existing_instance = true;
  gimp_set_data("gmic_instance",&is_existing_instance,sizeof(bool));

  // Read user-defined configuration files '.gmic_def' and '.gmic', when possible.
  unsigned size_update = 0, size_custom = 0, size_def = sizeof(data_gmic4gimp_def);
  char filename_update[1024] = { 0 }, filename_custom[1024] = { 0 };
  std::sprintf(filename_update,"%s/.gmic4gimp_def.%d",path_home,gmic_version);
  std::sprintf(filename_custom,"%s/.gmic4gimp",path_home);
  std::FILE
    *file_update = std::fopen(filename_update,"r"),
    *file_custom = std::fopen(filename_custom,"r");
  if (file_update) {
    std::fseek(file_update,0,SEEK_END);
    size_update = (unsigned int)std::ftell(file_update);
    std::rewind(file_update);
  }
  if (file_custom) {
    std::fseek(file_custom,0,SEEK_END);
    size_custom = (unsigned int)std::ftell(file_custom);
    std::rewind(file_custom);
  }
  const unsigned int size_final = size_update + size_custom + size_def + 1;
  char *ptrd = gmic_macros = new char[size_final];
  if (size_custom) { ptrd+=std::fread(ptrd,1,size_custom,file_custom); std::fclose(file_custom); }
  if (size_update) { ptrd+=std::fread(ptrd,1,size_update,file_update); std::fclose(file_update); }
  if (size_def)    { std::memcpy(ptrd,data_gmic4gimp_def,size_def); ptrd+=size_def; }
  *ptrd = 0;

  // Parse available G'MIC filters definitions.
  GtkTreeIter iter, parent[16];
  filter_store = gtk_tree_store_new(2,G_TYPE_UINT,G_TYPE_STRING);
  char line[256*1024] = { 0 }, entry[4096] = { 0 }, command[4096] = { 0 };
  char preview_command[4096] = { 0 }, arguments[4096] = { 0 };
  int level = 0;
  for (const char *data = gmic_macros; *data; ) {
    if (*data=='\n') ++data;
    else {
      if (std::sscanf(data,"%262143[^\n]\n",line)>0) data += cimg::strlen(line) + 1;
      arguments[0] = 0;
      if (line[0]=='#') {
        const int err = std::sscanf(line,"#@gimp %4095[^:]: %4095[^, ]%*c %4095[^, ]%*c %4095[^\n]",
                                    entry,command,preview_command,arguments);
        strparenthesis(entry);
        if (err==1) { // If entry is a menu folder.
          cimg::strclean(entry);
          char *nentry = entry;
          while (*nentry=='_') { ++nentry; --level; }
          if (level<0) level = 0;
          if (level>15) level = 15;
          cimg::strclean(nentry);
          if (*nentry) {
            gtk_tree_store_append(filter_store,&parent[level],level?&parent[level-1]:0);
            gtk_tree_store_set(filter_store,&parent[level],0,0,1,nentry,-1);
            ++level;
          }
        } else if (err>=2) { // If entry is a regular filter.
          cimg::strclean(entry);
          cimg::strclean(command);
          gmic_entries.insert(CImg<char>(entry,cimg::strlen(entry)+1));
          gmic_commands.insert(CImg<char>(command,cimg::strlen(command)+1));
          gmic_arguments.insert(CImg<char>(arguments,cimg::strlen(arguments)+1));
          if (err>=3) {
            cimg::strclean(preview_command);
            gmic_preview_commands.insert(CImg<char>(preview_command,cimg::strlen(preview_command)+1));
          }
          gtk_tree_store_append(filter_store,&iter,level?&parent[level-1]:0);
          gtk_tree_store_set(filter_store,&iter,0,gmic_entries.size,1,entry,-1);
        }
      }
    }
  }

  // Get currenty selected drawable and run image filter on it.
  GimpDrawable *drawable = gimp_drawable_get(param[2].data.d_drawable);
  gimp_tile_cache_ntiles(2*(drawable->width/gimp_tile_width()+1));
  if (run_mode==GIMP_RUN_INTERACTIVE) {
    if (create_dialog_gui(drawable)) {
      process_image(drawable,0);
      const char *commandline = get_commandline(false);
      if (commandline) { // Remember command line for the next use of the filter.
        char s_tmp[256] = { 0 };
        std::sprintf(s_tmp,"gmic_commandline%u",get_current_filter());
        gimp_set_data(s_tmp,commandline,cimg::strlen(commandline));
      }
    }
  } else if (run_mode==GIMP_RUN_WITH_LAST_VALS) {
    const unsigned int filter = get_current_filter();
    if (filter) {
      char s_tmp[256] = { 0 };
      std::sprintf(s_tmp,"gmic_commandline%u",filter);
      char commandline[4096] = { 0 };
      gimp_get_data(s_tmp,&commandline);
      process_image(drawable,commandline);
    }
  }

  // Free plug-in resources.
  delete[] gmic_macros;
  values[0].data.d_status = GIMP_PDB_SUCCESS;
  is_existing_instance = false;
  gimp_set_data("gmic_instance",&is_existing_instance,sizeof(bool));
}

// 'Query' function needed by GIMP plug-in API.
//---------------------------------------------
void gmic_query() {
  static const GimpParamDef args[] = {
    {GIMP_PDB_INT32,    "run_mode", "Interactive, non-interactive"},
    {GIMP_PDB_IMAGE,    "image", "(unused)"},
    {GIMP_PDB_DRAWABLE, "drawable", "Drawable to draw on"}
  };

  gimp_install_procedure("gmic",                     // name
                         "G'MIC Toolbox",            // blurb
                         "G'MIC Toolbox",            // help
                         "David Tschumperle",        // author
                         "David Tschumperle",        // copyright
                         "2008-12-02",               // date
                         "_G'MIC Toolbox...",        // menu_path
                         "RGB*, GRAY*",              // image_types
                         GIMP_PLUGIN,                // type
                         G_N_ELEMENTS(args),         // nparams
                         0,                          // nreturn_vals
                         args,                       // params
                         0);                         // return_vals

  gimp_plugin_menu_register("gmic", "<Image>/Filters");
}

GimpPlugInInfo PLUG_IN_INFO = { 0, 0, gmic_query, gmic_run };
MAIN();
