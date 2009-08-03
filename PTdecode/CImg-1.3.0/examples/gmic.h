/*
  #
  #  File        : gmic.h
  #                ( C++ header file )
  #
  #  Description : GREYC's Magic Image Converter
  #                ( http://gmic.sourceforge.net )
  #                This file is a part of the CImg Library project.
  #                ( http://cimg.sourceforge.net )
  #
  #  Note        : This file cannot be compiled on VC++ 6.
  #
  #  Copyright   : David Tschumperle
  #                ( http://www.greyc.ensicaen.fr/~dtschump/ )
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
  #  data to be ensured and,  more generally, to use and operate it in the
  #  same conditions as regards security.
  #
  #  The fact that you are presently reading this means that you have had
  #  knowledge of the CeCILL license and that you accept its terms.
  #
*/

#ifndef gmic_version
#include "CImg.h"
#define gmic_version 1304

// The lines below are necessary when using a non-standard compiler such as visualcpp6.
#ifdef cimg_use_visualcpp6
#define std
#endif
#ifdef min
#undef min
#undef max
#endif

// Define G'MIC Exception class.
//------------------------------
struct gmic_exception {
  char message[4096];
  gmic_exception();
  gmic_exception(const char *format, ...);
  gmic_exception(const char *format, std::va_list ap);
};

// Define G'MIC interpreter class.
//--------------------------------
struct gmic {

  // Internal variables.
  cimg_library::CImgList<char> command_line, filenames, macros, commands;
  cimg_library::CImgList<int> dowhile, repeatdone;
  bool is_released, is_debug, is_fullpath, is_begin, is_oriented3d;
  int verbosity_level, render3d, renderd3d;
  float focale3d, light3d_x, light3d_y, light3d_z, specular_light3d, specular_shine3d;
  unsigned char background3d[3];
  unsigned int position;

  // Constructors - Destructors.
  gmic();
  template<typename T> gmic(const int argc, const char *const *const argv, cimg_library::CImgList<T>& images,
                            const char *const custom_macros=0, const bool add_macros_at_start=true);
  template<typename T> gmic(const char *const command, cimg_library::CImgList<T>& images,
                            const char *const custom_macros=0, const bool add_macros_at_start=true);
  gmic& assign(const unsigned int size, const char *const custom_macros=0,
               const bool add_macros_at_start=true);

  // Messages procedures.
  const gmic& error(const char *format, ...) const;
  const gmic& warning(const char *format, ...) const;
  const gmic& debug(const char *format, ...) const;
  const gmic& print(const char *format, ...) const;

  // Add macros.
  gmic& add_macros(const char *const data_macros, const unsigned int data_size, const bool add_macros_at_start=true);
  gmic& add_macros(std::FILE *const file, const bool add_macros_at_start=true);

  // Return indices of the images from a string.
  cimg_library::CImg<unsigned int> indices2cimg(const char *const string, const unsigned int indice_max,
                                                const char *const command) const;

  // Return stringified version of indices or filenames.
  char* indices2string(const cimg_library::CImg<unsigned int>& indices, const bool display_indices) const;

  // Display image data.
  template<typename T>
  bool display_images(const cimg_library::CImgList<T>& images, const cimg_library::CImg<unsigned int>& indices,
                      const bool verbose) const;
  template<typename T>
  bool display_objects3d(const cimg_library::CImgList<T>& images, const cimg_library::CImg<unsigned int>& indices,
                         const bool verbose) const;
  template<typename T>
  bool display_plots(const cimg_library::CImgList<T>& images, const cimg_library::CImg<unsigned int>& indices,
                     const unsigned int plot_type, const unsigned int vertex_type,
                     const double xmin, const double xmax,
                     const double ymin, const double ymax,
                     const bool verbose) const;

  // Substitute '@' expressions.
  template<typename T>
  cimg_library::CImg<char> substitute_arobace(const char *const argument, const cimg_library::CImgList<T>& images) const;

  // Main parsing procedure.
  template<typename T>
  gmic& parse(cimg_library::CImgList<T> &images);
  gmic& parse_bool(cimg_library::CImgList<bool>& images);
  gmic& parse_uchar(cimg_library::CImgList<unsigned char>& images);
  gmic& parse_char(cimg_library::CImgList<char>& images);
  gmic& parse_ushort(cimg_library::CImgList<unsigned short>& images);
  gmic& parse_short(cimg_library::CImgList<short>& images);
  gmic& parse_uint(cimg_library::CImgList<unsigned int>& images);
  gmic& parse_int(cimg_library::CImgList<int>& images);
  gmic& parse_float(cimg_library::CImgList<float>& images);
  gmic& parse_double(cimg_library::CImgList<double>& images);

}; // End of the 'gmic' class.

#endif

// Local Variables:
// mode: c++
// End:
