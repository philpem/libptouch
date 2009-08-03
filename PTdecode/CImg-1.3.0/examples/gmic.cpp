/*
 #
 #  File        : gmic.cpp
 #                ( C++ source file )
 #
 #  Description : GREYC's Magic Image Converter (library and executable)
 #                ( http://gmic.sourceforge.net )
 #                This file is a part of the CImg Library project.
 #                ( http://cimg.sourceforge.net )
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

// Add specific G'MIC methods to the CImg<T> class.
//-------------------------------------------------
#ifdef cimg_plugin

template<typename t>
CImg<T> get_replace(const CImg<t>& img) const {
  return +img;
}

template<typename t>
CImg<T>& replace(CImg<t>& img) {
  return img.transfer_to(*this);
}

CImg<T>& gmic_set(const double value, const int x, const int y, const int z, const int v) {
  (*this).atXYZV(x,y,z,v,0) = (T)value;
  return *this;
}

CImg<T> get_gmic_set(const double value, const int x, const int y, const int z, const int v) const {
  return (+*this).gmic_set(value,x,y,z,v);
}

CImg<T> get_draw_point(const int x, const int y, const int z,
                       const CImg<T>& col, const float opacity) const {
  return (+*this).draw_point(x,y,z,col,opacity);
}

CImg<T> get_draw_line(const int x0, const int y0, const int x1, const int y1,
                      const CImg<T>& col, const float opacity) const {
  return (+*this).draw_line(x0,y0,x1,y1,col,opacity);
}

template<typename t>
CImg<T> get_draw_polygon(const CImg<t>& pts, const CImg<T>& col, const float opacity) const {
  return (+*this).draw_polygon(pts,col,opacity);
}

CImg<T> get_draw_ellipse(const int x, const int y, const float r0, const float r1,
                         const float ru, const float rv, const CImg<T>& col, const float opacity) const {
  return (+*this).draw_ellipse(x,y,r0,r1,ru,rv,col,opacity);
}

CImg<T> get_draw_text(const int x, const int y, const char *const text, const T *const col,
                      const int bg, const float opacity,const int siz) const {
  return (+*this).draw_text(x,y,text,col,bg,opacity,siz);
}

CImg<T> get_draw_image(const int x, const int y, const int z,
                       const CImg<T>& sprite, const CImg<T>& mask, const float opacity) const {
  return (+*this).draw_image(x,y,z,sprite,mask,opacity);
}

CImg<T> get_draw_image(const int x, const int y, const int z,
                       const CImg<T>& sprite, const float opacity) const {
  return (+*this).draw_image(x,y,z,sprite,opacity);
}

CImg<T> get_draw_plasma(const float alpha, const float beta, const float opacity) const {
  return (+*this).draw_plasma(alpha,beta,opacity);
}

CImg<T> get_draw_mandelbrot(const CImg<T>& color_palette, const float opacity,
                            const double z0r, const double z0i, const double z1r, const double z1i,
                            const unsigned int itermax, const bool normalized_iteration,
                            const bool julia_set, const double paramr, const double parami) const {
  return (+*this).draw_mandelbrot(color_palette,opacity,z0r,z0i,z1r,z1i,itermax,
                                  normalized_iteration,julia_set,paramr,parami);
}

CImg<T> get_draw_fill(const int x, const int y, const int z,
                      const CImg<T>& col, const float opacity, const float tolerance) const {
  return (+*this).draw_fill(x,y,z,col,opacity,tolerance);
}

bool is_CImg3d() const {
  const bool is_header = (width==1 && height>=8 && depth==1 && dim==1 &&
                          (*this)[0]=='C'+0.5f && (*this)[1]=='I'+0.5f &&
                          (*this)[2]=='m'+0.5f && (*this)[3]=='g'+0.5f &&
                          (*this)[4]=='3'+0.5f && (*this)[5]=='d'+0.5f);
  if (!is_header) return false;
  const int
    nbv = (int)(*this)[6],
    nbp = (int)(*this)[7];
  if (nbv<=0 || nbp<=0) return false;
  const T *ptrs = ptr() + 8 + 3*nbv, *const ptre = end();
  if (ptrs>=ptre) return false;
  for (int i = 0; i<nbp && ptrs<ptre; ++i) {
    const int N = (int)*(ptrs++);
    if (N<=0 || N>=8) return false;
    ptrs+=N;
  }
  ptrs+=4*nbp;
  if (ptrs>ptre) return false;
  return true;
}

template<typename tp, typename tf, typename tc, typename to>
CImg<T> get_draw_object3d(const float x0, const float y0, const float z0,
                          const CImg<tp>& points, const CImgList<tf>& primitives,
                          const CImgList<tc>& colors, const CImg<to>& opacities,
                          const unsigned int render_type, const bool double_sided,
                          const float focale, const float lightx, const float lighty,
                          const float lightz, const float specular_light, const float specular_shine,
                          float *const zbuffer) const {
  return (+*this).draw_object3d(x0,y0,z0,points,primitives,colors,opacities,render_type,double_sided,focale,
                                lightx,lighty,lightz,specular_light,specular_shine,zbuffer);
}

template<typename tp, typename tc, typename to>
CImg<T>& object3dtoCImg3d(CImgList<tp>& primitives, CImgList<tc>& colors, CImg<to>& opacities) {
  if (is_empty() || !primitives) { primitives.assign(); colors.assign(); opacities.assign(); return *this; }
  const unsigned int primitives_size = primitives.size;
  CImgList<floatT> res;
  res.insert(CImg<floatT>("CImg3d",1,6,1,1,false)+=0.5f);
  res.insert(CImg<floatT>::vector((float)width,(float)primitives.size));
  res.insert(1); resize(-100,3,1,1,0).transpose().unroll('y').transfer_to(res.last());
  cimglist_for(primitives,p) {
    res.insert(CImg<floatT>::vector((float)primitives[p].size())).insert(primitives[p]).last().unroll('y');
    primitives[p].assign();
  }
  primitives.assign();
  const unsigned int defined_colors = colors.size;
  cimglist_for(colors,c) { res.insert(colors[c]).last().resize(1,3,1,1,-1); colors[c].assign(); }
  colors.assign();
  if (defined_colors<primitives_size) res.insert(1).last().assign(1,3*(primitives_size-defined_colors),1,1,200);
  const unsigned int defined_opacities = opacities.size();
  res.insert(opacities).last().unroll('y');
  opacities.assign();
  if (defined_opacities<primitives.size) res.insert(1).last().assign(1,primitives_size-defined_opacities,1,1,1);
  return res.get_append('y').transfer_to(*this);
}

template<typename tp, typename tc, typename to>
CImg<T>& CImg3dtoobject3d(CImgList<tp>& primitives, CImgList<tc>& colors, CImg<to>& opacities) {
  const T *ptrs = ptr() + 6;
  const unsigned int
    nbv = (unsigned int)*(ptrs++),
    nbp = (unsigned int)*(ptrs++);
  CImg<T> points(nbv,3);
  primitives.assign(nbp);
  colors.assign(nbp,1,3,1,1);
  opacities.assign(nbp);
  cimg_forX(points,x) { points(x,0) = (T)*(ptrs++); points(x,1) = (T)*(ptrs++); points(x,2) = (T)*(ptrs++); }
  cimglist_for(primitives,p) {
    const unsigned int N = (unsigned int)*(ptrs++);
    primitives[p].assign(ptrs,1,N,1,1,false);
    ptrs+=N;
  }
  cimglist_for(colors,c) { colors(c,0) = (tc)*(ptrs++); colors(c,1) = (tc)*(ptrs++); colors(c,2) = (tc)*(ptrs++); }
  opacities.assign(ptrs,1,nbp,1,1,false);
  return assign(points);
}

CImg<T> get_appendCImg3d(const CImg<T>& img) const {
  CImg<T> res(1,img.size() + size() - 8);
  const T *ptrs = ptr() + 6, *ptrs0 = img.ptr() + 6;
  T *ptrd = res.ptr();
  *(ptrd++) = (T)('C' + 0.5f); *(ptrd++) = (T)('I' + 0.5f);
  *(ptrd++) = (T)('m' + 0.5f); *(ptrd++) = (T)('g' + 0.5f);
  *(ptrd++) = (T)('3' + 0.5f); *(ptrd++) = (T)('d' + 0.5f);
  const unsigned int
    nbv = (unsigned int)*(ptrs++),
    nbv0 = (unsigned int)*(ptrs0++),
    nbp = (unsigned int)*(ptrs++),
    nbp0 = (unsigned int)*(ptrs0++);
  *(ptrd++) = (T)(nbv + nbv0);
  *(ptrd++) = (T)(nbp + nbp0);
  std::memcpy(ptrd,ptrs,sizeof(T)*nbv*3);
  ptrd+=3*nbv; ptrs+=3*nbv;
  std::memcpy(ptrd,ptrs0,sizeof(T)*nbv0*3);
  ptrd+=3*nbv0; ptrs0+=3*nbv0;
  for (unsigned int i = 0; i<nbp; ++i) {
    const unsigned int N = (unsigned int)*(ptrs++);
    *(ptrd++) = (T)N;
    std::memcpy(ptrd,ptrs,sizeof(T)*N);
    ptrd+=N; ptrs+=N;
  }
  for (unsigned int i = 0; i<nbp0; ++i) {
    const unsigned int N = (unsigned int)*(ptrs0++);
    *(ptrd++) = (T)N;
    for (unsigned int j = 0; j<N; ++j) *(ptrd++) = (T)(*(ptrs0++) + nbv);
  }
  std::memcpy(ptrd,ptrs,sizeof(T)*nbp*3);
  ptrd+=3*nbp; ptrs+=3*nbp;
  std::memcpy(ptrd,ptrs0,sizeof(T)*nbp0*3);
  ptrd+=3*nbp0; ptrs0+=3*nbp0;
  std::memcpy(ptrd,ptrs,sizeof(T)*nbp);
  ptrd+=nbp;
  std::memcpy(ptrd,ptrs0,sizeof(T)*nbp0);
  return res;
}

CImg<T>& appendCImg3d(const CImg<T>& img) {
  return get_appendCImg3d(img).transfer_to(*this);
}

CImg<T>& centerCImg3d() {
  const unsigned int nbv = (unsigned int)(*this)[6];
  const T *ptrs = ptr() + 8;
  float xm = cimg::type<float>::max(), ym = xm, zm = xm, xM = cimg::type<float>::min(), yM = xM, zM = xM;
  for (unsigned int i = 0; i<nbv; ++i) {
    const float x = (float)*(ptrs++), y = (float)*(ptrs++), z = (float)*(ptrs++);
    if (x<xm) xm = x; if (x>xM) xM = x;
    if (y<ym) ym = y; if (y>yM) yM = y;
    if (z<zm) zm = z; if (z>zM) zM = z;
  }
  const float xc = (xm + xM)/2, yc = (ym + yM)/2, zc = (zm + zM)/2;
  T *ptrd = ptr() + 8;
  for (unsigned int i = 0; i<nbv; ++i) { *(ptrd++)-=(T)xc; *(ptrd++)-=(T)yc; *(ptrd++)-=(T)zc; }
  return *this;
}

CImg<T> get_centerCImg3d() const {
  return (+*this).centerCImg3d();
}

CImg<T>& normalizeCImg3d() {
  const unsigned int nbv = (unsigned int)(*this)[6];
  const T *ptrs = ptr() + 8;
  float xm = cimg::type<float>::max(), ym = xm, zm = xm, xM = cimg::type<float>::min(), yM = xM, zM = xM;
  for (unsigned int i = 0; i<nbv; ++i) {
    const float x = (float)*(ptrs++), y = (float)*(ptrs++), z = (float)*(ptrs++);
    if (x<xm) xm = x; if (x>xM) xM = x;
    if (y<ym) ym = y; if (y>yM) yM = y;
    if (z<zm) zm = z; if (z>zM) zM = z;
  }
  const float delta = cimg::max(xM-xm,yM-ym,zM-zm);
  if (delta>0) {
    T *ptrd = ptr() + 8;
    for (unsigned int i = 0; i<3*nbv; ++i) *(ptrd++)/=(T)delta;
  }
  return *this;
}

CImg<T> get_normalizeCImg3d() const {
  return (+*this).normalizeCImg3d();
}

template<typename t>
CImg<T>& rotateCImg3d(const CImg<t>& rot) {
  const unsigned int nbv = (unsigned int)(*this)[6];
  const T *ptrs = ptr() + 8;
  const float
    a = (float)rot(0,0), b = (float)rot(1,0), c = (float)rot(2,0),
    d = (float)rot(0,1), e = (float)rot(1,1), f = (float)rot(2,1),
    g = (float)rot(0,2), h = (float)rot(1,2), i = (float)rot(2,2);
  T *ptrd = ptr() + 8;
  for (unsigned int j = 0; j<nbv; ++j) {
    const float x = (float)*(ptrs++), y = (float)*(ptrs++), z = (float)*(ptrs++);
    *(ptrd++) = (T)(a*x + b*y + c*z);
    *(ptrd++) = (T)(d*x + e*y + f*z);
    *(ptrd++) = (T)(g*x + h*y + i*z);
  }
  return *this;
}

template<typename t>
CImg<T> get_rotateCImg3d(const CImg<t>& rot) const {
  return (+*this).rotateCImg3d(rot);
}

CImg<T>& translateCImg3d(const float tx, const float ty, const float tz) {
  const unsigned int nbv = (unsigned int)(*this)[6];
  T *ptrd = ptr() + 8;
  for (unsigned int j = 0; j<nbv; ++j) { *(ptrd++) += (T)tx; *(ptrd++) += (T)ty; *(ptrd++) += (T)tz; }
  return *this;
}

CImg<T> get_translateCImg3d(const float tx, const float ty, const float tz) const {
  return (+*this).translateCImg3d(tx,ty,tz);
}

CImg<T>& coloropacityCImg3d(const float R, const float G, const float B, const float opacity, const bool set_RGB, const bool set_opacity) {
  T *ptrd = ptr() + 6;
  const unsigned int
    nbv = (unsigned int)*(ptrd++),
    nbp = (unsigned int)*(ptrd++);
  ptrd+=3*nbv;
  for (unsigned int i = 0; i<nbp; ++i) { const unsigned int N = (unsigned int)*(ptrd++); ptrd+=N; }
  if (set_RGB) for (unsigned int c = 0; c<nbp; ++c) { *(ptrd++) = (T)R; *(ptrd++) = (T)G; *(ptrd++) = (T)B; } else ptrd+=3*nbp;
  if (set_opacity) for (unsigned int o = 0; o<nbp; ++o) *(ptrd++) = (T)opacity;
  return *this;
}

CImg<T> get_coloropacityCImg3d(const float R, const float G, const float B, const float opacity, const bool set_RGB, const bool set_opacity) const {
  return (+*this).coloropacityCImg3d(R,G,B,opacity,set_RGB,set_opacity);
}

#else  // eq. to #ifndef cimg_plugin

#define cimg_debug 1
#ifndef cimg_gmic_cpp
#define cimg_gmic_cpp "examples/gmic.cpp"
#define cimg_cimg_h "../CImg.h"
#endif
#define cimg_stdout stdout
#define cimg_plugin cimg_gmic_cpp
#include cimg_cimg_h
#include "gmic.h"
using namespace cimg_library;

// The lines below are necessary when using a non-standard compiler such as visualcpp6.
#ifdef cimg_use_visualcpp6
#define std
#endif
#ifdef min
#undef min
#undef max
#endif

#if !defined(gmic_main) || !defined(gmic_separate_compilation)

// Define some useful macros.
//---------------------------

// Code for validity checking of indices.
#define gmic_inds indices2string(indices,true)
#define gmic_check_indice(ind,funcname) { \
  const int indo = (int)ind; \
  if (ind<0) ind+=images.size; \
  if (ind<0 || ind>=(int)images.size) { \
    if (images.size) error(funcname " : Invalid indice '[%d]' (valid indice range is -%u...%u).",gmic_inds,indo,images.size,images.size-1); \
    else error(funcname " : Invalid indice '[%d]' (image list is empty).",gmic_inds,indo); \
  } \
}

// Code for having 'get' or 'non-get' versions of G'MIC commands.
#define gmic_apply(instance,function) { \
  if (get_version) { \
    unsigned int posi = 0; \
    if (images.contains(instance,posi)) filenames.insert(filenames[posi]); \
    else filenames.insert(CImg<char>("(gmic)",7,1,1,1,false)); \
    CImg<T> res = instance.get_##function; \
    images.insert(1); res.transfer_to(images.last()); \
  } else instance.function; \
}

// Code for simple commands that has no parameters and act on images.
#define gmic_simple_item(option,function,description) \
  if (!cimg::strcmp(option,item0)) { \
    print(description,gmic_inds); cimg_foroff(indices,l) gmic_apply(images[indices[l]],function()); \
    continue; \
}

// Code for the type cast command.
#define gmic_cast(pixel_type,st_type) \
  if (!cimg::strcmp(#pixel_type,argument)) { \
    print("Set pixel type to '%s'.",#pixel_type); ++position; \
    if (!cimg::strcmp(st_type,cimg::type<T>::string())) continue; \
    CImgList<pixel_type> casted_images; \
    while (images) { casted_images.insert(images[0]); images.remove(0); } \
    return parse_##pixel_type(casted_images); \
}

// Code for G'MIC arithmetic commands.
#define gmic_arithmetic_item(option1,option2,\
                             function1,description1,arg1_1,arg1_2,value_type1, \
                             function2,description2_1,description2_2,arg2_1,arg2_2,description3) \
 if (!cimg::strcmp(option1,item0) || !cimg::strcmp(option2,item0)) { \
   double value = 0; char inds[4096] = { 0 }, sep = 0, end = 0; \
    if (std::sscanf(argument,"%lf%c",&value,&end)==1) { \
      print(description1 ".",arg1_1,arg1_2); \
      cimg_foroff(indices,l) \
       if (get_version) { \
         images.insert(images[indices[l]]); images.last().function1((value_type1)value); \
         filenames.insert(filenames[indices[l]]); } \
       else images[indices[l]].function1((value_type1)value); \
      ++position; \
    } else if (std::sscanf(argument,"[%4095[0-9.eE%+-]%c%c",inds,&sep,&end)==2 && sep==']') { \
      const CImg<unsigned int> ind = indices2cimg(inds,images.size,option1); \
      if (ind.size()!=1) error(description2_1 " : Argument '[%s]' should contain one indice.",gmic_inds,inds); \
      print(description2_2 ".",arg2_1,arg2_2); \
      const CImg<T> img0 = images[ind[0]]; \
      cimg_foroff(indices,l) \
       if (get_version) { \
         images.insert(images[indices[l]]); images.last().function2(img0); \
         filenames.insert(filenames[indices[l]]); } \
       else images[indices[l]].function2(img0); \
      ++position; \
    } else { \
      print(description3 ".",gmic_inds); \
      if (images && indices) { \
        for (unsigned int siz = indices.size(), ind0 = indices[0], off = 0, l = 1; l<siz; ++l) { \
          const unsigned int ind = indices[l] - off; \
          images[ind0].function2(images[ind]); \
          images.remove(ind); filenames.remove(ind); \
          ++off; \
        }}} continue; \
}

// Constructors.
//--------------
#if defined(gmic_float) || !defined(gmic_separate_compilation)

#include "gmic_def.h"

gmic_exception::gmic_exception() {
  message[0] = '\0';
}

gmic_exception::gmic_exception(const char *format, ...) {
  std::va_list ap;
  va_start(ap,format);
  std::vsprintf(message,format,ap);
  va_end(ap);
}

gmic_exception::gmic_exception(const char *format, std::va_list ap) {
  std::vsprintf(message,format,ap);
}

gmic::gmic() {
  assign(0);
}

// Set default values of G'MIC parameters and macros.
//----------------------------------------------------
gmic& gmic::assign(const unsigned int size, const char *const custom_macros, const bool add_macros_start) {
  filenames.assign(size,CImg<char>("(gmic)",7,1,1,1,false));
  position = 0;
  verbosity_level = 0;
  is_released = true;
  is_debug = false;
  is_begin = true;
  background3d[0] = 120;
  background3d[1] = 120;
  background3d[2] = 140;
  render3d = 4;
  renderd3d = -1;
  is_oriented3d = false;
  focale3d = 500;
  light3d_x = 0;
  light3d_y = 0;
  light3d_z = -5000;
  specular_light3d = 0.15f;
  specular_shine3d = 0.8f;
  is_fullpath = false;
  add_macros(data_def,sizeof(data_def)-1,true);
  add_macros(custom_macros,cimg::strlen(custom_macros)-1,add_macros_start);
  return *this;
}

// Error procedure.
//-----------------
const gmic& gmic::error(const char *format, ...) const {
  va_list ap;
  va_start(ap,format);
  char message[1024] = { 0 };
  std::vsprintf(message,format,ap);
  va_end(ap);
  if (verbosity_level>=0) {
    std::fprintf(cimg_stdout,"\n<gmic-#%u> ** Error ** %s",filenames.size,message);
    std::fprintf(cimg_stdout,"\n<gmic-#%u> Abort G'MIC instance.\n",filenames.size);
    std::fflush(cimg_stdout);
  }
  throw gmic_exception(message);
  return *this;
}

// Warning procedure.
//-------------------
const gmic& gmic::warning(const char *format, ...) const {
  va_list ap;
  va_start(ap,format);
  if (verbosity_level>=0) {
    std::fprintf(cimg_stdout,"\n<gmic-#%u> ** Warning ** ",filenames.size);
    std::vfprintf(cimg_stdout,format,ap);
    std::fflush(cimg_stdout);
  }
  va_end(ap);
  return *this;
}

// Print debug messages.
//----------------------
const gmic& gmic::debug(const char *format, ...) const {
  const char t_normal[] = { 0x1b,'[','0',';','0',';','0','m','\0' };
  const char t_red[] = { 0x1b,'[','4',';','3','1',';','5','9','m','\0' };
  const char t_bold[] = { 0x1b,'[','1','m','\0' };
  if (is_debug) {
    va_list ap;
    va_start(ap,format);
    std::fprintf(cimg_stdout,"\n%s%s<gmic-debug-#%u>%s ",t_bold,t_red,filenames.size,t_normal);
    std::vfprintf(cimg_stdout,format,ap);
    va_end(ap);
    std::fflush(cimg_stdout);
  }
  return *this;
}

// Print status messages.
//-----------------------
const gmic& gmic::print(const char *format, ...) const {
  va_list ap;
  va_start(ap,format);
  if (verbosity_level>=0) {
    std::fprintf(cimg_stdout,"\n<gmic-#%u> ",filenames.size);
    std::vfprintf(cimg_stdout,format,ap);
    std::fflush(cimg_stdout);
  }
  va_end(ap);
  return *this;
}

// Add macros from a char* buffer.
//---------------------------------
gmic& gmic::add_macros(const char *const data_macros, const unsigned int data_size, const bool add_macros_at_start) {
  if (!data_macros || !data_size) return *this;
  char mac[4096] = { 0 }, com[256*1024] = { 0 }, line[256*1024] = { 0 }, sep = 0;
  const char *data = data_macros, *const data_end = data_macros + data_size;
  while (data<data_end) {
    if (*data=='\n') ++data;
    else {
      if (std::sscanf(data,"%262143[^\n]",line)>0) data += cimg::strlen(line) + 1;
      if (line[0]!='#') { // Useful line (not a comment)
        mac[0] = com[0] = 0;
        if (std::sscanf(line,"%4095[^: ] %c %262143[^\n]",mac,&sep,com)>=2 && sep==':' &&
            std::sscanf(mac,"%4095s",line)==1) { // Macro definition.
          macros.insert(CImg<char>(line,cimg::strlen(line)+1,1,1,1,false),add_macros_at_start?0:macros.size);
          commands.insert(CImg<char>(com,cimg::strlen(com)+1,1,1,1,false),add_macros_at_start?0:commands.size);
        } else { // Possible continuation of a previous macro definition.
          if (!macros) error("Fatal error : Invalid G'MIC macros data.");
          CImg<char> &last = commands[add_macros_at_start?0:commands.size-1];
          last[last.size()-1] = ' ';
          last.append(CImg<char>(line,cimg::strlen(line)+1,1,1,1,false),'x');
        }
      }
    }
  }
  return *this;
}

// Add macros from a macro file.
//------------------------------
gmic& gmic::add_macros(std::FILE *const file, const bool add_macros_at_start) {
  if (!file) return *this;
  char mac[4096] = { 0 }, com[256*1024] = { 0 }, line[256*1024] = { 0 }, sep = 0;
  int err = 0;
  while ((err=std::fscanf(file,"%262143[^\n] ",line)>=0)) {
    if (err) { // Non empty-line
      mac[0] = com[0] = 0;
      if (line[0]!='#') { // Useful line (not a comment).
        if (std::sscanf(line,"%4095[^: ] %c %262143[^\n]",mac,&sep,com)>=2 && sep==':' &&
            std::sscanf(mac,"%4095s",line)==1) { // Macro definition.
          macros.insert(CImg<char>(line,cimg::strlen(line)+1,1,1,1,false),add_macros_at_start?0:macros.size);
          commands.insert(CImg<char>(com,cimg::strlen(com)+1,1,1,1,false),add_macros_at_start?0:commands.size);
        } else { // Possible continuation of a previous macro definition.
          if (!macros) error("Fatal error : Invalid G'MIC macros data.");
          CImg<char> &last = commands[add_macros_at_start?0:commands.size-1];
          last[last.size()-1] = ' ';
          last.append(CImg<char>(line,cimg::strlen(line)+1,1,1,1,false),'x');
        }
      }
    }
  }
  return *this;
}

// Return indices of the images from a string.
//--------------------------------------------
CImg<unsigned int> gmic::indices2cimg(const char *const string, const unsigned int indice_max,
                                      const char *const command) const {
  if (!cimg::strlen(string)) return CImg<unsigned int>();
  CImgList<unsigned int> inds;
  const char *it = string;
  for (bool stopflag = false; !stopflag; ) {
    char sep = 0, end = 0, item0[4096] = { 0 }, item1[4096] = { 0 };
    float ind0 = 0, ind1 = 0, step = 1;
    if (std::sscanf(it,"%4095[^,]%c",item0,&end)!=2) stopflag = true;
    else it += 1 + cimg::strlen(item0);
    const int err = std::sscanf(item0,"%4095[^:]%c%f%c",item1,&sep,&step,&end);
    if (err!=1 && err!=3) error("Command '%s' : Invalid indice(s) '[%s]'.",command,string);
    if (std::sscanf(item1,"%f%%-%f%c%c",&ind0,&ind1,&sep,&end)==3 && sep=='%') {
      ind0 = (float)cimg::round(ind0*indice_max/100,1);
      ind1 = (float)cimg::round(ind1*indice_max/100,1);
    } else if (std::sscanf(item1,"%f%%-%f%c",&ind0,&ind1,&end)==2)
      ind0 = (float)cimg::round(ind0*indice_max/100,1);
    else if (std::sscanf(item1,"%f-%f%c%c",&ind0,&ind1,&sep,&end)==3 && sep=='%')
      ind1 = (float)cimg::round(ind1*indice_max/100,1);
    else if (std::sscanf(item1,"%f-%f%c",&ind0,&ind1,&end)==2) { }
    else if (std::sscanf(item1,"%f%c%c",&ind0,&sep,&end)==2 && sep=='%')
      ind1 = (ind0 = (float)cimg::round(ind0*indice_max/100,1));
    else if (std::sscanf(item1,"%f%c",&ind0,&end)==1)
      ind1 = ind0;
    else error("Command '%s' : Invalid indice(s) '[%s]'.",command,string);
    if (ind0<0) ind0+=indice_max;
    if (ind1<0) ind1+=indice_max;
    if (ind0<0 || ind0>=indice_max || ind1<0 || ind1>=indice_max || step<=0) {
      if (indice_max) error("Command '%s' : Invalid indice(s) '[%s]' (valid indice range is -%u...%u).",
                            command,string,indice_max,indice_max-1);
      else error("Command '%s' : Invalid indice(s) '[%s]' (image list is empty).",
                 command,string);
    }
    if (ind0>ind1) cimg::swap(ind0,ind1);
    const unsigned int
      iind0 = (unsigned int)ind0,
      _ind1 = (unsigned int)ind1,
      iind1 = (unsigned int)(_ind1 - cimg::mod((float)_ind1,step));
    if (iind0==iind1) inds.insert(CImg<unsigned int>::vector(iind0));
    else inds.insert(CImg<unsigned int>::sequence((unsigned int)(1+(iind1-iind0)/step),
                                                  (unsigned int)iind0,
                                                  (unsigned int)iind1).get_split('y'));
  }
  inds = inds.get_append('y').sort().get_split('y');
  cimglist_for(inds,l) if (l!=inds.size-1 && inds(l,0)==inds(l+1,0)) inds.remove(l--);
  if (is_debug) {
    debug("Indices : ");
    inds.get_append('y').print(); // List indices if debug mode is activated.
  }
  return inds.get_append('y').sort();
}

// Return stringified version of indices or filenames.
//----------------------------------------------------
char* gmic::indices2string(const CImg<unsigned int>& indices, const bool display_indices) const {
  static char res0[4096] = { 0 }, res1[4096] = { 0 };
  const unsigned int siz = indices.size();
  if (display_indices) {
    switch (siz) {
    case 0: std::sprintf(res0," []"); break;
    case 1: std::sprintf(res0," [%u]",indices[0]); break;
    case 2: std::sprintf(res0,"s [%u,%u]",indices[0],indices[1]); break;
    case 3: std::sprintf(res0,"s [%u,%u,%u]",indices[0],indices[1],indices[2]); break;
    case 4: std::sprintf(res0,"s [%u,%u,%u,%u]",indices[0],indices[1],indices[2],indices[3]); break;
    default: std::sprintf(res0,"s [%u,...,%u]",indices[0],indices[siz-1]);
    }
    return res0;
  }
  switch (siz) {
  case 0: std::sprintf(res1," "); break;
  case 1: std::sprintf(res1,"%s",filenames[indices[0]].ptr()); break;
  case 2: std::sprintf(res1,"%s, %s",filenames[indices[0]].ptr(),filenames[indices[1]].ptr()); break;
  case 3: std::sprintf(res1,"%s, %s, %s",filenames[indices[0]].ptr(),filenames[indices[1]].ptr(),
                       filenames[indices[2]].ptr()); break;
  case 4: std::sprintf(res1,"%s, %s, %s, %s",filenames[indices[0]].ptr(),filenames[indices[1]].ptr(),
                       filenames[indices[2]].ptr(), filenames[indices[3]].ptr()); break;
  default: std::sprintf(res1,"%s, ..., %s",filenames[indices[0]].ptr(),filenames[indices[siz-1]].ptr());
  }
  return res1;
}
#endif // #if defined(gmic_float) || !defined(gmic_separate_compilation)

// Template constructors.
//-----------------------
template<typename T>
gmic::gmic(const int argc, const char *const *const argv, CImgList<T>& images, const char *custom_macros, const bool add_macros_at_start) {
  assign(images.size,custom_macros,add_macros_at_start);
  for (int pos = 1; pos<argc; ++pos)
    command_line.insert(CImg<char>(argv[pos],cimg::strlen(argv[pos])+1,1,1,1,false));
  is_released = false;
  parse(images);
}

template<typename T>
gmic::gmic(const char *const command, CImgList<T>& images, const char *custom_macros, const bool add_macros_at_start) {
  assign(images.size,custom_macros,add_macros_at_start);
  char item[4096] = { 0 };
  for (const char *ncommand = command; *ncommand; ) {
    if (std::sscanf(ncommand,"%[^ ]",item)==1) {
      const int l = cimg::strlen(item);
      command_line.insert(CImg<char>(item,l+1,1,1,1,false));
      ncommand += l;
      while (*ncommand==' ') ++ncommand;
    } else break;
  }
  is_released = true;
  parse(images);
}

// Display specified image(s).
//-----------------------------
template<typename T>
bool gmic::display_images(const CImgList<T>& images, const CImg<unsigned int>& indices,
                          const bool verbose) const {
  if (!images || !indices) { print("Display image []."); return false; }
  CImgList<unsigned int> inds = indices.get_unroll('x').get_split('x');
  CImgList<T> visu;
  unsigned int max_height = 0;
  cimglist_for(inds,l) {
    const CImg<T>& img = images[inds(l,0)];
    if (img.height>max_height && !img.is_CImg3d()) max_height = img.height;
  }
  cimglist_for(inds,l) {
    const unsigned int ind = inds(l,0);
    const CImg<T> &img = images[ind];
    if (img) {
      if (!max_height || img.height<max_height) visu.insert(img,~0U,true);
      else visu.insert(img.get_lines(0,max_height-1));
    } else if (verbose) { warning("Display image : Image [%d] is empty.",ind); inds.remove(l--); }
  }
  const CImg<unsigned int> nindices = inds.get_append('x');
  const char *const fnames = indices2string(nindices,false);
  print("Display image%s = '%s'.\n\n",gmic_inds,fnames);
  if (visu.size) {
    if (visu.size!=1) visu.display(fnames,verbosity_level>=0,'x','p');
    else {
      const CImg<T> &img = visu[0];
      char title[4096] = { 0 };
      std::sprintf(title,"%s (%dx%dx%dx%d)",fnames,img.dimx(),img.dimy(),img.dimz(),img.dimv());
      img.display(title,verbosity_level>=0);
    }
  }
  return true;
}

// Display plots of specified image(s).
//--------------------------------------
template<typename T>
bool gmic::display_plots(const CImgList<T>& images, const CImg<unsigned int>& indices,
                         const unsigned int plot_type, const unsigned int vertex_type,
                         const double xmin, const double xmax,
                         const double ymin, const double ymax,
                         const bool verbose) const {
  if (!images || !indices) { print("Plot image []."); return false; }
  CImgDisplay disp(cimg_fitscreen(640,480,1),0,0);
  cimg_foroff(indices,l) {
    const unsigned int ind = indices[l];
    const CImg<T>& img = images[ind];
    if (img) {
      print("Plot image%s = '%s'.\n",gmic_inds,indices2string(indices,false));
      if (verbosity_level>=0) { std::fputc('\n',cimg_stdout); img.print(filenames[ind].ptr()); }
      char title[4096] = { 0 };
      std::sprintf(title,"%s (%dx%dx%dx%d)",
                   filenames[ind].ptr(),img.dimx(),img.dimy(),img.dimz(),img.dimv());
      img.display_graph(disp.set_title(title),plot_type,vertex_type,0,xmin,xmax,0,ymin,ymax);
    } else if (verbose) warning("Plot image : Image [%d] is empty.",ind);
  }
  return true;
}

// Display specified 3D object(s).
//--------------------------------
template<typename T>
bool gmic::display_objects3d(const CImgList<T>& images, const CImg<unsigned int>& indices,
                             const bool verbose) const {
  if (!indices) { print("Display 3D object []."); return false; }
  CImg<unsigned char> background;
  bool exist3d = false;
  CImgDisplay disp;
  cimg_foroff(indices,l) {
    const unsigned int ind = indices[l];
    const CImg<T> &img = images[ind];
    if (!img.is_CImg3d()) {
      if (verbose) warning("Display 3D object : Image [%d] is not a 3D object.",ind);
    } else {
      exist3d = true;
      if (!background || !disp) {
        background.assign(cimg_fitscreen(640,480,1),1,3);
        cimg_forV(background,k) background.get_shared_channel(k).fill(background3d[k]);
        disp.assign(background);
      }
      CImgList<unsigned int> primitives3d;
      CImgList<unsigned char> colors3d;
      CImg<float> opacities3d;
      CImg<float> points3d(img);
      points3d.CImg3dtoobject3d(primitives3d,colors3d,opacities3d);
      print("Display 3D object [%u] = '%s' (%d points, %u primitives).",
            ind,filenames[ind].ptr(),points3d.dimx(),primitives3d.size);
      disp.set_title("%s (%d points, %u primitives)",
                     filenames[ind].ptr(),points3d.dimx(),primitives3d.size);
      background.display_object3d(disp,points3d,primitives3d,colors3d,opacities3d,
                                  true,render3d,renderd3d,!is_oriented3d,focale3d,specular_light3d,specular_shine3d);
      if (disp.is_closed) break;
    }
  }
  return exist3d;
}

// Substitute '@' expressions.
//----------------------------
template<typename T>
CImg<char> gmic::substitute_arobace(const char *const argument, const CImgList<T>& images) const {
  if (!argument) return CImg<char>();
  CImgList<char> _largument;
  char range[4096] = { 0 };
  for (const char *nargument = argument; *nargument; ) {
    if (*nargument=='@') {
      char argx[4096] = { 0 }, argy[4096] = { 0 }, argz[4096] = { 0 }, argv[4096] = { 0 };
      int ind = 0, bcond = 0; *range = 0; char sepx = 0, sepy = 0, sepz = 0, sepv = 0, sep = 0, end = 0;
      float x = 0, y = 0, z = 0, v = 0, m = 0, M = 1;
      if (nargument[1]=='#' ||
          (std::sscanf(nargument,"@{#%c",&sep)==1 && sep=='}')) {
        std::sprintf(range,"%u",images.size);
        _largument.insert(CImg<char>(range,cimg::strlen(range),1,1,1,true));
        if (sep=='}') nargument+=4; else nargument+=2;
      } else if (std::sscanf(nargument,"@%d",&ind)==1 ||
                 (std::sscanf(nargument,"@{%d%c",&ind,&sep)==2 && sep=='}') ||
                 (std::sscanf(nargument,"@{%d%*c%4095[^}]%c",&ind,range,&sep)==3 && sep=='}')) {
        int nind = ind;
        if (nind<0) nind+=images.size;
        if (nind<0 || nind>=(int)images.size) {
          if (images.size) error("Invalid indice '%d' in '@argument' (valid indice range is -%u...%u).",
                                 ind,images.size,images.size-1);
          else error("Invalid indice '%d' in '@argument' (image list is empty).",ind);
        }
        const unsigned int sizrange = cimg::strlen(range);
        const CImg<T>& img = images[nind];
        CImg<T> values;
        if (sizrange) {
          const CImg<unsigned int> iinds = indices2cimg(range,img.size(),"parsing");
          values.assign(iinds.size());
          cimg_foroff(iinds,p) values[p] = img[iinds(p)];
        } else values = img.get_shared();
        const CImg<char> vs = values.value_string();
        const unsigned int vsl = vs.size();
        if (vsl>1) _largument.insert(CImg<char>(vs.ptr(),vsl-1,1,1,1,true));
        nargument+= 1 + (sep?2:0) + (sizrange?1:0) + sizrange + std::sprintf(range,"%d",ind);
      } else if (((std::sscanf(nargument,"@(%d%*c%4095[0-9.eE%+-]%c",&ind,argx,&sep)==3 && sep==')') ||
                  (std::sscanf(nargument,"@(%d%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",&ind,argx,argy,&sep)==4 && sep==')') ||
                  (std::sscanf(nargument,"@(%d%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",
                               &ind,argx,argy,argz,&sep)==5 && sep==')') ||
                  (std::sscanf(nargument,"@(%d%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",
                               &ind,argx,argy,argz,argv,&sep)==6 && sep==')') ||
                  (std::sscanf(nargument,"@(%d%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%d%c",
                               &ind,argx,argy,argz,argv,&bcond,&sep)==7 && sep==')')) &&
                 (!*argx ||
                  (std::sscanf(argx,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%') ||
                  std::sscanf(argx,"%f%c",&x,&end)==1) &&
                 (!*argy ||
                  (std::sscanf(argy,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%') ||
                  std::sscanf(argy,"%f%c",&y,&end)==1) &&
                 (!*argz ||
                  (std::sscanf(argz,"%f%c%c",&z,&sepz,&end)==2 && sepz=='%') ||
                  std::sscanf(argz,"%f%c",&z,&end)==1) &&
                 (!*argv ||
                  (std::sscanf(argv,"%f%c%c",&v,&sepv,&end)==2 && sepv=='%') ||
                  std::sscanf(argv,"%f%c",&v,&end)==1)) {
        int nind = ind;
        if (nind<0) nind+=images.size;
        if (nind<0 || nind>=(int)images.size) {
          if (images.size) error("Invalid indice '%d' in '@argument' (valid indice range is -%u...%u).",
                                 ind,images.size,images.size-1);
          else error("Invalid indice '%d' in '@argument' (image list is empty).",ind);
        }
        const CImg<T>& img = images[nind];
        const int
          nx = (int)cimg::round(sepx=='%'?x*(img.dimx()-1)/100:x,1),
          ny = (int)cimg::round(sepy=='%'?y*(img.dimy()-1)/100:y,1),
          nz = (int)cimg::round(sepz=='%'?z*(img.dimz()-1)/100:z,1),
          nv = (int)cimg::round(sepv=='%'?v*(img.dimv()-1)/100:v,1);
        std::sprintf(range,"%g",bcond?(double)img.atXYZV(nx,ny,nz,nv):(double)img.atXYZV(nx,ny,nz,nv,0));
        _largument.insert(CImg<char>(range,cimg::strlen(range)));
        while (*nargument!=')') ++nargument; ++nargument;
      } else if ((std::sscanf(nargument,"@%c",&sep)==1 && sep=='?') ||
                 (std::sscanf(nargument,"@{?%c",&sep)==1 && sep=='}') ||
                 (std::sscanf(nargument,"@[?%c",&sep)==1 && sep==']') ||
                 (std::sscanf(nargument,"@{?%*c%f%c",&M,&sep)==2 && sep=='}') ||
                 (std::sscanf(nargument,"@[?%*c%f%c",&M,&sep)==2 && sep==']') ||
                 (std::sscanf(nargument,"@{?%*c%f%*c%f%c",&m,&M,&sep)==3 && sep=='}') ||
                 (std::sscanf(nargument,"@[?%*c%f%*c%f%c",&m,&M,&sep)==3 && sep==']')) {
        const double r = m + (M-m)*cimg::rand();
        if (sep==']') std::sprintf(range,"%d",(int)cimg::round(r,1)); else std::sprintf(range,"%g",r);
        _largument.insert(CImg<char>(range,cimg::strlen(range),1,1,1,true));
        if (sep=='?') nargument+=2; else { while (*nargument!=sep) ++nargument; ++nargument; }
      } else _largument.insert(CImg<char>(nargument++,1,1,1,1,true));
    } else {
      std::sscanf(nargument,"%4095[^@]",range);
      const int l = cimg::strlen(range);
      _largument.insert(CImg<char>(range,l,1,1,1,true));
      nargument+=l;
    }
  }
  _largument.insert(CImg<char>::vector(0));
  return _largument.get_append('x');
}

// Main parsing procedure.
//------------------------
template<typename T>
gmic& gmic::parse(CImgList<T> &images) {
  const unsigned int command_line_maxsize = 65535;
  const int no_ind = (int)(~0U>>1);
  cimg::exception_mode() = 0;

  // Begin command line parsing.
  while (position<command_line.size && command_line.size<command_line_maxsize) {
    const char
      *const orig_item = command_line[position].ptr(),
      *const orig_argument = position+1<command_line.size?command_line[position+1].ptr():"";

    // Substitute '@' expressions in 'orig_item' and 'orig_argument' if necessary.
    CImg<char> _item, _argument, _argument_text;
    if (*orig_item=='-' || *orig_item=='[' || *orig_item=='(') {
      if (std::strchr(orig_item,'@')) {
        _item = substitute_arobace(orig_item,images);
        if (_item.size()>4095) error("Buffer overflow when substituting item '%s'.",orig_item);
      }
      if (*orig_item=='-' &&
          (*orig_argument!='-' ||
           (*orig_argument=='-' && (orig_argument[1]=='.' || orig_argument[1]=='@' ||
                                    (orig_argument[1]>='0' && orig_argument[1]<='9'))))
          && std::strchr(orig_argument,'@')) {
        _argument = substitute_arobace(orig_argument,images);
        if (_argument.size()>4095) error("Buffer overflow when substituting argument '%s'.",orig_argument);
      }
    }
    const char
      *item = _item?_item.ptr():orig_item,
      *argument = _argument?_argument.ptr():orig_argument;
    const char *argument_text = argument;
    if (cimg::strlen(argument)>=64) {
      _argument_text.assign(64,1,1,1,0);
      std::memcpy(_argument_text.ptr(),argument,60*sizeof(char));
      _argument_text[60] = _argument_text[61] = _argument_text[62] = '.'; _argument_text[63] = 0;
      argument_text = _argument_text.ptr();
    }

    // Get current item/command from the command line.
    char item0[4096] = { 0 }, item1[4096] = { 0 };
    bool get_version = false;
    CImg<unsigned int> indices;
    if (item[0]=='-' && item[1] && item[1]!='.') {
      char sep0 = 0, sep1 = 0, end = 0;
      if (item[1]=='-' && item[2] && item[2]!='[' && item[2]!='3' && item[3]!='d') { ++item; get_version = true; }
      const int err = std::sscanf(item,"%4095[^[]%c%4095[0-9.eE%,:+-]%c%c",item0,&sep0,item1,&sep1,&end);
      if (err==1) indices = CImg<unsigned int>::sequence(images.size,0,images.size-1);
      else if (err==4 && sep1==']')
        indices = indices2cimg(item1,(!strcmp(item0,"-i")||!strcmp(item0,"-input"))?images.size+1:images.size,item0);
      else { std::strcpy(item0,item); item1[0] = 0; }
    }
    ++position;

    // Check for verbosity commands.
    if (*item=='-') {
      if (!cimg::strcmp("-verbose+",item) || !cimg::strcmp("-v+",item)) ++verbosity_level;
      else if (!cimg::strcmp("-verbose-",item) || !cimg::strcmp("-v-",item)) --verbosity_level;
    }

    if (is_begin) { print("Start G'MIC instance."); is_begin = false; }
    debug("Item : '%s', Argument : '%s'.",item,argument);

    // Begin command interpretation.
    try {
      if (*item=='-') {

        //----------------
        // Global options
        //----------------

        // Verbosity (actually, just continue to next command since verbosity has been already processed above).
        if (!cimg::strcmp("-verbose+",item) || !cimg::strcmp("-v+",item) ||
            !cimg::strcmp("-verbose-",item) || !cimg::strcmp("-v-",item)) continue;

        // Load macro file.
        if (!cimg::strcmp("-macros",item) || !cimg::strcmp("-m",item)) {
          print("Load macro file '%s'",cimg::basename(argument));
          std::FILE *const file = cimg::fopen(argument,"r");
          if (file) {
            const unsigned int siz = macros.size;
            add_macros(file,true);
            cimg::fclose(file);
            if (verbosity_level>=0) std::fprintf(cimg_stdout," (%u macros added).",macros.size-siz);
          }
          else error("Load macro file '%s' : File not found.",argument);
          ++position; continue;
        }

        // Switch 'is_debug' flag.
        if (!cimg::strcmp("-debug",item)) {
          is_debug = !is_debug;
          continue;
        }

        // Switch 'is_fullpath' flag.
        if (!cimg::strcmp("-fullpath",item)) {
          is_fullpath = !is_fullpath;
          continue;
        }

        //----------------------
        // Arithmetic operators
        //----------------------

        // Common arithmetic operators.
        gmic_arithmetic_item("-add","-+",
                             operator+=,"Add %g to image%s",value,gmic_inds,T,
                             operator+=,"Add to image%s",
                             "Add image [%d] to image%s",ind[0],gmic_inds,
                             "Add image%s together");

        gmic_arithmetic_item("-sub","--",
                             operator-=,"Substract %g to image%s",value,gmic_inds,T,
                             operator-=,"Substract to image%s",
                             "Substract image [%d] to image%s",ind[0],gmic_inds,
                             "Substract image%s together");

        gmic_arithmetic_item("-mul","-*",
                             operator*=,"Multiply image%s by %g",gmic_inds,value,double,
                             mul,"Multiply image%s",
                             "Multiply image%s by image [%d]",gmic_inds,ind[0],
                             "Multiply image%s together");

        gmic_arithmetic_item("-div","-/",
                             operator/=,"Divide image%s by %g",gmic_inds,value,double,
                             div,"Divide image%s",
                             "Divide image%s by image [%d]",gmic_inds,ind[0],
                             "Divide image%s together");

        gmic_arithmetic_item("-pow","-^",
                             pow,"Compute image%s to the power of %g",gmic_inds,value,double,
                             pow,"Compute power of image%s",
                             "Compute image%s to the power of image [%d]",gmic_inds,ind[0],
                             "Compute the power of image%s together");

        gmic_arithmetic_item("-min","-min",
                             min,"Compute pointwise minimum between image%s and %g",gmic_inds,value,T,
                             min,"Compute pointwise minimum with image%s",
                             "Compute pointwise minimum between image%s and image [%d]",gmic_inds,ind[0],
                             "Compute pointwise minimum between image%s together");

        gmic_arithmetic_item("-max","-max",
                             max,"Compute pointwise maximum between image%s and %g",gmic_inds,value,T,
                             max,"Compute pointwise maximum with image%s",
                             "Compute pointwise maximum between image%s and image [%d]",gmic_inds,ind[0],
                             "Compute pointwise maximum between image%s together");

        gmic_arithmetic_item("-mod","-%",
                             operator%=,"Compute pointwise modulo between image%s and %g.",gmic_inds,value,T,
                             operator%=,"Compute pointwise modulo with image%s",
                             "Compute pointwise modulo between image%s and image [%d]",gmic_inds,ind[0],
                             "Compute pointwise modulo between image%s together");

        gmic_arithmetic_item("-and","-and",
                             operator&=,"Compute bitwise AND between image%s and %g.",gmic_inds,value,unsigned int,
                             operator&=,"Compute bitwise AND with image%s",
                             "Compute bitwise AND between image%s and image [%d]",gmic_inds,ind[0],
                             "Compute bitwise AND between image%s together");

        gmic_arithmetic_item("-or","-or",
                             operator|=,"Compute bitwise OR between image%s and %g.",gmic_inds,value,unsigned int,
                             operator|=,"Compute bitwise OR with image%s",
                             "Compute bitwise OR between image%s and image [%d]",gmic_inds,ind[0],
                             "Compute bitwise OR between image%s together");

        gmic_arithmetic_item("-xor","-xor",
                             operator^=,"Compute bitwise XOR between image%s and %g.",gmic_inds,value,unsigned int,
                             operator^=,"Compute bitwise XOR with image%s",
                             "Compute bitwise XOR between image%s and image [%d]",gmic_inds,ind[0],
                             "Compute bitwise XOR between image%s together");

        // Other arithmetic operators.
        gmic_simple_item("-cos",cos,"Compute cosine of image%s.");
        gmic_simple_item("-sin",sin,"Compute sine of image%s.");
        gmic_simple_item("-tan",tan,"Compute tangent of image%s.");
        gmic_simple_item("-acos",acos,"Compute arccosine of image%s.");
        gmic_simple_item("-asin",asin,"Compute arcsine of image%s.");
        gmic_simple_item("-atan",atan,"Compute arctangent of image%s.");
        gmic_simple_item("-abs",abs,"Compute absolute value of image%s.");
        gmic_simple_item("-sqr",sqr,"Compute square function of image%s.");
        gmic_simple_item("-sqrt",sqrt,"Compute square root of image%s.");
        gmic_simple_item("-exp",exp,"Compute exponential of image%s.");
        gmic_simple_item("-log",log,"Compute logarithm of image%s.");
        gmic_simple_item("-log10",log10,"Compute logarithm_10 of image%s.");

        //---------------------------------------
        // Pointwise operations on pixel values
        //---------------------------------------

        // Type cast.
        if (!cimg::strcmp("-type",item) || !cimg::strcmp("-t",item)) {
          typedef unsigned char uchar;
          typedef unsigned short ushort;
          typedef unsigned int uint;
#ifndef gmic_minimal
          gmic_cast(bool,"bool");
          gmic_cast(uchar,"unsigned char");
          gmic_cast(char,"char");
          gmic_cast(ushort,"unsigned short");
          gmic_cast(short,"short");
          gmic_cast(uint,"unsigned int");
          gmic_cast(int,"int");
          gmic_cast(double,"double");
#endif
          gmic_cast(float,"float");
          error("Set pixel type : Invalid argument '%s' "
                "(should be '{bool,uchar,char,ushort,short,uint,int,float,double}').",
                argument_text);
        }

        // Set scalar value.
        if (!cimg::strcmp("-set",item0) || !cimg::strcmp("-=",item0)) {
          double value = 0; float x = 0, y = 0, z = 0, v = 0; char sepx = 0, sepy = 0, sepz = 0, sepv = 0, end = 0;
          char argx[4096] = { 0 }, argy[4096] = { 0 }, argz[4096] = { 0 }, argv[4096] = { 0 };
          if ((std::sscanf(argument,"%lf%c",&value,&end)==1 ||
               std::sscanf(argument,"%lf%*c%4095[0-9.eE%+-]%c",&value,argx,&end)==2 ||
               std::sscanf(argument,"%lf%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",&value,argx,argy,&end)==3 ||
               std::sscanf(argument,"%lf%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",&value,argx,argy,argz,&end)==4 ||
               std::sscanf(argument,"%lf%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",&value,argx,argy,argz,argv,&end)==5) &&
              (!*argx ||
               (std::sscanf(argx,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%') ||
               std::sscanf(argx,"%f%c",&x,&end)==1) &&
              (!*argy ||
               (std::sscanf(argy,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%') ||
               std::sscanf(argy,"%f%c",&y,&end)==1) &&
              (!*argz ||
               (std::sscanf(argz,"%f%c%c",&z,&sepz,&end)==2 && sepz=='%') ||
               std::sscanf(argz,"%f%c",&z,&end)==1) &&
              (!*argv ||
               (std::sscanf(argv,"%f%c%c",&v,&sepv,&end)==2 && sepv=='%') ||
               std::sscanf(argv,"%f%c",&v,&end)==1)) {
            print("Set scalar value %g at position (%g%s,%g%s,%g%s,%g%s) in image%s",
                  value,x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",z,sepz=='%'?"%":"",v,sepv=='%'?"%":"",gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.dimx()-1)/100:x,1),
                ny = (int)cimg::round(sepy=='%'?y*(img.dimy()-1)/100:y,1),
                nz = (int)cimg::round(sepz=='%'?z*(img.dimz()-1)/100:z,1),
                nv = (int)cimg::round(sepv=='%'?v*(img.dimv()-1)/100:v,1);
              gmic_apply(images[indices[l]],gmic_set(value,nx,ny,nz,nv));
            }
          } else error("Set scalar value in image%s : Invalid argument '%s' "
                       "(should be 'value,x[,y[,z[,v]]]').",
                       argument_text);
          ++position; continue;
        }

        // Invert endianness.
        gmic_simple_item("-endian",invert_endianness,"Invert endianness of image%s.");

        // Fill.
        if (!cimg::strcmp("-fill",item0) || !cimg::strcmp("-f",item0)) {
          char sep = 0, end = 0; int ind0 = no_ind;
          if (std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']') {
            gmic_check_indice(ind0,"Fill image%s");
            print("Fill image%s with values of image [%d].",gmic_inds,ind0);
            const CImg<T> values = images[ind0];
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],fill(values));
          } else {
            print("Fill image%s with values '%s'.",gmic_inds,argument_text);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],fill(argument,true));
          }
          ++position; continue;
        }

        // Threshold.
        if (!cimg::strcmp("-threshold",item0)) {
          char sep = 0, end = 0; int soft = 0; double value = 0;
          if (std::sscanf(argument,"%lf%c",&value,&end)==1 ||
              (std::sscanf(argument,"%lf%c%c",&value,&sep,&end)==2 && sep=='%') ||
              std::sscanf(argument,"%lf%*c%d%c",&value,&soft,&end)==2 ||
              (std::sscanf(argument,"%lf%c%*c%d%c",&value,&sep,&soft,&end)==3 && sep=='%')) {
            print("Threshold image%s with %g%s (%s threshold).",gmic_inds,value,sep=='%'?"%":"",soft?"soft":"hard");
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              double vmin = 0, vmax = 0, nvalue = value;
              if (sep=='%') { vmin = img.minmax(vmax); nvalue = vmin + (vmax - vmin)*value/100; }
              gmic_apply(img,threshold((T)nvalue,soft?true:false));
            }
            ++position;
          } else {
            print("Threshold image%s : Interactive mode.",gmic_inds);
            CImgDisplay disp;
            char title[4096] = { 0 };
            cimg_foroff(indices,l) {
              CImg<T>
                &img = images[indices[l]],
                visu = img.depth>1?img.get_projections2d(img.dimx()/2,img.dimy()/2,img.dimz()/2).
                channels(0,cimg::min(3,img.dimv())-1):img.get_channels(0,cimg::min(3,img.dimv()-1));
              if (disp) disp.resize(cimg_fitscreen(visu.dimx(),visu.dimy(),1),false);
              else disp.assign(cimg_fitscreen(visu.dimx(),visu.dimy(),1),0,1);
              double
                vmin = 0, vmax = (double)img.maxmin(vmin),
                distmax = std::sqrt(cimg::sqr(disp.dimx()-1.0) + cimg::sqr(disp.dimy()-1.0)),
                amount = 50;
              bool stopflag = false, obutt = false;
              int omx = -1, omy = -1;
              CImg<T> res;
              for (disp.show().button = disp.key = 0; !stopflag; ) {
                const unsigned int key = disp.key;
                if (!res) {
                  std::sprintf(title,"%s : Interactive threshold %.3g%%",filenames[indices[l]].ptr(),amount);
                  disp.display(res=visu.get_threshold((T)(vmin + amount*(vmax-vmin)/100))).
                    set_title(title).wait();
                }
                const int mx = disp.mouse_x, my = disp.mouse_y;
                if (disp.button && mx>=0 && my>=0) {
                  if (omx==mx && omy==my && !obutt) break;
                  omx = mx; omy = my; obutt = true;
                  const double dist = std::sqrt((double)cimg::sqr(mx) + cimg::sqr(my));
                  amount = dist*100/distmax;
                  res.assign();
                } else if (!disp.button) obutt = false;
                if (disp.is_closed || (key && key!=cimg::keyCTRLLEFT)) stopflag = true;
                if (key==cimg::keyD && disp.is_keyCTRLLEFT &&
                    (disp.resize(cimg_fitscreen(3*disp.width/2,3*disp.height/2,1),stopflag=false).key=0)==0)
                  disp.is_resized = true;
                if (key==cimg::keyC && disp.is_keyCTRLLEFT &&
                    (disp.resize(cimg_fitscreen(2*disp.width/3,2*disp.height/3,1),stopflag=false).key=0)==0)
                  disp.is_resized = true;
                if (disp.is_resized) {
                  disp.resize(false).display(res);
                  distmax = std::sqrt(cimg::sqr(disp.dimx()-1.0) + cimg::sqr(disp.dimy()-1.0));
                }
              }
              gmic_apply(img,threshold((T)(vmin + amount*(vmax-vmin)/100)));
            }
          }
          continue;
        }

        // Cut.
        if (!cimg::strcmp("-cut",item0)) {
          char sep0 = 0, sep1 = 0, end = 0, arg0[4096] = { 0 }, arg1[4096] = { 0 };
          double value0 = 0, value1 = 0; int ind0 = no_ind, ind1 = no_ind;
          if (std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%c",arg0,arg1,&end)==2 &&
              ((std::sscanf(arg0,"[%d%c%c",&ind0,&sep0,&end)==2 && sep0==']') ||
               (std::sscanf(arg0,"%lf%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
               std::sscanf(arg0,"%lf%c",&value0,&end)==1) &&
              ((std::sscanf(arg1,"[%d%c%c",&ind1,&sep1,&end)==2 && sep1==']') ||
               (std::sscanf(arg1,"%lf%c%c",&value1,&sep1,&end)==2 && sep1=='%') ||
               std::sscanf(arg1,"%lf%c",&value1,&end)==1)) {
            if (ind0!=no_ind) { gmic_check_indice(ind0,"Cut image%s"); value0 = images[ind0].min(); sep0 = 0; }
            if (ind1!=no_ind) { gmic_check_indice(ind1,"Cut image%s"); value1 = images[ind1].max(); sep1 = 0; }
            print("Cut image%s in [%g%s,%g%s].",gmic_inds,value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"");
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              double vmin = 0, vmax = 0, nvalue0 = value0, nvalue1 = value1;
              if (sep0=='%') { vmin = img.minmax(vmax); nvalue0 = vmin + (vmax - vmin)*value0/100; }
              if (sep1=='%') { vmin = img.minmax(vmax); nvalue1 = vmin + (vmax - vmin)*value1/100; }
              gmic_apply(img,cut((T)nvalue0,(T)nvalue1));
            }
            ++position;
          } else if (std::sscanf(argument,"[%d%c%c",&(ind0=no_ind),&sep0,&end)==2) {
            if (ind0!=no_ind) gmic_check_indice(ind0,"Cut image%s");
            value0 = images[ind0].minmax(value1);
            print("Cut image%s in [%g,%g].",gmic_inds,value0,value1);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],cut((T)value0,(T)value1));
            ++position;
          } else {
            print("Cut image%s : Interactive mode.",gmic_inds);
            CImgDisplay disp;
            char title[4096] = { 0 };
            cimg_foroff(indices,l) {
              CImg<T>
                &img = images[indices[l]],
                visu = img.depth>1?img.get_projections2d(img.dimx()/2,img.dimy()/2,img.dimz()/2).
                channels(0,cimg::min(3,img.dimv())-1):img.get_channels(0,cimg::min(3,img.dimv()-1));
              if (disp) disp.resize(cimg_fitscreen(visu.dimx(),visu.dimy(),1),false);
              else disp.assign(cimg_fitscreen(visu.dimx(),visu.dimy(),1),0,1);
              double vmin = 0, vmax = (double)img.maxmin(vmin), amount0 = 0, amount1 = 100;
              bool stopflag = false, obutt = false;
              int omx = -1, omy = -1;
              CImg<T> res;
              for (disp.show().button = disp.key = 0; !stopflag; ) {
                const unsigned int key = disp.key;
                if (!res) {
                  std::sprintf(title,"%s : Interactive cut [%.3g%%,%.3g%%]",
                               filenames[indices[l]].ptr(),amount0,amount1);
                  disp.display(res = visu.get_cut((T)(vmin + amount0*(vmax-vmin)/100),
                                                  (T)(vmin + amount1*(vmax-vmin)/100))).
                    set_title(title).wait();
                }
                const int mx = disp.mouse_x, my = disp.mouse_y;
                if (disp.button && mx>=0 && my>=0) {
                  if (omx==mx && omy==my && !obutt) break;
                  omx = mx; omy = my; obutt = true;
                  amount0 = mx*100/disp.dimx(); amount1 = my*100/disp.dimy();
                  res.assign();
                } else if (!disp.button) obutt = false;
                if (disp.is_closed || (key && key!=cimg::keyCTRLLEFT)) stopflag = true;
                if (key==cimg::keyD && disp.is_keyCTRLLEFT &&
                    (disp.resize(cimg_fitscreen(3*disp.width/2,3*disp.height/2,1),stopflag=false).key=0)==0)
                  disp.is_resized = true;
                if (key==cimg::keyC && disp.is_keyCTRLLEFT &&
                    (disp.resize(cimg_fitscreen(2*disp.width/3,2*disp.height/3,1),stopflag=false).key=0)==0)
                  disp.is_resized = true;
                if (disp.is_resized) disp.resize(false).display(res);
              }
              gmic_apply(img,cut((T)(vmin + amount0*(vmax-vmin)/100),(T)(vmin + amount1*(vmax-vmin)/100)));
            }
          }
          continue;
        }

        // Normalize.
        if (!cimg::strcmp("-normalize",item0) || !cimg::strcmp("-n",item0)) {
          char sep0 = 0, sep1 = 0, end = 0, arg0[4096] = { 0 }, arg1[4096] = { 0 };
          double value0 = 0, value1 = 0; int ind0 = no_ind, ind1 = no_ind;
          if (std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%c",arg0,arg1,&end)==2 &&
              ((std::sscanf(arg0,"[%d%c%c",&ind0,&sep0,&end)==2 && sep0==']') ||
               (std::sscanf(arg0,"%lf%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
               std::sscanf(arg0,"%lf%c",&value0,&end)==1) &&
              ((std::sscanf(arg1,"[%d%c%c",&ind1,&sep1,&end)==2 && sep1==']') ||
               (std::sscanf(arg1,"%lf%c%c",&value1,&sep1,&end)==2 && sep1=='%') ||
               std::sscanf(arg1,"%lf%c",&value1,&end)==1)) {
            if (ind0!=no_ind) { gmic_check_indice(ind0,"Normalize image%s"); value0 = images[ind0].min(); sep0 = 0; }
            if (ind1!=no_ind) { gmic_check_indice(ind1,"Normalize image%s"); value1 = images[ind1].max(); sep1 = 0; }
            print("Normalize image%s in [%g%s,%g%s].",gmic_inds,value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"");
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              double vmin = 0, vmax = 0, nvalue0 = value0, nvalue1 = value1;
              if (sep0=='%') { vmin = img.minmax(vmax); nvalue0 = vmin + (vmax - vmin)*value0/100; }
              if (sep1=='%') { vmin = img.minmax(vmax); nvalue1 = vmin + (vmax - vmin)*value1/100; }
              gmic_apply(img,normalize((T)nvalue0,(T)nvalue1));
            }
          } else if (std::sscanf(argument,"[%d%c%c",&(ind0=no_ind),&sep0,&end)==2) {
            if (ind0!=no_ind) gmic_check_indice(ind0,"Normalize image%s");
            value0 = images[ind0].minmax(value1);
            print("Normalize image%s in [%g,%g].",gmic_inds,value0,value1);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],normalize((T)value0,(T)value1));
          } else error("Normalize image%s : Invalid argument '%s' "
                       "(should be '{value1[%%],[indice]},{value2[%%],[indice]}').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Round.
        if (!cimg::strcmp("-round",item0)) {
          char end = 0; double value = 0; int rtype = 0;
          if (std::sscanf(argument,"%lf%c",&value,&end)==1 ||
              std::sscanf(argument,"%lf%*c%d%c",&value,&rtype,&end)==2) {
            print("Round image%s with value %g (%s rounding).",
                  gmic_inds,value,rtype<0?"backward":rtype>0?"forward":"nearest");
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],round((float)value,rtype));
          } else error("Round image%s : Invalid argument '%s' "
                       "(should be 'round_value[,round_type]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Equalize histogram.
        if (!cimg::strcmp("-equalize",item0)) {
          float nb = 256; char sep = 0, end = 0;
          if (std::sscanf(argument,"%f%c",&nb,&end)==1 ||
              (std::sscanf(argument,"%f%c%c",&nb,&sep,&end)==2 && sep=='%')) {
            if (nb<=0) error("Equalize image%s : Invalid cluster number %g.",gmic_inds,nb);
            print("Equalize image%s with %g%s clusters.",gmic_inds,nb,sep=='%'?"%":"");
            cimg_foroff(indices,l) {
              CImg<T>& img = images[indices[l]];
              unsigned int N = (unsigned int)nb;
              if (sep=='%') { double m, M = img.maxmin(m); N = (unsigned int)cimg::round((M-m)*nb/100,1); }
              gmic_apply(img,equalize((int)nb));
            }
          } else error("Equalize image%s : Invalid argument '%s' "
                       "(should be 'nb_clusters[%%]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Quantize.
        if (!cimg::strcmp("-quantize",item0)) {
          int nb = 0; char end = 0;
          if (std::sscanf(argument,"%d%c",&nb,&end)==1) {
            print("Quantize image%s with %d levels.",gmic_inds,nb);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],quantize(nb));
          } else error("Quantize image%s : Invalid argument '%s' "
                       "(should be 'nb_levels').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Add noise.
        if (!cimg::strcmp("-noise",item0)) {
          float sigma = 0; char sep = 0, end = 0; int ntype = 0;
          if (std::sscanf(argument,"%f%c",&sigma,&end)==1 ||
              (std::sscanf(argument,"%f%c%c",&sigma,&sep,&end)==2 && sep=='%') ||
              std::sscanf(argument,"%f%*c%d%c",&sigma,&ntype,&end)==2 ||
              (std::sscanf(argument,"%f%c%*c%d%c",&sigma,&sep,&ntype,&end)==3 && sep=='%')) {
            const char *st_type = ntype==0?"gaussian":ntype==1?"uniform":ntype==2?"salt&pepper":"poisson";
            if (sep=='%') sigma = -sigma;
            print("Add %s noise with standard deviation %g%s to image%s.",
                  st_type,cimg::abs(sigma),sigma<0?"%":"",gmic_inds);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],noise(sigma,ntype));
          } else error("Add noise to image%s : Invalid argument '%s' "
                       "(should be 'std[%%][,noise_type]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Rand.
        if (!cimg::strcmp("-rand",item0)) {
          double value0 = 0, value1 = 0; char end = 0;
          if (std::sscanf(argument,"%lf%*c%lf%c",&value0,&value1,&end)==2) {
            print("Fill image%s with random values in [%g,%g].",gmic_inds,value0,value1);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],rand((T)value0,(T)value1));
          } else error("Fill image%s with random values : Invalid argument '%s' "
                       "(should be 'valmin,valmax').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Compute pointwise norms and orientations.
        gmic_simple_item("-norm",pointwise_norm,"Compute vector norm.");
        gmic_simple_item("-orientation",pointwise_orientation,"Compute vector orientation.");

        //------------------------
        // Color base conversions
        //------------------------
        gmic_simple_item("-rgb2hsv",RGBtoHSV,"Convert image%s from RGB to HSV colorbases.");
        gmic_simple_item("-rgb2hsl",RGBtoHSL,"Convert image%s from RGB to HSL colorbases.");
        gmic_simple_item("-rgb2hsi",RGBtoHSI,"Convert image%s from RGB to HSI colorbases.");
        gmic_simple_item("-rgb2yuv",RGBtoYUV,"Convert image%s from RGB to YUV colorbases.");
        gmic_simple_item("-rgb2ycbcr",RGBtoYCbCr,"Convert image%s from RGB to YCbCr colorbases.");
        gmic_simple_item("-rgb2xyz",RGBtoXYZ,"Convert image%s from RGB to XYZ colorbases.");
        gmic_simple_item("-rgb2lab",RGBtoLab,"Convert image%s from RGB to Lab colorbases.");
        gmic_simple_item("-rgb2cmy",RGBtoCMY,"Convert image%s from RGB to CMY colorbases.");
        gmic_simple_item("-rgb2cmyk",RGBtoCMYK,"Convert image%s from RGB to CMYK colorbases.");
        gmic_simple_item("-cmyk2rgb",CMYKtoRGB,"Convert image%s from CMYK to RGB colorbases.");
        gmic_simple_item("-cmy2rgb",CMYtoRGB,"Convert image%s from CMY to RGB colorbases.");
        gmic_simple_item("-lab2rgb",LabtoRGB,"Convert image%s from Lab to RGB colorbases.");
        gmic_simple_item("-xyz2rgb",XYZtoRGB,"Convert image%s from XYZ to RGB colorbases.");
        gmic_simple_item("-ycbcr2rgb",YCbCrtoRGB,"Convert image%s from YCbCr to RGB colorbases.");
        gmic_simple_item("-yuv2rgb",YUVtoRGB,"Convert image%s from YUV to RGB colorbases.");
        gmic_simple_item("-hsi2rgb",HSItoRGB,"Convert image%s from HSI to RGB colorbases.");
        gmic_simple_item("-hsl2rgb",HSLtoRGB,"Convert image%s from HSL to RGB colorbases.");
        gmic_simple_item("-hsv2rgb",HSVtoRGB,"Convert image%s from HSV to RGB colorbases.");

        // Apply LUT.
        if (!cimg::strcmp("-lut2rgb",item0)) {
          int nb = 0, ind0 = 0; char sep = 0, end = 0;
          if (std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']') {
            gmic_check_indice(ind0,"Map LUT on image%s");
            print("Map LUT [%d] on image%s.",ind0,gmic_inds);
            const CImg<T> palette = images[ind0];
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],pointwise_norm().LUTtoRGB(palette));
          } else if (std::sscanf(argument,"%d%c",&nb,&end)==1) {
            if (nb<0 || nb>2) error("Map LUT on image%s : Invalid LUT number %d.",gmic_inds,nb);
            print("Map %s LUT on image%s.",nb==0?"default":nb==1?"rainbow":"cluster",gmic_inds);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],pointwise_norm().
                                              LUTtoRGB(nb==0?CImg<T>::default_LUT8():nb==1?CImg<T>::rainbow_LUT8():CImg<T>::contrast_LUT8()));
          } else error("Map LUT on image%s : Invalid argument '%s' "
                       "(should be '[indice]' or '{0,1,2}').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Convert to LUT.
        if (!cimg::strcmp("-rgb2lut",item0)) {
          int nb = 0, ind0 = 0, dithering = 0; char sep = 0, end = 0;
          if (std::sscanf(argument,"[%d%c%*c%d",&ind0,&sep,&dithering)>=2 && sep==']') {
            gmic_check_indice(ind0,"Index image%s with LUT");
            print("Index image%s with LUT [%d] %s dithering.",gmic_inds,ind0,dithering?"with":"without");
            const CImg<T> palette = images[ind0];
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],RGBtoLUT(palette,dithering,true));
          } else if (std::sscanf(argument,"%d%*c%d%c",&nb,&dithering,&end)==2 ||
                     std::sscanf(argument,"%d%c",&nb,&end)==1) {
            if (nb<0 || nb>2) error("Index image%s with LUT : Invalid LUT number %d.",gmic_inds,nb);
            print("Index image%s with %s LUT.",gmic_inds,nb==0?"default":nb==1?"rainbow":"cluster");
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],
                                              RGBtoLUT(nb==0?CImg<T>::default_LUT8():nb==1?CImg<T>::rainbow_LUT8():CImg<T>::contrast_LUT8(),
                                                       dithering,true));
          } else error("Index image%s with LUT : Invalid argument '%s' "
                       "(should be '[indice][,dithering]', or '{0,1,2}[,dithering]').",gmic_inds,argument_text);
          ++position; continue;
        }

        //-----------------------
        // Geometric manipulation
        //-----------------------

        // Resize.
        if (!cimg::strcmp("-resize",item0) || !cimg::strcmp("-r",item0)) {
          char
            sepx = 0, sepy = 0, sepz = 0, sepv = 0, end = 0,
            argx[4096] = { 0 }, argy[4096] = { 0 }, argz[4096] = { 0 }, argv[4096] = { 0 };
            float valx = 0, valy = 0, valz = 0, valv = 0;
            int interpolation = 1, borders = -1, center = 0, indx = no_ind, indy = no_ind, indz = no_ind, indv = no_ind;
            if ((std::sscanf(argument,"[%d%c%c",&indx,&sepx,&end)==2 ||
                 std::sscanf(argument,"[%d%c%*c%d%c",&indx,&sepx,&interpolation,&end)==3 ||
                 std::sscanf(argument,"[%d%c%*c%d%*c%d%c",&indx,&sepx,&interpolation,&borders,&end)==4 ||
                 std::sscanf(argument,"[%d%c%*c%d%*c%d%*c%d%c",&indx,&sepx,&interpolation,&borders,&center,&end)==5)
                && sepx==']') {
              gmic_check_indice(indx,"Resize image%s");
              const int
                ivalx = images[indx].dimx(),
                ivaly = images[indx].dimy(),
                ivalz = images[indx].dimz(),
                ivalv = images[indx].dimv();
              print("Resize image%s to %dx%dx%dx%d with %s interpolation.",
                    gmic_inds,ivalx,ivaly,ivalz,ivalv,
                    interpolation==0?"no":interpolation==1?"nearest neighbor":
                    interpolation==2?"moving average":interpolation==3?"linear":
                    interpolation==4?"grid":"cubic");
              cimg_foroff(indices,l) gmic_apply(images[indices[l]],resize(ivalx,ivaly,ivalz,ivalv,interpolation,borders,center?true:false));
              ++position;
            } else if (((std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%d%*c%d%*c%d%c",
                                     argx,argy,argz,argv,&(interpolation=1),&(borders=-1),&(center=0),&end)==7 ||
                         std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%d%*c%d%c",
                                     argx,argy,argz,argv,&interpolation,&borders,&end)==6 ||
                         std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%d%c",
                                     argx,argy,argz,argv,&interpolation,&end)==5 ||
                         std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%c",
                                     argx,argy,argz,argv,&end)==4) &&
                        ((std::sscanf(argx,"[%d%c%c",&indx,&sepx,&end)==2 && sepx==']') ||
                         (std::sscanf(argx,"%f%c%c",&valx,&sepx,&end)==2 && sepx=='%') ||
                         std::sscanf(argx,"%f%c",&valx,&end)==1) &&
                        ((std::sscanf(argy,"[%d%c%c",&indy,&sepy,&end)==2 && sepy==']') ||
                         (std::sscanf(argy,"%f%c%c",&valy,&sepy,&end)==2 && sepy=='%') ||
                         std::sscanf(argy,"%f%c",&valy,&end)==1) &&
                        ((std::sscanf(argz,"[%d%c%c",&indz,&sepz,&end)==2 && sepz==']') ||
                         (std::sscanf(argz,"%f%c%c",&valz,&sepz,&end)==2 && sepz=='%') ||
                         std::sscanf(argz,"%f%c",&valz,&end)==1) &&
                        ((std::sscanf(argv,"[%d%c%c",&indv,&sepv,&end)==2 && sepv==']') ||
                         (std::sscanf(argv,"%f%c%c",&valv,&sepv,&end)==2 && sepv=='%') ||
                         std::sscanf(argv,"%f%c",&valv,&end)==1)) ||
                       ((std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%c",
                                     argx,argy,argz,&end)==3) &&
                        ((std::sscanf(argx,"[%d%c%c",&indx,&sepx,&end)==2 && sepx==']') ||
                         (std::sscanf(argx,"%f%c%c",&valx,&sepx,&end)==2 && sepx=='%') ||
                         std::sscanf(argx,"%f%c",&valx,&end)==1) &&
                        ((std::sscanf(argy,"[%d%c%c",&indy,&sepy,&end)==2 && sepy==']') ||
                         (std::sscanf(argy,"%f%c%c",&valy,&sepy,&end)==2 && sepy=='%') ||
                         std::sscanf(argy,"%f%c",&valy,&end)==1) &&
                        ((std::sscanf(argz,"[%d%c%c",&indz,&sepz,&end)==2 && sepz==']') ||
                         (std::sscanf(argz,"%f%c%c",&valz,&sepz,&end)==2 && sepz=='%') ||
                         std::sscanf(argz,"%f%c",&valz,&end)==1)) ||
                       ((std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%c",
                                     argx,argy,&end)==2) &&
                        ((std::sscanf(argx,"[%d%c%c",&indx,&sepx,&end)==2 && sepx==']') ||
                         (std::sscanf(argx,"%f%c%c",&valx,&sepx,&end)==2 && sepx=='%') ||
                         std::sscanf(argx,"%f%c",&valx,&end)==1) &&
                        ((std::sscanf(argy,"[%d%c%c",&indy,&sepy,&end)==2 && sepy==']') ||
                         (std::sscanf(argy,"%f%c%c",&valy,&sepy,&end)==2 && sepy=='%') ||
                         std::sscanf(argy,"%f%c",&valy,&end)==1)) ||
                       ((std::sscanf(argument,"%4095[][0-9.eE%+-]%c",
                                     argx,&end)==1) &&
                        ((std::sscanf(argx,"[%d%c%c",&indx,&sepx,&end)==2 && sepx==']') ||
                         (std::sscanf(argx,"%f%c%c",&valx,&sepx,&end)==2 && sepx=='%') ||
                         std::sscanf(argx,"%f%c",&valx,&end)==1))) {
              if (indx!=no_ind) { gmic_check_indice(indx,"Resize image%s"); valx = (float)images[indx].dimx(); sepx = 0; }
              if (indy!=no_ind) { gmic_check_indice(indy,"Resize image%s"); valy = (float)images[indy].dimy(); sepy = 0; }
              if (indz!=no_ind) { gmic_check_indice(indz,"Resize image%s"); valz = (float)images[indz].dimz(); sepz = 0; }
              if (indv!=no_ind) { gmic_check_indice(indv,"Resize image%s"); valv = (float)images[indv].dimv(); sepv = 0; }
              if (!valx) { valx = 100; sepx = '%'; }
              if (!valy) { valy = 100; sepy = '%'; }
              if (!valz) { valz = 100; sepz = '%'; }
              if (!valv) { valv = 100; sepv = '%'; }
              print("Resize image%s to %g%s%g%s%g%s%g%s with %s interpolation.",
                    gmic_inds,valx,sepx=='%'?"%x":"x",valy,sepy=='%'?"%x":"x",valz,
                    sepz=='%'?"%x":"x",valv,sepv=='%'?"% ":" ",
                    interpolation==0?"no":interpolation==1?"nearest neighbor":
                    interpolation==2?"moving average":interpolation==3?"linear":
                    interpolation==4?"grid":"cubic");

              cimg_foroff(indices,l) {
                CImg<T>& img = images[indices[l]];
                const int
                  ivalx0 = (int)cimg::round(sepx=='%'?valx*img.dimx()/100:valx,1),
                  ivaly0 = (int)cimg::round(sepy=='%'?valy*img.dimy()/100:valy,1),
                  ivalz0 = (int)cimg::round(sepz=='%'?valz*img.dimz()/100:valz,1),
                  ivalv0 = (int)cimg::round(sepv=='%'?valv*img.dimv()/100:valv,1),
                  ivalx = ivalx0?ivalx0:1,
                  ivaly = ivaly0?ivaly0:1,
                  ivalz = ivalz0?ivalz0:1,
                  ivalv = ivalv0?ivalv0:1;
                gmic_apply(img,resize(ivalx,ivaly,ivalz,ivalv,interpolation,borders,center?true:false));
              }
              ++position;
            } else {
              print("Resize image%s : Interactive mode.",gmic_inds);
              char title[4096] = { 0 };
              cimg_foroff(indices,l) {
                CImg<T>& img = images[indices[l]];
                CImgDisplay disp(img,0,1);
                std::sprintf(title,"%s : Interactive resize",filenames[indices[l]].ptr());
                disp.set_title(title);
                img.get_select(0,disp);
                print("Resize image [%d] to %dx%d.",indices[l],disp.dimx(),disp.dimy());
                gmic_apply(img,resize(disp));
              }
            }
            continue;
        }

        // Resize2x. and Resize3x.
        gmic_simple_item("-resize2x",resize_doubleXY,"Resize image%s using Scale2x algorithm.");
        gmic_simple_item("-resize3x",resize_doubleXY,"Resize image%s using Scale3x algorithm.");

        // Crop.
        if (!cimg::strcmp("-crop",item0) || !cimg::strcmp("-c",item0)) {
          char st0[4096] = { 0 }, st1[4096] = { 0 }, st2[4096] = { 0 }, st3[4096] = { 0 };
          char st4[4096] = { 0 }, st5[4096] = { 0 }, st6[4096] = { 0 }, st7[4096] = { 0 };
          char sep0 = 0, sep1 = 0, sep2 = 0, sep3 = 0, sep4 = 0, sep5 = 0, sep6 = 0, sep7 = 0, end = 0;
          float a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0, a6 = 0, a7 = 0; int borders = 0;
          if ((std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c"
                           "%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%d%c",
                           st0,st1,st2,st3,st4,st5,st6,st7,&borders,&end)==9 ||
               std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c"
                           "%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",
                           st0,st1,st2,st3,st4,st5,st6,st7,&end)==8) &&
              (std::sscanf(st0,"%f%c",&a0,&end)==1 || (std::sscanf(st0,"%f%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
              (std::sscanf(st1,"%f%c",&a1,&end)==1 || (std::sscanf(st1,"%f%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
              (std::sscanf(st2,"%f%c",&a2,&end)==1 || (std::sscanf(st2,"%f%c%c",&a2,&sep2,&end)==2 && sep2=='%')) &&
              (std::sscanf(st3,"%f%c",&a3,&end)==1 || (std::sscanf(st3,"%f%c%c",&a3,&sep3,&end)==2 && sep3=='%')) &&
              (std::sscanf(st4,"%f%c",&a4,&end)==1 || (std::sscanf(st4,"%f%c%c",&a4,&sep4,&end)==2 && sep4=='%')) &&
              (std::sscanf(st5,"%f%c",&a5,&end)==1 || (std::sscanf(st5,"%f%c%c",&a5,&sep5,&end)==2 && sep5=='%')) &&
              (std::sscanf(st6,"%f%c",&a6,&end)==1 || (std::sscanf(st6,"%f%c%c",&a6,&sep6,&end)==2 && sep6=='%')) &&
              (std::sscanf(st7,"%f%c",&a7,&end)==1 || (std::sscanf(st7,"%f%c%c",&a7,&sep7,&end)==2 && sep7=='%'))) {
            print("Crop image%s with (%g%s%g%s%g%s%g%s x (%g%s%g%s%g%s%g%s.",gmic_inds,
                  a0,sep0=='%'?"%,":",",a1,sep1=='%'?"%,":",",
                  a2,sep2=='%'?"%,":",",a3,sep3=='%'?"%)":")",
                  a4,sep4=='%'?"%,":",",a5,sep5=='%'?"%,":",",
                  a6,sep6=='%'?"%,":",",a7,sep7=='%'?"%)":")");
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int
                x0 = (int)cimg::round(sep0=='%'?a0*img.dimx()/100:a0,1),
                y0 = (int)cimg::round(sep1=='%'?a1*img.dimy()/100:a1,1),
                z0 = (int)cimg::round(sep2=='%'?a2*img.dimz()/100:a2,1),
                v0 = (int)cimg::round(sep3=='%'?a3*img.dimv()/100:a3,1),
                x1 = (int)cimg::round(sep4=='%'?a4*img.dimx()/100:a4,1),
                y1 = (int)cimg::round(sep5=='%'?a5*img.dimy()/100:a5,1),
                z1 = (int)cimg::round(sep6=='%'?a6*img.dimz()/100:a6,1),
                v1 = (int)cimg::round(sep7=='%'?a7*img.dimv()/100:a7,1);
              gmic_apply(img,crop(x0,y0,z0,v0,x1,y1,z1,v1,borders?true:false));
            }
            ++position;
          } else if ((std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c"
                                  "%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%d%c",
                                  st0,st1,st2,st3,st4,st5,&borders,&end)==7 ||
                      std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c"
                                  "%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",
                                  st0,st1,st2,st3,st4,st5,&end)==6) &&
                     (std::sscanf(st0,"%f%c",&a0,&end)==1 ||
                      (std::sscanf(st0,"%f%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
                     (std::sscanf(st1,"%f%c",&a1,&end)==1 ||
                      (std::sscanf(st1,"%f%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
                     (std::sscanf(st2,"%f%c",&a2,&end)==1 ||
                      (std::sscanf(st2,"%f%c%c",&a2,&sep2,&end)==2 && sep2=='%')) &&
                     (std::sscanf(st3,"%f%c",&a3,&end)==1 ||
                      (std::sscanf(st3,"%f%c%c",&a3,&sep3,&end)==2 && sep3=='%')) &&
                     (std::sscanf(st4,"%f%c",&a4,&end)==1 ||
                      (std::sscanf(st4,"%f%c%c",&a4,&sep4,&end)==2 && sep4=='%')) &&
                     (std::sscanf(st5,"%f%c",&a5,&end)==1 ||
                      (std::sscanf(st5,"%f%c%c",&a5,&sep5,&end)==2 && sep5=='%'))) {
            print("Crop image%s with (%g%s%g%s%g%s x (%g%s%g%s%g%s.",gmic_inds,
                  a0,sep0=='%'?"%,":",",a1,sep1=='%'?"%,":",",a2,sep2=='%'?"%)":")",
                  a3,sep3=='%'?"%,":",",a4,sep4=='%'?"%,":",",a5,sep5=='%'?"%)":")");
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int
                x0 = (int)cimg::round(sep0=='%'?a0*img.dimx()/100:a0,1),
                y0 = (int)cimg::round(sep1=='%'?a1*img.dimy()/100:a1,1),
                z0 = (int)cimg::round(sep2=='%'?a2*img.dimz()/100:a2,1),
                x1 = (int)cimg::round(sep3=='%'?a3*img.dimx()/100:a3,1),
                y1 = (int)cimg::round(sep4=='%'?a4*img.dimy()/100:a4,1),
                z1 = (int)cimg::round(sep5=='%'?a5*img.dimz()/100:a5,1);
              gmic_apply(img,crop(x0,y0,z0,x1,y1,z1,borders?true:false));
            }
            ++position;
          } else if ((std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c"
                                  "%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%d%c",
                                  st0,st1,st2,st3,&borders,&end)==5 ||
                      std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c"
                                  "%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",
                                  st0,st1,st2,st3,&end)==4) &&
                     (std::sscanf(st0,"%f%c",&a0,&end)==1 ||
                      (std::sscanf(st0,"%f%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
                     (std::sscanf(st1,"%f%c",&a1,&end)==1 ||
                      (std::sscanf(st1,"%f%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
                     (std::sscanf(st2,"%f%c",&a2,&end)==1 ||
                      (std::sscanf(st2,"%f%c%c",&a2,&sep2,&end)==2 && sep2=='%')) &&
                     (std::sscanf(st3,"%f%c",&a3,&end)==1 ||
                      (std::sscanf(st3,"%f%c%c",&a3,&sep3,&end)==2 && sep3=='%'))) {
            print("Crop image%s with (%g%s%g%s x (%g%s%g%s.",gmic_inds,
                  a0,sep0=='%'?"%,":",",a1,sep1=='%'?"%)":")",
                  a2,sep2=='%'?"%,":",",a3,sep3=='%'?"%)":")");
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int
                x0 = (int)cimg::round(sep0=='%'?a0*img.dimx()/100:a0,1),
                y0 = (int)cimg::round(sep1=='%'?a1*img.dimy()/100:a1,1),
                x1 = (int)cimg::round(sep2=='%'?a2*img.dimx()/100:a2,1),
                y1 = (int)cimg::round(sep3=='%'?a3*img.dimy()/100:a3,1);
              gmic_apply(img,crop(x0,y0,x1,y1,borders?true:false));
            }
            ++position;
          } else if ((std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%d%c",st0,st1,&borders,&end)==3 ||
                      std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",st0,st1,&end)==2) &&
                     (std::sscanf(st0,"%f%c",&a0,&end)==1 ||
                      (std::sscanf(st0,"%f%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
                     (std::sscanf(st1,"%f%c",&a1,&end)==1 ||
                      (std::sscanf(st1,"%f%c%c",&a1,&sep1,&end)==2 && sep1=='%'))) {
            print("Crop image%s with (%g%s x (%g%s.",gmic_inds,
                  a0,sep0=='%'?"%)":")",a1,sep1=='%'?"%)":")");
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int
                x0 = (int)cimg::round(sep0=='%'?a0*img.dimx()/100:a0,1),
                x1 = (int)cimg::round(sep1=='%'?a1*img.dimx()/100:a1,1);
              gmic_apply(img,crop(x0,x1,borders?true:false));
            }
            ++position;
          } else {
            print("Crop image%s : Interactive mode.",gmic_inds);
            char title[4096] = { 0 };
            cimg_foroff(indices,l) {
              CImg<T>& img = images[indices[l]];
              CImgDisplay disp(cimg_fitscreen(img.dimx(),img.dimy(),1),0,1);
              std::sprintf(title,"%s : Interactive crop",filenames[indices[l]].ptr());
              disp.set_title(title);
              const CImg<int> s = img.get_select(disp,2);
              print("Crop image [%d] with (%d,%d,%d) x (%d,%d,%d).",
                    indices[l],s[0],s[1],s[2],s[3],s[4],s[5]);
              gmic_apply(img,crop(s[0],s[1],s[2],s[3],s[4],s[5]));
            }
          }
          continue;
        }

        // Autocrop.
        if (!cimg::strcmp("-autocrop",item0)) {
          print("Autocrop image%s with color '%s'.",gmic_inds,argument_text);
          cimg_foroff(indices,l) {
            CImg<T>& img = images[indices[l]];
            const CImg<T> col = CImg<T>(img.dimv()).fill(argument,true);
            gmic_apply(img,autocrop(col));
          }
          ++position; continue;
        }

        // Select channels.
        if (!cimg::strcmp("-channels",item0)) {
          char sep0 = 0, sep1 = 0, end = 0, arg0[4096] = { 0 }, arg1[4096] = { 0 };
          float value0 = 0, value1 = 0; int ind0 = no_ind, ind1 = no_ind;
          if (std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%c",arg0,arg1,&end)==2 &&
              ((std::sscanf(arg0,"[%d%c%c]",&ind0,&sep0,&end)==2 && sep0==']') ||
               (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
               std::sscanf(arg0,"%f%c",&value0,&end)==1) &&
              ((std::sscanf(arg1,"[%d%c%c]",&ind1,&sep1,&end)==2 && sep1==']') ||
               (std::sscanf(arg1,"%f%c%c",&value1,&sep1,&end)==2 && sep1=='%') ||
               std::sscanf(arg1,"%f%c",&value1,&end)==1)) {
            if (ind0!=no_ind) { gmic_check_indice(ind0,"Keep channels of image%s"); value0 = images[ind0].dimv()-1.0f; sep0 = 0; }
            if (ind1!=no_ind) { gmic_check_indice(ind1,"Keep channels of image%s"); value1 = images[ind1].dimv()-1.0f; sep1 = 0; }
            print("Keep channels %g%s..%g%s of image%s.",value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int
                nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.dimv()-1)/100:value0,1),
                nvalue1 = (int)cimg::round(sep1=='%'?value1*(img.dimv()-1)/100:value1,1);
              gmic_apply(img,channels(nvalue0,nvalue1));
            }
          } else if (std::sscanf(argument,"%4095[][0-9.eE%+-]%c",arg0,&end)==1 &&
                     ((std::sscanf(arg0,"[%d%c%c]",&ind0,&sep0,&end)==2 && sep0==']') ||
                      (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
                      std::sscanf(arg0,"%f%c",&value0,&end)==1)) {
            if (ind0!=no_ind) { gmic_check_indice(ind0,"Keep channel of image%s"); value0 = images[ind0].dimv()-1.0f; sep0 = 0; }
            print("Keep channel %g%s of image%s.",value0,sep0=='%'?"%":"",gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.dimv()-1)/100:value0,1);
              gmic_apply(img,channel(nvalue0));
            }
          } else error("Keep channels of image%s : Invalid argument '%s' "
                       "(should be 'channel0[%%][,channel1[%%]]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Select slices.
        if (!cimg::strcmp("-slices",item0)) {
          char sep0 = 0, sep1 = 0, end = 0, arg0[4096] = { 0 }, arg1[4096] = { 0 };
          float value0 = 0, value1 = 0; int ind0 = no_ind, ind1 = no_ind;
          if (std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%c",arg0,arg1,&end)==2 &&
              ((std::sscanf(arg0,"[%d%c%c]",&ind0,&sep0,&end)==2 && sep0==']') ||
               (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
               std::sscanf(arg0,"%f%c",&value0,&end)==1) &&
              ((std::sscanf(arg1,"[%d%c%c]",&ind1,&sep1,&end)==2 && sep1==']') ||
               (std::sscanf(arg1,"%f%c%c",&value1,&sep1,&end)==2 && sep1=='%') ||
               std::sscanf(arg1,"%f%c",&value1,&end)==1)) {
            if (ind0!=no_ind) { gmic_check_indice(ind0,"Select slices of image%s"); value0 = images[ind0].dimz()-1.0f; sep0 = 0; }
            if (ind1!=no_ind) { gmic_check_indice(ind1,"Select slices of image%s"); value1 = images[ind1].dimz()-1.0f; sep1 = 0; }
            print("Select slices %g%s..%g%s of image%s.",value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int
                nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.dimv()-1)/100:value0,0),
                nvalue1 = (int)cimg::round(sep1=='%'?value1*(img.dimv()-1)/100:value1,0);
              gmic_apply(img,slices(nvalue0,nvalue1));
            }
          } else if (std::sscanf(argument,"%4095[][0-9.eE%+-]%c",arg0,&end)==1 &&
                     ((std::sscanf(arg0,"[%d%c%c]",&ind0,&sep0,&end)==2 && sep0==']') ||
                      (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
                      std::sscanf(arg0,"%f%c",&value0,&end)==1)) {
            if (ind0!=no_ind) { gmic_check_indice(ind0,"Select slice of image%s"); value0 = images[ind0].dimz()-1.0f; sep0 = 0; }
            print("Select slice %g%s of image%s.",value0,sep0=='%'?"%":"",gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.dimz()-1)/100:value0,1);
              gmic_apply(img,slice(nvalue0));
            }
          } else error("Select slices of image%s : Invalid argument '%s' "
                       "(should be 'slice0[%%][,slice1[%%]]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Select lines.
        if (!cimg::strcmp("-lines",item0) || !cimg::strcmp("-l",item0)) {
          char sep0 = 0, sep1 = 0, end = 0, arg0[4096] = { 0 }, arg1[4096] = { 0 };
          float value0 = 0, value1 = 0; int ind0 = no_ind, ind1 = no_ind;
          if (std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%c",arg0,arg1,&end)==2 &&
              ((std::sscanf(arg0,"[%d%c%c]",&ind0,&sep0,&end)==2 && sep0==']') ||
               (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
               std::sscanf(arg0,"%f%c",&value0,&end)==1) &&
              ((std::sscanf(arg1,"[%d%c%c]",&ind1,&sep1,&end)==2 && sep1==']') ||
               (std::sscanf(arg1,"%f%c%c",&value1,&sep1,&end)==2 && sep1=='%') ||
               std::sscanf(arg1,"%f%c",&value1,&end)==1)) {
            if (ind0!=no_ind) { gmic_check_indice(ind0,"Select lines of image%s"); value0 = images[ind0].dimy()-1.0f; sep0 = 0; }
            if (ind1!=no_ind) { gmic_check_indice(ind1,"Select lines of image%s"); value1 = images[ind1].dimy()-1.0f; sep1 = 0; }
            print("Select lines %g%s..%g%s of image%s.",value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int
                nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.dimy()-1)/100:value0,1),
                nvalue1 = (int)cimg::round(sep1=='%'?value1*(img.dimy()-1)/100:value1,1);
              gmic_apply(img,lines(nvalue0,nvalue1));
            }
          } else if (std::sscanf(argument,"%4095[][0-9.eE%+-]%c",arg0,&end)==1 &&
                     ((std::sscanf(arg0,"[%d%c%c]",&ind0,&sep0,&end)==2 && sep0==']') ||
                      (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
                      std::sscanf(arg0,"%f%c",&value0,&end)==1)) {
            if (ind0!=no_ind) { gmic_check_indice(ind0,"Select lines of image%s"); value0 = images[ind0].dimy()-1.0f; sep0 = 0; }
            print("Select lines %g%s of image%s.",value0,sep0=='%'?"%":"",gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.dimy()-1)/100:value0,1);
              gmic_apply(img,line(nvalue0));
            }
          } else error("Select lines of image%s : Invalid argument '%s' "
                       "(should be 'line0[%%][,line1[%%]]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Columns.
        if (!cimg::strcmp("-columns",item0)) {
          char sep0 = 0, sep1 = 0, end = 0, arg0[4096] = { 0 }, arg1[4096] = { 0 };
          float value0 = 0, value1 = 0; int ind0 = no_ind, ind1 = no_ind;
          if (std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%c",arg0,arg1,&end)==2 &&
              ((std::sscanf(arg0,"[%d%c%c]",&ind0,&sep0,&end)==2 && sep0==']') ||
               (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
               std::sscanf(arg0,"%f%c",&value0,&end)==1) &&
              ((std::sscanf(arg1,"[%d%c%c]",&ind1,&sep1,&end)==2 && sep1==']') ||
               (std::sscanf(arg1,"%f%c%c",&value1,&sep1,&end)==2 && sep1=='%') ||
               std::sscanf(arg1,"%f%c",&value1,&end)==1)) {
            if (ind0!=no_ind) { gmic_check_indice(ind0,"Select columns of image%s"); value0 = images[ind0].dimx()-1.0f; sep0 = 0; }
            if (ind1!=no_ind) { gmic_check_indice(ind1,"Select columns of image%s"); value1 = images[ind1].dimx()-1.0f; sep1 = 0; }
            print("Select columns %g%s..%g%s of image%s.",value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int
                nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.dimx()-1)/100:value0,1),
                nvalue1 = (int)cimg::round(sep1=='%'?value1*(img.dimx()-1)/100:value1,1);
              gmic_apply(img,lines(nvalue0,nvalue1));
            }
          } else if (std::sscanf(argument,"%4095[][0-9.eE%+-]%c",arg0,&end)==1 &&
                     ((std::sscanf(arg0,"[%d%c%c]",&ind0,&sep0,&end)==2 && sep0==']') ||
                      (std::sscanf(arg0,"%f%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
                      std::sscanf(arg0,"%f%c",&value0,&end)==1)) {
            if (ind0!=no_ind) { gmic_check_indice(ind0,"Select columns of image%s"); value0 = images[ind0].dimx()-1.0f; sep0 = 0; }
            print("Select columns %g%s of image%s.",value0,sep0=='%'?"%":"",gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int nvalue0 = (int)cimg::round(sep0=='%'?value0*(img.dimx()-1)/100:value0,1);
              gmic_apply(img,line(nvalue0));
            }
          } else error("Select columns of image%s : Invalid argument '%s' "
                       "(should be 'column0[%%][,column1[%%]]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Rotate.
        if (!cimg::strcmp("-rotate",item0)) {
          float angle = 0; int borders = 0, interpolation = 1; char end = 0;
          if (std::sscanf(argument,"%f%c",&angle,&end)==1 ||
              std::sscanf(argument,"%f%*c%d%c",&angle,&borders,&end)==2 ||
              std::sscanf(argument,"%f%*c%d%*c%d%c",&angle,&borders,&interpolation,&end)==3) {
            print("Rotate image%s with an angle of %g deg and %s interpolation.",
                  gmic_inds,angle,interpolation?"linear":"nearest-neighbor");
            if (borders>=0) { cimg_foroff(indices,l) gmic_apply(images[indices[l]],rotate(angle,borders,interpolation)); }
            else cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              gmic_apply(img,rotate(angle,img.dimx()/2.0f,img.dimy()/2.0f,1,-1-borders,interpolation));
            }
          } else error("Rotate image%s : Invalid argument '%s' "
                       "(should be 'angle[,border_conditions[,interpolation]]').",gmic_inds,argument_text);
          ++position;
          continue;
        }

        // Mirror.
        if (!cimg::strcmp("-mirror",item0)) {
          const char axis = cimg::uncase(*argument);
          if (cimg::strlen(argument)==1) {
            print("Mirror image%s along the %c-axis.",gmic_inds,axis);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],mirror(axis));
          } else error("Mirror image%s : Invalid argument '%s' "
                       "(should be '{x,y,z,v}').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Translate.
        if (!cimg::strcmp("-translate",item0)) {
          char stx[4096] = { 0 }, sty[4096] = { 0 }, stz[4096] = { 0 }, stv[4096] = { 0 };
          char sepx = 0, sepy = 0, sepz = 0, sepv = 0, end = 0;
          float dx = 0, dy = 0, dz = 0, dv = 0; int borders = 0;
          if (((std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%d%c",
                            stx,sty,stz,stv,&borders,&end)==5 ||
                std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",
                            stx,sty,stz,stv,&end)==4) &&
               (std::sscanf(stx,"%f%c",&dx,&end)==1 || (std::sscanf(stx,"%f%c%c",&dx,&sepx,&end)==2 && sepx=='%')) &&
               (std::sscanf(sty,"%f%c",&dy,&end)==1 || (std::sscanf(sty,"%f%c%c",&dy,&sepy,&end)==2 && sepy=='%')) &&
               (std::sscanf(stz,"%f%c",&dz,&end)==1 || (std::sscanf(stz,"%f%c%c",&dz,&sepz,&end)==2 && sepz=='%')) &&
               (std::sscanf(stv,"%f%c",&dv,&end)==1 || (std::sscanf(stv,"%f%c%c",&dv,&sepv,&end)==2 && sepv=='%'))) ||
              (std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",stx,sty,stz,&end)==3 &&
               (std::sscanf(stx,"%f%c",&dx,&end)==1 || (std::sscanf(stx,"%f%c%c",&dx,&sepx,&end)==2 && sepx=='%')) &&
               (std::sscanf(sty,"%f%c",&dy,&end)==1 || (std::sscanf(sty,"%f%c%c",&dy,&sepy,&end)==2 && sepy=='%')) &&
               (std::sscanf(stz,"%f%c",&dz,&end)==1 || (std::sscanf(stz,"%f%c%c",&dz,&sepz,&end)==2 && sepz=='%'))) ||
              (std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",stx,sty,&end)==2 &&
               (std::sscanf(stx,"%f%c",&dx,&end)==1 || (std::sscanf(stx,"%f%c%c",&dx,&sepx,&end)==2 && sepx=='%')) &&
               (std::sscanf(sty,"%f%c",&dy,&end)==1 || (std::sscanf(sty,"%f%c%c",&dy,&sepy,&end)==2 && sepy=='%'))) ||
              (std::sscanf(argument,"%4095[0-9.eE%+-]%c",stx,&end)==1 &&
               (std::sscanf(stx,"%f%c",&dx,&end)==1 || (std::sscanf(stx,"%f%c%c",&dx,&sepx,&end)==2 && sepx=='%')))) {
            print("Translate image%s with vector (%g%s,%g%s,%g%s,%g%s).",
                  gmic_inds,dx,sepx=='%'?"%":"",dy,sepy=='%'?"%":"",dz,sepz=='%'?"%":"",dv,sepv=='%'?"%":"");
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int
                ndx = (int)cimg::round(sepx=='%'?dx*img.dimx()/100:dx,1),
                ndy = (int)cimg::round(sepy=='%'?dy*img.dimy()/100:dy,1),
                ndz = (int)cimg::round(sepz=='%'?dz*img.dimz()/100:dz,1),
                ndv = (int)cimg::round(sepv=='%'?dv*img.dimv()/100:dv,1);
              gmic_apply(images[indices[l]],translate(ndx,ndy,ndz,ndv,borders));
            }
          } else error("Translate image%s : Invalid argument '%s' "
                       "(should be 'tx[%%][,ty[%%][,tz[%%][,tv[%%][,border_conditions]]]]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Transpose.
        gmic_simple_item("-transpose",transpose,"Transpose image%s.");

        // Invert.
        gmic_simple_item("-invert",invert,"Compute matrix inversion of image%s.");

        // Permute axes.
        if (!cimg::strcmp("-permute",item0)) {
          print("Permute axes of image%s with permutation '%s'.",gmic_inds,argument_text);
          cimg_foroff(indices,l) gmic_apply(images[indices[l]],permute_axes(argument));
          ++position; continue;
        }

        // Unroll.
        if (!cimg::strcmp("-unroll",item0)) {
          const char axis = cimg::uncase(*argument);
          if (cimg::strlen(argument)==1 && (axis=='x' || axis=='y' || axis=='z' || axis=='v')) {
            print("Unroll image%s along the %c-axis.",gmic_inds,axis);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],unroll(axis));
          } else error("Unroll image%s : Invalid argument '%s' "
                       "(should be '{x,y,z,v}').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Split image(s).
        if (!cimg::strcmp("-split",item0) || !cimg::strcmp("-s",item0)) {
          char axis = cimg::uncase(*argument), foo = 0, end = 0; int nb = 0, keep_value = 0; double value = 0;
          if ((std::sscanf(argument,"%c%c",&foo,&end)==1 ||
               std::sscanf(argument,"%c%*c%d%c",&foo,&nb,&end)==2) &
              (axis=='x' || axis=='y' || axis=='z' || axis=='v')) {
            if (nb<0) error("Split image%s along the %c-axis in %d part : Invalid number of parts.",
                            gmic_inds,axis,nb);
            if (nb>0) print("Split image%s along the %c-axis in %d parts.",gmic_inds,axis,nb);
            else print("Split image%s along the %c-axis.",gmic_inds,axis);
            unsigned int off = 0;
            cimg_foroff(indices,l) {
              const unsigned int ind = indices[l] + off;
              const CImg<T>& img = images[ind];
              const CImg<char> filename = filenames[ind];
              const CImgList<T> split = img.get_split(axis,nb);
              if (get_version) {
                images.insert(split);
                filenames.insert(split.size,filename);
              } else {
                images.remove(ind); images.insert(split,ind);
                filenames.remove(ind); filenames.insert(split.size,filename,ind);
                off+=split.size-1;
              }
            }
          } else if (std::sscanf(argument,"%lf%c",&value,&end)==1 ||
                     std::sscanf(argument,"%lf%*c%d%c",&value,&keep_value,&end)==2) {
            print("Split image%s according to value %g.",gmic_inds,value);
            unsigned int off = 0;
            cimg_foroff(indices,l) {
              const unsigned int ind = indices[l] + off;
              CImg<T>& img = images[ind];
              const CImg<char> filename = filenames[ind];
              const CImgList<T> split = img.get_split((T)value,keep_value,false);
              if (get_version) {
                images.insert(split);
                filenames.insert(split.size,filename);
              } else {
                images.remove(ind); images.insert(split,ind);
                filenames.remove(ind); filenames.insert(split.size,filename,ind);
                off+=split.size-1;
              }
            }
          } else error("Split image%s : Invalid argument '%s' "
                       "(should be 'axis[,nb_parts]' where 'axis' can be '{x,y,z,v}').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Append image(s).
        if (!cimg::strcmp("-append",item0) || !cimg::strcmp("-a",item0)) {
          char axis = 0, align='p', end = 0;
          if ((std::sscanf(argument,"%c%c",&axis,&end)==1 ||
               std::sscanf(argument,"%c%*c%c%c",&axis,&align,&end)==2)) {
            axis = cimg::uncase(axis);
            print("Append image%s along the %c-axis with %s alignment.",
                  gmic_inds,axis,align=='p'?"left":align=='c'?"center":"right");
            CImgList<T> subimages; cimg_foroff(indices,l) subimages.insert(images[indices[l]],~0U,true);
            if (get_version) {
              images.insert(subimages.get_append(axis,align));
              filenames.insert(filenames[indices[0]]);
            } else {
              images.insert(subimages.get_append(axis,align),indices[0]);
              filenames.insert(filenames[indices[0]],indices[0]);
              int off = 1;
              cimg_foroff(indices,l) {
                const int ind = indices[l] + off;
                images.remove(ind); filenames.remove(ind);
                --off;
              }
            }
          } else error("Append image%s : Invalid argument '%s' "
                       "(should be 'axis[,alignement]' where 'axis' can be '{x,y,z,v}' "
                       "and alignement '{p,c,n}').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Warp image(s).
        if (!cimg::strcmp("-warp",item0)) {
          int ind0 = no_ind, interpolation = 1, relative = 0, nb = 1, borders = 1; char end = 0, sep = 0;
          if ((std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']')||
              std::sscanf(argument,"[%d]%*c%d%c",&ind0,&relative,&end)==2 ||
              std::sscanf(argument,"[%d]%*c%d%*c%d%c",&ind0,&relative,&interpolation,&end)==3 ||
              std::sscanf(argument,"[%d]%*c%d%*c%d%*c%d%c",&ind0,&relative,&interpolation,&borders,&end)==4 ||
              std::sscanf(argument,"[%d]%*c%d%*c%d%*c%d%*c%d%c",&ind0,&relative,&interpolation,&borders,&nb,&end)==5) {
            gmic_check_indice(ind0,"Warp image%s");
            if (nb!=1) print("Warp image%s with %s field [%u] and %d frames.",
                             gmic_inds,relative?"relative":"absolute",ind0,nb);
            else print("Warp image%s with %s field [%u].",gmic_inds,relative?"relative":"absolute",ind0);
            if (nb>=1) {
              const CImg<T> warp = images[ind0];
              unsigned int off = 0;
              cimg_foroff(indices,l) {
                const unsigned int ind = indices[l] + off;
                CImg<T> &img = images[ind];
                CImgList<T> frames(nb);
                cimglist_for(frames,t) {
                  const CImg<T> nwarp = warp.get_resize(img.dimx(),img.dimy(),img.dimz(),warp.dimv(),3)*=(t+1.0f)/nb;
                  frames[t] = img.get_warp(nwarp,relative?true:false,interpolation?true:false,borders);
                }
                if (get_version) {
                  images.insert(frames);
                  filenames.insert(nb-1,filenames[ind]);
                } else {
                  images.remove(ind); images.insert(frames,ind);
                  filenames.insert(nb-1,filenames[ind],ind);
                  off+=nb-1;
                }
              }
            }
          } else error("Warp image%s : Invalid argument '%s' "
                       "(should be '[indice][,relative[,interpolation[,border_conditions[,nb_frames]]]]').",
                       gmic_inds,argument_text);
          ++position; continue;
        }

        //-----------------------
        // Image filtering
        //-----------------------

        // Gaussian blur.
        if (!cimg::strcmp("-blur",item0)) {
          float sigma = -1; int borders = 1; char end = 0;
          if ((std::sscanf(argument,"%f%c",&sigma,&end)==1 ||
               std::sscanf(argument,"%f%*c%d%c",&sigma,&borders,&end)==2)
              && sigma>=0) {
            print("Blur image%s with standard deviation %g.",gmic_inds,sigma);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],blur(sigma,borders?true:false));
          } else error("Blur image%s : Invalid argument '%s' "
                       "(should be 'stdev[,border_conditions]', with stdev>=0).",gmic_inds,argument_text);
          ++position; continue;
        }

        // Bilateral filter.
        if (!cimg::strcmp("-bilateral",item0)) {
          float sigmas = 0, sigmar = 0; char end = 0;
          if (std::sscanf(argument,"%f%*c%f%c",&sigmas,&sigmar,&end)==2) {
            print("Apply bilateral filter on image%s with standart deviations %g and %g.",
                  gmic_inds,sigmas,sigmar);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],blur_bilateral(sigmas,sigmar));
          } else error("Apply bilateral filter on image%s : Invalid argument '%s' "
                       "(should be 'stdevs,stdevr').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Smooth.
        if (!cimg::strcmp("-smooth",item0)) {
          float amplitude = 0, sharpness = 0.7f, anisotropy = 0.3f, alpha = 0.6f, sigma = 1.1f, dl =0.8f, da = 30.0f, gauss_prec = 2.0f;
          unsigned int interpolation_type = 0, fast_approx = 1;
          char end = 0;
          if (std::sscanf(argument,"%f%c",&amplitude,&end)==1 ||
              std::sscanf(argument,"%f%*c%f%c",&amplitude,&sharpness,&end)==2 ||
              std::sscanf(argument,"%f%*c%f%*c%f%c",&amplitude,&sharpness,&anisotropy,&end)==3 ||
              std::sscanf(argument,"%f%*c%f%*c%f%*c%f%c",&amplitude,&sharpness,&anisotropy,&alpha,&end)==4 ||
              std::sscanf(argument,"%f%*c%f%*c%f%*c%f%*c%f%c",&amplitude,&sharpness,&anisotropy,&alpha,&sigma,&end)==5 ||
              std::sscanf(argument,"%f%*c%f%*c%f%*c%f%*c%f%*c%f%c",&amplitude,&sharpness,&anisotropy,&alpha,&sigma,&dl,&end)==6 ||
              std::sscanf(argument,"%f%*c%f%*c%f%*c%f%*c%f%*c%f%*c%f%c",&amplitude,&sharpness,&anisotropy,&alpha,&sigma,&dl,&da,&end)==7 ||
              std::sscanf(argument,"%f%*c%f%*c%f%*c%f%*c%f%*c%f%*c%f%*c%f%c",
                          &amplitude,&sharpness,&anisotropy,&alpha,&sigma,&dl,&da,&gauss_prec,&end)==8 ||
              std::sscanf(argument,"%f%*c%f%*c%f%*c%f%*c%f%*c%f%*c%f%*c%f%*c%u%c",
                          &amplitude,&sharpness,&anisotropy,&alpha,&sigma,&dl,&da,&gauss_prec,&interpolation_type,&end)==9 ||
              std::sscanf(argument,"%f%*c%f%*c%f%*c%f%*c%f%*c%f%*c%f%*c%f%*c%u%*c%u%c",
                          &amplitude,&sharpness,&anisotropy,&alpha,&sigma,&dl,&da,&gauss_prec,&interpolation_type,&fast_approx,&end)==10) {
            print("Smooth image%s anisotropically with "
                  "amplitude %g, sharpness %g, anisotropy %g, alpha %g and sigma %g.",
                  gmic_inds,amplitude,sharpness,anisotropy,alpha,sigma);
            cimg_foroff(indices,l)
              gmic_apply(images[indices[l]],blur_anisotropic(amplitude,sharpness,anisotropy,alpha,sigma,
                                                             dl,da,gauss_prec,interpolation_type,fast_approx?true:false));
          } else error("Smooth image%s anisotropically : Invalid argument '%s' "
                       "(should be 'amplitude[,sharpness[,anisotropy[,alpha[,sigma[,dl[,da[,prec[,interp[,fast]]]]]]]]]').",
                       gmic_inds,argument_text);
          ++position; continue;
        }

        // Patch averaging.
        if (!cimg::strcmp("-denoise",item0)) {
          float sigmas = 10, sigmar = 10; int psize = 5, rsize = 6; char end = 0;
          if (std::sscanf(argument,"%f%c",&sigmas,&end)==1 ||
              std::sscanf(argument,"%f%*c%f%c",&sigmas,&sigmar,&end)==2 ||
              std::sscanf(argument,"%f%*c%f%*c%d%c",&sigmas,&sigmar,&psize,&end)==3 ||
              std::sscanf(argument,"%f%*c%f%*c%d%*c%d%c",&sigmas,&sigmar,&psize,&rsize,&end)==4) {
            if (sigmas<0 || sigmar<0 || psize<0 || rsize<0)
              error("Denoise image%s with %dx%d patches, standard deviations %lg,%g and lookup size %d : "
                    "Invalid parameters.",gmic_inds,psize,psize,sigmas,sigmar,rsize);
            print("Denoise image%s with %dx%d patches, standard deviations %lg,%g and lookup size %d.",
                  gmic_inds,psize,psize,sigmas,sigmar,rsize);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],blur_patch(psize,sigmas,sigmar,rsize));
          } else error("Denoise image%s : Invalid argument '%s' "
                       "(should be 'stdev_s[,stdev_p[,patch_size[,lookup_size]]]').",
                       gmic_inds,argument_text);
          ++position; continue;
        }

        // Median filter.
        if (!cimg::strcmp("-median",item0)) {
          int siz = 3; char end = 0;
          if (std::sscanf(argument,"%d%c",&siz,&end)==1) {
            if (siz<=0) error("Apply median filter on image%s : Invalid size %d.",gmic_inds,siz);
            print("Apply median filter of size %d on image%s.",siz,gmic_inds);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],blur_median(siz));
          } else error("Apply median filter on image%s : Invalid argument '%s' "
                       "(should be 'size').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Sharpen.
        if (!cimg::strcmp("-sharpen",item0)) {
          float amplitude = 0, edge = 1, alpha = 0, sigma = 0; int sharpen_type = 0; char end = 0;
          if (std::sscanf(argument,"%f%c",&amplitude,&end)==1 ||
              std::sscanf(argument,"%f%*c%d%c",&amplitude,&sharpen_type,&end)==2 ||
              std::sscanf(argument,"%f%*c%d%*c%f%c",&amplitude,&sharpen_type,&edge,&end)==3 ||
              std::sscanf(argument,"%f%*c%d%*c%f%*c%f%c",&amplitude,&sharpen_type,&edge,&alpha,&end)==4 ||
              std::sscanf(argument,"%f%*c%d%*c%f%*c%f%*c%f%c",&amplitude,&sharpen_type,&edge,&alpha,&sigma,&end)==5) {
            if (sharpen_type)
              print("Sharpen image%s with shock filters and amplitude %g, edge %g, alpha %g and sigma %g.",
                    gmic_inds,amplitude,edge,alpha,sigma);
            else
              print("Sharpen image%s with inverse diffusion and amplitude %g.",gmic_inds,amplitude);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],sharpen(amplitude,sharpen_type?true:false,edge,alpha,sigma));
          } else error("Sharpen image%s : Invalid argument '%s' "
                       "(should be 'amplitude[,sharpen_type[,edge[,alpha[,sigma]]]]', "
                       "where 'sharpen_type' can be '{0=inverse diffusion, 1=shock filters}').",
                       gmic_inds,argument_text);
          ++position; continue;
        }

        // Convolve.
        if (!cimg::strcmp("-convolve",item0)) {
          int ind0 = no_ind, borders = 1; char sep = 0, end = 0;
          if ((std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']') ||
              std::sscanf(argument,"[%d]%*c%d%c",&ind0,&borders,&end)==2) {
            gmic_check_indice(ind0,"Convolve image%s");
            print("Convolve image%s with mask [%d].",gmic_inds,ind0);
            const CImg<T> mask = images[ind0];
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],convolve(mask,borders));
          } else error("Convolve image%s : Invalid argument '%s' "
                       "(should be '[indice][,border_conditions]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Correlate.
        if (!cimg::strcmp("-correlate",item0)) {
          int ind0 = no_ind, borders = 1; char sep = 0, end = 0;
          if ((std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']') ||
              std::sscanf(argument,"[%d]%*c%d%c",&ind0,&borders,&end)==2) {
            gmic_check_indice(ind0,"Correlate image%s");
            print("Correlate image%s with mask [%d].",gmic_inds,ind0);
            const CImg<T> mask = images[ind0];
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],correlate(mask,borders));
          } else error("Correlate image%s : Invalid argument '%s' "
                       "(should be '[indice][,border_conditions]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Erode.
        if (!cimg::strcmp("-erode",item0)) {
          int siz = 3, ind0 = no_ind, borders = 1; char sep = 0, end = 0;
          if ((std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']') ||
              std::sscanf(argument,"[%d]%*c%d%c",&ind0,&borders,&end)==2) {
            gmic_check_indice(ind0,"Erode image%s");
            print("Erode image%s with mask [%d].",gmic_inds,ind0);
            const CImg<T> mask = images[ind0];
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],erode(mask,borders));
          } else if (std::sscanf(argument,"%d%c",&siz,&end)==1 ||
                     std::sscanf(argument,"%d%*c%d%c",&siz,&borders,&end)==2) {
            if (siz<=0) error("Erode image%s : Invalid size %d.",gmic_inds,siz);
            print("Erode image%s with size %d.",gmic_inds,siz);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],erode(siz,borders));
          } else error("Erode image%s : Invalid argument '%s' "
                       "(should be '[indice]' or 'size').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Dilate.
        if (!cimg::strcmp("-dilate",item0)) {
          int siz = 3, ind0 = no_ind, borders = 1; char sep = 0, end = 0;
          if ((std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']') ||
              std::sscanf(argument,"[%d]%*c%d%c",&ind0,&borders,&end)==2) {
            gmic_check_indice(ind0,"Dilate image%s");
            print("Dilate image%s with mask [%d].",gmic_inds,ind0);
            const CImg<T> mask = images[ind0];
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],dilate(mask,borders));
          } else if (std::sscanf(argument,"%d%c",&siz,&end)==1 ||
                     std::sscanf(argument,"%d%*c%d%c",&siz,&borders,&end)==2) {
            if (siz<=0) error("Dilate image%s : Invalid size %d.",gmic_inds,siz);
            print("Dilate image%s with size %d.",gmic_inds,siz);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],dilate(siz,borders));
          } else error("Dilate image%s : Invalid argument '%s' "
                       "(should be '[indice]' or 'size').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Compute gradient.
        if (!cimg::strcmp("-gradient",item0)) {
          char axes[4096] = { 0 }, *naxes = 0, end = 0; int scheme = 3;
          print("Compute gradient of image%s.",gmic_inds);
          if (std::sscanf(argument,"%4095[xyz]%c",axes,&end)==1 ||
              std::sscanf(argument,"%4095[xyz]%*c%d%c",axes,&scheme,&end)==2) { naxes = axes; ++position; }
          unsigned int off = 0;
          cimg_foroff(indices,l) {
            const unsigned int ind = indices[l] + off;
            CImg<T>& img = images[ind];
            const CImg<char> filename = filenames[ind];
            const CImgList<T> gradient = img.get_gradient(naxes,scheme);
            if (get_version) {
              images.insert(gradient);
              filenames.insert(gradient.size,filename);
            } else {
              images.remove(ind); images.insert(gradient,ind);
              filenames.remove(ind); filenames.insert(gradient.size,filename,ind);
              off+=gradient.size-1;
            }
          }
          continue;
        }

        // Compute Hessian.
        if (!cimg::strcmp("-hessian",item0)) {
          char axes[4096] = { 0 }, *naxes = 0, end = 0;
          print("Compute Hessian of image%s.",gmic_inds);
          if (std::sscanf(argument,"%4095[xyz]%c",axes,&end)==1) { naxes = axes; ++position; }
          unsigned int off = 0;
          cimg_foroff(indices,l) {
            const unsigned int ind = indices[l] + off;
            CImg<T>& img = images[ind];
            const CImg<char> filename = filenames[ind];
            const CImgList<T> hessian = img.get_hessian(naxes);
            if (get_version) {
              images.insert(hessian);
              filenames.insert(hessian.size,filename);
            } else {
              images.remove(ind); images.insert(hessian,ind);
              filenames.remove(ind); filenames.insert(hessian.size,filename,ind);
              off+=hessian.size-1;
            }
          }
          continue;
        }

        // Compute direct or inverse FFT.
        const bool inv_fft = !cimg::strcmp("-ifft",item0);
        if (!cimg::strcmp("-fft",item0) || inv_fft) {
          print("Compute %sFourier Transform of complex data",inv_fft?"inverse ":"");
          cimg_foroff(indices,l) {
            const unsigned int ind0 = indices[l], ind1 = l+1<_maxl?indices[l+1]:~0U;
            if (ind1!=~0U) {
              if (verbosity_level>=0) std::fprintf(cimg_stdout," ([%u],[%u])%c",ind0,ind1,l==_maxl-1?'.':',');
              CImgList<T> fft(images[ind0],images[ind1],!get_version);
              fft.FFT(inv_fft);
              if (get_version) {
                images.insert(2);
                fft[0].transfer_to(images[images.size-2]);
                fft[1].transfer_to(images[images.size-1]);
                filenames.insert(filenames[ind0]);
                filenames.insert(filenames[ind1]);
              } else {
                fft[0].transfer_to(images[ind0]);
                fft[1].transfer_to(images[ind1]);
              }
              ++l;
            } else {
              if (verbosity_level>=0) std::fprintf(cimg_stdout," ([%u],0)",ind0);
              CImgList<T> fft(images[ind0],!get_version);
              fft.insert(fft[0],~0U,false);
              fft[1].fill(0);
              fft.FFT(inv_fft);
              if (get_version) {
                images.insert(2);
                fft[0].transfer_to(images[images.size-2]);
                fft[1].transfer_to(images[images.size-1]);
                filenames.insert(2,filenames[ind0]);
              } else {
                fft[0].transfer_to(images[ind0]);
                images.insert(fft[1],1+ind0);
                filenames.insert(filenames[ind0],1+ind0);
              }
            }
          }
          continue;
        }

        //-----------------------------
        // Image creation and drawing
        //-----------------------------

        // Dimensions.
        if (!cimg::strcmp("-dimensions",item0)) {
          print("Get dimensions of image%s.",gmic_inds);
          cimg_foroff(indices,l) {
            CImg<T>& img = images[indices[l]];
            CImg<int> dims = CImg<int>::vector(img.dimx(),img.dimy(),img.dimz(),img.dimv());
            gmic_apply(img,replace(dims));
          }
          continue;
        }

        // Stats.
        if (!cimg::strcmp("-stats",item0)) {
          print("Get statistics of image%s.",gmic_inds);
          cimg_foroff(indices,l) gmic_apply(images[indices[l]],stats());
          continue;
        }

        // Histogram.
        if (!cimg::strcmp("-histogram",item0)) {
          int nb_levels = 256; char sep = 0, end = 0;
          if (std::sscanf(argument,"%d%c",&nb_levels,&end)==1 ||
              (std::sscanf(argument,"%d%c%c",&nb_levels,&sep,&end)==2 && sep=='%')) {
            print("Compute histogram of image%s using %d%s levels.",gmic_inds,nb_levels,sep=='%'?"%":"");
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              int nnb_levels = nb_levels;
              if (sep=='%') { double m, M = img.maxmin(m); nnb_levels = (int)cimg::round(nb_levels*(1+M-m)/100,1); }
              gmic_apply(images[indices[l]],histogram(nnb_levels));
            }
          } else error("Compute histogram of image%s : Invalid argument '%s' "
                       "(should be 'nb_levels[%%]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Distance function.
        if (!cimg::strcmp("-distance",item0)) {
          double value = 0; char sep = 0, end = 0;
          if (std::sscanf(argument,"%lf%c",&value,&end)==1 ||
              (std::sscanf(argument,"%lf%c%c",&value,&sep,&end)==2 && sep=='%')) {
            print("Compute distance map of image%s to isovalue %g%s.",gmic_inds,value,sep=='%'?"%":"");
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              double isovalue = value;
              if (sep=='%') { double m, M = img.maxmin(m); isovalue = m + value*(M - m)/100; }
              gmic_apply(img,distance((T)isovalue));
            }
          } else error("Compute distance function of image%s : Invalid argument '%s' "
                       "(should be 'value[%%]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Apply Hamilton-Jacobi PDE to compute distance to 0.
        if (!cimg::strcmp("-hamilton",item0)) {
          int nb_iter = 0; float band_size = 0; char end = 0;
          if (std::sscanf(argument,"%d%c",&nb_iter,&end)==1 ||
              std::sscanf(argument,"%d%*c%f%c",&nb_iter,&band_size,&end)==2) {
            print("Apply %d iterations of Hamilton-Jacobi PDE on image%s.",nb_iter,gmic_inds);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],distance_hamilton((unsigned int)nb_iter,band_size));
          } else error("Apply %d iterations of Hamilton-Jacobi PDE on image%s : Invalid argument '%s' "
                       "(should be 'nb_iter[,band_size]', with band_size>0).",nb_iter,gmic_inds,argument_text);
          ++position; continue;
        }

        // Label regions.
        gmic_simple_item("-label",label_regions,"Label regions on image%s.");

        // Displacement field.
        if (!cimg::strcmp("-displacement",item0)) {
          float smooth = 0.1f, precision = 0.1f; int ind0 = no_ind, nbscales = 0, itermax = 1000; char sep = 0, end = 0;
          if ((std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']') ||
              std::sscanf(argument,"[%d]%*c%f%c",&ind0,&smooth,&end)==2 ||
              std::sscanf(argument,"[%d]%*c%f%*c%f%c",&ind0,&smooth,&precision,&end)==3 ||
              std::sscanf(argument,"[%d]%*c%f%*c%f%*c%d%c",&ind0,&smooth,&precision,&nbscales,&end)==4 ||
              std::sscanf(argument,"[%d]%*c%f%*c%f%*c%d%*c%d%c",&ind0,&smooth,&precision,&nbscales,&itermax,&end)==5) {
            gmic_check_indice(ind0,"Compute displacement field of image%s");
            print("Compute displacement field of image%s with target [%u] and smoothness %g.",
                  gmic_inds,ind0,smooth);
            const CImg<T> target = images[ind0];
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],displacement_field(target,smooth,precision,nbscales,itermax));
          } else error("Compute displacement field of image%s : Invalid argument '%s' "
                       "(should be '[indice][,smoothness[,precision[,nbscales[,itermax]]]]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Sort.
        gmic_simple_item("-sort",sort,"Sort values in image%s.");

        // PSNR.
        if (!cimg::strcmp("-psnr",item0)) {
          double valmax = 255; char end = 0;
          if (std::sscanf(argument,"%lf%c",&valmax,&end)==1) ++position;
          if (images.size) {
            const unsigned int siz = indices.size();
            print("Compute a %ux%u matrix [%u] of PSNR values (max. pixel value is %g).",siz,siz,images.size,valmax);
            CImg<T> res(siz,siz,1,1,(T)-1);
            cimg_forXY(res,x,y) if (x>y) res(x,y) = res(y,x) = (T)images[indices[x]].PSNR(images[indices[y]],(float)valmax);
            images.insert(res);
            filenames.insert(CImg<char>("PSNR",5,1,1,1,false));
          } else error("Compute PSNR : image list is empty.");
          continue;
        }

        // Draw point.
        if (!cimg::strcmp("-point",item0)) {
          char arg0[4096] = { 0 }, arg1[4096] = { 0 }, arg2[4096] = { 0 }, color[4096] = { 0 };
          char sepx0 = 0, sepy0 = 0, sepz0 = 0, end = 0;
          float x0 = 0, y0 = 0, z0 = 0, opacity = 1;
          if (std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%f%*c%4095[0-9.eE,+-]",
                          arg0,arg1,arg2,&opacity,color)>=2 &&
              ((std::sscanf(arg0,"%f%c%c",&x0,&sepx0,&end)==2 && sepx0=='%') ||
               std::sscanf(arg0,"%f%c",&x0,&end)==1) &&
              ((std::sscanf(arg1,"%f%c%c",&y0,&sepy0,&end)==2 && sepy0=='%') ||
               std::sscanf(arg1,"%f%c",&y0,&end)==1) &&
              ((std::sscanf(arg2,"%f%c%c",&z0,&sepz0,&end)==2 && sepz0=='%') ||
               std::sscanf(arg2,"%f%c",&z0,&end)==1 || !arg2[0])) {
            print("Draw point (%g%s,%g%s,%g%s) with color '%s' and opacity %g on image%s.",
                  x0,sepx0=='%'?"%":"",y0,sepy0=='%'?"%":"",z0,sepz0=='%'?"%":"",
                  color[0]?color:"default",opacity,gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              CImg<T> col(img.dimv(),1,1,1,0);
              col.fill(color,true);
              const int
                nx0 = (int)cimg::round(sepx0=='%'?x0*(img.dimx()-1)/100:x0,1),
                ny0 = (int)cimg::round(sepy0=='%'?y0*(img.dimy()-1)/100:y0,1),
                nz0 = (int)cimg::round(sepz0=='%'?z0*(img.dimz()-1)/100:z0,1);
              gmic_apply(img,draw_point(nx0,ny0,nz0,col,opacity));
            }
          } else error("Draw point on image%s : Invalid argument '%s' "
                       "(should be 'x[%%],y[%%][,z[%%][,opacity[,color]]])",gmic_inds,argument_text);
          ++position; continue;
        }

        // Draw line.
        if (!cimg::strcmp("-line",item0)) {
          char arg0[4096] = { 0 }, arg1[4096] = { 0 }, arg2[4096] = { 0 }, arg3[4096] = { 0 }, color[4096] = { 0 };
          char sepx0 = 0, sepy0 = 0, sepx1 = 0, sepy1 = 0, end = 0;
          float x0 = 0, y0 = 0, x1 = 0, y1 = 0, opacity = 1;
          if (std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]"
                          "%*c%f%*c%4095[0-9.eE,+-]",
                          arg0,arg1,arg2,arg3,&opacity,color)>=4 &&
              ((std::sscanf(arg0,"%f%c%c",&x0,&sepx0,&end)==2 && sepx0=='%') ||
               std::sscanf(arg0,"%f%c",&x0,&end)==1) &&
              ((std::sscanf(arg1,"%f%c%c",&y0,&sepy0,&end)==2 && sepy0=='%') ||
               std::sscanf(arg1,"%f%c",&y0,&end)==1) &&
              ((std::sscanf(arg2,"%f%c%c",&x1,&sepx1,&end)==2 && sepx1=='%') ||
               std::sscanf(arg2,"%f%c",&x1,&end)==1) &&
              ((std::sscanf(arg3,"%f%c%c",&y1,&sepy1,&end)==2 && sepy1=='%') ||
               std::sscanf(arg3,"%f%c",&y1,&end)==1)) {
            print("Draw line (%g%s,%g%s) - (%g%s,%g%s) with color '%s' and opacity %g on image%s.",
                  x0,sepx0=='%'?"%":"",y0,sepy0=='%'?"%":"",x1,sepx1=='%'?"%":"",y1,sepy1=='%'?"%":"",
                  color[0]?color:"default",opacity,gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              CImg<T> col(img.dimv(),1,1,1,0);
              col.fill(color,true);
              const int
                nx0 = (int)cimg::round(sepx0=='%'?x0*(img.dimx()-1)/100:x0,1),
                ny0 = (int)cimg::round(sepy0=='%'?y0*(img.dimy()-1)/100:y0,1),
                nx1 = (int)cimg::round(sepx1=='%'?x1*(img.dimx()-1)/100:x1,1),
                ny1 = (int)cimg::round(sepy1=='%'?y1*(img.dimy()-1)/100:y1,1);
              gmic_apply(img,draw_line(nx0,ny0,nx1,ny1,col,opacity));
            }
          } else error("Draw line on image%s : Invalid argument '%s' "
                       "(should be 'x0[%%],y0[%%],x1[%%],y1[%%][,opacity[,color]]')",gmic_inds,argument_text);
          ++position; continue;
        }

        // Draw polygon.
        if (!cimg::strcmp("-polygon",item0)) {
          char arg0[4096] = { 0 }, arg1[4096] = { 0 }, tmp[4096] = { 0 }, sepx0 = 0, sepy0 = 0, end = 0;
          int N = 0; float x0 = 0, y0 = 0, opacity = 1;
          if (std::sscanf(argument,"%d%c",&N,&end)==2 && N>2) {
            const char
              *nargument = argument + std::sprintf(tmp,"%d",N) + 1,
              *const eargument = argument + cimg::strlen(argument);
            CImg<float> coords0(N,2,1,1,0);
            CImg<bool> percents(N,2,1,1,0);
            for (int n = 0; n<N; ++n) if (nargument<eargument) {
              if (std::sscanf(nargument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]",arg0,arg1)==2 &&
                  ((std::sscanf(arg0,"%f%c%c",&x0,&(sepx0=0),&end)==2 && sepx0=='%') ||
                   std::sscanf(arg0,"%f%c",&x0,&end)==1) &&
                  ((std::sscanf(arg1,"%f%c%c",&y0,&(sepy0=0),&end)==2 && sepy0=='%') ||
                   std::sscanf(arg1,"%f%c",&y0,&end)==1)) {
                coords0(n,0) = x0; percents(n,0) = (sepx0=='%');
                coords0(n,1) = y0; percents(n,1) = (sepy0=='%');
                nargument+=cimg::strlen(arg0) + cimg::strlen(arg1) + 2;
              } else error("Draw polygon on image%s : Invalid or incomplete argument '%s' "
                           "(should be 'N,x0[%%],y0[%%],x1[%%],y1[%%],..,xN[%%],yN[%%][,opacity[,color]]' with N>=3)",
                           gmic_inds,argument_text);
            } else error("Draw polygon on image%s : Incomplete argument '%s' "
                         "(%d xy-coordinates should be defined)",
                         gmic_inds,argument_text,N);
            if (nargument<eargument && std::sscanf(nargument,"%4095[0-9.eE+-]",arg0)==1 &&
                std::sscanf(arg0,"%f",&opacity)==1) nargument+=cimg::strlen(arg0)+1;
            const char *const color = nargument<eargument?nargument:&(end=0);
            print("Draw %d-vertices polygon with color '%s' and opacity %g on image%s.",
                  N,color[0]?color:"default",opacity,gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              CImg<int> coords(coords0);
              cimg_forX(coords,p) {
                if (percents(p,0)) coords(p,0) = (int)cimg::round(coords0(p,0)*(img.dimx()-1)/100,1);
                if (percents(p,1)) coords(p,1) = (int)cimg::round(coords0(p,1)*(img.dimy()-1)/100,1);
              }
              CImg<T> col(img.dimv(),1,1,1,0);
              col.fill(color,true);
              gmic_apply(img,draw_polygon(coords,col,opacity));
            }
          } else error("Draw polygon on image%s : Invalid argument '%s' "
                       "(should be 'N,x0[%%],y0[%%],x1[%%],y1[%%],..,xN[%%],yN[%%][,opacity[,color]]' with N>=3)",
                       gmic_inds,argument_text);
          ++position; continue;
        }

        // Draw ellipse.
        if (!cimg::strcmp("-ellipse",item0)) {
          char arg0[4096] = { 0 }, arg1[4096] = { 0 }, color[4096] = { 0 };
          char sepx0 = 0, sepy0 = 0, end = 0;
          float x0 = 0, y0 = 0, r0 = 0, r1 = 0, ru = 1, rv = 0, opacity = 1;
          if (std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%f%*c%f%*c%f%*c%f%*c%f%*c%4095[0-9.eE,+-]",
                          arg0,arg1,&r0,&r1,&ru,&rv,&opacity,color)>=4 &&
              ((std::sscanf(arg0,"%f%c%c",&x0,&sepx0,&end)==2 && sepx0=='%') ||
               std::sscanf(arg0,"%f%c",&x0,&end)==1) &&
              ((std::sscanf(arg1,"%f%c%c",&y0,&sepy0,&end)==2 && sepy0=='%') ||
               std::sscanf(arg1,"%f%c",&y0,&end)==1)) {
            print("Draw ellipse centered at (%g%s,%g%s) with radii (%g,%g), orientation (%g,%g), color '%s' "
                  "and opacity %g on image%s.",
                  x0,sepx0=='%'?"%":"",y0,sepy0=='%'?"%":"",
                  r0,r1,ru,rv,color[0]?color:"default",opacity,gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              CImg<T> col(img.dimv(),1,1,1,0);
              col.fill(color,true);
              const int
                nx0 = (int)cimg::round(sepx0=='%'?x0*(img.dimx()-1)/100:x0,1),
                ny0 = (int)cimg::round(sepy0=='%'?y0*(img.dimy()-1)/100:y0,1);
              gmic_apply(img,draw_ellipse(nx0,ny0,r0,r1,ru,rv,col,opacity));
            }
          } else error("Draw ellipse on image%s : Invalid argument '%s' "
                       "(should be 'x[%%],y[%%],r,R[,u,v[,opacity[,color]]])",
                       gmic_inds,argument_text);
          ++position; continue;
        }

        // Draw text.
        if (!cimg::strcmp("-text",item0)) {
          char arg0[4096] = { 0 }, arg1[4096] = { 0 }, color[4096] = { 0 }, text[4096] = { 0 };
          char sepx0 = 0, sepy0 = 0, end = 0;
          float x0 = 0, y0 = 0, opacity = 1; int siz = 11;
          if (std::sscanf(argument,"%4095[^,],%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%d%*c%f%*c%4095[0-9.eE,+-]",
                          text,arg0,arg1,&siz,&opacity,color)>=1 &&
              ((std::sscanf(arg0,"%f%c%c",&x0,&sepx0,&end)==2 && sepx0=='%') ||
               std::sscanf(arg0,"%f%c",&x0,&end)==1 || !arg0[0]) &&
              ((std::sscanf(arg1,"%f%c%c",&y0,&sepy0,&end)==2 && sepy0=='%') ||
               std::sscanf(arg1,"%f%c",&y0,&end)==1 || !arg1[0])) {
            cimg::strclean(text); cimg::strescape(text);
            print("Draw text \"%s\" at position (%g%s,%g%s) with font size %d, color '%s' "
                  "and opacity %f on image%s.",
                  text,x0,sepx0=='%'?"%":"",y0,sepy0=='%'?"%":"",siz,color[0]?color:"default",opacity,gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              CImg<T> col(img.dimv(),1,1,1,0);
              col.fill(color,true);
              const int
                nx0 = (int)cimg::round(sepx0=='%'?x0*(img.dimx()-1)/100:x0,1),
                ny0 = (int)cimg::round(sepy0=='%'?y0*(img.dimy()-1)/100:y0,1);
              gmic_apply(img,draw_text(nx0,ny0,text,col.ptr(),0,opacity,siz));
            }
          } else error("Draw text on image%s : Invalid argument '%s' "
                       "(should be 'text[,x[%%],y[%%][,size[,opacity[,color]]]]').",
                       gmic_inds,argument_text);
          ++position; continue;
        }

        // Draw image.
        if (!cimg::strcmp("-image",item0)) {
          char arg0[4096] = { 0 }, arg1[4096] = { 0 }, arg2[4096] = { 0 }, sep = 0, sepx = 0, sepy = 0, sepz = 0, end = 0;
          int ind0 = no_ind, indm0 = no_ind; float x = 0, y = 0, z = 0, opacity = 1;
          if (((std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==1 && sep==']') ||
               std::sscanf(argument,"[%d]%*c%4095[0-9.eE%+-]%c",&ind0,arg0,&end)==2 ||
               std::sscanf(argument,"[%d]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",&ind0,arg0,arg1,&end)==3 ||
               std::sscanf(argument,"[%d]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",&ind0,arg0,arg1,arg2,&end)==4 ||
               std::sscanf(argument,"[%d]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%f%c",&ind0,arg0,arg1,arg2,&opacity,&end)==5 ||
               std::sscanf(argument,"[%d]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%f%*c[%d%c%c",
                            &ind0,arg0,arg1,arg2,&opacity,&indm0,&sep,&end)==7) &&
              (!*arg0 ||
               std::sscanf(arg0,"%f%c",&x,&end)==1 ||
               (std::sscanf(arg0,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (!*arg1 ||
               std::sscanf(arg1,"%f%c",&y,&end)==1 ||
               (std::sscanf(arg1,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%')) &&
              (!*arg2 ||
               std::sscanf(arg2,"%f%c",&z,&end)==1 ||
               (std::sscanf(arg2,"%f%c%c",&z,&sepz,&end)==2 && sepz=='%'))) {
            gmic_check_indice(ind0,"Draw image on image%s");
            const CImg<T> sprite = images[ind0];
            CImg<T> mask;
            if (indm0!=no_ind) {
              gmic_check_indice(indm0,"Draw image on image%s");
              mask = images[indm0];
              print("Draw image [%d] at (%g%s,%g%s,%g%s), with mask [%d] and opacity %f on image%s.",
                    ind0,x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",z,sepz=='%'?"%":"",indm0,opacity,gmic_inds);
            } else print("Draw image [%d] at (%g%s,%g%s,%g%s) with opacity %f on image%s.",
                         ind0,x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",z,sepz=='%'?"%":"",opacity,gmic_inds);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.dimx()-1)/100:x,1),
                ny = (int)cimg::round(sepy=='%'?y*(img.dimy()-1)/100:y,1),
                nz = (int)cimg::round(sepz=='%'?z*(img.dimz()-1)/100:z,1);
              if (indm0!=no_ind) { gmic_apply(img,draw_image(nx,ny,nz,sprite,mask,opacity)); }
              else { gmic_apply(img,draw_image(nx,ny,nz,sprite,opacity)); }
            }
          } else error("Draw image on image%s : Invalid argument '%s' "
                       "(should be '[indice][,x[%%][,y[%%][,z[%%][,opacity[,indice_mask]]]]]').",
                       gmic_inds,argument_text);
          ++position; continue;
        }

        // Draw 3D object.
        if (!cimg::strcmp("-object3d",item0)) {
          char arg0[4096] = { 0 }, arg1[4096] = { 0 }, sep = 0, sepx = 0, sepy = 0, end = 0;
          float x = 0, y = 0, z = 0, opacity = 1;
          int ind0 = no_ind;
          if (((std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']') ||
               std::sscanf(argument,"[%d]%*c%4095[0-9.eE%+-]%c",&ind0,arg0,&end)==2 ||
               std::sscanf(argument,"[%d]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%c",&ind0,arg0,arg1,&end)==3 ||
               std::sscanf(argument,"[%d]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%f%c",&ind0,arg0,arg1,&z,&end)==4 ||
               std::sscanf(argument,"[%d]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%f%*c%f%c",&ind0,arg0,arg1,&z,&opacity,&end)==5) &&
              (!*arg0 ||
               std::sscanf(arg0,"%f%c",&x,&end)==1 ||
               (std::sscanf(arg0,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (!*arg1 ||
               std::sscanf(arg1,"%f%c",&y,&end)==1 ||
               (std::sscanf(arg1,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%'))) {
            gmic_check_indice(ind0,"Draw 3D object on image%s");
            if (!images[ind0].is_CImg3d())
              error("Draw 3D object on image%s : Image [%d] is not a 3D object.",gmic_inds,ind0);
            print("Draw 3D object [%d] at (%g%s,%g%s,%g) on image%s, with opacity %g.",
                  ind0,x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",z,gmic_inds,opacity);
            CImgList<unsigned int> primitives3d;
            CImgList<unsigned char> colors3d;
            CImg<float> opacities3d, points3d(images[ind0]);
            points3d.CImg3dtoobject3d(primitives3d,colors3d,opacities3d);
            opacities3d*=opacity;
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              const float
                nx = (float)cimg::round(sepx=='%'?x*(img.dimx()-1)/100:x,1),
                ny = (float)cimg::round(sepy=='%'?y*(img.dimy()-1)/100:y,1);
              gmic_apply(img,draw_object3d(nx,ny,z,points3d,primitives3d,colors3d,opacities3d,
                                           render3d,!is_oriented3d,focale3d,light3d_x,light3d_y,light3d_z,specular_light3d,
                                           specular_shine3d,0));
            }
          } else error("Draw 3D object on image%s : Invalid argument '%s' "
                       "(should be '[indice][,x[%%][,y[%%][,z[,opacity[,zoom[,u1,v1,w1,angle1[,...]]]]]]]').",
                       gmic_inds,argument_text);
          ++position; continue;
        }

        // Draw plasma fractal.
        if (!cimg::strcmp("-plasma",item0)) {
          float alpha = 1, beta = 1, opacity = 1; char end = 0;
          if (std::sscanf(argument,"%f%c",&alpha,&end)==1 ||
              std::sscanf(argument,"%f%*c%f%c",&alpha,&beta,&end)==2 ||
              std::sscanf(argument,"%f%*c%f%*c%f%c",&alpha,&beta,&opacity,&end)==3) {
            print("Draw plasma in image%s with alpha %g, beta %g and opacity %g.",gmic_inds,alpha,beta,opacity);
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],draw_plasma(alpha,beta,opacity));
          } else error("Draw plasma in image%d : Invalid argument '%s' "
                       "(should be 'alpha[,beta[,opacity]]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Draw Mandelbrot/Julia fractal.
        if (!cimg::strcmp("-mandelbrot",item0)) {
          double z0r = -2, z0i = -2, z1r = 2, z1i = 2, paramr = 0, parami = 0; char end = 0;
          float opacity = 1; int itermax = 100, julia = 0;
          if (std::sscanf(argument,"%lf%*c%lf%*c%lf%*c%lf%c",&z0r,&z0i,&z1r,&z1i,&end)==4 ||
              std::sscanf(argument,"%lf%*c%lf%*c%lf%*c%lf%*c%d%c",&z0r,&z0i,&z1r,&z1i,&itermax,&end)==5 ||
              std::sscanf(argument,"%lf%*c%lf%*c%lf%*c%lf%*c%d%*c%d%*c%lf%*c%lf%c",
                          &z0r,&z0i,&z1r,&z1i,&itermax,&julia,&paramr,&parami,&end)==8 ||
              std::sscanf(argument,"%lf%*c%lf%*c%lf%*c%lf%*c%d%*c%d%*c%lf%*c%lf%*c%f%c",
                          &z0r,&z0i,&z1r,&z1i,&itermax,&julia,&paramr,&parami,&opacity,&end)==9) {
            print("Draw %s fractal in image%s from complex area (%g,%g)-(%g,%g) with c0 = (%g,%g) (%d iterations).",
                  julia?"Julia":"Mandelbrot",gmic_inds,z0r,z0i,z1r,z1i,paramr,parami,itermax);
            cimg_foroff(indices,l)
              gmic_apply(images[indices[l]],draw_mandelbrot(CImg<T>(),opacity,z0r,z0i,z1r,z1i,itermax,true,
                                                            julia?true:false,paramr,parami));
          } else error("Draw fractal in image%s : Invalid argument '%s' "
                       "(should be 'z0r,z0i,z1r,z1i[,itermax[,julia,c0r,c0i[,opacity]]]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Flood fill.
        if (!cimg::strcmp("-flood",item0)) {
          char arg0[4096] = { 0 }, arg1[4096] = { 0 }, arg2[4096] = { 0 }, color[4096] = { 0 };
          char sepx = 0, sepy = 0, sepz = 0, end = 0;
          float x = 0, y = 0, z = 0, tolerance = 0, opacity = 1;
          if (std::sscanf(argument,"%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%4095[0-9.eE%+-]%*c%f%*c%f%*c%4095[0-9.eE,+-]",
                          arg0,arg1,arg2,&tolerance,&opacity,color)>=1 &&
              ((std::sscanf(arg0,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%') ||
               std::sscanf(arg0,"%f%c",&x,&end)==1) &&
              ((std::sscanf(arg1,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%') ||
               std::sscanf(arg1,"%f%c",&y,&end)==1 || !arg1[0]) &&
              ((std::sscanf(arg2,"%f%c%c",&z,&sepz,&end)==2 && sepz=='%') ||
               std::sscanf(arg2,"%f%c",&z,&end)==1 || !arg2[0])) {
            print("Flood fill image%s from (%g%s,%g%s,%g%s) with tolerance %g, opacity %g and color '%s'.",
                  gmic_inds,x,sepx=='%'?"%":"",y,sepy=='%'?"%":"",z,sepz=='%'?"%":"",tolerance,opacity,color);
            cimg_foroff(indices,l) {
              CImg<T> &img = images[indices[l]];
              CImg<T> col(img.dimv(),1,1,1,0);
              col.fill(color,true);
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.dimx()-1)/100:x,1),
                ny = (int)cimg::round(sepy=='%'?y*(img.dimy()-1)/100:y,1),
                nz = (int)cimg::round(sepz=='%'?z*(img.dimz()-1)/100:z,1);
              gmic_apply(img,draw_fill(nx,ny,nz,col,opacity,tolerance));
            }
          } else error("Flood fill image%s : Invalid argument '%s' "
                       "(should be 'x[,y[,z[,tolerance[,opacity[,color]]]]]').",gmic_inds,argument_text);
          ++position; continue;
        }

        //-------------------------
        // Image list manipulation
        //-------------------------

        // Remove specified image(s).
        if (!cimg::strcmp("-remove",item0) || !cimg::strcmp("-rm",item0)) {
          print("Remove image%s",gmic_inds);
          unsigned int off = 0;
          cimg_foroff(indices,l) {
            const unsigned int ind = indices[l] - off;
            images.remove(ind); filenames.remove(ind);
            ++off;
          }
          if (verbosity_level>=0) std::fprintf(cimg_stdout," (%u image%s left).",images.size,images.size==1?"":"s");
          continue;
        }

        // Keep specified image(s).
        if (!cimg::strcmp("-keep",item0) || !cimg::strcmp("-k",item0)) {
          print("Keep image%s",gmic_inds);
          CImgList<T> nimages(indices.size());
          cimg_foroff(indices,l) nimages[l].swap(images[indices[l]]);
          nimages.transfer_to(images);
          if (verbosity_level>=0) std::fprintf(cimg_stdout," (%u image%s left).",images.size,images.size==1?"":"s");
          continue;
        }

        // Move image(s) to specified position.
        if (!cimg::strcmp("-move",item0) || !cimg::strcmp("-mv",item0)) {
          int ind0 = no_ind; char end = 0;
          if (std::sscanf(argument,"%d%c",&ind0,&end)==1) {
            if (ind0<0) ind0+=images.size;
            if (ind0<0) ind0 = 0;
            if (ind0>(int)images.size) ind0 = images.size;
            print("Move image%s to position %d.",gmic_inds,ind0);
            CImgList<T> nimages;
            CImgList<char> nfilenames;
            cimg_foroff(indices,l) {
              const unsigned int ind = indices[l];
              nimages.insert(1); nimages.last().swap(images[ind]);
              nfilenames.insert(1); nfilenames.last().swap(filenames[ind]);
            }
            images.insert(nimages,ind0); filenames.insert(nfilenames,ind0);
            { cimglist_for(images,l) if (!images[l]) { images.remove(l); filenames.remove(l--); }}
          } else error("Move image%s : Invalid argument '%s' "
                       "(should be 'position').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Reverse images order.
        if (!cimg::strcmp("-reverse",item0)) {
          print("Reverse images order.");
          CImgList<T> nimages(indices.size());
          CImgList<char> nfilenames(indices.size());
          cimg_foroff(indices,l) { nimages[l].swap(images[indices[l]]); nfilenames[l].swap(filenames[indices[l]]); }
          nimages.reverse(); nfilenames.reverse();
          { cimg_foroff(indices,l) { nimages[l].swap(images[indices[l]]); nfilenames[l].swap(filenames[indices[l]]); }}
          continue;
        }

        // Set image name.
        if (!cimg::strcmp("-name",item0)) {
          cimg_foroff(indices,l) filenames[indices[l]].assign(argument,cimg::strlen(argument)+1,1,1,1,false);
          ++position; continue;
        }

        //-------------------------
        // 3D objects manipulation
        //-------------------------

        // Generate 3D cube.
        if (!cimg::strcmp("-cube3d",item)) {
          float size = 100; char end = 0;
          if (std::sscanf(argument,"%f%c",&size,&end)==1) {
            print("Generate 3D cube with size %g.",size);
            CImgList<unsigned int> primitives3d;
            CImg<float> points3d = CImg<T>::cube3d(primitives3d,size);
            CImgList<unsigned char> colors3d(primitives3d.size,1,3,1,1,200);
            CImg<float> opacities3d(1,primitives3d.size,1,1,1);
            points3d.object3dtoCImg3d(primitives3d,colors3d,opacities3d);
            images.insert(points3d);
            filenames.insert(CImg<char>("(gmic)",7,1,1,1,false));
          } else error("Generate 3D cube : Invalid argument '%s' "
                       "(should be 'size').",argument_text);
          ++position; continue;
        }

        // Generate 3D cone.
        if (!cimg::strcmp("-cone3d",item)) {
          float radius = 100, height = 200; char end = 0; unsigned int subdivisions = 24;
          if (std::sscanf(argument,"%f%c",&radius,&end)==1 ||
              std::sscanf(argument,"%f%*c%f%c",&radius,&height,&end)==2 ||
              std::sscanf(argument,"%f%*c%f%*c%u%c",&radius,&height,&subdivisions,&end)==3) {
            print("Generate 3D cone with radius %g, height %g and %u subdivisions.",radius,height,subdivisions);
            CImgList<unsigned int> primitives3d;
            CImg<float> points3d = CImg<T>::cone3d(primitives3d,radius,height,subdivisions);
            CImgList<unsigned char> colors3d(primitives3d.size,1,3,1,1,200);
            CImg<float> opacities3d(1,primitives3d.size,1,1,1);
            points3d.object3dtoCImg3d(primitives3d,colors3d,opacities3d);
            images.insert(points3d);
            filenames.insert(CImg<char>("(gmic)",7,1,1,1,false));
          } else error("Generate 3D cone : Invalid argument '%s' "
                       "(should be 'radius[,height[,subdivisions]]').",argument_text);
          ++position; continue;
        }

        // Generate 3D cylinder.
        if (!cimg::strcmp("-cylinder3d",item)) {
          float radius = 100, height = 200; char end = 0; unsigned int subdivisions = 24;
          if (std::sscanf(argument,"%f%c",&radius,&end)==1 ||
              std::sscanf(argument,"%f%*c%f%c",&radius,&height,&end)==2 ||
              std::sscanf(argument,"%f%*c%f%*c%u%c",&radius,&height,&subdivisions,&end)==3) {
            print("Generate 3D cylinder with radius %g, height %g and %u subdivisions.",radius,height,subdivisions);
            CImgList<unsigned int> primitives3d;
            CImg<float> points3d = CImg<T>::cylinder3d(primitives3d,radius,height,subdivisions);
            CImgList<unsigned char> colors3d(primitives3d.size,1,3,1,1,200);
            CImg<float> opacities3d(1,primitives3d.size,1,1,1);
            points3d.object3dtoCImg3d(primitives3d,colors3d,opacities3d);
            images.insert(points3d);
            filenames.insert(CImg<char>("(gmic)",7,1,1,1,false));
          } else error("Generate 3D cylinder : Invalid argument '%s' "
                       "(should be 'radius[,height[,subdivisions]]').",argument_text);
          ++position; continue;
        }

        // Generate 3D torus.
        if (!cimg::strcmp("-torus3d",item)) {
          float radius1 = 100, radius2 = 30; char end = 0; unsigned int subdivisions1 = 24, subdivisions2 = 12;
          if (std::sscanf(argument,"%f%*c%f%c",&radius1,&radius2,&end)==2 ||
              std::sscanf(argument,"%f%*c%f%*c%u%*c%u%c",&radius1,&radius2,&subdivisions1,&subdivisions2,&end)==4) {
            print("Generate 3D torus with radii %g and %g, and subdivisions %u and %u.",radius1,radius2,subdivisions1,subdivisions2);
            CImgList<unsigned int> primitives3d;
            CImg<float> points3d = CImg<T>::torus3d(primitives3d,radius1,radius2,subdivisions1,subdivisions2);
            CImgList<unsigned char> colors3d(primitives3d.size,1,3,1,1,200);
            CImg<float> opacities3d(1,primitives3d.size,1,1,1);
            points3d.object3dtoCImg3d(primitives3d,colors3d,opacities3d);
            images.insert(points3d);
            filenames.insert(CImg<char>("(gmic)",7,1,1,1,false));
          } else error("Generate 3D torus : Invalid argument '%s' "
                       "(should be 'radius1,radius2[,subdivisions1,subdivisions2]').",argument_text);
          ++position; continue;
        }

        // Generate 3D plane.
        if (!cimg::strcmp("-plane3d",item)) {
          float sizex = 100, sizey = 30; char end = 0; unsigned int subdivisionsx = 24, subdivisionsy = 12;
          if (std::sscanf(argument,"%f%*c%f%c",&sizex,&sizey,&end)==2 ||
              std::sscanf(argument,"%f%*c%f%*c%u%*c%u%c",&sizex,&sizey,&subdivisionsx,&subdivisionsy,&end)==4) {
            print("Generate 3D plane with dimensions %g and %g, and subdivisions %u and %u.",sizex,sizey,subdivisionsx,subdivisionsy);
            CImgList<unsigned int> primitives3d;
            CImg<float> points3d = CImg<T>::plane3d(primitives3d,sizex,sizey,subdivisionsx,subdivisionsy);
            CImgList<unsigned char> colors3d(primitives3d.size,1,3,1,1,200);
            CImg<float> opacities3d(1,primitives3d.size,1,1,1);
            points3d.object3dtoCImg3d(primitives3d,colors3d,opacities3d);
            images.insert(points3d);
            filenames.insert(CImg<char>("(gmic)",7,1,1,1,false));
          } else error("Generate 3D plane : Invalid argument '%s' "
                       "(should be 'sizex,sizey[,subdivisionsx,subdivisionsy]').",argument_text);
          ++position; continue;
        }

        // Generate 3D sphere.
        if (!cimg::strcmp("-sphere3d",item)) {
          float radius = 100; char end = 0; unsigned int subdivisions = 3;
          if (std::sscanf(argument,"%f%c",&radius,&end)==1 ||
              std::sscanf(argument,"%f%*c%u%c",&radius,&subdivisions,&end)==2) {
            print("Generate 3D sphere with radius %g and %u subdivisions.",radius,subdivisions);
            CImgList<unsigned int> primitives3d;
            CImg<float> points3d = CImg<T>::sphere3d(primitives3d,radius,subdivisions);
            CImgList<unsigned char> colors3d(primitives3d.size,1,3,1,1,200);
            CImg<float> opacities3d(1,primitives3d.size,1,1,1);
            points3d.object3dtoCImg3d(primitives3d,colors3d,opacities3d);
            images.insert(points3d);
            filenames.insert(CImg<char>("(gmic)",7,1,1,1,false));
          } else error("Generate 3D sphere : Invalid argument '%s' "
                       "(should be 'radius[,subdivisions]').",argument_text);
          ++position; continue;
        }

        // Build 3D elevation.
        if (!cimg::strcmp("-elevation3d",item0)) {
          float zfact = 0.2f; char end = 0, sep = 0; int ind0 = no_ind;
          if (std::sscanf(argument,"%f%c",&zfact,&end)==1 ||
              (std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']')) {
            CImg<typename CImg<T>::Tfloat> elev;
            if (ind0!=no_ind) {
              gmic_check_indice(ind0,"Build 3D elevation of image%s");
              print("Build 3D elevation of image%s with elevation map [%d].",gmic_inds,ind0);
              if (images[ind0].dimv()==1) elev = images[ind0];
              else elev = images[ind0].get_pointwise_norm();
            } else print("Build 3D elevation of image%s with z-factor %g.",gmic_inds,zfact);
            cimg_foroff(indices,l) {
              CImg<T>& img = images[indices[l]];
              CImgList<unsigned int> primitives3d;
              CImgList<unsigned char> colors3d;
              CImg<float> opacities3d, points3d;
              if (elev) points3d = img.get_elevation3d(primitives3d,colors3d,elev);
              else {
                if (img.dimv()==1) (elev = img)*=zfact; else (elev = img.get_pointwise_norm())*=zfact;
                points3d = img.get_elevation3d(primitives3d,colors3d,elev);
                elev.assign();
              }
              opacities3d.assign(1,primitives3d.size,1,1,1);
              points3d.object3dtoCImg3d(primitives3d,colors3d,opacities3d);
              gmic_apply(img,replace(points3d));
            }
          } else error("Build 3D elevation : invalid argument '%s' "
                       "(should be 'z-factor' or '[indice]').",argument_text);
          ++position; continue;
        }

        // Build 3D isovalue.
        if (!cimg::strcmp("-isovalue3d",item0)) {
          float value = 0; char end = 0;
          if (std::sscanf(argument,"%f%c",&value,&end)==1) {
            print("Build 3D isovalue %g of image%s.",value,gmic_inds);
            cimg_foroff(indices,l) {
              const unsigned int ind = indices[l];
              CImg<T>& img = images[ind];
              CImg<float> points3d;
              CImgList<unsigned int> primitives3d;
              CImgList<unsigned char> colors3d;
              CImg<float> opacities3d;
              CImg<unsigned char> palette;
              palette.assign(3,img.dim,1,1,220).noise(35,1);
              if (img.dim==1) palette(0) = palette(1) = palette(2) = 255;
              else {
                palette(0,0) = 255; palette(1,0) = 30; palette(2,0) = 30;
                palette(0,1) = 30; palette(1,1) = 255; palette(2,1) = 30;
                if (img.dim>=3) palette(0,2) = 30; palette(1,2) = 30; palette(2,2) = 255;
              }
              cimg_forV(img,k) {
                CImgList<unsigned int> prims;
                const CImg<float> pts = img.get_shared_channel(k).get_isovalue3d(prims,value);
                if (pts) {
                  points3d.append_object3d(primitives3d,pts,prims);
                  colors3d.insert(prims.size,
                                  CImg<unsigned char>::vector(palette(0,k),palette(1,k),palette(2,k)));
                }
              }
              opacities3d.assign(1,primitives3d.size,1,1,1);
              if (!points3d)
                warning("Build 3D isovalue of image [%u] : Isovalue %g not found.",ind,value);
              else points3d.object3dtoCImg3d(primitives3d,colors3d,opacities3d);
              gmic_apply(img,replace(points3d));
            }
          } else error("Build 3D isovalue of image%s : Invalid argument '%s' "
                       "(should be 'isovalue').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Center a 3D object.
        if (!cimg::strcmp("-center3d",item0) || !cimg::strcmp("-c3d",item0)) {
          print("Center 3D object%s.",gmic_inds);
          cimg_foroff(indices,l) {
            const unsigned int ind = indices[l];
            if (!images[ind].is_CImg3d())
              error("Center 3D object%s : Image [%d] is not a 3D object.",gmic_inds,ind);
            gmic_apply(images[ind],centerCImg3d());
          }
          continue;
        }

        // Normalize a 3D object.
        if (!cimg::strcmp("-normalize3d",item0) || !cimg::strcmp("-n3d",item0)) {
          print("Normalize 3D object%s.",gmic_inds);
          cimg_foroff(indices,l) {
            const unsigned int ind = indices[l];
            if (!images[ind].is_CImg3d())
              error("Normalize 3D object%s : Image [%d] is not a 3D object.",gmic_inds,ind);
            gmic_apply(images[ind],normalizeCImg3d());
          }
          continue;
        }

        // Rotate a 3D object.
        if (!cimg::strcmp("-rotate3d",item0) || !cimg::strcmp("-rot3d",item0)) {
          float u = 0, v = 0, w = 1, angle = 0; char end = 0;
          if (std::sscanf(argument,"%f%*c%f%*c%f%*c%f%c",&u,&v,&w,&angle,&end)==4) {
            print("Rotate 3D object%s around axis (%g,%g,%g) with angle %g.",gmic_inds,u,v,w,angle);
            const CImg<float> rot = CImg<float>::rotation_matrix(u,v,w,(float)(angle*cimg::valuePI/180));
            cimg_foroff(indices,l) {
              const unsigned int ind = indices[l];
              if (!images[ind].is_CImg3d())
                error("Rotate 3D object%s : Image [%d] is not a 3D object.",gmic_inds,ind);
              gmic_apply(images[ind],rotateCImg3d(rot));
            }
          } else error("Rotate 3D object%s : Invalid argument '%s' "
                       "(should be 'u,v,w,angle').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Add 3D objects together or translate a 3D object.
        if (!cimg::strcmp("-add3d",item0) || !cimg::strcmp("-+3d",item0)) {
          float tx = 0, ty = 0, tz = 0; int ind0 = no_ind; char sep = 0, end = 0;
          if (std::sscanf(argument,"%f%c",&tx,&end)==1 ||
              std::sscanf(argument,"%f%*c%f%c",&tx,&ty,&end)==2 ||
              std::sscanf(argument,"%f%*c%f%*c%f%c",&tx,&ty,&tz,&end)==3) {
            print("Translate 3D object%s with vector (%g,%g,%g).",gmic_inds,tx,ty,tz);
            cimg_foroff(indices,l) {
              const unsigned int ind = indices[l];
              if (!images[ind].is_CImg3d())
                error("Translate 3D object%s : Image [%d] is not a 3D object.",gmic_inds,ind);
              gmic_apply(images[ind],translateCImg3d(tx,ty,tz));
            }
            ++position;
          } else if (std::sscanf(argument,"[%d%c%c",&ind0,&sep,&end)==2 && sep==']') {
            gmic_check_indice(ind0,"Merge object with 3D object%s.");
            const CImg<T> img0 = images[ind0];
            if (!img0.is_CImg3d()) error("Merge object [%d] with 3D object%s : Image [%d] is not a 3D object.",ind0,gmic_inds,ind0);
            print("Merge object [%d] with 3D object%s.",ind0,gmic_inds);
            cimg_foroff(indices,l) {
              const unsigned int ind = indices[l];
              const CImg<T> &img = images[ind];
              if (!img.is_CImg3d())
                error("Merge object [%d] with 3D object%s : Image [%d] is not a 3D object.",ind0,gmic_inds,ind);
              gmic_apply(images[ind],appendCImg3d(img0));
            }
            ++position;
          } else {
            print("Merge 3D object%s together.",gmic_inds);
            if (indices) {
              const unsigned int ind0 = indices[0];
              if (!images[ind0].is_CImg3d())
                error("Merge 3D object%s together : Image [%d] is not a 3D object.",gmic_inds,ind0);
              for (unsigned int siz = indices.size(), off = 0, l = 1; l<siz; ++l) {
                const unsigned int ind = indices[l] - off;
                if (!images[ind].is_CImg3d())
                  error("Merge 3D object%s together : Image [%d] is not a 3D object.",gmic_inds,ind);
                images[ind0].appendCImg3d(images[ind]);
                images.remove(ind); filenames.remove(ind);
                ++off;
              }
            }
          }
          continue;
        }

        // Translate 3D object by the opposite vector.
        if (!cimg::strcmp("-sub3d",item0) || !cimg::strcmp("--3d",item0)) {
          float tx = 0, ty = 0, tz = 0; char end = 0;
          if (std::sscanf(argument,"%f%c",&tx,&end)==1 ||
              std::sscanf(argument,"%f%*c%f%c",&tx,&ty,&end)==2 ||
              std::sscanf(argument,"%f%*c%f%*c%f%c",&tx,&ty,&tz,&end)==3) {
            print("Translate 3D object%s with vector -(%g,%g,%g).",gmic_inds,tx,ty,tz);
            cimg_foroff(indices,l) {
              CImg<T>& img = images[indices[l]];
              CImgList<unsigned int> primitives3d;
              CImgList<unsigned char> colors3d;
              CImg<float> opacities3d;
              CImg<T> points3d;
              if (get_version) points3d.assign(img); else img.transfer_to(points3d);
              points3d.CImg3dtoobject3d(primitives3d,colors3d,opacities3d);
              points3d.get_shared_line(0)-=tx;
              points3d.get_shared_line(1)-=ty;
              points3d.get_shared_line(2)-=tz;
              points3d.object3dtoCImg3d(primitives3d,colors3d,opacities3d);
              if (get_version) {
                images.insert(1); points3d.transfer_to(images.last());
                filenames.insert(filenames[indices[l]]);
              } else points3d.transfer_to(images[indices[l]]);
            }
          } else error("Translate 3D object%s : Invalid argument '%s' "
                       "(should be 'tx,ty,tz').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Scale a 3D object.
        bool divide = false;
        if (!cimg::strcmp("-mul3d",item0) || !cimg::strcmp("-*3d",item0) ||
            ((divide=true)==true && (!cimg::strcmp("-div3d",item0) || !cimg::strcmp("-/3d",item0)))) {
          float sx = 0, sy = 1, sz = 1; char end = 0;
          if ((std::sscanf(argument,"%f%c",&sx,&end)==1 && (sy = sz = sx),1) ||
              std::sscanf(argument,"%f%*c%f%c",&sx,&sy,&end)==2 ||
              std::sscanf(argument,"%f%*c%f%*c%f%c",&sx,&sy,&sz,&end)==3) {
            if (divide) print("Scale 3D object%s with factors (1/%g,1/%g,1/%g).",gmic_inds,sx,sy,sz);
            else print("Scale 3D object%s with factors (%g,%g,%g).",gmic_inds,sx,sy,sz);
            cimg_foroff(indices,l) {
              CImg<T>& img = images[indices[l]];
              CImgList<unsigned int> primitives3d;
              CImgList<unsigned char> colors3d;
              CImg<float> opacities3d;
              CImg<T> points3d;
              if (get_version) points3d.assign(img); else img.transfer_to(points3d);
              points3d.CImg3dtoobject3d(primitives3d,colors3d,opacities3d);
              if (divide) {
                points3d.get_shared_line(0)/=sx;
                points3d.get_shared_line(1)/=sy;
                points3d.get_shared_line(2)/=sz;
              } else {
                points3d.get_shared_line(0)*=sx;
                points3d.get_shared_line(1)*=sy;
                points3d.get_shared_line(2)*=sz;
              }
              points3d.object3dtoCImg3d(primitives3d,colors3d,opacities3d);
              if (get_version) {
                images.insert(1); points3d.transfer_to(images.last());
                filenames.insert(filenames[indices[l]]);
              } else points3d.transfer_to(images[indices[l]]);
            }
          } else error("Scale 3D object%s : Invalid argument '%s' "
                       "(should be 'fact' or 'factx,facty[,factz]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Set color of 3D object(s).
        if (!cimg::strcmp("-color3d",item0) || !cimg::strcmp("-col3d",item0)) {
          float R = 200, G = 200, B = 200, opacity = -1; char end = 0;
          if (std::sscanf(argument,"%f%*c%f%*c%f%c",&R,&G,&B,&end)==3 ||
              std::sscanf(argument,"%f%*c%f%*c%f%*c%f%c",&R,&G,&B,&opacity,&end)==4) {
            const bool set_opacity = (opacity>=0);
            R = (float)cimg::round(R,1); G = (float)cimg::round(G,1); B = (float)cimg::round(B,1);
            if (R<0) R = 0; if (R>255) R = 255;
            if (G<0) G = 0; if (G>255) G = 255;
            if (B<0) B = 0; if (B>255) B = 255;
            if (set_opacity) print("Set colors of 3D object%s to (%g,%g,%g) and opacity to %g.",gmic_inds,R,G,B,opacity);
            else print("Set color of 3D object%s to (%g,%g,%g).",gmic_inds,R,G,B);
            cimg_foroff(indices,l) {
              const unsigned int ind = indices[l];
              if (!images[ind].is_CImg3d())
                error("Set color of 3D object%s : Image [%d] is not a 3D object.",gmic_inds,ind);
              gmic_apply(images[ind],coloropacityCImg3d(R,G,B,opacity,true,set_opacity));
            }
          } else error("Set color of 3D object%s : Invalid argument '%s' "
                       "(should be 'R,G,B[,opacity]').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Set opacity of 3D object(s).
        if (!cimg::strcmp("-opacity3d",item0) || !cimg::strcmp("-opac3d",item0)) {
          float opacity = 1; char end = 0;
          if (std::sscanf(argument,"%f%c",&opacity,&end)==1) {
            print("Set opacity of 3D object%s to %g.",gmic_inds,opacity);
            cimg_foroff(indices,l) {
              const unsigned int ind = indices[l];
              if (!images[ind].is_CImg3d())
                error("Set opacity of 3D object%s : Image [%d] is not a 3D object.",gmic_inds,ind);
              gmic_apply(images[ind],coloropacityCImg3d(0,0,0,opacity,false,true));
            }
          } else error("Set opacity of 3D object%s : Invalid argument '%s' "
                       "(should be 'opacity').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Invert 3D orientation.
        if (!cimg::strcmp("-invert3d",item0) || !cimg::strcmp("-i3d",item0)) {
          print("Invert orientation of 3D object%s.",gmic_inds);
          cimg_foroff(indices,l) {
            CImg<T> &img = images[indices[l]];
            CImgList<unsigned int> primitives3d;
            CImgList<unsigned char> colors3d;
            CImg<float> opacities3d;
            CImg<T> points3d;
            if (get_version) points3d.assign(img); else img.transfer_to(points3d);
            points3d.CImg3dtoobject3d(primitives3d,colors3d,opacities3d);
            if (primitives3d) primitives3d.invert_object3d();
            points3d.object3dtoCImg3d(primitives3d,colors3d,opacities3d);
            if (get_version) {
              images.insert(1); points3d.transfer_to(images.last());
              filenames.insert(filenames[indices[l]]);
            } else points3d.transfer_to(images[indices[l]]);
          }
          continue;
        }

        // Split 3D object(s) into 6 vector images {header,N,vertices,primitives,colors,opacities}
        if (!cimg::strcmp("-split3d",item0) || !cimg::strcmp("-s3d",item0)) {
          print("Split 3D object%s into its different characteristics.",gmic_inds);
          unsigned int off = 0;
          cimg_foroff(indices,l) {
            const unsigned int ind = indices[l] + off;
            CImg<T> &img = images[ind];
            const CImg<char> filename = filenames[ind];
            CImgList<unsigned int> primitives3d;
            CImgList<unsigned char> colors3d;
            CImg<float> opacities3d;
            CImg<T> points3d;
            if (get_version) points3d.assign(img); else img.transfer_to(points3d);
            points3d.CImg3dtoobject3d(primitives3d,colors3d,opacities3d);
            CImgList<T> split;
            split.insert(CImg<T>("CImg3d",1,6,1,1,false)+=0.5f);
            split.insert(CImg<T>::vector((T)points3d.dimx(),(T)primitives3d.size));
            split.insert(1); points3d.resize(-100,3,1,1,0).transpose().unroll('y').transfer_to(split.last());
            points3d.assign();
            CImgList<T> _prims;
            cimglist_for(primitives3d,p)
              _prims.insert(CImg<T>::vector((T)primitives3d[p].size())).insert(primitives3d[p]).last().unroll('y');
            primitives3d.assign();
            split.insert(_prims.get_append('y')); _prims.assign();
            split.insert(colors3d.get_append('x').transpose().unroll('y')); colors3d.assign();
            split.insert(1); opacities3d.transfer_to(split.last());
            if (get_version) {
              images.insert(split);
              filenames.insert(split.size,filename);
            } else {
              images.remove(ind); images.insert(split,ind);
              filenames.remove(ind); filenames.insert(split.size,filename,ind);
              off+=split.size-1;
            }
          }
          continue;
        }

        // Set 3D light position.
        if (!cimg::strcmp("-light3d",item) || !cimg::strcmp("-l3d",item)) {
          float lx = 0, ly = 0, lz = -5000; char end = 0;
          if (std::sscanf(argument,"%f%*c%f%*c%f%c",&lx,&ly,&lz,&end)==3) {
            print("Set 3D light position at (%g,%g,%g).",lx,ly,lz);
            light3d_x = lx;
            light3d_y = ly;
            light3d_z = lz;
          } else error("Set 3D light position : Invalid argument '%s' "
                       "(should be 'posx,posy,posz').",argument_text);
          ++position; continue;
        }

        // Set 3D focale.
        if (!cimg::strcmp("-focale3d",item) || !cimg::strcmp("-f3d",item)) {
          float focale = 500; char end = 0;
          if (std::sscanf(argument,"%f%c",&focale,&end)==1) {
            focale3d = focale;
            print("Set 3D focale to %g.",focale);
          } else error("Set 3D focale : Invalid argument '%s' "
                       "(should be 'value').");
          ++position; continue;
        }

        // Set 3D specular light parameters.
        if (!cimg::strcmp("-specl3d",item) || !cimg::strcmp("-sl3d",item)) {
          float value = 0; char end = 0;
          if (std::sscanf(argument,"%f%c",&value,&end)==1) {
            specular_light3d = value;
            print("Set amount of 3D specular light to %g.",specular_light3d);
          }
          else error("Set amount of 3D specular light : invalid argument '%s'"
                     "(should be 'value').",
                     argument_text);
          ++position; continue;
        }

        if (!cimg::strcmp("-specs3d",item) || !cimg::strcmp("-ss3d",item)) {
          float value = 0; char end = 0;
          if (std::sscanf(argument,"%f%c",&value,&end)==1) {
            specular_shine3d = value;
            print("Set shininess of 3D specular light to %g.",specular_shine3d);
          }
          else error("Set shininess of 3D specular light : invalid argument '%s'"
                     "(should be 'value').",
                     argument_text);
          ++position; continue;
        }

        // Switch double-sided mode for 3D rendering.
        if (!cimg::strcmp("-orient3d",item) || !cimg::strcmp("-o3d",item)) {
          is_oriented3d = !is_oriented3d;
          continue;
        }

        // Set 3D rendering mode.
        if (!cimg::strcmp("-render3d",item) || !cimg::strcmp("-r3d",item)) {
          unsigned int value = 0; char end = 0;
          if (std::sscanf(argument,"%u%c",&value,&end)==1) {
            render3d = value;
            print("Set static 3D render mode to %s.",
                  render3d==-1?"bounding-box":
                  render3d==0?"pointwise":render3d==1?"linear":render3d==2?"flat":
                  render3d==3?"flat-shaded":render3d==4?"Gouraud-shaded":
                  render3d==5?"Phong-shaded":"none");
          }
          else error("Set static 3D render mode : invalid argument '%s'"
                     "(should be '{0=pointwise, 1=linear, 2=flat, 3=flat shaded, 4=Gouraud shaded, 5=Phong-shaded}').",
                     argument_text);
          ++position; continue;
        }

        if (!cimg::strcmp("-renderd3d",item) || !cimg::strcmp("-rd3d",item)) {
          unsigned int value = 0; char end = 0;
          if (std::sscanf(argument,"%u%c",&value,&end)==1) {
            renderd3d = value;
            print("Set dynamic 3D render mode to %s.",
                  renderd3d==-1?"bounding-box":
                  renderd3d==0?"pointwise":renderd3d==1?"linear":renderd3d==2?"flat":
                  renderd3d==3?"flat-shaded":renderd3d==4?"Gouraud-shaded":
                  renderd3d==5?"Phong-shaded":"none");
          }
          else error("Set dynamic 3D render mode : invalid argument '%s'"
                     "(should be '{0=pointwise, 1=linear, 2=flat, 3=flat shaded, 4=Gouraud shaded, 5=Phong-shaded}').",
                     argument_text);
          ++position; continue;
        }

        // Set 3D background color.
        if (!cimg::strcmp("-background3d",item) || !cimg::strcmp("-b3d",item)) {
          int R = 0, G = 0, B = 0; char end = 0;
          const int nb = std::sscanf(argument,"%d%*c%d%*c%d%c",&R,&G,&B,&end);
          switch (nb) {
          case 1 : background3d[0] = background3d[1] = background3d[2] = R; break;
          case 2 : background3d[0] = R; background3d[1] = background3d[2] = G; break;
          case 3 : background3d[0] = R; background3d[1] = G; background3d[2] = B; break;
          default: error("Set 3D background color : Invalid argument '%s'.",argument_text);
          }
          print("Set 3D background color to (%d,%d,%d).",
                (int)background3d[0],(int)background3d[1],(int)background3d[2]);
          ++position; continue;
        }

        //----------------
        // Other commands.
        //----------------

        // No operations : do nothing
        if (!cimg::strcmp("-nop",item)) {
          continue;
        }

        // Skip next argument;
        if (!cimg::strcmp("-skip",item)) {
          ++position;
          continue;
        }

        // Echo.
        if (!cimg::strcmp("-echo",item) || !cimg::strcmp("-e",item)) {
          const int l = cimg::strlen(argument);
          if (l>=2 && argument[0]=='"' && argument[l-1]=='"') {
            if (l==2) print(""); else {
              CImg<char> nargument(argument+1,l-1,1,1,1,false);
              nargument(l-2)=0;
              print("%s",nargument.ptr());
            }
          } else print("%s",argument);
          ++position; continue;
        }

        // Print.
        if (!cimg::strcmp("-print",item0) || !cimg::strcmp("-p",item0)) {
          if (images.size) {
            print("Print image%s.\n\n",gmic_inds);
            char title[4096];
            if (verbosity_level>=0) cimg_foroff(indices,l) {
              const unsigned int ind = indices[l];
              std::sprintf(title,"image [%u] = '%s'",ind,filenames[ind].ptr());
              images[ind].print(title);
            }
            is_released = true;
          } else print("Print image[].");
          continue;
        }

        // Quit.
        if (!cimg::strcmp("-quit",item) || !cimg::strcmp("-q",item)) {
          print("Quit.");
          is_released = true;
          dowhile.assign();
          repeatdone.assign();
          position = command_line.size;
          continue;
        }

        // Do...while.
        if (!cimg::strcmp("-do",item)) {
          dowhile.insert(CImg<int>::vector((int)position));
          continue;
        }

        if (!cimg::strcmp("-while",item)) {
          double cond = 0; char end = 0;
          if (std::sscanf(argument,"%lf%c",&cond,&end)!=1) cond = 0;
          if (!dowhile) error("Directive '-while' is not associated with a '-do' command.");
          if (cond<=0) dowhile.remove();
          else { position = (unsigned int)dowhile.last()(0); continue; }
          ++position; continue;
        }

        // If..else..endif
        if (!cimg::strcmp("-if",item)) {
          double cond = 0; char end = 0;
          if (std::sscanf(argument,"%lf%c",&cond,&end)!=1) cond = 0;
          if (cond<=0) {
            for (int nbifs = 1; nbifs && position<command_line.size; ++position) {
              const char *it = command_line[position].ptr();
              if (!cimg::strcmp("-if",it)) ++nbifs;
              if (!cimg::strcmp("-endif",it)) --nbifs;
              if (!cimg::strcmp("-else",it) && nbifs==1) --nbifs;
            }
            continue;
          }
          ++position; continue;
        }
        if (!cimg::strcmp("-else",item)) {
          for (int nbifs = 1; nbifs && position<command_line.size; ++position) {
            if (!cimg::strcmp("-if",command_line[position].ptr())) ++nbifs;
            if (!cimg::strcmp("-endif",command_line[position].ptr())) --nbifs;
          }
          continue;
        }
        if (!cimg::strcmp("-endif",item)) continue;

        // Repeat...done
        if (!cimg::strcmp("-repeat",item)) {
          float fnb = 0; char end = 0;
          if (std::sscanf(argument,"%f%c",&fnb,&end)==1) {
            const int nb = (int)fnb;
            if (nb>0) repeatdone.insert(CImg<int>::vector((int)position+1,nb));
            else {
              int nbrepeats = 0;
              for (nbrepeats = 1; nbrepeats && position<command_line.size; ++position) {
                const char *it = command_line[position].ptr();
                if (!cimg::strcmp("-repeat",it)) ++nbrepeats;
                if (!cimg::strcmp("-done",it)) --nbrepeats;
              }
              if (nbrepeats && position>=command_line.size)
                error("Directive '-done' is missing after a '-repeat' command.");
              continue;
            }
          } else error("Repeat operation : Invalid argument '%s' "
                       "(should be a number).",argument_text);
          ++position; continue;
        }

        if (!cimg::strcmp("-done",item)) {
          if (!repeatdone) error("Directive '-done' is not associated with a '-repeat' command.");
          if (--repeatdone.last()(1))
            position = (unsigned int)repeatdone.last()(0);
          else repeatdone.remove();
          continue;
        }

        // Check argument type
        if (!cimg::strcmp("-int",item)) {
          char it[4096], end = 0, sep = 0; int value = 0;
          if (*argument) for (const char *nargument = argument; *nargument; ) {
            const int nb = std::sscanf(nargument,"%4095[^,]%c",it,&sep);
            if (nb) {
              if (std::sscanf(it,"%d%c",&value,&end)==1) nargument+=cimg::strlen(it) + nb -1;
              else error("Argument '%s' is not an integer value.",it);
            } else error("Argument '%s' is not an integer value.",argument_text);
          }
          ++position; continue;
        }

        if (!cimg::strcmp("-float",item)) {
          char it[4096], end = 0, sep = 0; double value = 0;
          if (*argument) for (const char *nargument = argument; *nargument; ) {
            const int nb = std::sscanf(nargument,"%4095[^,]%c",it,&sep);
            if (nb) {
              if (std::sscanf(it,"%lf%c",&value,&end)==1) nargument+=cimg::strlen(it) + nb -1;
              else error("Argument '%s' is not a float value.",it);
            } else error("Argument '%s' is not a float value.",argument_text);
          }
          ++position; continue;
        }

        //--------------------------
        // Input/Output and Display
        //--------------------------

        // Display.
        if (!cimg::strcmp("-display",item0) || !cimg::strcmp("-d",item0)) {
          if (display_images(images,indices,true)) is_released = true;
          continue;
        }

        // Display 3D object.
        if (!cimg::strcmp("-display3d",item0) || !cimg::strcmp("-d3d",item0)) {
          if (display_objects3d(images,indices,true)) is_released = true;
          continue;
        }

        // Display as a graph plot.
        if (!cimg::strcmp("-plot",item0)) {
          int plot_type = 1, vertex_type = 1; double ymin = 0, ymax = 0, xmin = 0, xmax = 0; char end = 0;
          const int nb = std::sscanf(argument,"%d%*c%d%*c%lf%*c%lf%*c%lf%*c%lf%c",&plot_type,&vertex_type,&xmin,&xmax,&ymin,&ymax,&end);
          if (nb==1 || nb==2 || nb==4 || nb==6) ++position;
          else { plot_type = 1; vertex_type = 0; ymin = ymax = xmin = xmax = 0; }
          is_released |= display_plots(images,indices,plot_type,vertex_type,xmin,xmax,ymin,ymax,true);
          continue;
        }

        // Select image feature.
        if (!cimg::strcmp("-select",item0)) {
          int select_type = 0; char end = 0;
          if (std::sscanf(argument,"%d%c",&select_type,&end)==1) {
            cimg_foroff(indices,l) gmic_apply(images[indices[l]],select(filenames[indices[l]].ptr(),select_type));
          } else error("Select image%s : Invalid argument '%s' "
                       "(should be 'select_type').",gmic_inds,argument_text);
          ++position; continue;
        }

        // Output.
        if (!cimg::strcmp("-output",item0) || !cimg::strcmp("-o",item0)) {
          char filename[4096] = { 0 }; char options[4096] = { 0 };
          if (std::sscanf(argument,"%4095[^,],%s",filename,options)!=2) std::strcpy(filename,argument);
          const char *const ext = cimg::split_filename(filename);
          if (!cimg::strcasecmp("off",ext)) {
            char nfilename[4096] = { 0 };
            std::strcpy(nfilename,filename);
            const unsigned int siz = indices.size();
            cimg_foroff(indices,l) {
              const unsigned int ind = indices[l];
              if (siz!=1) cimg::number_filename(filename,l,6,nfilename);
              if (!images[ind].is_CImg3d())
                error("Output 3D object [%u] as file '%s' : Image [%u] is not a 3D object.",ind,nfilename,ind);
              print("Output 3D object [%u] as file '%s'.",ind,nfilename);
              CImgList<unsigned int> primitives3d;
              CImgList<unsigned char> colors3d;
              CImg<float> opacities3d;
              CImg<float> points3d(images[ind]);
              points3d.CImg3dtoobject3d(primitives3d,colors3d,opacities3d).save_off(nfilename,primitives3d,colors3d);
            }
          } else if (!cimg::strcasecmp("jpeg",ext) || !cimg::strcasecmp("jpg",ext)) {
            int quality = 100; char end = 0;
            if (std::sscanf(options,"%d%c",&quality,&end)!=1) quality = 100;
            if (quality<0) quality = 0; else if (quality>100) quality = 100;
            CImgList<T> output_images;
            cimg_foroff(indices,l) output_images.insert(images[indices[l]],~0U,true);
            print("Output image%s as file '%s', with quality %u%%",gmic_inds,filename,quality);
            if (!output_images) throw CImgInstanceException("CImgList<%s>::save() : File '%s, instance list (%u,%p) is empty.",
                                                            output_images.pixel_type(),filename,
                                                            output_images.size,output_images.data);
            if (output_images.size==1) output_images[0].save_jpeg(filename,quality);
            else {
              char nfilename[1024];
              cimglist_for(output_images,l) {
                cimg::number_filename(filename,l,6,nfilename);
                output_images[l].save_jpeg(nfilename,quality);
              }
            }
          } else {
            CImgList<T> output_images;
            cimg_foroff(indices,l) output_images.insert(images[indices[l]],~0U,true);
            print("Output image%s as file '%s'.",gmic_inds,filename);
            output_images.save(filename);
          }
          is_released = true; ++position; continue;
        }

        // Substitute macros commands if necessary.
        if (cimg::strcmp("-i",item0) && cimg::strcmp("-input",item0)) {
          bool macro_found = false;
          cimglist_for(macros,l) {
            const char
              *const macro = macros[l].ptr(),
              *const command = commands[l].ptr();

            if (!cimg::strcmp(item+1,macro) && *command) {
              CImgList<char> arguments(256);
              unsigned int nb_arguments = 0;
              char s_argument[4096] = { 0 }, tmp[4096] = { 0 }, tmp2[4096] = { 0 };
              bool has_arguments = false;
              macro_found = true;
              debug("Found macro '%s', substituting by '%s'.",macro,command);

              // Get command-line values of macro arguments.
              if (argument)
                for (const char *nargument = argument; nb_arguments<256 && *nargument &&
                       std::sscanf(nargument,"%4095[^,]",s_argument)==1;) {
                  CImg<char>(s_argument,cimg::strlen(s_argument)+1,1,1,1,false).transfer_to(arguments[nb_arguments++]);
                  nargument+=cimg::strlen(s_argument);
                  if (*nargument) ++nargument;
                }

              // Substitute arguments in macro command expression.
              CImg<char> substituted_command;
              CImgList<char> lreplacement;
              for (const char *ncommand = command; *ncommand;) if (*ncommand=='$') {
                char *replace_text = 0, sep = 0;
                int ind = 0, ind1 = 0;

                // Replace $# and ${#}.
                if (ncommand[1]=='#' || (ncommand[1]=='{' && ncommand[2]=='#' && ncommand[3]=='}')) {
                  std::sprintf(replace_text=s_argument,"%u",nb_arguments);
                  ncommand+=(ncommand[1]=='#')?2:4;
                  has_arguments = true;

                  // Replace $* and ${*}.
                } else if (ncommand[1]=='*' || (ncommand[1]=='{' && ncommand[2]=='*' && ncommand[3]=='}')) {
                  replace_text = &(s_argument[0]=0);
                  for (unsigned int j = 1; j<=nb_arguments; ++j) {
                    replace_text+=std::sprintf(replace_text,"%s",arguments[j-1].ptr());
                    if (j<nb_arguments) *(replace_text++) = ',';
                  }
                  replace_text = s_argument;
                  ncommand+=(ncommand[1]=='*')?2:4;
                  has_arguments = true;

                  // Replace ${i*}.
                } else if (std::sscanf(ncommand,"${%d*%c",&ind,&sep)==2 &&
                           ind>0 && ind<256 && sep=='}') {
                  replace_text = &(s_argument[0]=0);
                  for (unsigned int j = ind; j<=nb_arguments; ++j) {
                    replace_text+=std::sprintf(replace_text,"%s",arguments[j-1].ptr());
                    if (j<nb_arguments) *(replace_text++) = ',';
                  }
                  replace_text = s_argument;
                  ncommand+=std::sprintf(tmp,"${%d*}",ind);
                  has_arguments = true;

                  // Replace $i and ${i}.
                } else if ((std::sscanf(ncommand,"$%d",&ind)==1 ||
                            (std::sscanf(ncommand,"${%d%c",&ind,&sep)==2 && sep=='}')) &&
                           ind>0 && ind<256) {
                  if (!arguments[ind-1]) {
                    if (sep=='}') error("Macro '%s' : Argument '$%d' is undefined (in expression '${%d}').",macro,ind,ind);
                    else error("Macro '%s' : Argument '$%d' is undefined (in expression '$%d').",macro,ind,ind);
                  }
                  replace_text = arguments[ind-1].ptr();
                  ncommand+=std::sprintf(tmp,"$%d",ind) + (sep=='}'?2:0);
                  has_arguments = true;

                  // Replace ${i=$#}.
                } else if (std::sscanf(ncommand,"${%d=$#%c",&ind,&sep)==2 &&
                           ind>0 && ind<256 && sep=='}') {
                  std::sprintf(replace_text=s_argument,"%g",(double)nb_arguments);
                  CImg<char>(s_argument,cimg::strlen(s_argument)+1,1,1,1,false).transfer_to(arguments[ind-1]);
                  ncommand+=std::sprintf(tmp,"${%d=$#}",ind);
                  has_arguments = true;

                  // Replace ${i=$j}.
                } else if (std::sscanf(ncommand,"${%d=$%d%c",&ind,&ind1,&sep)==3 && sep=='}' &&
                           ind>0 && ind<256 && ind1>0 && ind1<256) {
                  if (!arguments[ind1-1])
                    error("Macro '%s' : Argument '$%d' is undefined (in expression '${%d=$%d}').",macro,ind1,ind,ind1);
                  if (!arguments[ind-1]) arguments[ind-1] = arguments[ind1-1];
                  replace_text = arguments[ind-1].ptr();
                  ncommand+=std::sprintf(tmp,"${%d=$%d}",ind,ind1);
                  has_arguments = true;

                  // Replace ${i=default}.
                } else if (std::sscanf(ncommand,"${%d=%4095[^}]%c",&ind,tmp,&sep)==3 && sep=='}' &&
                           ind>0 && ind<256) {
                  if (!arguments[ind-1]) CImg<char>(tmp,cimg::strlen(tmp)+1,1,1,1,false).transfer_to(arguments[ind-1]);
                  replace_text = arguments[ind-1].ptr();
                  ncommand+=cimg::strlen(tmp) + 4 + std::sprintf(tmp2,"%d",ind);
                  has_arguments = true;

                  // Any other expression starting by '$'.
                } else {
                  replace_text = &(s_argument[0]='$');
                  if (std::sscanf(ncommand,"%4095[^$]",s_argument+1)!=1) { s_argument[1] = 0; ++ncommand; }
                  else ncommand+=cimg::strlen(s_argument);
                }

                const int replace_length = cimg::strlen(replace_text);
                if (replace_length) {
                  lreplacement.insert(1);
                  CImg<char>(replace_text,replace_length,1,1,1,false).transfer_to(lreplacement.last());
                }

              } else {
                std::sscanf(ncommand,"%4095[^$]",s_argument);
                const int replace_length = cimg::strlen(s_argument);
                if (replace_length) {
                  lreplacement.insert(1);
                  CImg<char>(s_argument,replace_length,1,1,1,false).transfer_to(lreplacement.last());
                  ncommand+=cimg::strlen(s_argument);
                }
              }
              const CImg<char> zero(1,1,1,1,0);
              lreplacement.insert(zero).get_append('x').transfer_to(substituted_command);

              // Substitute macro expression in command line.
              bool is_dquote = false;
              cimg_foroff(substituted_command,k)
                if (substituted_command[k]=='"') is_dquote = !is_dquote;
                else if (is_dquote && substituted_command[k]==' ') substituted_command[k] = 30;
              CImgList<char> command_items = substituted_command.get_split(' ',false,false);
              cimglist_for(command_items,k) {
                CImg<char> &item = command_items[k];
                cimg_foroff(item,l) if (item[l]==30) item[l]=' ';
              }
              cimglist_for(command_items,k) command_items[k].append(zero,'y');
              if (position<command_line.size && has_arguments) command_line.remove(position);
              command_line.remove(--position);
              command_line.insert(command_items,position);
              break;
            }
          }
          if (macro_found) continue;
        }
      }

      // Input.
      if (!cimg::strcmp("-i",item0) || !cimg::strcmp("-input",item0)) ++position;
      else { if (get_version) --item; argument = item; item1[0] = 0; }
      if (!cimg::strlen(item1)) indices.assign(1,1,1,1,images.size);
      CImgList<T> input_images;
      CImgList<char> input_filenames;
      bool obj3d = false;
      char st_inds[4096] = { 0 }, stx[4096] = { 0 }, sty[4096] = { 0 }, stz[4096] = { 0 }, stv[4096] = { 0 };
      char end = 0, sep = 0, sepx = 0, sepy = 0, sepz = 0, sepv = 0;
      int nb = 1, indx = no_ind, indy = no_ind, indz = no_ind, indv = no_ind;
      float dx = 0, dy = 1, dz = 1, dv = 1;

      if (std::sscanf(argument,"[%4095[0-9%,:-]]%*c%d%c",st_inds,&nb,&end)==2 ||
          (std::sscanf(argument,"[%4095[0-9%,:-]%c%c",st_inds,&sep,&end)==2 && sep==']')) {

        // nb copies of existing sub-images.
        const CImg<unsigned int> indices0 = indices2cimg(st_inds,images.size,"-input");
        char st_tmp[4096] = { 0 }; std::strcpy(st_tmp,indices2string(indices0,true));
        if (nb<=0) error("Input %d copies of image%s : Invalid argument '%s'.",
                         nb,st_tmp,argument_text);
        if (nb!=1) print("Input %d copies of image%s at position%s",nb,st_tmp,gmic_inds);
        else print("Input copy of image%s at position%s",st_tmp,gmic_inds);
        for (int i = 0; i<nb; ++i) cimg_foroff(indices0,l) {
          input_images.insert(images[indices0[l]]);
          input_filenames.insert(filenames[indices0[l]]);
        }
      } else if (((std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%d%c",
                               stx,sty,stz,stv,&(nb=1),&end)==5 ||
                   std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%c",
                               stx,sty,stz,stv,&end)==4) &&
                  (std::sscanf(stx,"%f%c",&dx,&end)==1 ||
                   (std::sscanf(stx,"%f%c%c",&dx,&sepx,&end)==2 && sepx=='%') ||
                   (std::sscanf(stx,"[%d%c%c",&indx,&sepx,&end)==2 && sepx==']')) &&
                  (std::sscanf(sty,"%f%c",&dy,&end)==1 ||
                   (std::sscanf(sty,"%f%c%c",&dy,&sepy,&end)==2 && sepy=='%') ||
                   (std::sscanf(sty,"[%d%c%c",&indy,&sepy,&end)==2 && sepy==']')) &&
                  (std::sscanf(stz,"%f%c",&dz,&end)==1 ||
                   (std::sscanf(stz,"%f%c%c",&dz,&sepz,&end)==2 && sepz=='%') ||
                   (std::sscanf(stz,"[%d%c%c",&indz,&sepz,&end)==2 && sepz==']')) &&
                  (std::sscanf(stv,"%f%c",&dv,&end)==1 ||
                   (std::sscanf(stv,"%f%c%c",&dv,&sepv,&end)==2 && sepv=='%') ||
                   (std::sscanf(stv,"[%d%c%c",&indv,&sepv,&end)==2 && sepv==']'))) ||
                 (std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%c",stx,sty,stz,&end)==3 &&
                  (std::sscanf(stx,"%f%c",&dx,&end)==1 ||
                   (std::sscanf(stx,"%f%c%c",&dx,&sepx,&end)==2 && sepx=='%') ||
                   (std::sscanf(stx,"[%d%c%c",&indx,&sepx,&end)==2 && sepx==']')) &&
                  (std::sscanf(sty,"%f%c",&dy,&end)==1 ||
                   (std::sscanf(sty,"%f%c%c",&dy,&sepy,&end)==2 && sepy=='%') ||
                   (std::sscanf(sty,"[%d%c%c",&indy,&sepy,&end)==2 && sepy==']')) &&
                  (std::sscanf(stz,"%f%c",&dz,&end)==1 ||
                   (std::sscanf(stz,"%f%c%c",&dz,&sepz,&end)==2 && sepz=='%') ||
                   (std::sscanf(stz,"[%d%c%c",&indz,&sepz,&end)==2 && sepz==']'))) ||
                 (std::sscanf(argument,"%4095[][0-9.eE%+-]%*c%4095[][0-9.eE%+-]%c",stx,sty,&end)==2 &&
                  (std::sscanf(stx,"%f%c",&dx,&end)==1 ||
                   (std::sscanf(stx,"%f%c%c",&dx,&sepx,&end)==2 && sepx=='%') ||
                   (std::sscanf(stx,"[%d%c%c",&indx,&sepx,&end)==2 && sepx==']')) &&
                  (std::sscanf(sty,"%f%c",&dy,&end)==1 ||
                   (std::sscanf(sty,"%f%c%c",&dy,&sepy,&end)==2 && sepy=='%') ||
                   (std::sscanf(sty,"[%d%c%c",&indy,&sepy,&end)==2 && sepy==']'))) ||
                 (std::sscanf(argument,"%4095[][0-9.eE%+-]%c",stx,&end)==1 &&
                  (std::sscanf(stx,"%f%c",&dx,&end)==1 ||
                   (std::sscanf(stx,"%f%c%c",&dx,&sepx,&end)==2 && sepx=='%') ||
                   (std::sscanf(stx,"[%d%c%c",&indx,&sepx,&end)==2 && sepx==']')))) {

        // nb new black image.
        if (indx!=no_ind) { gmic_check_indice(indx,"Input black image%s"); dx = (float)images[indx].dimx(); sepx = 0; }
        if (indy!=no_ind) { gmic_check_indice(indy,"Input black image%s"); dy = (float)images[indy].dimy(); sepy = 0; }
        if (indz!=no_ind) { gmic_check_indice(indz,"Input black image%s"); dz = (float)images[indz].dimz(); sepz = 0; }
        if (indv!=no_ind) { gmic_check_indice(indv,"Input black image%s"); dv = (float)images[indv].dimv(); sepv = 0; }
        if (sepx=='%') { dx = images.size?dx*images.last().dimx()/100:0; if (!(int)dx) ++dx; }
        if (sepy=='%') { dy = images.size?dy*images.last().dimy()/100:0; if (!(int)dy) ++dy; }
        if (sepz=='%') { dz = images.size?dz*images.last().dimz()/100:0; if (!(int)dz) ++dz; }
        if (sepv=='%') { dv = images.size?dv*images.last().dimv()/100:0; if (!(int)dv) ++dv; }

        if (nb<=0) error("Input %d black image%s : Invalid number of copies.",nb,gmic_inds);
        if (dx<=0 || dy<=0 || dz<=0 || dv<=0)
          error("Input %d black image%s : Invalid image dimensions %gx%gx%gx%g.",
                nb,gmic_inds,dx,dy,dz,dv);
        if (nb!=1) print("Input %d black images at position%s",nb,gmic_inds);
        else print("Input black image at position%s",gmic_inds);
        CImg<T> empty((int)dx,(int)dy,(int)dz,(int)dv,0);
        input_images.insert(nb-1,empty); input_images.insert(1);
        input_images.last().swap(empty);
        filenames.insert(input_images.size,CImg<char>("(gmic)",7,1,1,1,false));
      } else if (std::sscanf(argument,"(%4095[^)])x%d%c",stx,&(nb=1),&end)==2 ||
                 (std::sscanf(argument,"(%4095[^)]%c%c",stx,&sep,&end)==2 && sep==')')) {

        // Insert nb IxJxKxL image(s) with specified values.
        if (nb<=0) error("Input %d images : Invalid number of copies.",nb);
        unsigned int cx = 0, cy = 0, cz = 0, cv = 0, maxcx = 0, maxcy = 0, maxcz = 0;
        const char *nargument = 0;
        for (nargument = argument+1; *nargument!=')'; ) {
          char s_value[256] = { 0 }, separator = 0; double value = 0;
          if (std::sscanf(nargument,"%255[0-9.eE+-]%c",s_value,&separator)==2 && std::sscanf(s_value,"%lf",&value)==1) {
            if (cx>maxcx) maxcx = cx;
            if (cy>maxcy) maxcy = cy;
            if (cz>maxcz) maxcz = cz;
            switch (separator) {
            case '^' : cx = cy = cz = 0; ++cv; break;
            case '/' : cx = cy = 0; ++cz; break;
            case ';' : cx = 0; ++cy; break;
            default : ++cx;
            }
            nargument+=cimg::strlen(s_value) + (separator==')'?0:1);
          } else break;
        }
        if (*nargument!=')') error("Input %d images : Invalid input string '%s'.",nb,argument);

        CImg<T> img(maxcx+1,maxcy+1,maxcz+1,cv+1,0);
        cx = cy = cz = cv = 0;
        for (nargument = argument+1; *nargument; ) {
          char s_value[256] = { 0 }, separator = 0; double value = 0;
          if (std::sscanf(nargument,"%255[0-9.eE+-]%c",s_value,&separator)==2 && std::sscanf(s_value,"%lf",&value)==1) {
            img(cx,cy,cz,cv) = (T)value;
            switch (separator) {
            case '^' : cx = cy = cz = 0; ++cv; break;
            case '/' : cx = cy = 0; ++cz; break;
            case ';' : cx = 0; ++cy; break;
            default : ++cx;
            }
            nargument+=cimg::strlen(s_value) + (separator==')'?0:1);
          } else break;
        }
        if (nb==1) print("Input image %dx%dx%dx%d",img.dimx(),img.dimy(),img.dimz(),img.dimv());
        else print("Input %d images %dx%d",nb,img.dimx(),img.dimy(),img.dimz(),img.dimv());
        input_images.insert(nb,img); filenames.insert(nb,CImg<char>("(gmic)",7,1,1,1,false));
      } else {

        // Insert image as a loaded filename.
        char filename[4096] = { 0 }, options[4096] = { 0 }; const char *ext = 0, *basename = 0;
        if (argument[0]!='-' || (argument[1] && argument[1]!='.')) {
          std::FILE *file = std::fopen(argument,"r");
          if (file) { std::fclose(file); std::strcpy(filename,argument); }
          else {
            std::sscanf(argument,"%4095[^,],%s",filename,options);
            if (!(file=std::fopen(filename,"r"))) {
              if (filename[0]=='-') error("Input '%s' : Command not found.",filename);
              else error("Input '%s' : File not found.",filename);
            }
            std::fclose(file);
          }
        } else std::strcpy(filename,argument);
        basename = cimg::basename(filename);
        ext = cimg::split_filename(filename);

        if (!cimg::strcasecmp("off",ext)) {  // 3D object file.
          print("Input 3D object '%s'",filename);
          CImgList<unsigned int> primitives3d;
          CImgList<unsigned char> colors3d;
          CImg<float> opacities3d, points3d = CImg<float>::get_load_off(filename,primitives3d,colors3d);
          opacities3d.assign(1,primitives3d.size,1,1,1);
          points3d.object3dtoCImg3d(primitives3d,colors3d,opacities3d);
          input_images.insert(1); points3d.transfer_to(input_images[0]);
          input_filenames.insert(CImg<char>(is_fullpath?filename:basename,
                                            cimg::strlen(is_fullpath?filename:basename)+1,1,1,1,false));
          obj3d = true;
        } else if (!cimg::strcasecmp(ext,"avi") ||
                   !cimg::strcasecmp(ext,"mov") ||
                   !cimg::strcasecmp(ext,"asf") ||
                   !cimg::strcasecmp(ext,"divx") ||
                   !cimg::strcasecmp(ext,"flv") ||
                   !cimg::strcasecmp(ext,"mpg") ||
                   !cimg::strcasecmp(ext,"m1v") ||
                   !cimg::strcasecmp(ext,"m2v") ||
                   !cimg::strcasecmp(ext,"m4v") ||
                   !cimg::strcasecmp(ext,"mjp") ||
                   !cimg::strcasecmp(ext,"mkv") ||
                   !cimg::strcasecmp(ext,"mpe") ||
                   !cimg::strcasecmp(ext,"movie") ||
                   !cimg::strcasecmp(ext,"ogm") ||
                   !cimg::strcasecmp(ext,"qt") ||
                   !cimg::strcasecmp(ext,"rm") ||
                   !cimg::strcasecmp(ext,"vob") ||
                   !cimg::strcasecmp(ext,"wmv") ||
                   !cimg::strcasecmp(ext,"xvid") ||
                   !cimg::strcasecmp(ext,"mpeg")) {
          unsigned int value0 = 0, value1 = 0, step = 1; char sep0 = 0, sep1 = 0, end = 0;
          if ((std::sscanf(options,"%u%c%*c%u%c%*c%u%c",&value0,&sep0,&value1,&sep1,&step,&end)==5 && sep0=='%' && sep1=='%') ||
              (std::sscanf(options,"%u%c%*c%u%*c%u%c",&value0,&sep0,&value1,&step,&end)==4 && sep0=='%') ||
              (std::sscanf(options,"%u%*c%u%c%*c%u%c",&value0,&value1,&sep1,&step,&end)==4 && sep1=='%') ||
              (std::sscanf(options,"%u%*c%u%*c%u%c",&value0,&value1,&step,&end)==3) ||
              (std::sscanf(options,"%u%c%*c%u%c%c",&value0,&sep0,&value1,&sep1,&end)==4 && sep0=='%' && sep1=='%') ||
              (std::sscanf(options,"%u%c%*c%u%c",&value0,&sep0,&value1,&end)==3 && sep0=='%') ||
              (std::sscanf(options,"%u%*c%u%c%c",&value0,&value1,&sep1,&end)==3 && sep1=='%') ||
              (std::sscanf(options,"%u%*c%u%c",&value0,&value1,&end)==2)) { // Read several frames
            print("Input frames %u%s...%u%s with step %u of file '%s'",
                  value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",step,filename);
            if (sep0=='%' || sep1=='%') {
              const unsigned int nb_frames = CImg<unsigned int>::get_load_ffmpeg(filename,0,0,0)[0];
              if (sep0=='%') value0 = value0*nb_frames/100;
              if (sep1=='%') value1 = value1*nb_frames/100;
            }
          } else if ((std::sscanf(options,"%u%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
                     (std::sscanf(options,"%u%c",&value0,&end)==1)) { // Read one frame
            print("Input frame %u%s of file '%s'",value0,sep0=='%'?"%":"",filename);
            if (sep0=='%') {
              const unsigned int nb_frames = CImg<unsigned int>::get_load_ffmpeg(filename,0,0,0)[0];
              value0 = value0*nb_frames/100;
            }
            value1 = value0; step = 1;
          } else { // Read all frames
            print("Input all frames of file '%s'",filename);
            value0 = 0; value1 = ~0U; sep0 = sep1 = 0; step = 1;
          }
          input_images.load_ffmpeg(filename,value0,value1,step);
          if (input_images)
            input_filenames.insert(input_images.size,CImg<char>(is_fullpath?filename:basename,
                                                                cimg::strlen(is_fullpath?filename:basename)+1,1,1,1,false));
        } else if (!cimg::strcasecmp("raw",ext)) { // Raw file.
          int dx = 0, dy = 1, dz = 1, dv = 1;
          if (std::sscanf(options,"%d%*c%d%*c%d%*c%d",&dx,&dy,&dz,&dv)>0) {
            if (dx<=0 || dy<=0 || dz<=0 || dv<=0)
              error("Input raw file '%s' : Invalid specified dimensions %dx%dx%dx%d.",filename,dx,dy,dz,dv);
            print("Input raw file '%s'",filename);
            input_images.insert(1); input_images[0].load_raw(filename,dx,dy,dz,dv);
            input_filenames.insert(CImg<char>(is_fullpath?filename:basename,
                                              cimg::strlen(is_fullpath?filename:basename)+1,1,1,1,false));
          } else error("Input raw file '%s' : Image dimensions must be specified.",filename);
        } else if (!cimg::strcasecmp("yuv",ext)) { // YUV file.
          int dx = 0, dy = 0; unsigned int first = 0, last = ~0U, step = 1;
          if (std::sscanf(options,"%d%*c%d%*c%u%*c%u%*c%u",&dx,&dy,&first,&last,&step)>0) {
            if (dx<=0 || dy<=0)
              error("Input yuv file '%s' : Invalid specified dimensions %dx%d.",filename,dx,dy);
            print("Input yuv file '%s'",filename);
            input_images.load_yuv(filename,dx,dy,first,last,step);
            input_filenames.insert(input_images.size,CImg<char>(is_fullpath?filename:basename,
                                                                cimg::strlen(is_fullpath?filename:basename)+1,1,1,1,false));
          } else error("Input yuv file '%s' : Image dimensions must be specified.",filename);
        } else { // Other file type.
          print("Input file '%s'",filename);
          input_images.load(filename);
          input_filenames.insert(input_images.size,
                                 CImg<char>(is_fullpath?filename:basename,
                                            cimg::strlen(is_fullpath?filename:basename)+1,1,1,1,false));
        }
      }

      if (verbosity_level>=0) {
        if (input_images) {
          const unsigned int last = input_images.size-1;
          if (obj3d)
            std::fprintf(cimg_stdout," (%d points, %u primitives, %u colors).",
                         (unsigned int)input_images(0,6),
                         (unsigned int)input_images(0,7),
                         (unsigned int)input_images(0,8));
          else if (input_images.size==1)
            std::fprintf(cimg_stdout," (1 image %ux%ux%ux%u).",
                         input_images[0].width,input_images[0].height,input_images[0].depth,
                         input_images[0].dim);
          else std::fprintf(cimg_stdout," (%u images [0] = %ux%ux%ux%u, %s[%u] = %ux%ux%ux%u).",
                            input_images.size,
                            input_images[0].width,input_images[0].height,input_images[0].depth,
                            input_images[0].dim,
                            last==1?"":"...,",
                            last,
                            input_images[last].width,input_images[last].height,input_images[last].depth,
                            input_images[last].dim);
        } else std::fprintf(cimg_stdout," (no available data).");
      }

      for (unsigned int l = 0, siz = indices.size()-1, off = 0; l<=siz; ++l) {
        const unsigned int ind = indices[l] + off;
        if (l!=siz) images.insert(input_images,ind);
        else {
          images.insert(input_images.size,ind);
          cimglist_for(input_images,k) images[ind+k].swap(input_images[k]);
        }
        filenames.insert(input_filenames,ind);
        off+=input_images.size;
      }

    } catch (CImgException &e) {
      const char *error_message = e.message;
      char tmp[4096] = { 0 }, sep = 0;
      if (std::sscanf(error_message,"%4095[^>]>:%c",tmp,&sep)==2 && sep==':') error_message+=cimg::strlen(tmp)+3;
      error(error_message);
    }
  }

  // Check if command line has grown too much (possible recursive macro calls).
  if (command_line.size>=command_line_maxsize)
    error("Command line overflow : There are too much commands specified (possible recursive macro substitution).");

  // Check if some loops have not been terminated.
  if (dowhile) warning("A '-while' directive is missing somewhere.");
  if (repeatdone) warning("A '-done' directive is missing somewhere.");

  // Display final result if necessary (not 'released' before).
  if (images.size && !is_released) {
    if (!display_objects3d(images,CImg<unsigned int>::sequence(images.size,0,images.size-1),false))
      display_images(images,CImg<unsigned int>::sequence(images.size,0,images.size-1),true);
  }

  print("End G'MIC instance.\n");
  return *this;
}

// Small hack to separate the compilation of G'MIC in different pixel types.
// (only intended to save computer memory when compiling !)
//--------------------------------------------------------------------------
#ifdef gmic_minimal
gmic& gmic::parse_float(CImgList<float>& images) { return parse(images); }
template gmic::gmic(const int, const char *const *const, CImgList<float>&, const char *const custom_macros, const bool add_macros_at_start);
template gmic::gmic(const char* const, CImgList<float>&, const char *const custom_macros, const bool add_macros_at_start);
#else
#if defined(gmic_bool) || !defined(gmic_separate_compilation)
gmic& gmic::parse_bool(CImgList<bool>& images) { return parse(images); }
template gmic::gmic(const int, const char *const *const, CImgList<bool>&, const char *const custom_macros, const bool add_macros_at_start);
template gmic::gmic(const char* const, CImgList<bool>&, const char *const custom_macros, const bool add_macros_at_start);
#endif
#if defined(gmic_uchar) || !defined(gmic_separate_compilation)
gmic& gmic::parse_uchar(CImgList<unsigned char>& images) { return parse(images); }
template gmic::gmic(const int, const char *const *const, CImgList<unsigned char>&, const char *const custom_macros, const bool add_macros_at_start);
template gmic::gmic(const char* const, CImgList<unsigned char>&, const char *const custom_macros, const bool add_macros_at_start);
#endif
#if defined(gmic_char) || !defined(gmic_separate_compilation)
gmic& gmic::parse_char(CImgList<char>& images) { return parse(images); }
template gmic::gmic(const int, const char *const *const, CImgList<char>&, const char *const custom_macros, const bool add_macros_at_start);
template gmic::gmic(const char* const, CImgList<char>&, const char *const custom_macros, const bool add_macros_at_start);
#endif
#if defined(gmic_ushort) || !defined(gmic_separate_compilation)
gmic& gmic::parse_ushort(CImgList<unsigned short>& images) { return parse(images); }
template gmic::gmic(const int, const char *const *const, CImgList<unsigned short>&, const char *const custom_macros, const bool add_macros_at_start);
template gmic::gmic(const char* const, CImgList<unsigned short>&, const char *const custom_macros, const bool add_macros_at_start);
#endif
#if defined(gmic_short) || !defined(gmic_separate_compilation)
gmic& gmic::parse_short(CImgList<short>& images) { return parse(images); }
template gmic::gmic(const int, const char *const *const, CImgList<short>&, const char *const custom_macros, const bool add_macros_at_start);
template gmic::gmic(const char* const, CImgList<short>&, const char *const custom_macros, const bool add_macros_at_start);
#endif
#if defined(gmic_uint) || !defined(gmic_separate_compilation)
gmic& gmic::parse_uint(CImgList<unsigned int>& images) { return parse(images); }
template gmic::gmic(const int, const char *const *const, CImgList<unsigned int>&, const char *const custom_macros, const bool add_macros_at_start);
template gmic::gmic(const char* const, CImgList<unsigned int>&, const char *const custom_macros, const bool add_macros_at_start);
#endif
#if defined(gmic_int) || !defined(gmic_separate_compilation)
gmic& gmic::parse_int(CImgList<int>& images) { return parse(images); }
template gmic::gmic(const int, const char *const *const, CImgList<int>&, const char *const custom_macros, const bool add_macros_at_start);
template gmic::gmic(const char* const, CImgList<int>&, const char *const custom_macros, const bool add_macros_at_start);
#endif
#if defined(gmic_float) || !defined(gmic_separate_compilation)
gmic& gmic::parse_float(CImgList<float>& images) { return parse(images); }
template gmic::gmic(const int, const char *const *const, CImgList<float>&, const char *const custom_macros, const bool add_macros_at_start);
template gmic::gmic(const char* const, CImgList<float>&, const char *const custom_macros, const bool add_macros_at_start);
#endif
#if defined(gmic_double) || !defined(gmic_separate_compilation)
gmic& gmic::parse_double(CImgList<double>& images) { return parse(images); }
template gmic::gmic(const int, const char *const *const, CImgList<double>&, const char *const custom_macros, const bool add_macros_at_start);
template gmic::gmic(const char* const, CImgList<double>&, const char *const custom_macros, const bool add_macros_at_start);
#endif
#endif
#endif

//-----------------------
// Start main procedure.
//-----------------------
#if defined(gmic_main) || (!defined(gmic_separate_compilation) && !defined(gmic_minimal))
extern char data_def[];

int main(int argc, char **argv) {

  // Display help if necessary.
  //---------------------------
  if (argc==1) {
    std::fprintf(cimg_stdout,"<gmic> No options or data provided. Try '%s -h' for help.\n",cimg::basename(argv[0]));
    std::exit(0);
  }

  if (cimg_option("-h",false,0) || cimg_option("-help",false,0) || cimg_option("--help",false,0)) {
    cimg_usage("GREYC's Magic Image Converter");

    char version[1024] = { 0 };
    std::sprintf(version,"        Version %d.%d.%d.%d, Copyright (C) 2008-2009, David Tschumperle (http://gmic.sourceforge.net)",
                 gmic_version/1000,(gmic_version/100)%10,(gmic_version/10)%10,gmic_version%10);
    cimg_help(version);

    cimg_help("\n  Usage\n"
              "  -----");
    cimg_help("  gmic [file_1] [file_2] .. [file_n] [command_1] .. [command_n] [file_n+1] ...\n");
    cimg_help("  G'MIC is an interpreter of image processing macros whose goal is to convert, manipulate and");
    cimg_help("  visualize generic 1D/2D/3D multi-spectral image and video files. It follows these simple rules :\n");
    cimg_help("    - G'MIC handles a numbered list of images which are all stored in computer memory.");
    cimg_help("    - The first image of the list has indice '[0]'.");
    cimg_help("    - Negative indices are treated in a cyclic way (i.e. '[-1]' is the last image,");
    cimg_help("      '[-2]' the penultimate one, and so on...).");
    cimg_help("    - Command line items tell how to add/remove/manipulate/display images of the list.");
    cimg_help("    - Items are read and executed in the order they appear on the command line, from the left to the right.");
    cimg_help("    - Items can thus appear more than one time on the command line.");
    cimg_help("    - An item starting by '-' is a G'MIC instruction.");
    cimg_help("    - One instruction may have two equivalent names (regular and short).");
    cimg_help("    - A G'MIC instruction may have mandatory or optional arguments.");
    cimg_help("    - When multiple arguments are needed, they are separated by commas ','.");
    cimg_help("    - Items that are not instructions are considered either as input filenames or input strings.");
    cimg_help("    - When an input filename is encountered, the corresponding image data are loaded");
    cimg_help("      and added to the end of the image list.");
    cimg_help("      (see section 'Filename options' below for more informations on file input/output).");
    cimg_help("    - Special filenames '-' or '-.ext' mean 'standard input/output' (optionally. in 'ext' format).");
    cimg_help("    - Special input strings can be used to insert new images to the list. They can be :");
    cimg_help("        - 'width[%][xheight[%][xdepth[%][xdim[%][xN]]]]' : Insert 'N' black images with specified size.");
    cimg_help("          (adding '%' to a dimension means 'percentage of to the same dimension in the last image'),");
    cimg_help("        - '[indice]' or '[indice]xN' : Insert 1 or N copies of the existing image [indice].");
    cimg_help("        - '(v1,v2,...)' or '(v1,v2,...)xN' : Insert 1 or N copies of the specified IxJxKxL image.");
    cimg_help("          Separators inside '(..)' can be ',' (column), ';' (line), '/' (slice) or '^' (channel).");
    cimg_help("    - A G'MIC instruction may be restricted to a specific subset of the list, by adding '[subset]' to");
    cimg_help("      the instruction name. Several usual expressions are possible for 'subset', for instance : ");
    cimg_help("        '-command[0,1,3]' : Apply instruction on images 0,1 and 3.");
    cimg_help("        '-command[3-5]' : Apply instruction on images 3 to 5.");
    cimg_help("        '-command[50%-100%] : Apply instruction on the second half of the image list.");
    cimg_help("        '-command[0,-2,-1]' : Apply instruction on the first and two latest images.");
    cimg_help("        '-command[0-9:3]' : Apply instruction on images 0 to 9, with a step of 3 (i.e. images 0,3,6,9).");
    cimg_help("        '-command[0,2-4,50%--1]' : Apply instruction on images 0,2,3,4 and the second half of the list.");
    cimg_help("    - When no image subset is specified, a G'MIC instruction is applied on all images of the list.");
    cimg_help("    - Native (non-macro) instructions starting with '--' instead of '-' do not act 'in-place' but");
    cimg_help("      insert their result as a new image at the end of the list.");
    cimg_help("    - On the command line, any item of the form '@indice' or '@{indice}' is replaced");
    cimg_help("      by the values of the image '[indice]' separated by commas.");
    cimg_help("    - Items '@?' (or '@{?}'), '@{?,max}' or '@{?,min,max}' are replaced by a float random value between [0,1], [0,max] or [min,max].");
    cimg_help("    - Items '@[?]', '@[?,max]' or '@[?,min,max]' do the same but return integer random values.");
    cimg_help("    - Restrictions to a subset of image values can be specified with '@{indice,subset}' (as in '@{2,0-50%}').");
    cimg_help("    - Restriction to a particular pixel coordinate can be specified with '@(indice,x[,y[,z[,v[,borders]]]])'.");
    cimg_help("    - On the command line, the item '@#' is replaced by the number of images in the list.");
    cimg_help("    - Input filenames or commands may result to the generation of 3D objects.");
    cimg_help("    - A 3D object viewer, based on software rendering, is included in G'MIC.");
    cimg_help("    - A 3D object is stored as a single-column image containing all object data, in the following order :");
    cimg_help("        { magic header, vertices, faces, colors, opacities }.");
    cimg_help("    - Custom user-defined G'MIC instructions can be defined with the use of a macro file.");
    cimg_help("    - A macro file is a simple ASCII text file, each line being of the form");
    cimg_help("        'instruction_name : substitution' or 'substitution (continuation)' or '# comment'.");
    cimg_help("    - Each invoked macro instruction is substituted as its defined content on the command line.");
    cimg_help("    - A macro file 'gmic_def.raw' is distributed within the G'MIC package.");
    cimg_help("    - The macros defined in 'gmic_def.raw' are included by default in G'MIC.");
    cimg_help("    - Macros arguments, separated by commas, can be added after the invokation of a macro instruction.");
    cimg_help("    - In macros definitions, expressions starting with '$' are used to reference macro arguments :");
    cimg_help("        $i and ${i} are replaced by the value of the i-th macro argument.");
    cimg_help("        $# and ${#} are replaced by the number of macro arguments.");
    cimg_help("        $* and ${*} are replaced by the entire string of macro arguments.");
    cimg_help("        ${i*} is replaced by all macro arguments following the i-th argument (included).");
    cimg_help("        ${i=$#} is replaced by $i if defined, or by its new value $# else.");
    cimg_help("        ${i=$j} is replaced by $i if defined, or by its new value $j else.");
    cimg_help("        ${i=default} is replaced by $i if defined, or by its new value 'default' else.");
    cimg_help("\n  A list of available native and macro instructions is available below.");
    cimg_help("  A parameter specified in '[]' is optional, except when standing for '[indices]' where it");
    cimg_help("  corresponds to one or several indices of the image list, as described above. In this case, the '[' and ']'");
    cimg_help("  characters must explicitly appear when writting the item.");

    cimg_help("\n  Global options\n"
              "  --------------");
    cimg_option("-help","(no args)","Display this help (eq. to '-h').");
    cimg_option("-verbose+","(no args)","Increment verbosity level (eq. to '-v+').");
    cimg_option("-verbose-","(no args)","Decrement verbosity level (eq. to '-v-').");
    cimg_option("-macros","'filename'","Load macro file from specified filename (eq. to '-m').");
    cimg_option("-debug","(no args)","Switch debug flag (when on, displays internal infos for debugging).");
    cimg_option("-fullpath","(no args)","Switch full path flag (when on, displays full filename paths).");

    cimg_help("\n  Arithmetic operators\n"
              "  --------------------");
    cimg_option("-add","'value', '[indice]' or (no args)","Add 'value' or '[indice]' to image(s)");
    cimg_help("                                              "
              "or add image(s) together (eq. to '-+').");
    cimg_option("-sub","'value', '[indice]' or (no args)","Substract 'value' or '[indice]' to image(s)");
    cimg_help("                                              "
              "or substract image(s) together (eq. to '--').");
    cimg_option("-mul","'value', '[indice]' or (no args)","Multiply image(s) by 'value' or '[indice]'");
    cimg_help("                                              "
              "or multiply image(s) together (eq. to '-*').");
    cimg_option("-div","'value', '[indice]' or (no args)","Divide image(s) by 'value' or '[indice]'");
    cimg_help("                                              "
              "or divide image(s) together (eq. to '-/').");
    cimg_option("-pow","'value', '[indice]' or (no args)","Compute image(s) to the power of 'value' or '[indice]'");
    cimg_help("                                              "
              "or power of the image(s) together (eq. to '-^').");
    cimg_option("-min","'value', '[indice]' or (no args)","Compute minimum between image(s) and 'value' or '[indice]'");
    cimg_help("                                              "
              "or minimum of image(s) together.");
    cimg_option("-max","'value', '[indice]' or (no args)","Compute maximum between image(s) and 'value' or '[indice]'");
    cimg_help("                                              "
              "or maximum of image(s) together.");
    cimg_option("-mod","'value', '[indice]' or (no args)","Compute modulo of image(s) with 'value' or '[indice]'");
    cimg_help("                                              "
              "or modulo with image(s) together.");
    cimg_option("-and","'value', '[indice]' or (no args)","Compute bitwise AND of image(s) with 'value' or '[indice]'");
    cimg_help("                                              "
              "or bitwise AND of image(s) together.");
    cimg_option("-or","'value', '[indice]' or (no args)","Compute bitwise OR of image(s) with 'value' or '[indice]'");
    cimg_help("                                              "
              "or bitwise OR of image(s) together.");
    cimg_option("-xor","'value', '[indice]' or (no args)","Compute bitwise XOR of image(s) with 'value' '[indice]'");
    cimg_help("                                              "
              "or bitwise XOR of image(s) together.");
    cimg_option("-cos","(no args)","Compute cosine of image(s) values.");
    cimg_option("-sin","(no args)","Compute sine of image(s) values.");
    cimg_option("-tan","(no args)","Compute tangent of image(s) values.");
    cimg_option("-acos","(no args)","Compute arccosine of image(s) values.");
    cimg_option("-asin","(no args)","Compute arcsine of image(s) values.");
    cimg_option("-atan","(no args)","Compute arctangent of image(s) values.");
    cimg_option("-abs","(no args)","Compute absolute value of image(s) values.");
    cimg_option("-sqr","(no args)","Compute square of image(s) values.");
    cimg_option("-sqrt","(no args)","Compute square root of image(s) values.");
    cimg_option("-exp","(no args)","Compute exponential of image(s) values.");
    cimg_option("-log","(no args)","Compute logarithm of image(s) values.");
    cimg_option("-log10","(no args)","Compute logarithm_10 of image(s) values.");

    cimg_help("\n  Pointwise pixel manipulation\n"
              "  ----------------------------");
    cimg_option("-type","'value_type'","Cast all images into specified value type (eq. to '-t').");
    cimg_help("                                              "
              "('value_type' can be 'bool','uchar','char','ushort','short',");
    cimg_help("                                              "
              "'uint','int','float','double').");
    cimg_option("-set","'value[,x[,y[,z[,v]]]]'","Set scalar value in image(s) at specified position (eq. to '-=').");
    cimg_option("-endian","(no args)","Invert endianness of the image(s) buffers.");
    cimg_option("-fill","'value1,value2,...'","Fill image(s) with scalar values in a repetitive way (eq. to '-f').");
    cimg_option("-threshold","'value[%][,soft]' or (noargs)","Threshold pixel values ((noargs) for interactive mode).");
    cimg_option("-cut","'{value1[%],[indice]},{value2[%],[indice]}' or (noargs)","Cut pixel values in specified range ");
    cimg_help("                                              "
              "((noargs) for interactive mode).");
    cimg_option("-normalize","'{value1[%],[indice]},{value2[%],[indice]}'",
                "Normalize pixel values in specified range (eq. to '-n').");
    cimg_option("-round","'round_value[,round_type]'","Round pixel values.");
    cimg_option("-equalize","'nb_levels'","Equalize image(s) histogram(s) using specified number of levels.");
    cimg_option("-quantize","'nb_levels'","Quantize image(s) with 'nb_levels' levels.");
    cimg_option("-noise","'std[%][,noise_type]'","Add noise with specified standard deviation");
    cimg_help("                                              "
              "('noise_type' can be '{0=gaussian, 1=uniform, 2=salt&pepper, 3=poisson}'.");
    cimg_option("-norm","(no args)","Compute pointwise L2-norm of vector-valued image(s).");
    cimg_option("-orientation","(no args)","Compute pointwise orientation of vector-valued image(s).");

    cimg_help("\n  Color bases conversions\n"
              "  -----------------------");
    cimg_option("-rgb2hsv","(no args)","Convert image(s) from RGB to HSV colorbases.");
    cimg_option("-rgb2hsl","(no args)","Convert image(s) from RGB to HSL colorbases.");
    cimg_option("-rgb2hsi","(no args)","Convert image(s) from RGB to HSI colorbases.");
    cimg_option("-rgb2yuv","(no args)","Convert image(s) from RGB to YUV colorbases.");
    cimg_option("-rgb2ycbcr","(no args)","Convert image(s) from RGB to YCbCr colorbases.");
    cimg_option("-rgb2xyz","(no args)","Convert image(s) from RGB to XYZ colorbases.");
    cimg_option("-rgb2lab","(no args)","Convert image(s) from RGB to Lab colorbases.");
    cimg_option("-rgb2cmy","(no args)","Convert image(s) from RGB to CMY colorbases.");
    cimg_option("-rgb2cmyk","(no args)","Convert image(s) from RGB to CMYK colorbases.");
    cimg_option("-rgb2lut","'[indice][,dithering]' or 'LUT_type[,dithering]'","Index image(s) with color palette ");
    cimg_help("                                              "
              "('LUT_type' can be '{0=default, 1=rainbow, 2=contrast}'");
    cimg_help("                                              "
              "'dithering' can be '{0=off, 1=on}').");
    cimg_option("-hsv2rgb","(no args)","Convert image(s) from HSV to RGB colorbases.");
    cimg_option("-hsl2rgb","(no args)","Convert image(s) from HSL to RGB colorbases.");
    cimg_option("-hsi2rgb","(no args)","Convert image(s) from HSI to RGB colorbases.");
    cimg_option("-yuv2rgb","(no args)","Convert image(s) from YUV to RGB colorbases.");
    cimg_option("-ycbcr2rgb","(no args)","Convert image(s) from YCbCr to RGB colorbases.");
    cimg_option("-xyz2rgb","(no args)","Convert image(s) from XYZ to RGB colorbases.");
    cimg_option("-lab2rgb","(no args)","Convert image(s) from Lab to RGB colorbases.");
    cimg_option("-cmy2rgb","(no args)","Convert image(s) from CMY to RGB colorbases.");
    cimg_option("-cmyk2rgb","(no args)","Convert image(s) from CMYK to RGB colorbases.");
    cimg_option("-lut2rgb","'[indice]' or 'LUT_type'","Map color palette to image(s) ");
    cimg_help("                                              "
              "('LUT_type' can be '{0=default, 1=rainbow, 2=contrast}'.");

    cimg_help("\n  Geometric manipulation\n"
              "  ----------------------");
    cimg_option("-resize","'[indice][,interpolation[,borders[,center]]]' or ","");
    cimg_help("                     "
              "'{[indice],width[%]}[x{[indice],height[%]}[x{[indice],depth[%]}[x{[indice],dim[%]}[,interpolation[,borders[,center]]]]]]'");
    cimg_help("                                              "
              "or (noargs)");
    cimg_help("                                              "
              "Resize image(s) to specified geometry ((noargs) for interactive mode) (eq. to '-r')");
    cimg_help("                                              "
              "('interpolation' can be '{0=none, 1=nearest, 2=average, 3=linear, 4=grid, 5=cubic}').");
    cimg_option("-resize2x","(no args)","Resize image(s) with Scale2x.");
    cimg_option("-resize3x","(no args)","Resize image(s) with Scale3x.");
    cimg_option("-crop","'x0[%],x1[%][,border_conditions]' or 'x0[%],y0[%],x1[%],y1[%][,border_conditions]' or ","");
    cimg_help("                                              "
              "'x0[%],y0[%],z0[%],x1[%],y1[%],z1[%][,border_conditions]' or ");
    cimg_help("                                              "
              "'x0[%],y0[%],z0[%],v0[%],x1[%],y1[%],z1[%],v1[%][,border_conditions]' or (noargs).");
    cimg_help("                                              "
              "Crop image(s) using specified geometry ((noargs) for interactive mode) (eq. to '-c') ");
    cimg_help("                                              "
              "('border_conditions' can be '{0=zero, 1=nearest}')");
    cimg_help("                                              "
              "((no args) for interactive mode).");
    cimg_option("-autocrop","'value1,value2,...'","Autocrop image(s) using specified background color.");
    cimg_option("-channels","'{[ind0],v0[%]}[,{[ind1],v1[%]}]'","Select channels v0..v1 of multi-spectral image(s).");
    cimg_option("-slices","'{[ind0],z0[%]}[,{[ind1],z1[%]}]'","Select slices z0..z1 of volumetric image(s).");
    cimg_option("-lines","'{[ind0],y0[%]}[,{[ind1],y1[%]}]'","Select lines y0..y1 of image(s).");
    cimg_option("-columns","'{[ind0],x0[%]}[,{[ind1],x1[%]}]'","Select columns x0..x1 of image(s).");
    cimg_option("-rotate","'angle[,border_conditions]'","Rotate image(s) with a given angle ");
    cimg_help("                                              "
              "('border_conditions' can be '{-3=cyclic (in-place), -2=nearest(ip), -1=zero(ip), 0=zero, 1=nearest, 2=cyclic}'");
    cimg_help("                                              "
              "and 'interpolation' can be '{0=none, 1=linear, 2=cubic}').");
    cimg_option("-mirror","'axis'",
                "Mirror image(s) along specified axis ('axis' can be '{x,y,z,v}').");
    cimg_option("-translate","'tx[%][,ty[%][,tz[%][,tv[%][,border_conditions]]]]'",
                "Translate image(s) by vector (dx,dy,dz,dv)");
    cimg_help("                                              "
              "('border_conditions' can be '{0=zero, 1=nearest, 2=cyclic}').");
    cimg_option("-transpose","(no args)","Transpose image(s).");
    cimg_option("-invert","(no args)","Compute inverse matrix.");
    cimg_option("-permute","'permutation'","Permute image axes by the specified permutation "
                "('permutation' can be 'yxzv',...).");
    cimg_option("-unroll","'axis'",
                "Unroll image(s) along specified axis ('axis' can be '{x,y,z,v}').");
    cimg_option("-split","'axis[,nb_parts]' or 'value[,keep]'",
                "Split image(s) along specified axis or value ('axis' can be '{x,y,z,v}') (eq. to '-s').");
    cimg_option("-append","'axis,[alignement]'","Append image(s) along specified axis and alignement (eq. to '-a')");
    cimg_help("                                              "
              "('axis' can be '{x,y,z,v}' and 'alignement' can be '{p=left, c=center, n=right)'.");
    cimg_option("-warp","'[indice][,relative[,interpolation[,border_conditions[,nb_frames]]]]'",
                "Warp image(s) in 'nb_frames' with field '[indice]' ");
    cimg_help("                                              "
              "('relative' can be '{0,1}', 'interpolation' can be '{0,1}', "
              "'border_conditions' can be '{0=zero, 1=nearest}').");

    cimg_help("\n  Image filtering\n"
              "  ---------------");
    cimg_option("-blur","'std[,border_conditions]'",
                "Apply gaussian blur of specified standard deviation");
    cimg_help("                                              "
              "('border_conditions' can be '{0=zero, 1=nearest}').");
    cimg_option("-bilateral","'stdevs,stdevr'",
                "Apply bilateral filtering with specified standard deviations 'stdevs' and 'stdevr'.");
    cimg_option("-smooth","'amplitude[,sharpness[,anisotropy[,alpha[,sigma[,dl[,da[,prec[,interp[,fast]]]]]]]]]'","");
    cimg_help("                                              "
              "Smooth image(s) anisotropically with specified parameters.");
    cimg_option("-denoise","'patch_size[,stdev_p[,stdev_s[,lookup_size]]]'",
                "Denoise image(s) with a patch-averaging procedure.");
    cimg_option("-median","'size'","Apply median filter with specified size.");
    cimg_option("-sharpen","'amplitude[,0]' or 'amplitude,1[,edge[,alpha[,sigma]]]'",
                "Sharpen image(s) with inverse diffusion or shock filters.");
    cimg_option("-convolve","'[indice][,border_conditions]'",
                "Convolve image(s) by the specified mask");
    cimg_help("                                              "
              "('border_conditions' can be '{0=zero, 1=nearest}').");
    cimg_option("-correlate","'[indice][,border_conditions]'",
                "Correlate image(s) by the specified mask (same parameters as above).");
    cimg_option("-erode","'size[,border_conditions]' or '[indice][,border_conditions]'","");
    cimg_help("                                              "
              "Erode image(s) by the specified mask (same parameters as above)').");
    cimg_option("-dilate","'size[,border_conditions]' or '[indice][,border_conditions]'","");
    cimg_help("                                              "
              "Dilate image(s) by the specified mask (same parameters as above).");
    cimg_option("-gradient","'x', 'xy', 'xyz' or (no args)","Compute image gradient.");
    cimg_option("-hessian","'{xx,xy,xz,yy,yz,zz}' or (no args)","Compute image Hessian.");
    cimg_option("-fft","(no args)","Compute direct Fourier transform.");
    cimg_option("-ifft","(no args)","Compute inverse Fourier transform.");

    cimg_help("\n  Image creation and drawing\n"
              "  --------------------------");
    cimg_option("-dimensions","(no args)","Get dimensions of the image(s) as a 1x4 vector.");
    cimg_option("-stats","(no args)","Get statistics of the image(s) as a 1x6 vector.");
    cimg_option("-histogram","'nb_values[%]'","Compute histogram of image(s) with 'nb_values' values.");
    cimg_option("-distance","'isovalue'","Compute distance function(s) to specified isovalue.");
    cimg_option("-hamilton","'nb_iter[,band_size]'","Apply Hamilton-Jacobi PDE to compute distance to 0.");
    cimg_option("-label","(no args)","Label connected components of image(s).");
    cimg_option("-displacement","'[indice][,smoothness[,precision[,nbscales[,itermax]]]]","Estimate smooth displacement field between image(s) "
                "and specified target '[indice]'.");
    cimg_option("-sort","(no args)","Sort values of image(s) in increasing order.");
    cimg_option("-psnr","'max_value' or (noargs)","Compute PSNR between specified image(s).");
    cimg_option("-point","'x[%],y[%][,z[%][,opacity[,color]]]'","Draw 3D colored point on specified image(s).");
    cimg_option("-line","'x0[%],y0[%],x1[%],y1[%][,opacity[,color]]'","Draw 2D colored line on specified image(s).");
    cimg_option("-polygon","'N,x0[%],y0[%],..,xN[%],yN[%][,opacity[,color]]'","Draw a 2D colored N-vertices polygon on specified image(s).");
    cimg_option("-ellipse","'x[%],y[%],r,R[,u,v[,opacity[,color]]]'","Draw 2D colored ellipse on specified image(s).");
    cimg_option("-text","text,x[%],y[%],size[,opacity[,color]]'",
                "Draw specified text at position (x,y) and with specified font size.");
    cimg_option("-image","'[indice][,x[%][,y[%][,z[%][,opacity[,ind_mask]]]]]'","Draw sprite image on specified image(s).");
    cimg_option("-object3d","'[indice][,x[%][,y[%][,z[,opacity]]]]'","Draw 3D object on specified image(s).");
    cimg_option("-plasma","'alpha[,beta[,opacity]]'","Draw plasma on specified image(s).");
    cimg_option("-mandelbrot","'z0r,z0i,z1r,z1i[,itermax[,julia,c0r,c0i[,opacity]]]'","Draw Mandelbrot/Julia fractals on specified image(s).");
    cimg_option("-flood","'x[%][,y[%][,z[%][,tolerance[,opacity[,color]]]]]'",
                "Flood-fill image(s) starting from (x,y,z) with specified tolerance.");

    cimg_help("\n  List manipulation\n"
              "  -----------------");
    cimg_option("-remove","(no args)","Remove image(s) from list (eq. to '-rm').");
    cimg_option("-keep","(no args)","Keep only specified image(s) (eq. to '-k').");
    cimg_option("-move","'position'","Move image(s) at specified position (eq. to '-mv').");
    cimg_option("-reverse","(no args)","Reverse image(s) order.");
    cimg_option("-name","\"name\"","Set display name of image(s).");

    cimg_help("\n  3D Rendering\n"
              "  ------------");
    cimg_option("-cube3d","'size'","Insert a 3D cube at the end of the list.");
    cimg_option("-cone3d","'radius[,height[,subdivisions]]'","Insert a 3D cube at the end of the list.");
    cimg_option("-cylinder3d","'radius[,height[,subdivisions]]'","Insert a 3D cylinder at the end of the list.");
    cimg_option("-torus3d","'radius1,radius2[,subdivisions1,subdivisions2]'","Insert a 3D torus at the end of the list.");
    cimg_option("-plane3d","'sizex,sizey[,subdivisionsx,subdisivionsy]'","Insert a 3D plane at the end of the list.");
    cimg_option("-sphere3d","'radius[,subdivisions]'","Insert a 3D sphere at the end of the list.");
    cimg_option("-elevation3d","'z-factor' or '[indice]'",
                "Generate 3D elevation(s) of image(s) using specified z-factor or elevation map.");
    cimg_option("-isovalue3d","'value'","Generate 3D object(s) by retrieving isophote or isosurface of image(s).");
    cimg_option("-center3d","(no args)","Center 3D object(s) (eq. to '-c3d').");
    cimg_option("-normalize3d","(no args)","Normalize 3D object(s) to a unit size (eq. to '-n3d').");
    cimg_option("-rotate3d","'u,v,w,angle'","Rotate 3D object(s) around axis (u,v,w) with specified angle (eq. to '-rot3d').");
    cimg_option("-add3d","'[indice]' or 'tx,ty,tz' or (noargs)","Append or translate 3D object(s), or append 3D object(s) together (eq. to '-+3d').");
    cimg_option("-sub3d","'tx,ty,tz'","Translate 3D object(s) with the opposite of the specified vector (eq. to '--3d').");
    cimg_option("-mul3d","'fact' or 'factx,facty[,factz]'","Scale 3D object(s) with specified factor (eq. to '-*3d').");
    cimg_option("-div3d","'fact' or 'factx,facty[,factz]'","Scale 3D object(s) with specified inverse factor (eq. to '-/3d').");
    cimg_option("-color3d","'R,G,B[,opacity]'","Set color of 3D object(s) (eq. to '-col3d').");
    cimg_option("-opacity3d","'opacity'","Set opacity of 3D object(s) (eq. to '-opac3d').");
    cimg_option("-invert3d","(no args)","Invert primitive orientations of 3D object(s) (eq. to '-i3d').");
    cimg_option("-split3d","(no args)","Split 3D object data into 6 data vectors 'header,N,vertices,primitives,colors,opacities' (eq. to '-s3d').");
    cimg_option("-light3d","'posx,posy,posz'","Set the 3D position of the light for 3D rendering (eq. to '-l3d').");
    cimg_option("-focale3d","'value'","Set focale value for 3D rendering (eq. to '-f3d').");
    cimg_option("-specl3d","'value'","Set amount of specular light for 3D rendering (eq. to '-sl3d').");
    cimg_option("-specs3d","'value'","Set shininess of specular light for 3D rendering (eq. to '-ss3d').");
    cimg_option("-orient3d","(no args)","Switch double-sided mode for 3D rendering (eq. to '-o3d').");
    cimg_option("-render3d","'mode'","Set 3D rendering mode");
    cimg_help("                                              "
              "(can be '{-1=bounding-box, 0=pointwise, 1=linear, 2=flat, 3=flat-shaded, 4=Gouraud-shaded, 5=Phong-shaded}') (eq. to '-r3d').");
    cimg_option("-renderd3d","'mode'","Set dynamic rendering mode in 3D viewer (same values as above) (eq. to '-rd3d').");
    cimg_option("-background3d","'R,G,B'","Define background color in 3D viewer (eq. to '-b3d').");

    cimg_help("\n  Program controls\n"
              "  ----------------");
    cimg_option("-nop","(no args)","Do nothing.");
    cimg_option("-skip","(any args)","Do nothing but skip the next argument.");
    cimg_option("-echo","'text'","Output specified message (eq. to '-e').");
    cimg_option("-print","(no args)","Print image(s) informations (eq. to '-p').");
    cimg_option("-quit","(no args)","Force interpreter to quit (eq. to '-q').");
    cimg_option("-do","(no args)","Start a 'do-while' code bloc.");
    cimg_option("-while","'cond'","End a 'do-while' code bloc and go back to associated '-do' if 'cond' is a strictly positive value.");
    cimg_option("-if","'cond'","Start a 'if-then-else' code bloc and test if 'cond' is a strictly positive value.");
    cimg_option("-else","(no args)","Execute following commands if previous '-if' condition failed.");
    cimg_option("-endif","(no args)","End a 'if-then-else' code bloc");
    cimg_option("-repeat","'N'","Start a 'repeat-done' code bloc.");
    cimg_option("-done","(no args)","End a 'repeat-done' code bloc, and go to associated '-repeat' if iterations remain.");
    cimg_option("-int","'arg1,...,argN'","Check if all specified arguments are integer. If not, print an error message and exit.");
    cimg_option("-float","'arg1,...,argN","Check if all specified arguments are float values. If not, print an error message and exit.");

    cimg_help("\n  Input/output\n"
              "  ------------");
    cimg_option("-input","'filename' or 'width[%][xheight[%][xdepth[%][xdim[%][xN]]]]'","");
    cimg_help("                     "
              "or '[indice][xN]' or '(v11{,;/^}v21...vLM)[xN]'");
    cimg_help("                                              "
              "Input filename, empty image, image copy, or image with specified values (eq. to '-i' or (no args)).");
    cimg_option("-output","'filename'","Output image(s) in specified filename (eq. to '-o').");
    cimg_option("-display","(no args)","Display image(s) (eq. to '-d').");
    cimg_option("-display3d","(no args)","Display 3D object(s) (eq. to '-d3d').");
    cimg_option("-plot","'[plot_type[,vertex_type[,xmin[,xmax[,ymin[,ymax]]]]]]'",
                "Display image(s) as 1D plot(s)");
    cimg_help("                                              "
              "('plot_type' can be '{0=none, 1=lines, 2=splines, 3=bar}').");
    cimg_help("                                              "
              "('vertex_type' can be in '[0-7]').");
    cimg_option("-select","'select_type'","Select feature from image(s) in an interactive way");
    cimg_help("                                              "
              "('select_type' can be in '{0=point, 1=line, 2=rectangle, 3=circle').");

    // Print descriptions of default macros.
    char line[256*1024] = { 0 }, name[4096] = { 0 }, args[4096] = { 0 }, desc[4096] = { 0 };
    bool first_description = true;
    for (const char *data = data_def; *data; ) {
      if (*data=='\n') ++data;
      else {
        if (std::sscanf(data,"%262143[^\n]",line)>0) data += cimg::strlen(line);
        if (line[0]=='#' && std::sscanf(line,"#@gmic %4095[^:] : %4095[^:] : %4095[^:]",name,args,desc)>0) {
          if (first_description) cimg_help("\n  Commands : Default macros\n"
                                           "  -------------------------");
          std::fprintf(cimg_stdout,"%s    %s-%-15s%s %-24s %s%s%s",
                       first_description?"":"\n",
                       cimg::t_bold,name,cimg::t_normal,args,cimg::t_green,desc,cimg::t_normal);
          first_description = false;
        }
      }
    }

    // Print descriptions of use-defined macros.
    first_description = true;
    for (int i = 1; i<argc-1; ++i) if (!cimg::strcmp("-m",argv[i]) || !cimg::strcmp("-macros",argv[i])) {
      std::FILE *file = cimg::fopen(argv[i+1],"r");
      if (file) {
        int err = 0;
        while ((err=std::fscanf(file,"%262143[^\n] ",line)>=0)) {
          if (err) { // Non empty-line
            name[0] = args[0] = desc[0] = 0;
            if (line[0]=='#' && std::sscanf(line,"#@gmic %4095[^:] : %4095[^:] : %4095[^:]",name,args,desc)>0) {
              if (first_description) cimg_help("\n\n  Commands : User-defined macros\n"
                                               "  ------------------------------");
              std::fprintf(cimg_stdout,"%s    %s-%-15s%s %-24s %s%s%s",
                           first_description?"":"\n",
                           cimg::t_bold,name,cimg::t_normal,args,cimg::t_green,desc,cimg::t_normal);
              first_description = false;
            }
          }
        }
      }
      cimg::fclose(file);
    }

    cimg_help("\n\n  Viewers shortcuts\n"
              "  -----------------");
    cimg_help("  When displaying image(s) or 3D object(s) with G'MIC, you can use these shortcuts in viewers :");
    cimg_help("   - CTRL+D : Increase window size.");
    cimg_help("   - CTRL+C : Decrease window size.");
    cimg_help("   - CTRL+R : Reset window size.");
    cimg_help("   - CTRL+F : Toggle fullscreen mode.");
    cimg_help("   - CTRL+S : Save current window snapshot.");
    cimg_help("   - CTRL+O : Save current instance of viewed image or 3D object.\n");
    cimg_help("  Specific options for the viewer of image(s) :");
    cimg_help("   - CTRL+P             : Play stack of frames as a movie.");
    cimg_help("   - CTRL+(mousewheel)  : Zoom in/out.");
    cimg_help("   - SHIFT+(mousewheel) : Go left/right.");
    cimg_help("   - ALT+(mousewheel)   : Go up/down.");
    cimg_help("   - Numeric PAD        : Zoom in/out (+/-) and move zoomed region (numbers).");
    cimg_help("   - BACKSPACE          : Reset zoom.\n");
    cimg_help("  Specific options for the viewer of 3D object(s) :");
    cimg_help("   - (mouse) + (left mouse button)   : Rotate object.");
    cimg_help("   - (mouse) + (right mouse button)  : Zoom object.");
    cimg_help("   - (mouse) + (middle mouse button) : Translate object.");
    cimg_help("   - (mousewheel)                    : Zoom in/out.");
    cimg_help("   - CTRL + Z : Enable/disable Z-buffer rendering");

    cimg_help("\n  File options\n"
              "  ------------");
    cimg_help("  G'MIC is able to read/write most of the classical image file formats, including :");
    cimg_help("   - 2D grayscale/color images : PNG,JPEG,GIF,PNM,TIFF,BMP,....");
    cimg_help("   - 3D volumetric images : DICOM,HDR,NII,PAN,CIMG,INR,....");
    cimg_help("   - Video files : MPEG,AVI,MOV,OGG,FLV,...");
    cimg_help("   - Generic data files : DLM,ASC,RAW,TXT.");
    cimg_help("   - 3D objects : OFF.\n");
    cimg_help("  Specific options :");
    cimg_help("   - For video files : you can read only sub-frames of the image sequence (recommended) with the expression");
    cimg_help("     'video.ext,[first_frame[%][,last_frame[%][,step]]]'.");
    cimg_help("   - For RAW files : you must specify the image dimensions with the expression");
    cimg_help("     'file.raw,width[,height[,depth[,dim]]]]'.");
    cimg_help("   - For YUV files : you must specify the image dimensions and can read only sub-frames of the image sequence with the expression");
    cimg_help("     'file.yuv,width,height[,first_frame[,last_frame[,step]]]'.");
    cimg_help("   - For JPEG files : you can specify the quality (in %) of an output jpeg file format with the expression");
    cimg_help("     'file.jpg,30%'.");

    cimg_help("\n  Examples of use\n"
              "  ---------------");
    cimg_help("  G'MIC is a simple but quite complete interpreter of image processing instructions, and can be used for a wide variety of");
    cimg_help("  image processing tasks. Here are (very few) examples of how the command line tool G'MIC can be used :\n");
    cimg_help("   - View image data : ");
    cimg_help("     gmic file1.bmp file2.jpeg");
    cimg_help("   - Convert image files : ");
    cimg_help("     gmic input.bmp -o output.jpg");
    cimg_help("   - Create volumetric image(s) from movie sequence : ");
    cimg_help("     gmic input.mpg -a z -o output.hdr");
    cimg_help("   - Compute image gradient norm : ");
    cimg_help("     gmic input.bmp -gradient_norm");
    cimg_help("   - Create G'MIC 3D logo : ");
    cimg_help("     gmic 180x70x1x3 -text G\\'MIC,30,5,50,1,1 -blur 2 -n 0,100 [0] -plasma[1] \\");
    cimg_help("     -+ -blur 1 -elevation -0.1 -rd3d 4");
    cimg_help("\n  See also the macros defined in the provided macro file 'gmic_def.raw' for other examples.");

    cimg_help("\n  ** G'MIC comes with ABSOLUTELY NO WARRANTY; "
              "for details visit http://gmic.sourceforge.net **");
    std::exit(0);
  }

  // Launch G'MIC instance.
  //-----------------------
  CImgList<float> images;
  try {
    gmic(argc,argv,images);
  } catch (gmic_exception &e) {
    std::fprintf(cimg_stdout,"\n** Error occurred : %s **\n",e.message);
  }
  return 0;
}
#endif

#endif // #ifdef cimg_plugin ... #else ...
