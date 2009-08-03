#
#  File        : CMakeLists.txt
#                ( Cmake configuration file )
#
#  Description : CMake can use this to generate a working Makefile.
#                This file is a part of the CImg Library project.
#                ( http://cimg.sourceforge.net )
#
#  Copyright   : Baptiste Mougel
#                ( http://www.mougel.org )
#                David Tschumperle
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

#-------------------
# IMPORTANT NOTICE : This file is intended to be copied in the 'examples/' directory before use.
#-------------------
project(CIMG)

IF(UNIX)
   CMAKE_MINIMUM_REQUIRED (VERSION 2.4.0 FATAL_ERROR)
   MESSAGE(STATUS "Compilation mode : UNIX")
ELSE(UNIX)
  IF(WIN32)
    CMAKE_MINIMUM_REQUIRED(VERSION 2.4.6 FATAL_ERROR)
    MESSAGE(STATUS "Compilation mode : Win32")
  ELSE(WIN32)
    IF(APPLE)
       CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0 FATAL_ERROR)
       IF(${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} EQUAL 4 AND ${CMAKE_PATCH_VERSION} LESS 7)
         MESSAGE("Warning: A critical CMake bug exists in 2.4.6 and below. Trying to build Universal Binaries will result in a compile error that seems unrelated. Either avoid building Universal Binaries by changing the CMAKE_OSX_ARCHITECTURES field to list only your architecture, or upgrade to the current CVS version of CMake or a newer stable version if it exists.")
       ENDIF(${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} EQUAL 4 AND ${CMAKE_PATCH_VERSION} LESS 7)
       MESSAGE(STATUS "Compilation mode : MAC OSX")
    ELSE(APPLE)
        CMAKE_MINIMUM_REQUIRED(VERSION 2.4.0 FATAL_ERROR)
    ENDIF(APPLE)
  ENDIF(WIN32)
ENDIF(UNIX)

#-----------------------------------------------------------------------------
# CImg version number.  An even minor number corresponds to releases.
SET(CIMG_VERSION "1.3.0")
SET(CIMG_H "CImg.h")

MESSAGE(STATUS "CImg version : ${CIMG_VERSION}")

#-----------------------------------------------------------------------------
# Output directories.
INCLUDE_DIRECTORIES (
  ${CIMG_BINARY_DIR}
  ${CIMG_SOURCE_DIR}
)

IF(NOT LIBRARY_OUTPUT_PATH)
   SET (LIBRARY_OUTPUT_PATH ${CIMG_BINARY_DIR}/ CACHE INTERNAL "Single output directory for building all libraries.")
ENDIF(NOT LIBRARY_OUTPUT_PATH)
IF(NOT EXECUTABLE_OUTPUT_PATH)
   SET (EXECUTABLE_OUTPUT_PATH ${CIMG_BINARY_DIR}/ CACHE INTERNAL "Single output directory for building all executables.")
ENDIF(NOT EXECUTABLE_OUTPUT_PATH)

#MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

#-----------------------------------------------------------------------------
# Compilation.

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -W -ansi -pedantic")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -W -ansi -pedantic")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -lpthread")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lm -lpthread")
SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -lm -lpthread")

OPTION(BUILD_EXAMPLES "Build the Examples directory." ON)
IF (BUILD_EXAMPLES)
	ADD_EXECUTABLE(CImg_demo CImg_demo.cpp)
	ADD_EXECUTABLE(captcha captcha.cpp)
	ADD_EXECUTABLE(curve_editor curve_editor.cpp)
	ADD_EXECUTABLE(dtmri_view dtmri_view.cpp)
	ADD_EXECUTABLE(edge_explorer edge_explorer.cpp)
	ADD_EXECUTABLE(fade_images fade_images.cpp)
	ADD_EXECUTABLE(generate_loop_macros generate_loop_macros.cpp)
        ADD_EXECUTABLE(gmic gmic.cpp)
	ADD_EXECUTABLE(greycstoration greycstoration.cpp)
	ADD_EXECUTABLE(hough_transform hough_transform.cpp)
	ADD_EXECUTABLE(image2ascii image2ascii.cpp)
	ADD_EXECUTABLE(image_registration image_registration.cpp)
	ADD_EXECUTABLE(image_surface image_surface.cpp)
	ADD_EXECUTABLE(jawbreaker jawbreaker.cpp)
	ADD_EXECUTABLE(mcf_levelsets3d mcf_levelsets3d.cpp)
	ADD_EXECUTABLE(mcf_levelsets mcf_levelsets.cpp)
	ADD_EXECUTABLE(odykill odykill.cpp)
	ADD_EXECUTABLE(pde_heatflow2d pde_heatflow2d.cpp)
	ADD_EXECUTABLE(pde_TschumperleDeriche2d pde_TschumperleDeriche2d.cpp)
	ADD_EXECUTABLE(radon_transform radon_transform.cpp)
	ADD_EXECUTABLE(scene3d scene3d.cpp)
	ADD_EXECUTABLE(tetris tetris.cpp)
	ADD_EXECUTABLE(tron tron.cpp)
	ADD_EXECUTABLE(tutorial tutorial.cpp)
	ADD_EXECUTABLE(wavelet_atrous wavelet_atrous.cpp)
	ADD_EXECUTABLE(use_draw_gradient use_draw_gradient.cpp)
	ADD_EXECUTABLE(use_greycstoration use_greycstoration.cpp)
	ADD_EXECUTABLE(use_nlmeans use_nlmeans.cpp)
	ADD_EXECUTABLE(use_RGBclass use_RGBclass.cpp)
	ADD_EXECUTABLE(use_skeleton use_skeleton.cpp)
ENDIF (BUILD_EXAMPLES)

# Flags to enable code debugging.
OPTION(CIMG_DEBUG_FLAGS "Flags to enable code debugging." OFF)
IF (CIMG_DEBUG_FLAGS)
  ADD_DEFINITIONS(-Dcimg_debug=3)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
ENDIF(CIMG_DEBUG_FLAGS)

# Flags to enable color output messages.
# (requires a VT100 compatible terminal)
OPTION(CIMG_VT100_FLAGS "Flags to enable color output messages. (requires a VT100 compatible terminal)" ON)
IF (CIMG_VT100_FLAGS)
  ADD_DEFINITIONS(-Dcimg_use_vt100)
ENDIF(CIMG_VT100_FLAGS)

# Flags to enable code optimization by the compiler.
OPTION(CIMG_OPT_FLAGS "Flags to enable code optimization by the compiler." ON)
IF (CIMG_OPT_FLAGS)
   SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ffast-math")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math")
ENDIF(CIMG_OPT_FLAGS)

# Flags to enable OpenMP support.
OPTION(CIMG_OPENMP_FLAGS "Flags to enable code optimization by the compiler." OFF)
IF (CIMG_OPENMP_FLAGS)
   SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -fopenmp")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fopenmp")
   ADD_DEFINITIONS(-Dcimg_use_openmp)
ENDIF(CIMG_OPENMP_FLAGS)

# Flags to enable the use of the X11 library.
# (X11 is used by CImg to handle display windows)
# !!! For 64bits systems : replace -L$(X11PATH)/lib by -L$(X11PATH)/lib64 !!!
OPTION(CIMG_X11_FLAGS "Flags to enable fast image display, using the XSHM library." ON)
IF (CIMG_X11_FLAGS)

  FIND_PATH(X11_PATH /bin/xterm /usr/X11R6)

  OPTION(CIMG_SYSTEMS_64BITS "For 64bits systems." OFF)
  IF (CIMG_SYSTEMS_64BITS)
    INCLUDE_DIRECTORIES(${X11_PATH}/include)
    LINK_DIRECTORIES(${X11_PATH}/lib64)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${X11_PATH}/lib64 -lX11")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -L${X11_PATH}/lib64 -lX11")
    SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -L${X11_PATH}/lib64 -lX11")
  ELSE (CIMG_SYSTEMS_64BITS)
    INCLUDE_DIRECTORIES(${X11_PATH}/include)
    LINK_DIRECTORIES(${X11_PATH}/lib)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${X11_PATH}/lib -lX11")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -L${X11_PATH}/lib -lX11")
    SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -L${X11_PATH}/lib -lX11")
  ENDIF(CIMG_SYSTEMS_64BITS)

  # Flags to enable fast image display, using the XSHM library.
  # (CIMG_X11_FLAGS must be defined also)
  OPTION(CIMG_XSHM_FLAGS "Flags to enable fast image display, using the XSHM library." ON)
  IF (CIMG_XSHM_FLAGS)
    ADD_DEFINITIONS(-Dcimg_use_xshm)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lXext")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lXext")
    SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -lXext")
  ENDIF(CIMG_XSHM_FLAGS)

  # Flags to enable screen mode switching, using the XRandr library.
  # (CIMG_X11_FLAGS must be defined also)
  OPTION(CIMG_XRANDR_FLAGS "Flags to enable screen mode switching, using the XRandr library." ON)
  IF (CIMG_XRANDR_FLAGS)
    ADD_DEFINITIONS(-Dcimg_use_xrandr)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lXrandr")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lXrandr")
    SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -lXrandr")
  ENDIF(CIMG_XRANDR_FLAGS)

  # Flags to compile on Solaris
  #CIMG_SOLARIS_FLAGS = -R$(X11PATH)/lib -lrt -lnsl -lsocket
  OPTION(CIMG_SOLARIS_FLAGS "Flags to compile on Solaris." OFF)
  IF (CIMG_SOLARIS_FLAGS)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -R$(X11_PATH)/lib -lrt -lnsl -lsocket")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -R$(X11_PATH)/lib -lrt -lnsl -lsocket")
    SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -R$(X11_PATH)/lib -lrt -lnsl -lsocket")
  ENDIF(CIMG_SOLARIS_FLAGS)

ENDIF(CIMG_X11_FLAGS)


# Flags to enable native support for PNG image files, using the PNG library.
#CIMG_PNG_FLAGS = -Dcimg_use_png -lpng -lz
OPTION(CIMG_PNG_FLAGS "Flags to enable native support for PNG image files, using the PNG library." OFF)
IF (CIMG_PNG_FLAGS)
  ADD_DEFINITIONS(-Dcimg_use_png)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpng -lz")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lpng -lz")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -lpng -lz")
ENDIF(CIMG_PNG_FLAGS)

# Flags to enable native support for JPEG image files, using the JPEG library.
#CIMG_JPEG_FLAGS = -Dcimg_use_jpeg -ljpeg
OPTION(CIMG_JPEG_FLAGS "Flags to enable native support for JPEG image files, using the JPEG library." OFF)
IF (CIMG_JPEG_FLAGS)
  ADD_DEFINITIONS(-Dcimg_use_jpeg)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ljpeg")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -ljpeg")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -ljpeg")
ENDIF(CIMG_JPEG_FLAGS)

# Flags to enable native support for TIFF image files, using the TIFF library.
#CIMG_TIFF_FLAGS = -Dcimg_use_tiff -ltiff
OPTION(CIMG_TIFF_FLAGS "Flags to enable native support for TIFF image files, using the TIFF library." OFF)
IF (CIMG_TIFF_FLAGS)
  ADD_DEFINITIONS(-Dcimg_use_tiff)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltiff")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -ltiff")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -ltiff")
ENDIF(CIMG_TIFF_FLAGS)

# Flags to enable native support for various video files, using the FFMPEG library.
#CIMG_FFMPEG_FLAGS = -Dcimg_use_ffmpeg -lavcodec -lavformat
OPTION(CIMG_FFMPEG_FLAGS "Flags to enable native support for various video files, using the FFMPEG library." OFF)
IF (CIMG_FFMPEG_FLAGS)
  ADD_DEFINITIONS(-Dcimg_use_ffmpeg)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lavcodec -lavformat")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lavcodec -lavformat")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -lavcodec -lavformat")
ENDIF(CIMG_FFMPEG_FLAGS)

# Flags to enable native support for compressed .cimgz files, using the Zlib library.
#CIMG_ZLIB_FLAGS = -Dcimg_use_zlib -lz
OPTION(CIMG_ZLIB_FLAGS "Flags to enable native support for compressed .cimgz files, using the Zlib library." OFF)
IF (CIMG_ZLIB_FLAGS)
  ADD_DEFINITIONS(-Dcimg_use_zlib)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lz")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lz")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -lz")
ENDIF(CIMG_ZLIB_FLAGS)

# Flags to enable native support of most classical image file formats, using the Magick++ library.
#CIMG_MAGICK_FLAGS = -Dcimg_use_magick `Magick++-config --cppflags` `Magick++-config --cxxflags` `Magick++-config --ldflags` `Magick++-config --libs`
OPTION(CIMG_MAGICK_FLAGS "Flags to enable native support of most classical image file formats, using the Magick++ library." OFF)
IF (CIMG_MAGICK_FLAGS)
  ADD_DEFINITIONS(-Dcimg_use_magick)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} `Magick++-config --cppflags` `Magick++-config --cxxflags`")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} `Magick++-config --cppflags` `Magick++-config --cxxflags`")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} `Magick++-config --ldflags` `Magick++-config --libs`")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} `Magick++-config --ldflags` `Magick++-config --libs`")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} `Magick++-config --ldflags` `Magick++-config --libs`")
ENDIF(CIMG_MAGICK_FLAGS)

# Flags to enable faster Discrete Fourier Transform computation, using the FFTW3 library
#CIMG_FFTW3_FLAGS = -Dcimg_use_fftw3 -lfftw3
OPTION(CIMG_FFTW3_FLAGS "Flags to enable faster Discrete Fourier Transform computation, using the FFTW3 library." OFF)
IF (CIMG_FFTW3_FLAGS)
  ADD_DEFINITIONS(-Dcimg_use_fftw3)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lfftw3")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lfftw3")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -lfftw3")
ENDIF(CIMG_FFTW3_FLAGS)

# Flags to enable the use of LAPACK routines for matrix computation
#CIMG_LAPACK_FLAGS = -Dcimg_use_lapack -lblas -lg2c -llapack
OPTION(CIMG_LAPACK_FLAGS "Flags to enable the use of LAPACK routines for matrix computation" OFF)
IF (CIMG_LAPACK_FLAGS)
  ADD_DEFINITIONS(-Dcimg_use_lapack)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lblas -lg2c -llapack")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lblas -lg2c -llapack")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -lblas -lg2c -llapack")
ENDIF(CIMG_LAPACK_FLAGS)

# Flags to enable the use of the Board library
#CIMG_BOARD_FLAGS = -Dcimg_use_board -I/usr/include/board -lboard
OPTION(CIMG_BOARD_FLAGS "Flags to enable the use of the Board library." OFF)
IF (CIMG_BOARD_FLAGS)
  ADD_DEFINITIONS(-Dcimg_use_board)
  INCLUDE_DIRECTORIES(/usr/include/board)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lboard")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lboard")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -lboard")
ENDIF(CIMG_BOARD_FLAGS)

# Flags to compile on MacOSX with Carbon-based display support
#CIMG_CARBON_FLAGS = -Dcimg_display=3 -framework Carbon
IF(APPLE)
  OPTION(CIMG_CARBON_FLAGS "Flags to compile on MacOSX with Carbon-based display support." ON)
  IF (CIMG_CARBON_FLAGS)
    ADD_DEFINITIONS(-Dcimg_display=3)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Carbon")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -framework Carbon")
    SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -framework Carbon")
  ENDIF(CIMG_CARBON_FLAGS)
ENDIF(APPLE)
