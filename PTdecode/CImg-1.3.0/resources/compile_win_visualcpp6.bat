REM A dirty batch file to compile all CImg examples
REM Using the microsoft's CL compiler.
REM -----------------------------------------------
set SDKPATH=C:\Program Files\Microsoft Platform SDK for Windows Server 2003 R2
set CPP=cl /W3 /Zm800 /Ox /Ob2 /Oi /Ot /c /EHsc /I"%SDKPATH%\Include" /I"..\\"
set LDD=link /LIBPATH:"%SDKPATH%\Lib"
set LDDLIBS=user32.lib gdi32.lib shell32.lib

cd ..\examples\
set CPPFILE=CImg_demo captcha curve_editor dtmri_view edge_explorer fade_images generate_loop_macros greycstoration hough_transform image2ascii image_registration image_surface jawbreaker mcf_levelsets3D mcf_levelsets odykill pde_heatflow2D pde_TschumperleDeriche2D radon_transform scene3d tetris tron tutorial wavelet_atrous use_draw_gradient use_greycstoration use_nlmeans use_RGBclass use_skeleton gmic
FOR %%F IN (%CPPFILE%) DO (
  %CPP% %%F.cpp
  %LDD% %%F.obj %LDDLIBS%
)





