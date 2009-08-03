@echo off
Rem Intel(R) C++ Compiler Build Environment for 32-bit applications

echo.
echo Intel(R) C++ Compiler 8.0 Build Environment for 32-bit applications
echo Copyright (C) 1985-2003 Intel Corporation. All rights reserved.
echo.

@call "C:\Program Files\Microsoft Visual Studio .NET 2003\Vc7\Bin\Vcvars32.bat"

echo.

SET INTEL_COMPILER80=C:\Program Files\Intel\CPP\Compiler80
SET INTEL_SHARED=C:\Program Files\Fichiers communs\Intel\Shared Files
SET INTEL_LICENSE_FILE=C:\Program Files\Fichiers communs\Intel\Licenses
SET PATH=%INTEL_COMPILER80%\Ia32\Bin;%INTEL_SHARED%\Ia32\Bin;%PATH%
SET LIB=%INTEL_COMPILER80%\Ia32\Lib;%INTEL_SHARED%\Ia32\Lib;%LIB%
SET INCLUDE=%INTEL_COMPILER80%\Ia32\Include;%INCLUDE%

cd ..\examples\
echo.
echo ** Compiling 'CImg_demo'
echo.
icl /GX /Ox CImg_demo.cpp gdi32.lib user32.lib
