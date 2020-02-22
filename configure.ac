#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([CASM Utilities], [0.0.0], [jgg@ucsb.edu])
AC_CONFIG_SRCDIR([src/primify/primify.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([build-aux/m4])
AM_INIT_AUTOMAKE([subdir-objects])

#We'll be making shared libraries that call shared libraries    #You sure about that chief? I think we're doing static. Try deleting this
LT_INIT

# Checks for programs.
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_17

#Is this how you get it to find the python related things for compiling the wrapper?
# AX_PYTHON_DEVEL(=='3.5')
AX_PYTHON_DEVEL
AX_PYTHON
AM_PATH_PYTHON

# Checks for libraries.
BOOSTV=1.52
AX_BOOST_BASE([$BOOSTV],[],[AC_MSG_ERROR(You need Boost $BOOSTV libraries or higher.)])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_PROGRAM_OPTIONS

#If you're compiling casm for them, you need these also
AX_BOOST_CHRONO
AX_BOOST_REGEX

AC_CHECK_LIB([dl], [dlopen], [], AC_MSG_ERROR(dlopen from dl library not found!))
AX_CHECK_ZLIB(,[AC_MSG_ERROR([Could not find zlib])])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
