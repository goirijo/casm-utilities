#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([CASM Utilities], [0.0.0], [jgg@umail.ucsb.edu])
AC_CONFIG_SRCDIR([src/structure.cxx])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([build-aux/m4])
AM_INIT_AUTOMAKE([subdir-objects])

#We'll be making shared libraries that call shared libraries
LT_INIT

# Checks for programs.
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11

#Is this how you get it to find the python related things for compiling the wrapper?
AM_PATH_PYTHON
AX_PYTHON
# AX_PYTHON_DEVEL(=='3.5')
AX_PYTHON_DEVEL

# Checks for libraries.
BOOSTV=1.52
AX_BOOST_BASE([$BOOSTV],[],[AC_MSG_ERROR(You need Boost $BOOSTV libraries or higher.)])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_PYTHON


AC_CHECK_LIB([casm],[main], [], [AC_MSG_ERROR(The CASM libraries are required!)])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
