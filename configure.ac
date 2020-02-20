#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([CASM Utilities], [0.0.0], [jgg@ucsb.edu])
AC_CONFIG_SRCDIR([src/casmutils/structure.cxx])
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

# AX_BOOST_PYTHON

#You will need pybind11 for this operation
#Why won't it find the libraries???
# AC_CHECK_HEADER([pybind11/pybind11.h],[],[AC_MSG_ERROR(Pybind11 is required for this operation)])
AC_CHECK_LIB([dl], [dlopen], [], AC_MSG_ERROR(dlopen from dl library not found!))
AX_CHECK_ZLIB(,[AC_MSG_ERROR([Could not find zlib])])


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
