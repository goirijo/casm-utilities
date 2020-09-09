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
AX_CHECK_ZLIB([],AC_MSG_ERROR(zlib library not found!))
AX_PTHREAD([],AC_MSG_WARN(pthread is required to run tests!))

AC_SEARCH_LIBS(dlopen, dl, [], AC_MSG_ERROR(dl library not found!))
# AC_SEARCH_LIBS(std::filesystem::status_known,stdc++fs,[], AC_MSG_ERROR(stdc++fs not found!),[])

# Why is this such a pain in the butt? I just want to -lstdc++fs, I could be doing this in the
# Makefile, but instead I'm following this example:
# https://github.com/LibreOffice/online/blob/master/configure.ac
AC_MSG_CHECKING(for std::filesystem features in -lstdc++fs)
STD_FILESYSTEM=
    save_CXXFLAGS=$CXXFLAGS
    CXXFLAGS="$CXXFLAGS -Werror"
    save_LIBS=$LIBS
    LIBS="$save_LIBS -lstdc++fs"
    AC_LANG_PUSH([C++])
    AC_LINK_IFELSE([AC_LANG_SOURCE([[
            #include <filesystem>
            int main()
            {
                if (!std::filesystem::temp_directory_path().empty())
                    return 0;
                return 1;
            }
            ]])],[STD_FILESYSTEM=TRUE], AC_MSG_ERROR(stdc++fs not found!))
    AC_LANG_POP([C++])
    CXXFLAGS=$save_CXXFLAGS
    LIBS=$save_LIBS

if test -n "$STD_FILESYSTEM" ; then
    LIBS="$LIBS -lstdc++fs"
    AC_MSG_RESULT(yes)
else
    AC_MSG_ERROR("libstdc++fs library not found!)
fi

AC_CONFIG_FILES([Makefile])

##############Prepares all tests scripts related to python modules################# 
AC_CONFIG_FILES([tests/py/stage.py],[chmod +x tests/py/stage.py])
AC_CONFIG_FILES([tests/py/casmutils/xtal/coordinate.py],[chmod +x tests/py/casmutils/xtal/coordinate.py])
AC_CONFIG_FILES([tests/py/casmutils/xtal/site.py],[chmod +x tests/py/casmutils/xtal/site.py])
AC_CONFIG_FILES([tests/py/casmutils/xtal/lattice.py],[chmod +x tests/py/casmutils/xtal/lattice.py])
AC_CONFIG_FILES([tests/py/casmutils/xtal/structure.py],[chmod +x tests/py/casmutils/xtal/structure.py])
AC_CONFIG_FILES([tests/py/casmutils/xtal/structure_tools.py],[chmod +x tests/py/casmutils/xtal/structure_tools.py])
AC_CONFIG_FILES([tests/py/casmutils/mapping/mapping.py],[chmod +x tests/py/casmutils/mapping/mapping.py])
AC_CONFIG_FILES([tests/py/casmutils/mapping/structure.py],[chmod +x tests/py/casmutils/mapping/structure.py])
AC_CONFIG_FILES([tests/py/casmutils/sym/cart.py],[chmod +x tests/py/casmutils/sym/cart.py])
AC_CONFIG_FILES([tests/py/casmutils/xtal/symmetry.py],[chmod +x tests/py/casmutils/xtal/symmetry.py])

AC_CONFIG_FILES([tests/integration/phony_module_installation.sh])
AC_CONFIG_FILES([tests/integration/cu-import.py],[chmod +x tests/integration/cu-import.py])
####################################################################################


AC_OUTPUT
