#Store macros locally
CLOCAL_AMFLAGS = -I build-aux/m4

#===========================================================================#
#Values are filled by Makemodue.am files
lib_LTLIBRARIES=
pyexec_LTLIBRARIES=
noinst_LIBRARIES=
bin_PROGRAMS=
check_PROGRAMS=
man1_MANS=
dist_bin_SCRIPTS=
TESTS=

AM_CXXFLAGS=

AM_CPPFLAGS= -I$(srcdir)/include\
			 $(BOOST_CPPFLAGS)\
			 $(PYTHON_CPPFLAGS)

AM_LDFLAGS=$(BOOST_LDFLAGS)


#===========================================================================#
#Now do the programs
include $(srcdir)/src/Makemodule.am
include $(srcdir)/utilipy/Makemodule.am
include $(srcdir)/utilities/Makemodule.am


#===========================================================================#
EXTRA_DIST=

var-check:
	echo $(LIBS)
	echo $(PYTHON_LIBS)
	echo $(BOOST_SYSTEM_LIB)
	echo $(BOOST_PYTHON_LIB)
	echo $(BOOST_PYTHON_LIBS)
	echo $(BOOST_LDFLAGS)
	echo $(BOOST_CPPFLAGS)
	echo $(PYTHON)
	echo $(PYTHON_CPPFLAGS)
	echo "lol"
