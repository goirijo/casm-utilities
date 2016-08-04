#Store macros locally
CLOCAL_AMFLAGS = -I build-aux/m4

#===========================================================================#
#Values are filled by Makemodue.am files
bin_PROGRAMS=
noinst_LIBRARIES=

AM_CXXFLAGS=

AM_CPPFLAGS=

AM_LDFLAGS=$(BOOST_LDFLAGS)


#===========================================================================#
#First do all the library stuff
include $(srcdir)/lib/splash/Makemodule.am


#===========================================================================#
#Now do the programs
include $(srcdir)/src/utilities/Makemodule.am


#===========================================================================#
EXTRA_DIST=
