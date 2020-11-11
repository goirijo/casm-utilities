#Store macros locally
CLOCAL_AMFLAGS = -I build-aux/m4

#===========================================================================#
#Values are filled by Makemodue.am files
lib_LTLIBRARIES=
noinst_LIBRARIES=
noinst_LTLIBRARIES=
bin_PROGRAMS=
noinst_PROGRAMS=
check_PROGRAMS=
man1_MANS=
dist_bin_SCRIPTS=
nobase_include_HEADERS=
TESTS=

BUILT_SOURCES=

#Completely overkill, but we want the copyright files too
EXTRA_DIST=submodules

AM_CXXFLAGS=\
			  -DEIGEN_DEFAULT_DENSE_INDEX_TYPE=long\
			  -DGZSTREAM_NAMESPACE=gz

AM_CPPFLAGS=\
			 -I$(srcdir)/include\
			 -I$(srcdir)/submodules/pybind11/include\
			 -I$(srcdir)/submodules/nlohmann-json/include\
			 -I$(srcdir)/submodules/CLI11/include\
			 -I$(srcdir)/submodules/CASMcode/include\
			 -I$(srcdir)/submodules/googletest/googletest\
			 -I$(srcdir)/submodules/googletest/googletest/include\
			 -I$(srcdir)/submodules/CASMcode/include/casm/external\
			 -I$(srcdir)/submodules/CASMcode/include/casm/external/qhull/libqhullcpp\
			 -I$(srcdir)/submodules/CASMcode/include/casm/external/qhull/libqhull_r\
			 -I$(srcdir)/submodules/CASMcode/include/casm/external/gzstream\
			 $(PYTHON_CPPFLAGS)\
			 -DABS_SRCDIR=\"$(abs_srcdir)\"\
			 -DABS_TOP_BUILDDIR=\"$(abs_top_builddir)\"

#===========================================================================#
#Now do the programs
#You may comment out any of these lines to suppress parts of the installation process
include $(srcdir)/lib/Makemodule.am			#utilities library

if CASMUTILS_LIBRARY
include $(srcdir)/include/Makemodule.am		#utilities headers
endif

if CASMUTILS_PYTHON
pyexec_LTLIBRARIES=
include $(srcdir)/lib-py/Makemodule.am 		#python library wrappers to utilities library
endif

if WITH_PLUGINS
# include $(srcdir)/src/Makemodule.am			#utilities script for calling plugins
include $(srcdir)/plugins/Makesocket.am		#CLIe utilities, or anything else to plug in via submodule
endif

include $(srcdir)/tests/Makemodule.am 		#All types of tests

if FULL_CASM
include $(srcdir)/lib/casmblob/Makefullcasm.am
else
include $(srcdir)/lib/casmblob/Makepartialcasm.am
endif

#===========================================================================#
var-check:
	echo ${prefix}
	echo ${srcdir}
	echo ${abs_srcdir}
	echo ${abs_top_builddir}
	echo $(LIBS)
	echo $(PYTHON_LIBS)
	echo $(BOOST_SYSTEM_LIB)
	echo $(BOOST_PYTHON_LIB)
	echo $(BOOST_PYTHON_LIBS)
	echo $(BOOST_LDFLAGS)
	echo $(BOOST_CPPFLAGS)
	echo $(PYTHON)
	echo $(PYTHON_CPPFLAGS)
	echo $(AM_V_CXX)
	echo $(CXXCOMPILE)
	echo $(DEPDIR)
	echo $(OUTSIDE)
	echo "lol"


soft-clean:
	make clean-libtool
	make clean-pyexecLTLIBRARIES
	make clean-frankensteinpyLTLIBRARIES
	make clean-structurepyLTLIBRARIES

