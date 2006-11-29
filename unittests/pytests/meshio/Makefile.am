# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

subpackage = meshio
include $(top_srcdir)/subpackage.am

TESTS = testdriver.py

check_SCRIPTS = testdriver.py

noinst_PYTHON = \
	TestMeshIO.py \
	TestMeshIOAscii.py


# module
subpkgpyexec_LTLIBRARIES = testmeshiomodule.la

testmeshiomodule_la_LDFLAGS = -module

testmeshiomodule_la_SOURCES = testmeshio.pyxe

nodist_testmeshiomodule_la_SOURCES = \
	testmeshio.c testmeshio_embed.cpp testmeshio_embed.h \
	$(top_srcdir)/unittests/libtests/meshio/TestMeshIOAscii.cc

testmeshiomodule_la_LIBADD = \
	$(PETSC_LIB)

INCLUDES += \
	-I$(PYTHON_INCDIR) \
	$(PETSC_INCLUDE) \
	-I$(top_srcdir)/unittests/libtests/meshio

testmeshio.pyx testmeshio_embed.cpp  testmeshio_embed.h: testmeshio.pyxe
	cp $< . && pyrexembed testmeshio.pyxe && rm -f testmeshio.pyxe
testmeshio_embed.cpp: testmeshio_embed.h
testmeshio_embed.h: testmeshio.pyx

.pyx.c:
	pyrexc $<

CLEANFILES = \
	testmeshio.pyx testmeshio.c *_embed.* \
	$(am__installdirs)/testmeshiomodule.*

# End of file 
