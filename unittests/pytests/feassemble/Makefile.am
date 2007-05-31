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

subpackage = feassemble
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testfeassemble.py

check_SCRIPTS = testfeassemble.py

noinst_PYTHON = \
	TestFIATSimplex.py \
	TestQuadrature.py \
	TestIntegrator.py


# module
subpkgpyexec_LTLIBRARIES = testfeassemblemodule.la

testfeassemblemodule_la_LDFLAGS = -module

testfeassemblemodule_la_SOURCES = testfeassemble.pyxe

nodist_testfeassemblemodule_la_SOURCES = \
	testfeassemble.c testfeassemble_embed.cpp testfeassemble_embed.h

testfeassemblemodule_la_LIBADD = \
	$(PETSC_LIB)

INCLUDES += -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

testfeassemble.pyx testfeassemble_embed.cpp  testfeassemble_embed.h: testfeassemble.pyxe
	pyrexembed testfeassemble.pyxe
testfeassemble.pyxe: $(srcdir)/testfeassemble.pyxe.src
	cp $(srcdir)/testfeassemble.pyxe.src $@
testfeassemble_embed.cpp: testfeassemble_embed.h
testfeassemble_embed.h: testfeassemble.pyx

.pyx.c:
	pyrexc $<

CLEANFILES = \
	testfeassemble.pyxe testfeassemble.pyx testfeassemble.c *_embed.* \
	$(am__installdirs)/testfeassemblemodule.*

TESTS_ENVIRONMENT = $(PYTHON)


# End of file 
