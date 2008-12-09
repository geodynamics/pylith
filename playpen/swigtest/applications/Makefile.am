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

bin_PROGRAMS = testswig

bin_SCRIPTS = testswig.py

INCLUDES = -I$(top_srcdir)/libsrc

testswig_SOURCES = \
	testswig.cc
testswig_LDADD = \
	$(top_builddir)/libsrc/libswigtest.la

CLEANFILES = $(bin_SCRIPTS)


# End of file 
