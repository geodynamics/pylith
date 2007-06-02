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

ACLOCAL_AMFLAGS = -I ./m4

SUBDIRS = \
	libsrc \
	modulesrc \
	pylith \
	applications

if ENABLE_TESTING
  SUBDIRS += \
	unittests \
	tests
endif

if ENABLE_DOCUMENTATION
  SUBDIRS += doc
endif

install-data-hook:
	top_builddir=`$(am__cd) $(top_builddir) && pwd` && \
		cd $(top_srcdir) && \
		$(PYTHON) setup.py install_deps -H None -f $$top_builddir/deps --prefix=$(prefix) && \
		$(PYTHON) setup.py egg_info -e $(pythondir)

# version
# $Id$

# End of file 
