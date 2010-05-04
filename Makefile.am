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

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = \
	libsrc \
	modulesrc \
	pylith \
	applications \
	templates \
	share

if ENABLE_TESTING
  SUBDIRS += \
	unittests

# Full-scale tests in progress.
#	tests_auto
endif

DIST_SUBDIRS = $(SUBDIRS) examples tests doc

if ENABLE_DOCUMENTATION
  SUBDIRS += doc
endif

EXTRA_DIST = \
	CHANGES \
	DEPENDENCIES \
	setup.py \
	archimedes/__init__.py

install-data-hook:
	top_builddir=`$(am__cd) $(top_builddir) && pwd` && \
		cd $(top_srcdir) && \
		$(PYTHON) setup.py install_deps -H None -f $$top_builddir/deps --prefix=$(prefix) && \
		$(PYTHON) setup.py egg_info -e $(pythondir)


# End of file 
