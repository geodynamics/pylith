## Process this file with automake to produce Makefile.in
# $Id: Makefile.am 2370 2005-09-09 19:49:10Z leif $

ACLOCAL_AMFLAGS = -I ./m4

SUBDIRS = pylith3d

install-data-hook:
	top_builddir=`$(am__cd) $(top_builddir) && pwd` && \
		cd $(top_srcdir) && \
		$(PYTHON) setup.py install_deps -H None -f $$top_builddir/deps --prefix=$(prefix) && \
		$(PYTHON) setup.py egg_info -e $(pythondir)

## end of Makefile.am
