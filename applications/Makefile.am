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

bin_SCRIPTS = pylithic.py
do_build = sed -e s%[@]PYTHONPATH[@]%$(PYTHONPATH)%g
pylithic.py: $(srcdir)/pylithic.py.in Makefile
	$(do_build) <  $(srcdir)/pylithic.py.in > $@ || (rm -f $@ && exit 1)
	chmod +x $@
EXTRA_DIST = pylithic.py.in
CLEANFILES = $(bin_SCRIPTS)

# version
# $Id$

# End of file 
