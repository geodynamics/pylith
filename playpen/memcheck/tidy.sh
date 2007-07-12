#!/bin/bash
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

# Remove autoconf stuff
autoconf_files="aclocal.m4 configure portinfo.in"
autoconf_dirs="autom4te.cache aux-config"
rm -r $autoconf_dirs
rm $autoconf_files
rm `find . -name Makefile.in`

# Remove emacs backup stuff
rm `find . -name "*~"`

# Remove compiled python
rm `find . -name "*.pyc"`

# version
# $Id: tidy.sh,v 1.1 2005/10/10 17:21:19 brad Exp $

# End of file
