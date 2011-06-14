#!/usr/bin/env nemesis
# -*- Python -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

__requires__ = "PyLith"

if __name__ == "__main__":

    from MesherApp import MesherApp
    from pyre.applications import start
    start(applicationClass=MesherApp)

# End of file 
