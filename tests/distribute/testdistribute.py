#!/usr/bin/env nemesis
# -*- Python -*-
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

__requires__ = "PyLith"

if __name__ == "__main__":

    from MesherApp import MesherApp
    from pyre.applications import start
    start(applicationClass=MesherApp)

# End of file 
