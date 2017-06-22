#!/usr/bin/env pvpython
# -*- Python -*- (syntax highlighting)
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------

# Plot the domain, colored by materials.


# User-specified parameters.
#
# These are used if running from within the ParaView GUI via the
# Python shell or as defaults if running outside the ParaView GUI via
# pvpython.

EXODUS_FILE = "mesh/mesh_tet.exo"

# ----------------------------------------------------------------------
from paraview.simple import *
import os

def visualize(filename):

    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()

    if not os.path.isfile(filename):
        raise IOError("Exodus file '%s' does not exist." % filename)
    dataDomain = ExodusIIReader(FileName=[filename])
    RenameSource("domain", dataDomain)

    view = GetActiveViewOrCreate('RenderView')

    # Show domain, colored by block.
    domainDisplay = Show(dataDomain, view)
    ColorBy(domainDisplay, ("FIELD", "vtkBlockColors"))
    domainDisplay.RescaleTransferFunctionToDataRange(True, False)
    domainDisplay.SetScalarBarVisibility(view, False)
    domainDisplay.SetRepresentationType("Surface")
    domainDisplay.Opacity = 0.5

    # Add coordinate axes
    axes = Axes()
    axes.ScaleFactor = 1.0e+5

    axesDisplay = Show(axes, view)
    axesDisplay.SetRepresentationType('Wireframe')
    axesDisplay.LineWidth = 4.0
    axesDisplay.SetScalarBarVisibility(view, False)
    axesDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
    
    view.ResetCamera()
    Render()

# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Running from outside the ParaView GUI via pvpython
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", action="store", dest="filename", default=EXODUS_FILE)
    args = parser.parse_args()

    visualize(args.filename)

    view = GetRenderView()
    #view.CameraPosition = [-80160, -1130898, 133977]
    #view.CameraFocalPoint = [-55107, 446810, 283137]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.CameraViewAngle = 45.0
    view.ViewSize = [960, 540]
    view.Update()

    Interact()

else:
    # Running inside the ParaView GUI

    visualize(EXODUS_FILE)


# End of file
