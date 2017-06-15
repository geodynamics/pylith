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

# Plot the undeformed domain as a gray wireframe and then the deformed
# domain, colored by the value of the x-displacemenet.


# User-specified parameters.
#
# These are used if running from within the ParaView GUI via the
# Python shell or as defaults if running outside the ParaView GUI via
# pvpython.

# Root name for simulation.
SIM_NAME = "step01"

# Scale used to exaggerate deformation.
DISPLACEMENT_SCALE = 10.0e+3

# ----------------------------------------------------------------------
from paraview.simple import *

def visualize(sim, exaggeration):
    
    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()

    # Read data
    dataDomain = XDMFReader(FileNames=["output/%s-domain.xmf" % sim])
    RenameSource("%s-domain" % sim, dataDomain)

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    view = GetActiveViewOrCreate('RenderView')

    # Gray wireframe for undeformed domain.
    domainDisplay = Show(dataDomain, view)
    domainDisplay.Representation = 'Wireframe'
    domainDisplay.AmbientColor = [0.5, 0.5, 0.5]

    # Warp domain to show deformation
    warp = WarpByVector(Input=dataDomain)
    warp.Vectors = ['POINTS', 'displacement']
    warp.ScaleFactor = exaggeration

    warpDisplay = Show(warp, view)
    ColorBy(warpDisplay, ('POINTS', 'displacement', 'X'))
    warpDisplay.RescaleTransferFunctionToDataRange(True)
    warpDisplay.SetScalarBarVisibility(view, True)
    warpDisplay.SetRepresentationType('Surface With Edges')

    # Rescale color and/or opacity maps used to exactly fit the current data range
    displacementLUT = GetColorTransferFunction('displacement')
    warpDisplay.RescaleTransferFunctionToDataRange(False, False)
    # Update a scalar bar component title.
    UpdateScalarBarsComponentTitle(displacementLUT, warpDisplay)

    view.ResetCamera()
    view.Update()
    Render()

# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Running from outside the ParaView GUI via pvpython
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", action="store", dest="sim")
    parser.add_argument("--exaggeration", action="store", type=float, dest="exaggeration")
    args = parser.parse_args()

    sim = args.sim
    exaggeration = args.exaggeration
    if sim is None:
        sim = SIM_NAME
    if exaggeration is None:
        exaggeration = DISPLACEMENT_SCALE
    visualize(sim, exaggeration)
    Interact()

else:
    # Running inside the ParaView GUI

    visualize(SIM_NAME, DISPLACEMENT_SCALE)


# End of file
