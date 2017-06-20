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

# Plot the domain, colored by the magnitude of the displacement
# vector, with white arrows showing the displacement vectors.
#
# This Python script runs using pvpython or within the ParaView Python
# shell.


# User-specified parameters
#
# These are used if running from within the ParaView GUI via the
# Python shell or as defaults if running outside the ParaView GUI via
# pvpython.

# Root name for simulation.
SIM_NAME = "step01"

# Scale used for displacement vectors.
DISPLACEMENT_SCALE = 10.0e+3


# ----------------------------------------------------------------------
from paraview.simple import *
import os

def visualize(sim, dispScale):
    
    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()

    # Read data
    filename = "output/%s-domain.xmf" % sim
    if not os.path.isfile(filename):
        raise IOError("File '%s' does not exist." % filename)
    dataDomain = XDMFReader(FileNames=[filename])
    RenameSource("%s-domain" % sim, dataDomain)

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    view = GetActiveViewOrCreate('RenderView')

    # Show undeformed domain, colored by magnitude of displacement vector.
    domainDisplay = Show(dataDomain, view)
    ColorBy(domainDisplay, ("POINTS", "displacement", "Magnitude"))
    domainDisplay.RescaleTransferFunctionToDataRange(True)
    domainDisplay.SetScalarBarVisibility(view, True)
    domainDisplay.SetRepresentationType("Surface With Edges")
    
    # Rescale color and/or opacity maps used to exactly fit the current data range
    displacementLUT = GetColorTransferFunction('displacement')
    domainDisplay.RescaleTransferFunctionToDataRange(False, False)
    # Update scalar bar component title.
    UpdateScalarBarsComponentTitle(displacementLUT, domainDisplay)


    # Add arrows to show displacement vectors.
    glyph = Glyph(Input=dataDomain, GlyphType="Arrow")
    glyph.Vectors = ["POINTS", "displacement"]
    glyph.ScaleFactor = dispScale
    glyph.ScaleMode = "vector"
    glyph.GlyphMode = "All Points"

    glyphDisplay = Show(glyph, view)
    glyphDisplay.Representation = "Surface"
    
    view.ResetCamera()
    
    Render()

# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Running from outside the ParaView GUI via pvpython
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", action="store", dest="sim")
    parser.add_argument("--vector-scale", action="store", type=float, dest="scale")
    args = parser.parse_args()

    sim = args.sim
    scale = args.scale
    if sim is None:
        sim = SIM_NAME
    if scale is None:
        scale = DISPLACEMENT_SCALE
    visualize(sim, scale)
    Interact()

else:
    # Running inside the ParaView GUI

    visualize(SIM_NAME, DISPLACEMENT_SCALE)


# End of file
