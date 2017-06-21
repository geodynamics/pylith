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

def visualize(sim, dispScale, showFinalTimeStep=False):
    
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
    if showFinalTimeStep:
        scene.GoToLast()
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
    glyph.GlyphType.TipRadius = 0.2
    glyph.GlyphType.ShaftRadius = 0.05

    glyphDisplay = Show(glyph, view)
    glyphDisplay.Representation = "Surface"

    # Annotate time
    tstamp = AnnotateTimeFilter(dataDomain)
    tstamp.Format = 'Time: %2.0f yr'
    tstamp.Scale = 3.168808781402895e-08 # seconds to years

    tstampDisplay = Show(tstamp, view)
    tstampDisplay.FontFamily = "Courier"
    tstampDisplay.FontSize = 14
    
    view.ResetCamera()
    view.Update()
    Render()

# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Running from outside the ParaView GUI via pvpython
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", action="store", dest="sim", default=SIM_NAME)
    parser.add_argument("--vector-scale", action="store", type=float, dest="scale", default=DISPLACEMENT_SCALE)
    parser.add_argument("--screenshot", action="store", dest="screenshot")
    args = parser.parse_args()

    visualize(args.sim, args.scale, showFinalTimeStep=True)

    view = GetRenderView()
    view.CameraPosition = [78002.89373974672, -1531813.1739094853, 595774.2094961794]
    view.CameraFocalPoint = [-45014.6313325238, 149523.68421156122, -335271.271063906]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.ViewSize = [960, 540]
    view.Update()
    
    if args.screenshot:
        WriteImage(args.screenshot)

    Interact()

else:
    # Running inside the ParaView GUI

    visualize(SIM_NAME, DISPLACEMENT_SCALE)

    

# End of file
