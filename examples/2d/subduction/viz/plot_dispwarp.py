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
SIM_NAME = "step05"

# Scale used to exaggerate deformation.
DISPLACEMENT_SCALE = 10.0e+3
FIELD = "displacement"
FIELD_COMPONENT = "X"

# ----------------------------------------------------------------------
from paraview.simple import *
import os

def visualize(sim, exaggeration, field, component, showFinalTimeStep=False):
    
    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()

    # Read data
    filename = "output/%s.xmf" % sim
    if not os.path.isfile(filename):
        raise IOError("File '%s' does not exist." % filename)
    dataDomain = XDMFReader(FileNames=[filename])
    RenameSource("%s-domain" % sim, dataDomain)

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    if showFinalTimeStep:
        scene.GoToLast()

    view = GetActiveViewOrCreate('RenderView')

    # Gray wireframe for undeformed domain.
    domainDisplay = Show(dataDomain, view)
    domainDisplay.Representation = 'Wireframe'
    domainDisplay.AmbientColor = [0.5, 0.5, 0.5]
    Hide(dataDomain, view)

    # Warp domain to show deformation
    warp = WarpByVector(Input=dataDomain)
    warp.Vectors = ['POINTS', 'displacement']
    warp.ScaleFactor = exaggeration

    warpDisplay = Show(warp, view)
    ColorBy(warpDisplay, ('POINTS', field, component))
    warpDisplay.RescaleTransferFunctionToDataRange(True)
    warpDisplay.SetScalarBarVisibility(view, True)
    warpDisplay.SetRepresentationType('Surface With Edges')
    # Rescale color bar to exactly fit the current data range
    warpDisplay.RescaleTransferFunctionToDataRange(False, False)

    # Customize colorbar
    displacementLUT = GetColorTransferFunction('displacement')
    colorbar = GetScalarBar(displacementLUT, view)
    if component.lower() == "magnitude":
        colorbar.Title = "Displacement Mag. (m)"
    else:
        colorbar.Title = "%s-displacement (m)" % component.lower()
    colorbar.ComponentTitle = ""

    # Annotate time
    tstamp = AnnotateTimeFilter(warp)
    tstamp.Format = 'Time: %5.1f yr'
    tstamp.Scale = 3.168808781402895e-08 # seconds to years

    tstampDisplay = Show(tstamp, view)
    tstampDisplay.FontFamily = "Courier"
    tstampDisplay.FontSize = 12
    
    view.ResetCamera()
    view.Update()
    Render()

# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Running from outside the ParaView GUI via pvpython
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", action="store", dest="sim", default=SIM_NAME)
    parser.add_argument("--exaggeration", action="store", type=float, dest="exaggeration", default=DISPLACEMENT_SCALE)
    parser.add_argument("--field", action="store", dest="field", default=FIELD)
    parser.add_argument("--component", action="store", dest="component", default=FIELD_COMPONENT)
    parser.add_argument("--screenshot", action="store", dest="screenshot")
    args = parser.parse_args()

    visualize(args.sim, args.exaggeration, args.field, args.component, showFinalTimeStep=True)

    view = GetRenderView()
    view.ViewSize = [1024, 540]
    view.CameraPosition = [68527.89880980579, -152111.39463431376, 1405120.1034155919]
    view.CameraFocalPoint = [68527.89880980579, -152111.39463431376, -1004573.5784798338]
    view.Update()

    if args.screenshot:
        WriteImage(args.screenshot)    

    Interact()

else:
    # Running inside the ParaView GUI

    visualize(SIM_NAME, DISPLACEMENT_SCALE, FIELD, FIELD_COMPONENT)


# End of file
