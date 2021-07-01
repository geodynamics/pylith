#!/usr/bin/env pvpython
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
import os
from paraview.simple import *

"""
Plot the undeformed domain as a gray wireframe and then the deformed
domain, colored by the value of the y-displacemenet.


User-specified parameters.

Default values for parameters. To use different values, overwrite
them in the ParaView Python shell or on the command line. For
example, set OUTPUT_DIR to the absolute path if not starting
ParaView from the terminal shell where you ran PyLith:

import os
OUTPUT_DIR = os.path.join(os.environ["HOME"], "src", "pylith", "examples", "2d", "box", "output")
"""

DEFAULTS = {
    "OUTPUT_DIR": "output",
    "SIM": "step01_slip",
    "WARP_SCALE": 2.0e+3,
    "FIELD": "displacement",
    "FIELD_COMPONENT": "Y",
    "TIMESTEP": 0,  # Use 0 for first, -1 for last.
}


# ----------------------------------------------------------------------
def visualize(parameters):

    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()

    # Read data
    filename = os.path.join(parameters.output_dir, "%s-domain.xmf" % parameters.sim)
    if not os.path.isfile(filename):
        raise IOError("File '%s' does not exist." % filename)
    dataDomain = XDMFReader(FileNames=[filename])
    RenameSource("%s-domain" % parameters.sim, dataDomain)

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    if parameters.timestep == -1:
        scene.GoToLast()

    view = GetActiveViewOrCreate('RenderView')

    # Gray wireframe for undeformed domain.
    domainDisplay = Show(dataDomain, view)
    domainDisplay.Representation = 'Wireframe'
    domainDisplay.AmbientColor = [0.5, 0.5, 0.5]

    # Warp domain to show deformation
    warp = WarpByVector(Input=dataDomain)
    warp.Vectors = ['POINTS', 'displacement']
    warp.ScaleFactor = parameters.warp_scale

    warpDisplay = Show(warp, view)
    ColorBy(warpDisplay, ('POINTS', parameters.field, parameters.field_component))
    warpDisplay.RescaleTransferFunctionToDataRange(True)
    warpDisplay.SetScalarBarVisibility(view, True)
    warpDisplay.SetRepresentationType('Surface With Edges')
    # Rescale color bar to exactly fit the current data range
    warpDisplay.RescaleTransferFunctionToDataRange(False, False)

    # Customize colorbar
    displacementLUT = GetColorTransferFunction(parameters.field)
    colorbar = GetScalarBar(displacementLUT, view)
    if parameters.field_component.lower() == "magnitude":
        colorbar.Title = "Displacement Mag. (m)"
    else:
        colorbar.Title = "%s-displacement (m)" % parameters.field_component.lower()
    colorbar.ComponentTitle = ""

    # Annotate time
    tstamp = AnnotateTimeFilter(warp)
    tstamp.Format = 'Time: %5.1f yr'
    tstamp.Scale = 3.168808781402895e-08  # seconds to years

    tstampDisplay = Show(tstamp, view)
    tstampDisplay.FontFamily = "Courier"
    tstampDisplay.FontSize = 12

    view.ResetCamera()
    view.Update()
    Render()
    return


class Parameters(object):
    keys = ("OUTPUT_DIR", "SIM", "WARP_SCALE", "FIELD", "FIELD_COMPONENT", "TIMESTEP")

    def __init__(self):
        globalVars = globals()
        for key in Parameters.keys:
            if key in globalVars.keys():
                setattr(self, key.lower(), globalVars[key])
            else:
                setattr(self, key.lower(), DEFAULTS[key])
        return


# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Running from outside the ParaView GUI via pvpython

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", action="store", dest="output_dir", default=DEFAULTS["OUTPUT_DIR"])
    parser.add_argument("--sim", action="store", dest="sim", default=DEFAULTS["SIM"])
    parser.add_argument("--warp_scale", action="store", type=float, dest="warp_scale", default=DEFAULTS["WARP_SCALE"])
    parser.add_argument("--field", action="store", dest="field", default=DEFAULTS["FIELD"])
    parser.add_argument("--component", action="store", dest="field_component", default=DEFAULTS["FIELD_COMPONENT"])
    parser.add_argument("--timestep", action="store", type=int, dest="timestep", default=-1)
    parser.add_argument("--screenshot", action="store", dest="screenshot")
    args = parser.parse_args()

    visualize(args)

    view = GetRenderView()
    view.ViewSize = [1024, 540]
    view.CameraPosition = [68527.89880980579, -152111.39463431376, 1405120.1034155919]
    view.CameraFocalPoint = [68527.89880980579, -152111.39463431376, -1004573.5784798338]
    view.Update()

    if args.screenshot:
        WriteImage(args.screenshot)

    Interact()

else:
    visualize(Parameters())


# End of file
