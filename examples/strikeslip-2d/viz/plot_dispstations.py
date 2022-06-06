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
# Copyright (c) 2010-2022 University of California, Davis
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
OUTPUT_DIR = os.path.join(os.environ["HOME"], "src", "pylith", "examples", "strikeslip-2d", "output")
"""

DEFAULTS = {
    "OUTPUT_DIR": "output",
    "SIM": "step01_slip",
    "GLYPH_SCALE": 50.0e+3,
    "FIELD": "displacement",
    "TIMESTEP": 0,  # Use 0 for first, -1 for last.
}


# ----------------------------------------------------------------------
def visualize(parameters):

    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()

    # Read data
    filename = os.path.join(parameters.output_dir, "%s-gps_stations.xmf" % parameters.sim)
    if not os.path.isfile(filename):
        raise IOError("File '%s' does not exist." % filename)
    dataStations = XDMFReader(FileNames=[filename])
    RenameSource("%s-gps_stations" % parameters.sim, dataStations)

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    if parameters.timestep == -1:
        scene.GoToLast()

    view = GetActiveViewOrCreate('RenderView')

    # Add displacement vectors
    glyph = Glyph(registrationName='Glyph', Input=dataStations, GlyphType='Arrow')
    glyph.OrientationArray = ['POINTS', parameters.field]
    glyph.ScaleArray = ['POINTS', parameters.field]
    glyph.ScaleFactor = parameters.glyph_scale
    glyph.GlyphTransform = 'Transform2'

    glyphDisplay = Show(glyph, view, 'GeometryRepresentation')
    glyphDisplay.AmbientColor = [1.0, 1.0, 0.0]
    glyphDisplay.DiffuseColor = [1.0, 1.0, 0.0]
    glyphDisplay.Representation = 'Surface'

    # Annotate time
    tstamp = AnnotateTimeFilter(glyph)
    tstamp.Format = 'Time: {time:5.1f} yr'
    tstamp.Scale = 3.168808781402895e-08  # seconds to years

    tstampDisplay = Show(tstamp, view)
    tstampDisplay.FontFamily = "Courier"
    tstampDisplay.FontSize = 14

    view.ResetCamera()
    view.Update()
    Render()
    return


class Parameters(object):
    keys = ("OUTPUT_DIR", "SIM", "GLYPH_SCALE", "FIELD", "TIMESTEP")

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
    parser.add_argument("--glyph-scale", action="store", type=float, dest="glyph_scale", default=DEFAULTS["GLYPH_SCALE"])
    parser.add_argument("--field", action="store", dest="field", default=DEFAULTS["FIELD"])
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
