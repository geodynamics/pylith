#!/usr/bin/env pvpython
# -*- Python -*- (syntax highlighting)
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
# See LICENSE.md.md for license information.
#
# ----------------------------------------------------------------------

# Plot the domain, colored by the magnitude of the displacement
# vector, with white arrows showing the displacement vectors.
#
# This Python script runs using pvpython or within the ParaView Python
# shell.


# User-specified parameters
#
# Default values for parameters. To use different values, overwrite
# them in the ParaView Python shell or on the command line. For
# example, set OUTPUT_DIR to the absolute path if not starting
# ParaView from the terminal shell where you ran PyLith:
#
# import os
# OUTPUT_DIR = os.path.join(os.environ["HOME"], "src", "pylith", "examples", "2d", "subduction", "output")

DEFAULTS = {
    "OUTPUT_DIR": "output",
    "SIM": "step02",
    "VECTOR_SCALE": 10.0e+3,
    "FIELD": "displacement",
    "FIELD_COMPONENT": "Magnitude",
    "TIMESTEP": 0, # Use 0 for first, -1 for last.
    }

# ----------------------------------------------------------------------
from paraview.simple import *
import os

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

    # Show undeformed domain, colored by magnitude of displacement vector.
    domainDisplay = Show(dataDomain, view)
    ColorBy(domainDisplay, ("POINTS", parameters.field, parameters.field_component))
    domainDisplay.RescaleTransferFunctionToDataRange(True)
    domainDisplay.SetScalarBarVisibility(view, True)
    domainDisplay.SetRepresentationType("Surface With Edges")
    # Rescale color maps to exactly fit the current data range
    domainDisplay.RescaleTransferFunctionToDataRange(False, False)

    # Customize colorbar
    displacementLUT = GetColorTransferFunction(parameters.field)
    colorbar = GetScalarBar(displacementLUT, view)
    if parameters.field_component.lower() == "magnitude":
        colorbar.Title = "Displacement Mag. (m)"
    else:
        colorbar.Title = "%s-displacement (m)" % parameters.field_component.lower()
    colorbar.ComponentTitle = ""

    # Add arrows to show vectors.
    glyph = Glyph(Input=dataDomain, GlyphType="Arrow")
    glyph.Vectors = ["POINTS", parameters.field]
    glyph.ScaleFactor = parameters.vector_scale
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

class Parameters(object):
    keys = ("OUTPUT_DIR", "SIM", "VECTOR_SCALE", "FIELD", "FIELD_COMPONENT", "TIMESTEP")
    
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
    parser.add_argument("--vector-scale", action="store", type=float, dest="vector_scale", default=DEFAULTS["VECTOR_SCALE"])
    parser.add_argument("--field", action="store", dest="field", default=DEFAULTS["FIELD"])
    parser.add_argument("--component", action="store", dest="field_component", default=DEFAULTS["FIELD_COMPONENT"])
    parser.add_argument("--timestep", action="store", type=int, dest="timestep", default=-1)
    parser.add_argument("--screenshot", action="store", dest="screenshot")
    args = parser.parse_args()

    visualize(args)

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

    visualize(Parameters())

    

# End of file
