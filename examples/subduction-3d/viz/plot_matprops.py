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
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

# Plot the domain, colored by material properties.


# User-specified parameters.
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
    "FIELD": "mu",
    "MATERIALS": ["crust", "mantle", "wedge", "slab"],
    }

# ----------------------------------------------------------------------
from paraview.simple import *
import os

def visualize(parameters):

    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()


    dataAll = []
    # Read data
    for material in parameters.materials:
        filename = os.path.join(parameters.output_dir, "%s-%s_info.xmf" % (parameters.sim, material))
        if not os.path.isfile(filename):
            raise IOError("File '%s' does not exist." % filename)
        dataMaterial = XDMFReader(FileNames=[filename])
        RenameSource("%s-%s" % (parameters.sim, material), dataMaterial)
        dataAll.append(dataMaterial)
    groupMaterials = GroupDatasets(Input=dataAll)

    view = GetActiveViewOrCreate('RenderView')

    # Show domain, colored by magnitude of displacement vector.
    materialDisplay = Show(groupMaterials, view)
    ColorBy(materialDisplay, ("CELLS", parameters.field))
    materialDisplay.RescaleTransferFunctionToDataRange(True)
    materialDisplay.SetScalarBarVisibility(view, True)
    materialDisplay.SetRepresentationType("Surface With Edges")

    # Rescale color and/or opacity maps used to exactly fit the current data range
    materialLUT = GetColorTransferFunction(parameters.field)
    materialDisplay.RescaleTransferFunctionToDataRange(False, False)
    # Update scalar bar component title.
    UpdateScalarBarsComponentTitle(materialLUT, materialDisplay)

    view.ResetCamera()
    Render()

class Parameters(object):
    keys = ("OUTPUT_DIR", "SIM", "FIELD", "MATERIALS")
    
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
    parser.add_argument("--field", action="store", dest="field", default=DEFAULTS["FIELD"])
    parser.add_argument("--materials", action="store", dest="materials")
    args = parser.parse_args()

    if args.materials:
        args.materials = args.materials.split(",")
    else:
        args.materials = DEFAULTS["MATERIALS"]
    
    visualize(args)
    Interact()

else:
    # Running inside the ParaView GUI

    visualize(Parameters())


# End of file
