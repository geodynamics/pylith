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
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md.md for license information.
#
# ----------------------------------------------------------------------

# Plot the undeformed domain as a gray wireframe and then the fault
# surfaces, colored by the magnitude of fault slip.
#
# This Python script runs using pvpython or within the ParaView Python
# shell.

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
    "SIM": "step05",
    "FAULTS": ["fault-slabtop"],
    "TIMESTEP": 0, # Use 0 for first, -1 for last.
    }

# ----------------------------------------------------------------------
from paraview.simple import *
import os
import numpy

def visualize(parameters):
    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()

    
    # Read fault data
    dataFaults = []
    for fault in parameters.faults:
        filename = os.path.join(parameters.output_dir, "%s-%s.xmf" % (parameters.sim, fault))
        if not os.path.isfile(filename):
            raise IOError("File '%s' does not exist." % filename)
        data = XDMFReader(FileNames=[filename])
        RenameSource("%s-%s" % (parameters.sim, fault), data)
        dataFaults.append(data)

    groupFaults = GroupDatasets(Input=dataFaults)

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    view = GetActiveViewOrCreate('RenderView')

    # Ratio of shear to normal traction
    calculatorRatio = Calculator(Input=groupFaults)
    calculatorRatio.Function = '-abs(traction_X)/traction_Y'
    calculatorRatio.ResultArrayName = 'shearRatio'

    ratioDisplay = Show(calculatorRatio, view)
    ColorBy(ratioDisplay, ('POINTS', 'shearRatio'))
    ratioDisplay.RescaleTransferFunctionToDataRange(True)
    ratioDisplay.SetScalarBarVisibility(view, True)
    ratioDisplay.SetRepresentationType('Wireframe')
    ratioDisplay.LineWidth = 8.0
    
    # Rescale color and/or opacity maps used to exactly fit the current data range
    ratioLUT = GetColorTransferFunction('shearDivNormal')
    ratioDisplay.RescaleTransferFunctionToDataRange(False, False)
    # Update a scalar bar component title.
    UpdateScalarBarsComponentTitle(ratioLUT, ratioDisplay)

    # Annotate time
    tstamp = AnnotateTimeFilter(groupFaults)
    tstamp.Format = 'Time: %2.0f yr'
    tstamp.Scale = 3.168808781402895e-08 # seconds to years

    tstampDisplay = Show(tstamp, view)
    tstampDisplay.FontFamily = "Courier"
    tstampDisplay.FontSize = 14
    
    view.ResetCamera()
    view.Update()
    Render()

class Parameters(object):
    keys = ("OUTPUT_DIR", "SIM", "FAULTS")
    
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
    parser.add_argument("--sim", action="store", dest="sim", default=DEFAULTS["SIM"])
    parser.add_argument("--faults", action="store", dest="faults")
    args = parser.parse_args()

    if args.faults:
        args.faults = args.faults.split(",")
    else:
        args.faults = DEFAULT["FAULTS"]
    visualize(args.sim)

    view = GetRenderView()
    view.ViewSize = [960, 540]
    view.Update()

    Interact()

else:
    # Running inside the ParaView GUI

    visualize(Parameters())


# End of file
