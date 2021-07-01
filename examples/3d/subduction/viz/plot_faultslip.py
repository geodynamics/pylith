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
    "SIM": "step02",
    "FAULTS": ["fault-slab"],
    "TIMESTEP": 0, # Use 0 for first, -1 for last.
    }

# ----------------------------------------------------------------------
from paraview.simple import *
import os

def visualize(parameters):
    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()

    # Read domain data
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

    faultDisplay = Show(groupFaults, view)
    ColorBy(faultDisplay, ('POINTS', 'slip', 'Magnitude'))
    faultDisplay.RescaleTransferFunctionToDataRange(True)
    faultDisplay.SetScalarBarVisibility(view, True)
    faultDisplay.SetRepresentationType('Surface With Edges')

    # Rescale color and/or opacity maps used to exactly fit the current data range
    slipLUT = GetColorTransferFunction('slip')
    faultDisplay.RescaleTransferFunctionToDataRange(False, False)
    # Update a scalar bar component title.
    UpdateScalarBarsComponentTitle(slipLUT, faultDisplay)

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
    keys = ("OUTPUT_DIR", "SIM", "FAULTS", "TIMESTEP")
    
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
    parser.add_argument("--faults", action="store", dest="faults")
    parser.add_argument("--timestep", action="store", dest="timestep", default=-1)
    parser.add_argument("--screenshot", action="store", dest="screenshot")
    args = parser.parse_args()

    if args.faults:
        args.faults = args.faults.split(",")
    else:
        args.faults = DEFAULTS["FAULTS"]
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
