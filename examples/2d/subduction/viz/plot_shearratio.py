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

# Plot the undeformed domain as a gray wireframe and then the fault
# surfaces, colored by the magnitude of fault slip.
#
# This Python script runs using pvpython or within the ParaView Python
# shell.

# User-specified parameters.
#
# These are used if running from within the ParaView GUI via the
# Python shell or as defaults if running outside the ParaView GUI via
# pvpython.

# Root name for simulation.
SIM_NAME = "stepXX"

# Names of faults (with spontaneous rupture) for output files.
FAULTS = ["fault-slabtop"]

# ----------------------------------------------------------------------
from paraview.simple import *
import os
import numpy

def visualize(sim, faults):
    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()

    
    # Read fault data
    dataFaults = []
    for fault in faults:
        filename = "output/%s-%s.xmf" % (sim, fault)
        if not os.path.isfile(filename):
            raise IOError("File '%s' does not exist." % filename)
        data = XDMFReader(FileNames=[filename])
        RenameSource("%s-%s" % (sim, fault), data)
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

# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Running from outside the ParaView GUI via pvpython
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", action="store", dest="sim", default=SIM_NAME)
    parser.add_argument("--faults", action="store", dest="faults")
    args = parser.parse_args()

    if args.faults:
        faults = args.faults.split(",")
    else:
        faults = FAULTS
    visualize(args.sim, faults)

    view = GetRenderView()
    view.ViewSize = [960, 540]
    view.Update()

    Interact()

else:
    # Running inside the ParaView GUI

    visualize(SIM_NAME, FAULTS)


# End of file
