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
SIM_NAME = "step05"

# Names of faults for output files.
FAULTS = ["fault-slab"]

# ----------------------------------------------------------------------
from paraview.simple import *
import os
import numpy

def visualize(sim, faults):
    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()

    # Read domain data
    filename = "output/%s-domain.xmf" % sim
    if not os.path.isfile(filename):
        raise IOError("File '%s' does not exist." % filename)
    dataDomain = XDMFReader(FileNames=[filename])
    RenameSource("%s-domain" % sim, dataDomain)

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    view = GetActiveViewOrCreate('RenderView')

    # Gray wireframe for undeformed domain.
    domainDisplay = Show(dataDomain, view)
    domainDisplay.Representation = 'Wireframe'
    domainDisplay.AmbientColor = [0.5, 0.5, 0.5]

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

    # Ratio of shear to normal traction
    calculator = Calculator(Input=groupFaults)
    calculator.Function = 'mag(traction_X*iHat+traction_Y*jHat)/abs(traction_Z)'
    calculator.ResultArrayName = 'shearDivNormal'

    ratioDisplay = Show(groupFaults, view)
    ColorBy(ratioDisplay, ('POINTS', 'shearDivNormal'))
    ratioDisplay.RescaleTransferFunctionToDataRange(True)
    ratioDisplay.SetScalarBarVisibility(view, True)
    ratioDisplay.SetRepresentationType('Surface With Edges')

    # Rescale color and/or opacity maps used to exactly fit the current data range
    ratioLUT = GetColorTransferFunction('shearDivNormal')
    ratioDisplay.RescaleTransferFunctionToDataRange(False, False)
    # Update a scalar bar component title.
    UpdateScalarBarsComponentTitle(ratioLUT, ratioDisplay)

    calculator = Calculator(Input=groupFaults)
    calculator.Function = 'mag(slip)'
    calculator.ResultArrayName = 'slipMag'
    
    contour = Contour(Input=groupFaults)
    contour.ContourBy = ['POINTS', 'slipMag']
    contour.Isosurfaces = numpy.arange(0.0, 4.01, 0.5)
    contour.PointMergeMethod = 'Uniform Binning'

    contourDisplay = Show(contour, view)


    
    view.ResetCamera()
    view.Update()
    Render()

# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Running from outside the ParaView GUI via pvpython
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", action="store", dest="sim")
    parser.add_argument("--faults", action="store", dest="faults")
    args = parser.parse_args()

    sim = args.sim
    if args.faults:
        faults = args.faults.split(",")
    else:
        faults = None

    if sim is None:
        sim = SIM_NAME
    if faults is None:
        faults = FAULTS
        
    visualize(sim, faults)
    Interact()

else:
    # Running inside the ParaView GUI

    visualize(SIM_NAME, FAULTS)


# End of file
