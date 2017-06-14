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

# User-specified parameters

# Root name for simulation.
SIM_NAME = "step02"

# Names of faults for output files.
FAULTS = ["fault-slab"]

# ----------------------------------------------------------------------
from paraview.simple import *
# Disable automatic camera reset on "Show"
paraview.simple._DisableFirstRenderCameraReset()

# Read domain data
dataDomain = XDMFReader(FileNames=["output/%s-domain.xmf" % SIM_NAME])
RenameSource("%s-domain" % SIM_NAME, dataDomain)

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
view = GetActiveViewOrCreate('RenderView')

# Gray wireframe for undeformed domain.
domainDisplay = Show(dataDomain, view)
domainDisplay.Representation = 'Wireframe'
domainDisplay.AmbientColor = [0.5, 0.5, 0.5]

# Read fault data
dataFaults = []
for fault in FAULTS:
    data = XDMFReader(FileNames=["output/%s-%s.xmf" % (SIM_NAME, fault)])
    RenameSource("%s-%s" % (SIM_NAME, fault), data)
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

view.ResetCamera()
view.Update()

Render()

# Uncomment if running from shell outside ParaView.
#Interact()


# End of file
