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
#
# This Python script runs using pvpython or within the ParaView Python
# shell.

# User-specified parameters

# Root name for simulation.
SIM_NAME = "step01"

# Scale used to exaggerate deformation.
DISPLACEMENT_SCALE = 10.0e+3

# ----------------------------------------------------------------------
from paraview.simple import *
# Disable automatic camera reset on "Show"
paraview.simple._DisableFirstRenderCameraReset()

# Read data
dataDomain = XDMFReader(FileNames=["output/%s-domain.xmf" % SIM_NAME])
RenameSource("%s-domain" % SIM_NAME, dataDomain)

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
view = GetActiveViewOrCreate('RenderView')

# Gray wireframe for undeformed domain.
domainDisplay = Show(dataDomain, view)
domainDisplay.Representation = 'Wireframe'
domainDisplay.AmbientColor = [0.5, 0.5, 0.5]

# Warp domain to show deformation
warp = WarpByVector(Input=dataDomain)
warp.Vectors = ['POINTS', 'displacement']
warp.ScaleFactor = DISPLACEMENT_SCALE

warpDisplay = Show(warp, view)
ColorBy(warpDisplay, ('POINTS', 'displacement', 'X'))
warpDisplay.RescaleTransferFunctionToDataRange(True)
warpDisplay.SetScalarBarVisibility(view, True)
warpDisplay.SetRepresentationType('Surface With Edges')

# Rescale color and/or opacity maps used to exactly fit the current data range
displacementLUT = GetColorTransferFunction('displacement')
warpDisplay.RescaleTransferFunctionToDataRange(False, False)
# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(displacementLUT, warpDisplay)

view.ResetCamera()
view.Update()

Render()

# Uncomment if running from shell outside ParaView.
#Interact()


# End of file
