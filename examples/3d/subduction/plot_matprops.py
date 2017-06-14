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

# Plot the domain, colored by material properties.
#
# This Python script runs using pvpython or within the ParaView Python
# shell.

# User-specified parameters

# Root name for simulation.
SIM_NAME = "step01"

# Material property and materials to plot.
INFO_FIELD = "mu"
MATERIALS = ["crust", "mantle", "wedge", "slab"]

# ----------------------------------------------------------------------
from paraview.simple import *
# Disable automatic camera reset on "Show"
paraview.simple._DisableFirstRenderCameraReset()


dataAll = []
# Read data
for material in MATERIALS:
    dataMaterial = XDMFReader(FileNames=["output/%s-%s_info.xmf" % (SIM_NAME, material)])
    RenameSource("%s-%s" % (SIM_NAME, material), dataMaterial)
    dataAll.append(dataMaterial)
groupMaterials = GroupDatasets(Input=dataAll)

view = GetActiveViewOrCreate('RenderView')

# Show domain, colored by magnitude of displacement vector.
materialDisplay = Show(groupMaterials, view)
ColorBy(materialDisplay, ("CELLS", INFO_FIELD))
materialDisplay.RescaleTransferFunctionToDataRange(True)
materialDisplay.SetScalarBarVisibility(view, True)
materialDisplay.SetRepresentationType("Surface With Edges")

# Rescale color and/or opacity maps used to exactly fit the current data range
materialLUT = GetColorTransferFunction(INFO_FIELD)
materialDisplay.RescaleTransferFunctionToDataRange(False, False)
# Update scalar bar component title.
UpdateScalarBarsComponentTitle(materialLUT, materialDisplay)


view.ResetCamera()

Render()

# Uncomment if running from shell outside ParaView.
#Interact()


# End of file
