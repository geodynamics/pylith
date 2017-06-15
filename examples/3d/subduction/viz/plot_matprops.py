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


# User-specified parameters.
#
# These are used if running from within the ParaView GUI via the
# Python shell or as defaults if running outside the ParaView GUI via
# pvpython.

# Root name for simulation.
SIM_NAME = "step01"

# Material property and materials to plot.
INFO_FIELD = "mu"
MATERIALS = ["crust", "mantle", "wedge", "slab"]

# ----------------------------------------------------------------------
from paraview.simple import *

def visualize(sim, field, materials):

    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()


    dataAll = []
    # Read data
    for material in materials:
        dataMaterial = XDMFReader(FileNames=["output/%s-%s_info.xmf" % (sim, material)])
        RenameSource("%s-%s" % (sim, material), dataMaterial)
        dataAll.append(dataMaterial)
    groupMaterials = GroupDatasets(Input=dataAll)

    view = GetActiveViewOrCreate('RenderView')

    # Show domain, colored by magnitude of displacement vector.
    materialDisplay = Show(groupMaterials, view)
    ColorBy(materialDisplay, ("CELLS", field))
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

# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Running from outside the ParaView GUI via pvpython
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", action="store", dest="sim")
    parser.add_argument("--field", action="store", dest="field")
    parser.add_argument("--materials", action="store", dest="materials")
    args = parser.parse_args()

    sim = args.sim
    field = args.field
    if args.materials:
        materials = args.materials.split(",")
    else:
        materials = None
    
    if sim is None:
        sim = SIM_NAME
    if field is None:
        field = INFO_FIELD
    if materials is None:
        materials = MATERIALS
    visualize(sim, field, materials)
    Interact()

else:
    # Running inside the ParaView GUI

    visualize(SIM_NAME, INFO_FIELD, MATERIALS)


# End of file
