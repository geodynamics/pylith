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
SIM_NAME = "step02"

# Names of faults for output files.
FAULTS = ["fault-slab"]
FIELD = "strike_dir"

# ----------------------------------------------------------------------
from paraview.simple import *
import os

def visualize(sim, faults, direction):
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
        filename = "output/%s-%s_info.xmf" % (sim, fault)
        if not os.path.isfile(filename):
            raise IOError("File '%s' does not exist." % filename)
        data = XDMFReader(FileNames=[filename])
        RenameSource("%s-%s" % (sim, fault), data)
        dataFaults.append(data)

    groupFaults = GroupDatasets(Input=dataFaults)

    faultDisplay = Show(groupFaults, view)
    faultDisplay.SetRepresentationType('Surface With Edges')
    faultDisplayProperties = GetDisplayProperties(groupFaults, view=view)
    faultDisplayProperties.DiffuseColor = [0.25, 0.25, 1.0]
    
    # Add arrows to show displacement vectors.
    glyph = Glyph(Input=groupFaults, GlyphType="Arrow")
    glyph.Vectors = ["POINTS", direction]
    glyph.GlyphMode = "All Points"

    glyphDisplay = Show(glyph, view)
    glyphDisplay.Representation = "Surface"
    
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
    parser.add_argument("--direction", action="store", dest="direction", default=FIELD)
    args = parser.parse_args()

    if args.faults:
        faults = args.faults.split(",")
    else:
        faults = FAULTS
    visualize(args.sim, faults, args.direction)
    Interact()

else:
    # Running inside the ParaView GUI

    visualize(SIM_NAME, FAULTS, FIELD)


# End of file
