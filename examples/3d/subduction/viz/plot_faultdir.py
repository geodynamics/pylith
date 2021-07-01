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
    "FIELD": "normal_dir",
    "FAULTS": ["fault-slab"],
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
    view = GetActiveViewOrCreate('RenderView')

    # Gray wireframe for undeformed domain.
    domainDisplay = Show(dataDomain, view)
    domainDisplay.Representation = 'Wireframe'
    domainDisplay.AmbientColor = [0.5, 0.5, 0.5]

    # Read fault data
    dataFaults = []
    for fault in parameters.faults:
        filename = os.path.join(parameters.output_dir, "%s-%s_info.xmf" % (parameters.sim, fault))
        if not os.path.isfile(filename):
            raise IOError("File '%s' does not exist." % filename)
        data = XDMFReader(FileNames=[filename])
        RenameSource("%s-%s" % (parameters.sim, fault), data)
        dataFaults.append(data)

    groupFaults = GroupDatasets(Input=dataFaults)

    faultDisplay = Show(groupFaults, view)
    faultDisplay.SetRepresentationType('Surface With Edges')
    faultDisplayProperties = GetDisplayProperties(groupFaults, view=view)
    faultDisplayProperties.DiffuseColor = [0.25, 0.25, 1.0]
    
    # Add arrows to show displacement vectors.
    glyph = Glyph(Input=groupFaults, GlyphType="Arrow")
    glyph.Vectors = ["POINTS", parameters.field]
    glyph.GlyphMode = "All Points"

    glyphDisplay = Show(glyph, view)
    glyphDisplay.Representation = "Surface"
    
    view.ResetCamera()
    view.Update()
    Render()

class Parameters(object):
    keys = ("OUTPUT_DIR", "SIM", "FIELD", "FAULTS")
    
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
    parser.add_argument("--field", action="store", dest="field", default=DEFAULTS["FIELD"])
    args = parser.parse_args()

    if args.faults:
        args.faults = args.faults.split(",")
    else:
        args.faults = DEFAULTS["FAULTS"]
    visualize(args)
    Interact()

else:
    # Running inside the ParaView GUI

    visualize(Parameters())


# End of file
