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
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

# Plot the domain, colored by materials. Optionally, plot nodesets and
# mesh quality.
#
# User-specified parameters.
#
# To use different values, overwrite them in the ParaView Python shell
# or on the command line. For example, set EXODUS_FILE to the absolute
# path for the Exodus file if not starting ParaView from the terminal
# shell where you ran PyLith:
#
# import os
# EXODUS_FILE = os.path.join(os.environ["HOME"], "src", "pylith", "examples", "reverse-2d", "mesh", "mesh_tri.exo")

DEFAULTS = {
    "EXODUS_FILE": "mesh_tri.exo",
    "SHOW_NODESETS": True,
    }

# ----------------------------------------------------------------------
from paraview.simple import *
import os

def plot_nodesets(dataDomain):
    """Create animation with each nodeset displayed as a different
    frame. We use a Python script to set the nodeset to display.
    """
    # Nodeset information
    nodeSetInfo = dataDomain.GetProperty("NodeSetInfo").GetData()
    nodeSets = nodeSetInfo[::2]
    dataDomain.NodeSetArrayStatus = nodeSets[0]
    numNodeSets = len(nodeSets)

    view = GetActiveViewOrCreate('RenderView')

    nsLabel = Text()
    nsLabel.Text = "Nodeset: %s" % nodeSets[0]
    RenameSource("nodeset-label")
    labelDisplay = Show(nsLabel, view)
    labelDisplay.FontSize = 10

    scene = GetAnimationScene()
    scene.NumberOfFrames = numNodeSets
    scene.StartTime = 0
    scene.EndTime = float(numNodeSets-1)
    
    cue = PythonAnimationCue()
    cue.Script = """
from paraview.simple import *

def tick(self):
    scene = GetAnimationScene()
    i = int(scene.TimeKeeper.Time)
    
    domain = FindSource("domain")
    nodeSetInfo = domain.GetProperty("NodeSetInfo").GetData()
    nodeSets = nodeSetInfo[::2]
    nodeSet = nodeSets[i]
    domain.NodeSetArrayStatus = nodeSet

    label = FindSource("nodeset-label")
    label.Text = "Nodeset: %s" % nodeSet
"""
    scene.Cues.append(cue)
    
    return


def visualize(parameters):

    # Disable automatic camera reset on "Show"
    paraview.simple._DisableFirstRenderCameraReset()

    if not os.path.isfile(parameters.exodus_file):
        raise IOError("Exodus file '%s' does not exist." % parameters.exodus_file)
    dataDomain = ExodusIIReader(FileName=[parameters.exodus_file])
    RenameSource("domain", dataDomain)

    
    view = GetActiveViewOrCreate('RenderView')

    # Show domain, colored by block.
    domainDisplay = Show(dataDomain, view)
    ColorBy(domainDisplay, ("FIELD", "vtkBlockColors"))
    domainDisplay.RescaleTransferFunctionToDataRange(True, False)
    domainDisplay.SetScalarBarVisibility(view, False)
    domainDisplay.SetRepresentationType("Surface With Edges")
    domainDisplay.PointSize = 6.0
    domainDisplay.Opacity = 0.5

    # Add coordinate axes
    axes = Axes()
    axes.ScaleFactor = 1.0e+4

    axesDisplay = Show(axes, view)
    axesDisplay.SetRepresentationType('Wireframe')
    axesDisplay.LineWidth = 4.0
    axesDisplay.SetScalarBarVisibility(view, False)
    axesDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]

    if parameters.show_nodesets:
        plot_nodesets(dataDomain)
    
    view.ResetCamera()
    Render()

class Parameters(object):
    """Object for managing default values and overriding them from the
    current Python shell.
    """
    keys = ("EXODUS_FILE", "SHOW_NODESETS")
    
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
    parser.add_argument("--filename", action="store", dest="exodus_file", default=DEFAULTS["EXODUS_FILE"])
    parser.add_argument("--hide-nodesets", action="store_false", dest="show_nodesets", default=DEFAULTS["SHOW_NODESETS"])
    parser.add_argument("--screenshot", action="store", dest="screenshot")
    args = parser.parse_args()

    visualize(args)

    view = GetRenderView()
    view.CameraPosition = [0, 0, 282700]
    view.CameraFocalPoint = [0, 0, -229160]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.ViewSize = [1000, 1500]
    view.Update()

    if args.screenshot:
        WriteImage(args.screenshot)

    Interact()

else:
    # Running inside the ParaView GUI

    visualize(Parameters())


# End of file
