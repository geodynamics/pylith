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
# EXODUS_FILE = os.path.join(os.environ["HOME"], "src", "pylith", "examples", "3d", "subduction", "mesh", "mesh_tet.exo")

DEFAULTS = {
    "EXODUS_FILE": "mesh_tri.exo",
    "SHOW_NODESETS": True,
    "SHOW_QUALITY": False,
    "QUALITY_METRIC": "Condition",
    "QUALITY_THRESHOLD": 2.0,
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

def plot_quality(dataDomain, parameters):
    """Plot the mesh quality. Only cells with a quality greater (poorer
    quality) than the threshold are displayed.
    """

    quality = MeshQuality(Input=dataDomain)
    quality.TriangleQualityMeasure = parameters.quality_metric
    quality.QuadQualityMeasure = parameters.quality_metric
    quality.TetQualityMeasure = parameters.quality_metric


    view = GetActiveViewOrCreate('RenderView')

    # Threshold
    threshold = Threshold(Input=quality)
    threshold.Scalars = ['CELLS', 'Quality']
    threshold.ThresholdRange = [parameters.quality_threshold, 1.0e+3]

    thresholdDisplay = Show(threshold, view)
    thresholdDisplay.Representation = 'Surface'
    thresholdDisplay.ColorArrayName = ['CELLS', 'Quality']
    thresholdDisplay.SetScalarBarVisibility(view, True)

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
    axes.ScaleFactor = 1.0e+5

    axesDisplay = Show(axes, view)
    axesDisplay.SetRepresentationType('Wireframe')
    axesDisplay.LineWidth = 4.0
    axesDisplay.SetScalarBarVisibility(view, False)
    axesDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]

    if parameters.show_nodesets:
        plot_nodesets(dataDomain)
    
    if parameters.show_quality:
        plot_quality(dataDomain, parameters)
    
    view.ResetCamera()
    Render()

class Parameters(object):
    """Object for managing default values and overriding them from the
    current Python shell.
    """
    keys = ("EXODUS_FILE", "SHOW_NODESETS", "SHOW_QUALITY", "QUALITY_METRIC", "QUALITY_THRESHOLD",)
    
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
    parser.add_argument("--hide-quality", action="store_false", dest="show_quality", default=DEFAULTS["SHOW_QUALITY"])
    parser.add_argument("--quality-metric", action="store", dest="quality_metric", default=DEFAULTS["QUALITY_METRIC"])
    parser.add_argument("--quality-threshold", action="store", type=float, dest="quality_threshold", default=DEFAULTS["QUALITY_THRESHOLD"])
    parser.add_argument("--screenshot", action="store", dest="screenshot")
    args = parser.parse_args()

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
