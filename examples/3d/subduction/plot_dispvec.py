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

# Plot the domain, colored by the magnitude of the displacement
# vector, with white arrows showing the displacement vectors.
#
# This Python script runs using pvpython or within the ParaView Python
# shell.

# User-specified parameters

# Root name for simulation.
SIM_NAME = "step01"

# Format of simulation output (choices=["vtk", "hdf5"], case insentitive)
DATA_FORMAT = "hdf5"

# Scale used to exaggerate deformation.
DISPLACEMENT_SCALE = 10.0e+3


# ----------------------------------------------------------------------
# We create sources, filters, etc using the servermanager so that we can set the
# name of the proxy shown in the pipeline using Register(registrationName=NAME).

from paraview.simple import *
# Disable automatic camera reset on "Show"
paraview.simple._DisableFirstRenderCameraReset()

# Read data
if DATA_FORMAT.lower() == "vtk":
    dataDomain = servermanager.sources.LegacyVTKReader(FileNames=['output/%s-domain_t0000000.vtk' % SIM_NAME])
elif DATA_FORMAT.lower() == "hdf5":
    dataDomain = servermanager.sources.XDMFReader(FileNames=["output/%s-domain.xmf" % SIM_NAME])
else:
     raise ValueError("Unknown file format '%s' when choosing reader in Python script." % DATA_FORMAT)
servermanager.Register(dataDomain, registrationName="%s-domain" % SIM_NAME)
view = GetActiveViewOrCreate('RenderView')

# Show undeformed domain, colored by magnitude of displacement vector.
domainDisplay = Show(dataDomain, view)
ColorBy(domainDisplay, ("POINTS", "displacement", "Magnitude"))
domainDisplay.RescaleTransferFunctionToDataRange(True)
domainDisplay.SetScalarBarVisibility(view, True)
domainDisplay.SetRepresentationType("Surface With Edges")

# Rescale color and/or opacity maps used to exactly fit the current data range
displacementLUT = GetColorTransferFunction('displacement')
domainDisplay.RescaleTransferFunctionToDataRange(False, False)
# Update scalar bar component title.
UpdateScalarBarsComponentTitle(displacementLUT, domainDisplay)


# Add arrows to show displacement vectors.
glyph = servermanager.filters.Glyph(Input=dataDomain, GlyphType="Arrow")
servermanager.Register(glyph, registrationName="%s-domain-glyph" % SIM_NAME)
glyph.Vectors = ["POINTS", "displacement"]
glyph.ScaleFactor = DISPLACEMENT_SCALE
glyph.ScaleMode = "vector"
glyph.GlyphMode = "All Points"

glyphDisplay = Show(glyph, view)
glyphDisplay.Representation = "Surface"

view.ResetCamera()

Render()

# Uncomment if running from shell outside ParaView.
#Interact()


# End of file
