#!/usr/bin/env nemesis
"""
This script generates a plot showing slip or fault tractions.
"""

# The code requires the numpy, h5py, and matplotlib packages.
import numpy
import h5py
import matplotlib.pyplot as pyplot

# ----------------------------------------------------------------------
import sys

plot = sys.argv[1]
if not plot in ['step01_slip', 
                'step01_stress',
                'step04_bg',
                'step04_initial',
                ]:
    raise ValueError("Unknown plot '%s'." % plot)

# ----------------------------------------------------------------------
def getStep01():
  """Function to get slip, tractions, and fault coordinates from step01.
  """

  # Open solution file and get slip and coordinates.
  h5 = h5py.File("output/step01-fault.h5", "r")
  vertices = h5['geometry/vertices'][:]
  slip = h5['vertex_fields/slip'][0,:, 0].squeeze()
  traction_change = h5['vertex_fields/traction_change'][0,:,:].squeeze()
  h5.close()

  # Sort by y-coordinate (elevation).
  ids = numpy.argsort(vertices[:, 1])
  vertices = vertices[ids,:]
  slip = slip[ids,:]
  traction_change = traction_change[ids]

  return (vertices, slip, traction_change)

# ======================================================================
(vertices, slip, traction_change) = getStep01()

# Background stress field (from afterslip_tractions.py)
density = 2900.0
gacc = 9.80665
mu_s = 0.6

# Background normal tractions are overburden and compressive
# (negative, y is negative)
traction_bg_normal = density*gacc*(vertices[:, 1])

# Background shear tractions are reverse (in 2-D right-lateral is negative)
# because the normal tractions are negative.
traction_bg_shear = mu_s*traction_bg_normal

# ----------------------------------------------------------------------
figure = pyplot.figure(figsize=(5.0, 5.0), facecolor='white', dpi=150)
figure.set_facecolor('white')

axes = figure.add_axes([0.15, 0.1, 0.80, 0.87])

if plot == "step01_slip":
  axes.plot(slip, vertices[:, 1]/1.0e+3)
  axes.set_xlabel("Slip (m)")
  axes.set_ylabel("Elevation (km)")
  axes.set_ylim((-60.0, 0.0))

elif plot == "step01_stress":
  axes.plot(traction_change[:, 0]/1.0e+6, vertices[:, 1]/1.0e+3)
  axes.set_xlabel("Traction Change (MPa)")
  axes.set_ylabel("Elevation (km)")
  axes.plot([0, 0], [0, -600], linestyle='--', color='gray', alpha=0.5)
  axes.set_ylim((-60.0, 0.0))

elif plot == "step04_bg":
  axes.plot(traction_bg_normal/1.0e+6, vertices[:, 1]/1.0e+3, 
            color='green', linestyle='--')
  axes.plot(traction_bg_shear/1.0e+6, vertices[:, 1]/1.0e+3, color='blue')
  axes.set_xlabel("Traction (MPa)")
  axes.set_ylabel("Elevation (km)")
  axes.set_ylim((-60.0, 0.0))
  axes.set_xlim((-2000.0, 0.0))
  axes.legend(('Normal', 'Shear'), loc='upper left')


elif plot == "step04_initial":
  traction_initial_normal = traction_bg_normal + traction_change[:, 1]
  traction_initial_shear = traction_bg_shear + traction_change[:, 0]
  traction_friction = -2.0e+6 + mu_s*traction_initial_normal

  axes.plot(traction_initial_normal/1.0e+6, vertices[:, 1]/1.0e+3, 
            color='green', linestyle='--')
  axes.plot(traction_initial_shear/1.0e+6, vertices[:, 1]/1.0e+3, color='blue')
  axes.plot(traction_friction/1.0e+6, vertices[:, 1]/1.0e+3, 
            color='red', linestyle='-.')
  axes.set_xlabel("Traction (MPa)")
  axes.set_ylabel("Elevation (km)")
  axes.set_ylim((-60.0, 0.0))
  axes.set_xlim((-2000.0, 0.0))
  axes.legend(('Normal', 'Shear', "Failure"), loc='upper left')


pyplot.show()
pyplot.savefig(plot)
