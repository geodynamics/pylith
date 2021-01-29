#!/usr/bin/env python
import sys

if len(sys.argv) != 2:
    raise ValueError("usage: plot_data.py SIM")

sim = sys.argv[1]

if not sim in ['sin']:
    raise ValueError("Unknown sim '%s'." % sim)


# ======================================================================
import tables
import pylab
import numpy

# ----------------------------------------------------------------------
def get_data(sim, dim, form):
    h5 = tables.openFile("output/%s%dd_%s-statevars.h5" % (sim, dim, form), "r")
    plasticStrain = h5.root.cell_fields.plastic_strain[:]
    totalStrain = h5.root.cell_fields.total_strain[:]
    stress = h5.root.cell_fields.stress[:]
    timeStamps =  h5.root.time[:].ravel()

    # Scales for time and stress
    from pythia.pyre.units.time import year
    timescale = 1.0*year
    stressscale = 1.0e+6

    h5.close()
    return {'time': timeStamps/timescale.value,
            'plastic_strain': plasticStrain,
            'total_strain': totalStrain,
            'stress': stress/stressscale}



# ----------------------------------------------------------------------
quasi2D = get_data(sim, dim=2, form='quasistatic')
dyn2D = get_data(sim, dim=2, form='dynamic')
quasi3D = get_data(sim, dim=3, form='quasistatic')
dyn3D = get_data(sim, dim=3, form='dynamic')

# Plot
fig = pylab.Figure()

# Total Strain
ax = pylab.subplot(3, 3, 1)
pylab.plot(quasi2D['time'], quasi2D['total_strain'][:, 0, 0], 'b-',
           dyn2D['time'], dyn2D['total_strain'][:, 0, 0], 'r--',
           quasi3D['time'], quasi3D['total_strain'][:, 0, 0], 'g-',
           dyn3D['time'], dyn3D['total_strain'][:, 0, 0], 'm--',
    )
ax.set_title("Total Strain")

ax = pylab.subplot(3, 3, 4)
pylab.plot(quasi2D['time'], quasi2D['total_strain'][:, 0, 1], 'b-',
           dyn2D['time'], dyn2D['total_strain'][:, 0, 1], 'r--',
           quasi3D['time'], quasi3D['total_strain'][:, 0, 1], 'g-',
           dyn3D['time'], dyn3D['total_strain'][:, 0, 1], 'm--',
    )

ax = pylab.subplot(3, 3, 7)
pylab.plot(quasi2D['time'], quasi2D['total_strain'][:, 0, 2], 'b-',
           dyn2D['time'], dyn2D['total_strain'][:, 0, 2], 'r--',
           quasi3D['time'], quasi3D['total_strain'][:, 0, 3], 'g-',
           dyn3D['time'], dyn3D['total_strain'][:, 0, 3], 'm--',
    )
ax.set_xlabel("Time (yr)")

# Plastic Strain
ax = pylab.subplot(3, 3, 2)
pylab.plot(quasi2D['time'], quasi2D['plastic_strain'][:, 0, 0], 'b-',
           dyn2D['time'], dyn2D['plastic_strain'][:, 0, 0], 'r--',
           quasi3D['time'], quasi3D['plastic_strain'][:, 0, 0], 'g-',
           dyn3D['time'], dyn3D['plastic_strain'][:, 0, 0], 'm--',
    )
ax.set_title("Plastic Strain")

ax = pylab.subplot(3, 3, 5)
pylab.plot(quasi2D['time'], quasi2D['plastic_strain'][:, 0, 1], 'b-',
           dyn2D['time'], dyn2D['plastic_strain'][:, 0, 1], 'r--',
           quasi3D['time'], quasi3D['plastic_strain'][:, 0, 1], 'g-',
           dyn3D['time'], dyn3D['plastic_strain'][:, 0, 1], 'm--',
    )

ax = pylab.subplot(3, 3, 8)
pylab.plot(quasi2D['time'], quasi2D['plastic_strain'][:, 0, 2], 'b-',
           dyn2D['time'], dyn2D['plastic_strain'][:, 0, 2], 'r--',
           quasi3D['time'], quasi3D['plastic_strain'][:, 0, 2], 'g-',
           dyn3D['time'], dyn3D['plastic_strain'][:, 0, 2], 'm--',
    )
ax.set_xlabel("Time (yr)")


# Stress
ax = pylab.subplot(3, 3, 3)
pylab.plot(quasi2D['time'], quasi2D['stress'][:, 0, 0], 'b-',
           dyn2D['time'], dyn2D['stress'][:, 0, 0], 'r--',
           quasi3D['time'], quasi3D['stress'][:, 0, 0], 'g-',
           dyn3D['time'], dyn3D['stress'][:, 0, 0], 'm--',
    )
ax.set_title("Stress (MPa)")

ax = pylab.subplot(3, 3, 6)
pylab.plot(quasi2D['time'], quasi2D['stress'][:, 0, 1], 'b-',
           dyn2D['time'], dyn2D['stress'][:, 0, 1], 'r--',
           quasi3D['time'], quasi3D['stress'][:, 0, 1], 'b-',
           dyn3D['time'], dyn3D['stress'][:, 0, 1], 'r--',
    )

ax = pylab.subplot(3, 3, 9)
pylab.plot(quasi2D['time'], quasi2D['stress'][:, 0, 2], 'b-',
           dyn2D['time'], dyn2D['stress'][:, 0, 2], 'r--',
           quasi3D['time'], quasi3D['stress'][:, 0, 3], 'b-',
           dyn3D['time'], dyn3D['stress'][:, 0, 3], 'r--',
    )
ax.set_xlabel("Time (yr)")



pylab.show()
