#!/usr/bin/env nemesis
#
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
"""Create spatial database for power law bulk rheology for axial traction
problem with n=1 that is consistent with the analytical solution for the
linear Maxwell bulk rheology.
"""

import numpy

from spatialdata.spatialdb.SimpleIOAscii import createWriter
from spatialdata.geocoords.CSCart import CSCart

import axialtraction_maxwell_soln as soln


def run():
    """Generate databases for solution at t=0.
    """
    time = numpy.array([0.0], dtype=numpy.float64)
    xyz = numpy.array([0.0, 0.0, 0.0], dtype=numpy.float64).reshape(1,3)
    cs = CSCart()
    cs.inventory.spaceDim = 3
    cs._configure()
    dataPowerLaw = {'points': xyz,
                    'coordsys': cs,
                    'data_dim': 0,
                    'values': [{'name': "density",
                                'units': "kg/m**3",
                                'data': numpy.array(soln.p_density).reshape(1)},
                               {'name': "vs",
                                'units': "m/s",
                                'data': numpy.array(soln.p_vs).reshape(1)},
                               {'name': "vp",
                                'units': "m/s",
                                'data': numpy.array(soln.p_vp).reshape(1)},
                               {'name': "power_law_reference_strain_rate",
                                'units': "1/s",
                                'data': numpy.array(soln.p_power_law_reference_strain_rate).reshape(1)},
                               {'name': "power_law_reference_stress",
                                'units': "Pa",
                                'data': numpy.array(soln.p_power_law_reference_stress).reshape(1)},
                               {'name': "power_law_exponent",
                                'units': "None",
                                'data': numpy.array(soln.p_power_law_exponent).reshape(1)},
                               {'name': "viscous_strain_xx",
                                'units': "None",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "viscous_strain_yy",
                                'units': "None",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "viscous_strain_zz",
                                'units': "None",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "viscous_strain_xy",
                                'units': "None",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "viscous_strain_yz",
                                'units': "None",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "viscous_strain_xz",
                                'units': "None",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "deviatoric_stress_xx",
                                'units': "Pa",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "deviatoric_stress_yy",
                                'units': "Pa",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "deviatoric_stress_zz",
                                'units': "Pa",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "deviatoric_stress_xy",
                                'units': "Pa",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "deviatoric_stress_yz",
                                'units': "Pa",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "deviatoric_stress_xz",
                                'units': "Pa",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "reference_stress_xx",
                                'units': "Pa",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "reference_stress_yy",
                                'units': "Pa",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "reference_stress_zz",
                                'units': "Pa",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "reference_stress_xy",
                                'units': "Pa",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "reference_stress_yz",
                                'units': "Pa",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "reference_stress_xz",
                                'units': "Pa",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "reference_strain_xx",
                                'units': "None",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "reference_strain_yy",
                                'units': "None",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "reference_strain_zz",
                                'units': "None",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "reference_strain_xy",
                                'units': "None",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "reference_strain_yz",
                                'units': "None",
                                'data': numpy.array(0.0).reshape(1)},
                               {'name': "reference_strain_xz",
                                'units': "None",
                                'data': numpy.array(0.0).reshape(1)},
                               ]}
    (dispz, stress, devStress, strain, devStrain, maxwellVisStrain, powerLawVisStrain) = soln.AnalyticalSoln(time, xyz)
    dataPowerLaw['values'][6]['data'] = numpy.array(powerLawVisStrain[0,0]).reshape(1)
    dataPowerLaw['values'][7]['data'] = numpy.array(powerLawVisStrain[0,1]).reshape(1)
    dataPowerLaw['values'][8]['data'] = numpy.array(powerLawVisStrain[0,2]).reshape(1)
    dataPowerLaw['values'][12]['data'] = numpy.array(devStress[0,0]).reshape(1)
    dataPowerLaw['values'][13]['data'] = numpy.array(devStress[0,1]).reshape(1)
    dataPowerLaw['values'][14]['data'] = numpy.array(devStress[0,2]).reshape(1)

    ioPowerLaw = createWriter('mat_powerlaw_n1.spatialdb')
    ioPowerLaw.write(dataPowerLaw)


# End of file
