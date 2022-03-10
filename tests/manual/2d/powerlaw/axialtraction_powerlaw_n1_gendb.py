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
#
# @file tests/manual/2d/powerlaw/axialtraction_powerlaw_n1_gendb.py
#
# @brief Python script to generate spatial database of material properties
# and state variables for the axial traction test.

import numpy

fileNamePowerLaw = 'mat_powerlaw_n1' + '.spatialdb'

# Initial state from elastic solution.
time = numpy.array([0.0], dtype=numpy.float64)

def run(params):
    """Generate databases for solution at t=0.
    """
    from axialtraction_maxwell_soln import AnalyticalSoln
    from spatialdata.spatialdb.SimpleIOAscii import createWriter

    xy = numpy.array([0.0, 0.0], dtype=numpy.float64).reshape(1,2)
    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.inventory.spaceDim = 2
    cs._configure()
    dataPowerLaw = {'points': xy,
                    'coordsys': cs,
                    'data_dim': 0,
                    'values': [{'name': "density",
                                'units': "kg/m**3",
                                'data': numpy.array(params['p_density']).reshape(1)},
                               {'name': "vs",
                                'units': "m/s",
                                'data': numpy.array(params['p_vs']).reshape(1)},
                               {'name': "vp",
                                'units': "m/s",
                                'data': numpy.array(params['p_vp']).reshape(1)},
                               {'name': "power_law_reference_strain_rate",
                                'units': "1/s",
                                'data': numpy.array(params['p_power_law_reference_strain_rate']).reshape(1)},
                               {'name': "power_law_reference_stress",
                                'units': "Pa",
                                'data': numpy.array(params['p_power_law_reference_stress']).reshape(1)},
                               {'name': "power_law_exponent",
                                'units': "None",
                                'data': numpy.array(params['p_power_law_exponent']).reshape(1)},
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
                               ]}
    (dispy, stress, devStress, strain, devStrain, maxwellVisStrain, powerLawVisStrain) = AnalyticalSoln(time, xy, params)
    dataPowerLaw['values'][6]['data'] = numpy.array(powerLawVisStrain[0,0]).reshape(1)
    dataPowerLaw['values'][7]['data'] = numpy.array(powerLawVisStrain[0,1]).reshape(1)
    dataPowerLaw['values'][8]['data'] = numpy.array(powerLawVisStrain[0,2]).reshape(1)
    dataPowerLaw['values'][12]['data'] = numpy.array(devStress[0,0]).reshape(1)
    dataPowerLaw['values'][13]['data'] = numpy.array(devStress[0,1]).reshape(1)
    dataPowerLaw['values'][14]['data'] = numpy.array(devStress[0,2]).reshape(1)

    ioPowerLaw = createWriter(fileNamePowerLaw)
    ioPowerLaw.write(dataPowerLaw)
        
    return


# End of file
