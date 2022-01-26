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
# @file tests/manual/3d/powerlaw/axialtraction_powerlaw_n1_gendb.py
#
# @brief Python script to generate spatial database of material properties
# and state variables for the axial traction test.

import numpy

def run():
    """Generate databases for different time step sizes.
    """
    from axialtraction_maxwell_soln import AnalyticalSoln
    from axialtraction_maxwell_soln import p_density, p_vs, p_vp, p_viscosity
    from axialtraction_maxwell_soln import p_power_law_exponent, p_power_law_reference_strain_rate, p_power_law_reference_stress
    from spatialdata.spatialdb.SimpleIOAscii import createWriter

    # Time steps to use.
    dt = numpy.array([0.01, 0.02, 0.05, 0.1, 0.2], dtype=numpy.float64)
    dtStr = ['0.01', '0.02', '0.05', '0.1', '0.2']
    numSteps = dt.shape[0]
    tensorSize = 6

    xyz = numpy.array([0.0, 0.0, 0.0], dtype=numpy.float64).reshape(1,3)
    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.inventory.spaceDim = 3
    cs._configure()
    dataPowerLaw = {'points': xyz,
                    'coordsys': cs,
                    'data_dim': 0,
                    'values': [{'name': "density",
                                'units': "kg/m**3",
                                'data': numpy.array(p_density).reshape(1)},
                               {'name': "vs",
                                'units': "m/s",
                                'data': numpy.array(p_vs).reshape(1)},
                               {'name': "vp",
                                'units': "m/s",
                                'data': numpy.array(p_vp).reshape(1)},
                               {'name': "power_law_reference_strain_rate",
                                'units': "1/s",
                                'data': numpy.array(p_power_law_reference_strain_rate).reshape(1)},
                               {'name': "power_law_reference_stress",
                                'units': "Pa",
                                'data': numpy.array(p_power_law_reference_stress).reshape(1)},
                               {'name': "power_law_exponent",
                                'units': "None",
                                'data': numpy.array(p_power_law_exponent).reshape(1)},
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
    for stepNum in range(numSteps):
        fileNamePowerLaw = 'mat_powerlaw_n1_norefstate_dt' + dtStr[stepNum] + '.spatialdb'
        time = numpy.array(dt[stepNum], dtype=numpy.float64).reshape(1)
        (dispz, stress, devStress, strain, devStrain, maxwellVisStrain, powerLawVisStrain) = AnalyticalSoln(time, xyz)
        dataPowerLaw['values'][6]['data'] = numpy.array(powerLawVisStrain[0,0]).reshape(1)
        dataPowerLaw['values'][7]['data'] = numpy.array(powerLawVisStrain[0,1]).reshape(1)
        dataPowerLaw['values'][8]['data'] = numpy.array(powerLawVisStrain[0,2]).reshape(1)
        dataPowerLaw['values'][12]['data'] = numpy.array(devStress[0,0]).reshape(1)
        dataPowerLaw['values'][13]['data'] = numpy.array(devStress[0,1]).reshape(1)
        dataPowerLaw['values'][14]['data'] = numpy.array(devStress[0,2]).reshape(1)

        ioPowerLaw = createWriter(fileNamePowerLaw)
        ioPowerLaw.write(dataPowerLaw)
        
                               
    return


# ======================================================================
if __name__ == "__main__":
    run()


# End of file
