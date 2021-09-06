#!/usr/bin/env nemesis
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/manual/3d/powerlaw/cylinderpres_gendb.py
#
# @brief Python script to generate spatial database with auxiliary
# fields for test with pressurized power-law viscoelastic cylinder and initial
# stress.

import numpy

def generate_refstate_db(vertices, spatialdbFile):
    """Generate the database.
        """

    npts = vertices.shape[0]
    from cylinderpres_soln import AnalyticalSoln
    from cylinderpres_soln import p_density, p_vs, p_vp, \
        p_powerLawExponent, p_powerLawReferenceStrainRate, p_powerLawReferenceStress
    soln = AnalyticalSoln()
    stress = soln.stress(vertices)
    zeroes = numpy.zeros(npts, dtype=numpy.float64)

    from spatialdata.geocoords.CSCart import CSCart
    cs = CSCart()
    cs.inventory.spaceDim = 3
    cs._configure()
    data = {
        'points': vertices,
        'coordsys': cs,
        'data_dim': 3,
        'values': [
            {
                'name': "density",
                'units': "kg/m**3",
                'data': p_density * numpy.ones((npts,)),
            }, {
                'name': "vs",
                'units': "m/s",
                'data': p_vs * numpy.ones((npts,)),
            }, {
                'name': "vp",
                'units': "m/s",
                'data': p_vp * numpy.ones((npts,)),
            }, {
                'name': "power_law_reference_strain_rate",
                'units': "1/s",
                'data': p_powerLawReferenceStrainRate * numpy.ones((npts,)),
            }, {
                'name': "power_law_reference_stress",
                'units': "Pa",
                'data': p_powerLawReferenceStress * numpy.ones((npts,)),
            }, {
                'name': "power_law_exponent",
                'units': "None",
                'data': p_powerLawExponent * numpy.ones((npts,)),
            }, {
                'name': "viscous_strain_xx",
                'units': "None",
                'data': zeroes,
            }, {
                'name': "viscous_strain_yy",
                'units': "None",
                'data': zeroes,
            }, {
                'name': "viscous_strain_zz",
                'units': "None",
                'data': zeroes,
            }, {
                'name': "viscous_strain_xy",
                'units': "None",
                'data': zeroes,
            }, {
                'name': "viscous_strain_yz",
                'units': "None",
                'data': zeroes,
            }, {
                'name': "viscous_strain_xz",
                'units': "None",
                'data': zeroes,
            }, {
                'name': "deviatoric_stress_xx",
                'units': "Pa",
                'data': zeroes,
            }, {
                'name': "deviatoric_stress_yy",
                'units': "Pa",
                'data': zeroes,
            }, {
                'name': "deviatoric_stress_zz",
                'units': "Pa",
                'data': zeroes,
            }, {
                'name': "deviatoric_stress_xy",
                'units': "Pa",
                'data': zeroes,
            }, {
                'name': "deviatoric_stress_yz",
                'units': "Pa",
                'data': zeroes,
            }, {
                'name': "deviatoric_stress_xz",
                'units': "Pa",
                'data': zeroes,
            }, {
                'name': "reference_stress_xx",
                'units': "Pa",
                'data': numpy.ravel(stress[0, :, 0]),
            }, {
                'name': "reference_stress_yy",
                'units': "Pa",
                'data': numpy.ravel(stress[0, :, 1]),
            }, {
                'name': "reference_stress_zz",
                'units': "Pa",
                'data': numpy.ravel(stress[0, :, 2]),
            }, {
                'name': "reference_stress_yz",
                'units': "Pa",
                'data': numpy.ravel(stress[0, :, 4]),
            }, {
                'name': "reference_stress_xz",
                'units': "Pa",
                'data': numpy.ravel(stress[0, :, 5]),
            }, {
                'name': "reference_stress_xy",
                'units': "Pa",
                'data': numpy.ravel(stress[0, :, 3]),
            }, {
                'name': "reference_strain_xx",
                'units': "none",
                'data': zeroes,
            }, {
                'name': "reference_strain_yy",
                'units': "none",
                'data': zeroes,
            }, {
                'name': "reference_strain_zz",
                'units': "none",
                'data': zeroes,
            }, {
                'name': "reference_strain_xy",
                'units': "none",
                'data': zeroes,
            }, {
                'name': "reference_strain_yz",
                'units': "none",
                'data': zeroes,
            }, {
                'name': "reference_strain_xz",
                'units': "none",
                'data': zeroes,
            }
        ]
    }

    from spatialdata.spatialdb.SimpleIOAscii import createWriter
    io = createWriter(spatialdbFile)
    io.write(data)
    return


# End of file
