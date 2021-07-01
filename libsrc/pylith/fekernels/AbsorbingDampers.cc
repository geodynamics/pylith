/* -*- C -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University at Buffalo
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2021 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/AbsorbingDampers.hh"

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for computing value from parameters for time-dependent boundary conditions.
 * ======================================================================
 */

/* ----------------------------------------------------------------------
 * g0 function for absorbing dampers boundary condition.
 *
 * g_0(x)
 */
void
pylith::fekernels::AbsorbingDampers::g0(const PylithInt dim,
                                        const PylithInt numS,
                                        const PylithInt numA,
                                        const PylithInt sOff[],
                                        const PylithInt sOff_x[],
                                        const PylithScalar s[],
                                        const PylithScalar s_t[],
                                        const PylithScalar s_x[],
                                        const PylithInt aOff[],
                                        const PylithInt aOff_x[],
                                        const PylithScalar a[],
                                        const PylithScalar a_t[],
                                        const PylithScalar a_x[],
                                        const PylithReal t,
                                        const PylithReal x[],
                                        const PylithReal n[],
                                        const PylithInt numConstants,
                                        const PylithScalar constants[],
                                        PylithScalar g0[]) {
    assert(2 == dim || 3 == dim);

    const PylithInt _numA = 3;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    const PylithInt i_density = aOff[0];
    const PylithInt i_vp = aOff[1];
    const PylithInt i_vs = aOff[2];

    const PylithInt _numS = 2;
    assert(sOff);
    assert(s);
    assert(numS >= _numS);
    const PylithInt i_vel = sOff[1];

    const PylithScalar density = a[i_density];
    const PylithScalar vp = a[i_vp];
    const PylithScalar vs = a[i_vs];

    PylithScalar velN[3];
    PylithScalar velT[3];
    PylithScalar velNMag = 0;
    for (PylithInt i = 0; i < dim; ++i) {
        velNMag += s[i_vel+i] * n[i];
    } // for
    for (PylithInt i = 0; i < dim; ++i) {
        velN[i] = velNMag * n[i];
    } // for
    for (PylithInt i = 0; i < dim; ++i) {
        velT[i] = s[i_vel+i] - velN[i];
    } // for

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] -= density * (vs * velT[i] + vp * velN[i]);
    } // for
} // g0


// End of file
