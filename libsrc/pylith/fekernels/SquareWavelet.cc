/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University of Chicago
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2015 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/SquareWavelet.hh"

#include <cassert> // USES assert()
#include <iostream> // debugging.

// =====================================================================================================================
// Kernels for the Square Source Time Function in 2D.
// =====================================================================================================================

// ----------------------------------------------------------------------
// g1 entry function for velocity equation
void
pylith::fekernels::SquareWaveletPlaneStrain::g1v(const PylithInt dim,
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
                                                 const PylithScalar x[],
                                                 const PylithInt numConstants,
                                                 const PylithScalar constants[],
                                                 PylithScalar g1[]) {
    assert(sOff);
    assert(s);
    assert(g1);

    PylithInt _dim = 2;

    // Incoming re-packed solution field.

    // Incoming re-packed auxiliary field.
    const PylithInt i_momentTensor = 0;
    const PylithInt i_timeDelay = 1;
    const PylithInt i_squarewaveletCenterFrequency = numA - 1;

    const PylithScalar* momentTensor = &a[aOff[i_momentTensor]];
    const PylithScalar timeDelay = a[aOff[i_timeDelay]];
    const PylithScalar squarewaveletCenterFrequency = a[aOff[i_squarewaveletCenterFrequency]];

    // SquareWavelet source time function (time domain)

    PylithScalar rt = t - timeDelay;
    PylithScalar squarewavelet = (rt >= 0.0) ? 1.0 : 0.0;

    // PetscPrintf(PETSC_COMM_WORLD, "timeDelay %f\n", (double)timeDelay);
    // PetscPrintf(PETSC_COMM_WORLD, "t %f\n", (double)t);
    // PetscPrintf(PETSC_COMM_WORLD, "Center Freq %f\n", (double)squarewaveletCenterFrequency);
    // PetscPrintf(PETSC_COMM_WORLD, "squareWavelet %f\n", (double)squarewavelet);

    for (PylithInt i = 0; i < dim*dim; ++i) {
        g1[i] -= momentTensor[i] * squarewavelet;
        // PetscPrintf(PETSC_COMM_WORLD, "g1[%i]: %f - square\n", (int)i, (double)g1[i]);
    } // for
} // g1v


// =====================================================================================================================
// Kernels for the Square Source Time Function in 3D.
// =====================================================================================================================

// ----------------------------------------------------------------------
// g1 entry function for velocity equation
void
pylith::fekernels::SquareWavelet3D::g1v(const PylithInt dim,
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
                                        const PylithScalar x[],
                                        const PylithInt numConstants,
                                        const PylithScalar constants[],
                                        PylithScalar g1[]) {
    assert(sOff);
    assert(s);
    assert(g1);

    PylithInt _dim = 3;

    // Incoming re-packed solution field.

    // Incoming re-packed auxiliary field.
    const PylithInt i_momentTensor = 0;
    const PylithInt i_timeDelay = 1;
    const PylithInt i_squarewaveletCenterFrequency = numA - 1;

    const PylithScalar* momentTensor = &a[aOff[i_momentTensor]];
    const PylithScalar timeDelay = a[aOff[i_timeDelay]];
    const PylithScalar squarewaveletCenterFrequency = a[aOff[i_squarewaveletCenterFrequency]];

    // SquareWavelet source time function (time domain)

    PylithScalar rt = t - timeDelay;
    PylithScalar squarewavelet = (rt >= 0.0) ? 1.0 : 0.0;

    for (PylithInt i = 0; i < dim*dim; ++i) {
        g1[i] -= momentTensor[i] * squarewavelet;
    } // for
} // g1v


// End of file
