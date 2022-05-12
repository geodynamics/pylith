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

#include "pylith/fekernels/GaussianWavelet.hh"

#include <cassert> // USES assert()
#include <iostream> // debugging.

// =====================================================================================================================
// Kernels for the Ricker Source Time Function in 2D.
// =====================================================================================================================

// ----------------------------------------------------------------------
// g1 entry function for velocity equation
void
pylith::fekernels::GaussianWaveletPlaneStrain::g1v(const PylithInt dim,
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
    const PylithInt i_gaussianwaveletCenterFrequency = numA - 1;
    
    const PylithScalar* momentTensor = &a[aOff[i_momentTensor]];
    const PylithScalar timeDelay = a[aOff[i_timeDelay]];
    const PylithScalar f0 = a[aOff[i_gaussianwaveletCenterFrequency]];

    // GaussianWavelet source time function (time domain)

    PylithScalar rt = t - timeDelay;
    PylithScalar gaussianwavelet = PetscExpReal( (PETSC_PI*PETSC_PI * f0*f0) * rt*rt) / (2.0 * (PETSC_PI*PETSC_PI * f0*f0) );

    for (PylithInt i = 0; i < dim*dim; ++i) {
        g1[i] -= momentTensor[i] * gaussianwavelet;
    } // for
} // g1v


// =====================================================================================================================
// Kernels for the Ricker Source Time Function in 3D.
// =====================================================================================================================

// ----------------------------------------------------------------------
// g1 entry function for velocity equation
void
pylith::fekernels::GaussianWavelet3D::g1v(const PylithInt dim,
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
    const PylithInt i_gaussianwaveletCenterFrequency = numA - 1;
    
    const PylithScalar* momentTensor = &a[aOff[i_momentTensor]];
    const PylithScalar timeDelay = a[aOff[i_timeDelay]];
    const PylithScalar f0 = a[aOff[i_gaussianwaveletCenterFrequency]];

    // GaussianWavelet source time function (time domain)

    PylithScalar rt = t - timeDelay;
    PylithScalar gaussianwavelet = PetscExpReal( (PETSC_PI*PETSC_PI * f0*f0) * rt*rt) / (2.0 * (PETSC_PI*PETSC_PI * f0*f0) );

    for (PylithInt i = 0; i < dim*dim; ++i) {
        g1[i] -= momentTensor[i] * gaussianwavelet;
    } // for
} // g1v


// End of file
