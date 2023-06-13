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

#include "pylith/fekernels/RickerWavelet.hh"

#include <cassert> // USES assert()
#include <iostream> // debugging.

// =====================================================================================================================
// Kernels for the Ricker Source Time Function in 2D.
// =====================================================================================================================

// ----------------------------------------------------------------------
// g1 entry function for velocity equation
void
pylith::fekernels::RickerWaveletPlaneStrain::g1v(const PylithInt dim,
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
    const PylithInt i_rickerwaveletCenterFrequency = numA - 1;
    
    const PylithScalar* momentTensor = &a[aOff[i_momentTensor]];
    const PylithScalar timeDelay = a[aOff[i_timeDelay]];
    const PylithScalar rickerwaveletCenterFrequency = a[aOff[i_rickerwaveletCenterFrequency]];

    // RickerWavelet source time function (time domain)

    PylithScalar rt = t - timeDelay;
    PylithScalar rickerwavelet = (1.0 - 2.0*PETSC_PI*PETSC_PI*rickerwaveletCenterFrequency*rickerwaveletCenterFrequency*rt*rt) * 
                           PetscExpReal(-PETSC_PI*PETSC_PI*rickerwaveletCenterFrequency*rickerwaveletCenterFrequency*rt*rt);

    for (PylithInt i = 0; i < dim*dim; ++i) {
        g1[i] -= momentTensor[i] * rickerwavelet;
    } // for
} // g1v


// =====================================================================================================================
// Kernels for the Ricker Source Time Function in 3D.
// =====================================================================================================================

// ----------------------------------------------------------------------
// g1 entry function for velocity equation
void
pylith::fekernels::RickerWavelet3D::g1v(const PylithInt dim,
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
    const PylithInt i_rickerwaveletCenterFrequency = numA - 1;
    
    const PylithScalar* momentTensor = &a[aOff[i_momentTensor]];
    const PylithScalar timeDelay = a[aOff[i_timeDelay]];
    const PylithScalar rickerwaveletCenterFrequency = a[aOff[i_rickerwaveletCenterFrequency]];

    // RickerWavelet source time function (time domain)

    PylithScalar rt = t - timeDelay;
    PylithScalar rickerwavelet = (1.0 - 2.0*PETSC_PI*PETSC_PI*rickerwaveletCenterFrequency*rickerwaveletCenterFrequency*rt*rt) * 
                           PetscExpReal(-PETSC_PI*PETSC_PI*rickerwaveletCenterFrequency*rickerwaveletCenterFrequency*rt*rt);

    for (PylithInt i = 0; i < dim*dim; ++i) {
        g1[i] -= momentTensor[i] * rickerwavelet;
    } // for
} // g1v


// End of file
