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

#include "pylith/fekernels/PointForce.hh"

#include <cassert> // USES assert()
#include <iostream> // debugging.

// ----------------------------------------------------------------------
// g1 entry function for velocity equation
void
pylith::fekernels::PointForce::g1v(const PylithInt dim,
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

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;

    // Incoming re-packed auxiliary field.
    const PylithInt i_momentTensor = 0;
    const PylithInt i_rickerTimeDelay = 1;
    const PylithInt i_rickerCenterFrequency = 2;
    
    const PylithScalar* momentTensor = &a[aOff[i_momentTensor]];
    const PylithScalar rickerTimeDelay = a[aOff[i_rickerTimeDelay]];
    const PylithScalar rickerCenterFrequency = a[aOff[i_rickerCenterFrequency]];

    // Ricker source  time function (time domain)

    PylithScalar rt = t - rickerTimeDelay;
    PylithScalar ricker = (1.0 - 2.0*PETSC_PI*PETSC_PI*rickerCenterFrequency*rickerCenterFrequency*rt*rt) * 
                           PetscExpReal(-PETSC_PI*PETSC_PI*rickerCenterFrequency*rickerCenterFrequency*rt*rt);

    for (PylithInt i = 0; i < dim*dim; ++i) {
        g1[i] -= momentTensor[i] * ricker;
    } // for
} // g1v


// End of file
