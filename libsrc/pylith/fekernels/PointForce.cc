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
// f0 entry function for displacement equation: f0u = \dot{u}.
void
pylith::fekernels::PointForce::f0u(const PylithInt dim,
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
                                   PylithScalar f0[]) {
    assert(sOff);
    assert(s);
    assert(s_t);
    assert(f0);

    // Incoming re-packed solution field.
    const PylithInt i_displacement = 0;

    // Incoming re-packed auxiliary field.
    const PylithInt i_pointForce = 0;

    const PylithScalar *pointForce = &a[aOff[i_pointForce]];

    for (PylithInt i = 0; i < dim; ++i) {
        f0[i] += pointForce[i];
    } // for
}


// End of file
