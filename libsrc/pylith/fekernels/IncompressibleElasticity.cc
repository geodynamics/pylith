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

#include "pylith/fekernels/IncompressibleElasticity.hh"

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for pressure volume integral.
 *
 * Solution fields = [disp(dim), pres]
 * Auxiliary fields = bulkModulus
 *
 * 0 = \int_V \phi_p \cdot
 *  \left( \vec {\nabla} \cdot \vec{u} + \frac{p}{\kappa} \right) \, dV.
 *
 * ======================================================================
 */

// ----------------------------------------------------------------------
/* g0 function for pressure equation.
 */
void
pylith::fekernels::IncompressibleElasticity::g0p(const PylithInt dim,
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
                                                 PylithScalar g0[]) {
    const PylithInt _numS = 2;

    const PylithInt i_disp = 0;
    const PylithInt i_pres = 1;

    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar pres = s[sOff[i_pres]];

    const PylithInt i_bulkModulus = 2;

    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    PylithInt i;

    assert(_numS == numS);
    assert(3 <= numA);
    assert(sOff);
    assert(s);
    assert(aOff);
    assert(a);
    assert(g0);

    PylithScalar strainTrace = 0;

    for (i = 0; i < dim; ++i) {
        strainTrace += disp_x[i];
    } // for
    g0[0] += strainTrace + pres/bulkModulus;
} // g0p


// ----------------------------------------------------------------------
/* Jg0 function for pressure equation.
 */
void
pylith::fekernels::IncompressibleElasticity::Jg0pp(const PylithInt dim,
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
                                                   const PylithReal utshift,
                                                   const PylithScalar x[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar Jg0[]) {
    const PylithInt _numS = 2;

    const PylithInt i_bulkModulus = 1;
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    assert(_numS == numS);
    assert(3 <= numA);
    assert(aOff);
    assert(a);
    assert(Jg0);


    Jg0[0] += 1.0 / bulkModulus;
} // Jg0pp_implicit

// End of file
