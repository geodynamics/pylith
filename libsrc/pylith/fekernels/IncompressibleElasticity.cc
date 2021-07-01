/* -*- C++ -*-
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

#include "pylith/fekernels/IncompressibleElasticity.hh"

#include <cassert> // USES assert()

// =====================================================================================================================
// Generic incompressible elasticity kernels.
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
// Jf1pu function for pressure equation for incompressible elasticity.
void
pylith::fekernels::IncompressibleElasticity::Jf1pu(const PylithInt dim,
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
                                                   const PylithReal s_tshift,
                                                   const PylithScalar x[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar Jf1[]) {
    assert(Jf1);

    /* j(f,g,dg), f=0, g=0..dim, dg=0..dim
     *
     * j == 1 if g==dg, otherwise 0.
     *
     * 3-D
     * 0: j000 = 1
     * 1: j001 = 0
     * 2: j002 = 0
     * 3: j010 = 0
     * 4: j011 = 1
     * 5: j012 = 0
     * 6: j020 = 0
     * 7: j021 = 0
     * 8: j022 = 1
     */

    for (PylithInt i = 0; i < dim; ++i) {
        Jf1[i*dim+i] += 1.0;
    } // for
} // Jf1pu


// ---------------------------------------------------------------------------------------------------------------------
// Jf2up function for elasticity equation.
void
pylith::fekernels::IncompressibleElasticity::Jf2up(const PylithInt dim,
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
                                                   const PylithReal s_tshift,
                                                   const PylithScalar x[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar Jf2[]) {
    assert(Jf2);

    /* j(f,g,df), f=0..dim, df=0..dim, g=0
     *
     * j == 1 if f==df, otherwise 0.
     *
     * 3-D
     * 0: j000 = 1
     * 1: j001 = 0
     * 2: j002 = 0
     * 3: j100 = 0
     * 4: j101 = 1
     * 5: j102 = 0
     * 6: j200 = 0
     * 7: j201 = 0
     * 8: j202 = 1
     */

    for (PylithInt i = 0; i < dim; ++i) {
        Jf2[i*dim+i] += 1.0;
    } // for
} // Jf2up


// End of file
