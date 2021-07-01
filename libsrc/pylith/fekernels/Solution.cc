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

#include "pylith/fekernels/Solution.hh"

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
// Identify function kernel.
void
pylith::fekernels::Solution::passThruSubfield(const PylithInt dim,
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
                                              PylithScalar field[]) {
    assert(s);
    assert(sOff);
    assert(field);
    const PetscInt subfieldIndex = PetscInt(t); // :KLUDGE: Easiest way to get subfield to extract into fn.

    const PylithInt sEnd = sOff[subfieldIndex+1];
    for (PylithInt iS = sOff[subfieldIndex], iF = 0; iS < sEnd; ++iS, ++iF) {
        field[iF] = s[iS];
    } // for
} // passThruSubfield


// End of file
