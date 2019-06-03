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

#include "pylith/fekernels/FaultCohesiveKin.hh"

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for time evolution equation with displacement and velocity
 * solution fields.
 *
 * Solution fields = [disp(dim), vel(dim), ...]
 * Auxiliary fields = None
 *
 * \int_V \vec{\phi}_v \cdot \left( \frac{\partial \vec{u}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_v \cdot \vec{v}(t) \, dV.
 *
 * ======================================================================
 */

// ----------------------------------------------------------------------
// g0 function for integration of the elasticity equation: g0u = -\lambda.
void
pylith::fekernels::FaultCohesiveKin::g0u(const PylithInt dim,
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
                                         const PylithReal n[],
                                         const PylithInt numConstants,
                                         const PylithScalar constants[],
                                         PylithScalar g0[]) {
    assert(sOff);
    assert(s);
    assert(g0);

    assert(numS >= 2);

    const PylithInt spaceDim = dim+1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt i_disp = 0;
    const PylithInt i_lagrange = numS-1;
    const PylithInt o_dispN = sOff[i_disp];
    const PylithInt o_dispP = sOff[i_disp]+spaceDim;

    // :KLUDGE: sOff doesn't account for DOF on two sides of the fault.
    //
    // Offset in s is off by number of subfields on two sides of the fault (numS-1)*dim.
    const PylithInt o_lagrange = sOff[i_lagrange]+(numS-1)*spaceDim;

    const PylithScalar* lagrange = &s[o_lagrange];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[o_dispN+i] += +lagrange[i];
        g0[o_dispP+i] += -lagrange[i];
    } // for
} // g0u


// ----------------------------------------------------------------------
// g0 function for integration of the slip constraint equation: g0\lambda = d - u^+ + u^-.
void
pylith::fekernels::FaultCohesiveKin::g0l(const PylithInt dim,
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
                                         const PylithReal n[],
                                         const PylithInt numConstants,
                                         const PylithScalar constants[],
                                         PylithScalar g0[]) {
    assert(sOff);
    assert(aOff);
    assert(s);
    assert(a);
    assert(g0);

    assert(numS >= 2);

    assert(numA >= 1);
    const PylithInt i_slip = numA-1;
    const PylithScalar* slip = &a[aOff[i_slip]];

    const PylithInt spaceDim = dim+1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt i_disp = 0;
    const PylithInt i_lagrange = numS-1;
    const PylithInt o_dispN = sOff[i_disp];
    const PylithInt o_dispP = sOff[i_disp]+spaceDim;

    // :KLUDGE: sOff doesn't account for DOF on two sides of the fault.
    //
    // Offset in s is off by number of subfields on two sides of the fault (numS-1)*dim.
    const PylithInt o_lagrange = sOff[i_lagrange]+(numS-1)*spaceDim;

    const PylithScalar* dispN = &s[o_dispN];
    const PylithScalar* dispP = &s[o_dispP];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[o_lagrange+i] += slip[i] - dispP[i] + dispN[i];
    } // for
} // g0l


// ----------------------------------------------------------------------
/* Jg0 function for integration of the displacement equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKin::Jg0ul(const PylithInt dim,
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
                                           const PylithReal n[],
                                           const PylithInt numConstants,
                                           const PylithScalar constants[],
                                           PylithScalar Jg0[]) {
    assert(numS >= 2);

    const PylithInt spaceDim = dim+1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt i_disp = 0;
    const PylithInt i_lagrange = numS-1;
    const PylithInt o_dispN = sOff[i_disp];
    const PylithInt o_dispP = sOff[i_disp]+spaceDim;

    // :KLUDGE: sOff doesn't account for DOF on two sides of the fault.
    //
    // Offset in s is off by number of subfields on two sides of the fault (numS-1)*dim.
    const PylithInt o_lagrange = sOff[i_lagrange]+(numS-1)*spaceDim;

    for (PylithInt i = 0; i < dim; ++i) {
        Jg0[(o_dispN+i)*dim+o_lagrange+i] += +1.0;
        Jg0[(o_dispP+i)*dim+o_lagrange+i] += -1.0;
    } // for
} // Jg0ul


// ----------------------------------------------------------------------
/* Jg0 function for integration of the slip constraint equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKin::Jg0lu(const PylithInt dim,
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
                                           const PylithReal n[],
                                           const PylithInt numConstants,
                                           const PylithScalar constants[],
                                           PylithScalar Jg0[]) {
    assert(numS >= 2);

    const PylithInt spaceDim = dim+1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt i_disp = 0;
    const PylithInt i_lagrange = numS-1;
    const PylithInt o_dispN = sOff[i_disp];
    const PylithInt o_dispP = sOff[i_disp]+spaceDim;

    // :KLUDGE: sOff doesn't account for DOF on two sides of the fault.
    //
    // Offset in s is off by number of subfields on two sides of the fault (numS-1)*dim.
    const PylithInt o_lagrange = sOff[i_lagrange]+(numS-1)*spaceDim;

    for (PylithInt i = 0; i < dim; ++i) {
        Jg0[(o_lagrange+i)*dim+o_dispN+i] += +1.0;
        Jg0[(o_lagrange+i)*dim+o_dispP+i] += -1.0;
    } // for
} // Jg0lu


// End of file
