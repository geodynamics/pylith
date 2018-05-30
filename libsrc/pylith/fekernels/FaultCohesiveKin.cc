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

#include <stdexcept> // USES std::logic_error
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
// g0 function for elasticity equation: g0u = +-\lambda.
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
                                         const PylithInt numConstants,
                                         const PylithScalar constants[],
                                         PylithScalar g0[]) {
    assert(sOff);
    assert(s);
    assert(g0);

    assert(numS >= 2);

    const PylithInt i_lagrange = numS-1;
    const PylithScalar* lagrange = &s[sOff[i_lagrange]];

    throw std::logic_error(":TODO: @matt @brad How to assemble for each side of the fault? -lambda for u+, +lambda for u-");

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += lagrange[i];
    } // for
} // f0u


// ----------------------------------------------------------------------
// g0 function for slip constraint equation: g0\lambda = d - u^+ + u^-.
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

    const PylithInt i_disp = 0;
    const PylithScalar* dispPos = &s[sOff[i_disp]];
    const PylithScalar* dispNeg = dispPos;

    throw std::logic_error(":TODO: @matt @brad How to get displacements for each side of the fault?");

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += slip[i] - dispPos[i] + dispNeg[i];
    } // for
} // f0v


// ----------------------------------------------------------------------
/* Jg0 function for displacement equation.
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
                                           const PylithReal utshift,
                                           const PylithScalar x[],
                                           const PylithInt numConstants,
                                           const PylithScalar constants[],
                                           PylithScalar Jg0[]) {
    assert(numS >= 2);

    throw std::logic_error(":TODO: @matt @brad How to assemble for each side of the fault? -1 for u+, +1 for u-");

    for (PylithInt i = 0; i < dim; ++i) {
        Jg0[i*dim+i] += 1.0;
    } // for
} // Jg0ul


// ----------------------------------------------------------------------
/* Jg0 function for slip constraint equation.
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
                                           const PylithReal utshift,
                                           const PylithScalar x[],
                                           const PylithInt numConstants,
                                           const PylithScalar constants[],
                                           PylithScalar Jg0[]) {
    assert(numS >= 2);

    throw std::logic_error(":TODO: @matt @brad How to assemble for each side of the fault? -1 for u+, +1 for u-");

    for (PylithInt i = 0; i < dim; ++i) {
        Jg0[i*dim+i] += 1.0;
    } // for
} // Jg0ul


// End of file
