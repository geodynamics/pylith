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
// g0 function for integration of the elasticity equation on positive side of the fault: g0u = -\lambda.
void
pylith::fekernels::FaultCohesiveKin::g0u_pos(const PylithInt dim,
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

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += -lagrange[i];
    } // for
} // g0u_pos


// ----------------------------------------------------------------------
// g0 function for integration of the elasticity equation on the negative side of the fault: g0u = +\lambda.
void
pylith::fekernels::FaultCohesiveKin::g0u_neg(const PylithInt dim,
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

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += +lagrange[i];
    } // for
} // g0u_neg


// ----------------------------------------------------------------------
// g0 function for integration of the slip constraint equation on the positive side of the fault:
// g0\lambda = 0.5*d - u^+.
void
pylith::fekernels::FaultCohesiveKin::g0l_pos(const PylithInt dim,
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
    const PylithScalar* disp = &s[sOff[i_disp]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += 0.5*slip[i] - disp[i];
    } // for
} // g0v_pos


// ----------------------------------------------------------------------
// g0 function for integration of the slip constraint equation on the positive side of the fault:
// g0\lambda = 0.5*d+- u^-.
void
pylith::fekernels::FaultCohesiveKin::g0l_neg(const PylithInt dim,
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
    const PylithScalar* disp = &s[sOff[i_disp]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += 0.5*slip[i] + disp[i];
    } // for
} // g0v_neg


// ----------------------------------------------------------------------
/* Jg0 function for integration of the displacement equation on the positive side of the fault.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKin::Jg0ul_pos(const PylithInt dim,
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

    for (PylithInt i = 0; i < dim; ++i) {
        Jg0[i*dim+i] += -1.0;
    } // for
} // Jg0ul_pos


// ----------------------------------------------------------------------
/* Jg0 function for integration of the displacement equation on the negative side of the fault.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKin::Jg0ul_neg(const PylithInt dim,
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

    for (PylithInt i = 0; i < dim; ++i) {
        Jg0[i*dim+i] += +1.0;
    } // for
} // Jg0ul_neg


// ----------------------------------------------------------------------
/* Jg0 function for integration of the slip constraint equation on the positive side of the fault.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKin::Jg0lu_pos(const PylithInt dim,
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

    for (PylithInt i = 0; i < dim; ++i) {
        Jg0[i*dim+i] += -1.0;
    } // for
} // Jg0lu_pos


// ----------------------------------------------------------------------
/* Jg0 function for integration of the slip constraint equation on the negative side of the fault.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKin::Jg0lu_neg(const PylithInt dim,
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

    for (PylithInt i = 0; i < dim; ++i) {
        Jg0[i*dim+i] += +1.0;
    } // for
} // Jg0lu_neg


// End of file
