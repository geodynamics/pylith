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

#include "pylith/fekernels/DispVel.hh"

#include <cassert> // USES assert()
#include <iostream> // debugging.

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
// f0 entry function for displacement equation: f0u = \dot{u}.
void
pylith::fekernels::DispVel::f0u(const PylithInt dim,
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

    const PylithInt _numS = 2;
    assert(_numS == numS);

    const PylithInt i_disp = 0;
    const PylithScalar* disp_t = &s_t[sOff[i_disp]];

    for (PylithInt i = 0; i < dim; ++i) {
        f0[i] += disp_t[i];
    } // for
} // f0u


// ----------------------------------------------------------------------
// f0 entry function for displacement equation: f0u = \dot{v}.
void
pylith::fekernels::DispVel::f0v(const PylithInt dim,
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

    const PylithInt _numS = 2;
    assert(_numS == numS);

    const PylithInt i_vel = 1;
    const PylithScalar* vel_t = &s_t[sOff[i_vel]];

    for (PylithInt i = 0; i < dim; ++i) {
        f0[i] += vel_t[i];
    } // for
} // f0v


// ----------------------------------------------------------------------
// g0 function for displacement equation: g0u = v.
void
pylith::fekernels::DispVel::g0u(const PylithInt dim,
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

    const PylithInt _numS = 2;
    assert(_numS == numS);

    const PylithInt i_vel = 1;
    const PylithScalar* vel = &s[sOff[i_vel]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += vel[i];
    } // for
} // g0u


// ----------------------------------------------------------------------
// Jf0 function for displacement equation with zero values on diagonal.
void
pylith::fekernels::DispVel::Jf0uu_zero(const PylithInt dim,
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
                                       PylithScalar Jf0[]) {
    // No work to do for zero values.
} //


// ----------------------------------------------------------------------
// Jf0 function for displacement equation: Jf0uu = s_tshift.
void
pylith::fekernels::DispVel::Jf0uu_stshift(const PylithInt dim,
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
                                          PylithScalar Jf0[]) {
    const PylithInt _numS = 2;
    assert(_numS == numS);
    assert(s_tshift > 0);

    for (PylithInt i = 0; i < dim; ++i) {
        Jf0[i*dim+i] += s_tshift;
    } // for
} // Jf0uu


// ----------------------------------------------------------------------
/* Jg0 function for disp/velocity equation.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::DispVel::Jg0uv(const PylithInt dim,
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
                                  PylithScalar Jg0[]) {
    const PylithInt _numS = 2;
    assert(_numS == numS);

    for (PylithInt i = 0; i < dim; ++i) {
        Jg0[i*dim+i] += 1.0;
    } // for
} // Jg0uv


// End of file
