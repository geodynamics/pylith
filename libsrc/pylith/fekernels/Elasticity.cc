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

#include "pylith/fekernels/Elasticity.hh"

#include <cassert> // USES assert()

namespace pylith {
    namespace fekernels {
        namespace _Elasticity {
            const PylithScalar tolerance = 1.0e-30;

            PylithInt lagrange_sOff(const PylithInt sOff[],
                                    const PylithInt numS);

        }
    }
}

// ================================================================================================
// Generic elasticity kernels for inertia and body forces.
// ================================================================================================

// ------------------------------------------------------------------------------------------------
// f0 function for elasticity equation.
void
pylith::fekernels::Elasticity::f0v(const PylithInt dim,
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
    const PylithInt _numS = 2;
    const PylithInt _numA = 1;

    // Incoming solution fields.
    const PylithInt i_vel = 1;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_numS == numS);
    assert(_numA <= numA);
    assert(sOff);
    assert(sOff[i_vel] >= 0);
    assert(s_t);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(a);

    const PylithScalar* vel_t = &s_t[sOff[i_vel]]; // acceleration
    const PylithScalar density = a[aOff[i_density]];

    for (PylithInt i = 0; i < dim; ++i) {
        f0[i] += vel_t[i] * density;
    } // for
} // f0v


// ------------------------------------------------------------------------------------------------
// Jf0 function for elasticity equation.
void
pylith::fekernels::Elasticity::Jf0vv(const PylithInt dim,
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
    const PylithInt _numA = 1;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_numA <= numA);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(a);

    const PylithScalar density = a[aOff[i_density]];

    for (PetscInt i = 0; i < dim; ++i) {
        Jf0[i*dim+i] += s_tshift * density;
    } // for
} // Jf0vv


// ------------------------------------------------------------------------------------------------
// g0 function for elasticity equation with gravitational body forces.
void
pylith::fekernels::Elasticity::g0v_grav(const PylithInt dim,
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
    const PylithInt _numA = 2;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_gravityField = 1;

    assert(_numA <= numA);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(a);

    const PylithScalar density = a[aOff[i_density]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += density*gravityField[i];
    } // for
} // g0v_grav


// ------------------------------------------------------------------------------------------------
// g0 function for elasticity equation with body forces.
void
pylith::fekernels::Elasticity::g0v_bodyforce(const PylithInt dim,
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
    const PylithInt _numA = 2;

    // Incoming auxiliary fields.
    const PylithInt i_bodyForce = 1;

    assert(_numA <= numA);
    assert(aOff);
    assert(aOff[i_bodyForce] >= 0);
    assert(a);
    assert(g0);

    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += bodyForce[i];
    } // for
} // g0v_bodyforce


// ------------------------------------------------------------------------------------------------
// g0 function for elasticity with both gravitational and body forces.
void
pylith::fekernels::Elasticity::g0v_gravbodyforce(const PylithInt dim,
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
    const PylithInt _numA = 3;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_bodyForce = 1;
    const PylithInt i_gravityField = 2;

    assert(_numA <= numA);
    assert(aOff);

    const PylithInt numSGrav = 0; // Number passed on to g0_grav.
    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    g0v_grav(dim, numSGrav, numAGrav, NULL, NULL, NULL, NULL, NULL, aOffGrav, NULL, a, a_t, NULL,
             t, x, numConstants, constants, g0);

    const PylithInt numSBody = 0; // Number passed on to g0_bodyforce.
    const PylithInt numABody = 2; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[2] = { aOff[i_density], aOff[i_bodyForce] };
    g0v_bodyforce(dim, numSBody, numABody, NULL, NULL, NULL, NULL, NULL, aOffBody, NULL, a, a_t, NULL,
                  t, x, numConstants, constants, g0);
} // g0v_gravbodyforce


// ------------------------------------------------------------------------------------------------
/** f0 function for negative fault face (+lambda).
 *
 * Solution fields: [disp(dim), vel(dim), lagrange_multiplier(dim)]
 * Auxiliary fields: [density(1), ...]
 */
void
pylith::fekernels::Elasticity::f0l_neg(const PylithInt dim,
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
                                       PylithScalar f0[]) {
    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(numS >= 3);
    assert(numA >= 1);

    assert(sOff);
    assert(s);
    assert(f0);

    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(a);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt fOffN = 0;
    const PylithInt sOffLagrange = pylith::fekernels::_Elasticity::lagrange_sOff(sOff, numS);
    assert(sOffLagrange >= 0);

    const PylithScalar density = a[aOff[i_density]];assert(density > _Elasticity::tolerance);
    const PylithScalar* lagrange = &s[sOffLagrange];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        f0[fOffN+i] += +lagrange[i] / density;
    } // for
} // f0l_neg


// ------------------------------------------------------------------------------------------------
/** f0 function for positive fault face (+lambda).
 *
 * Solution fields: [disp(dim), vel(dim), lagrange_multiplier(dim)]
 * Auxiliary fields: [density(1), ...]
 */
void
pylith::fekernels::Elasticity::f0l_pos(const PylithInt dim,
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
                                       PylithScalar f0[]) {
    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(numS >= 3);
    assert(numA > 1);

    assert(sOff);
    assert(s);
    assert(f0);

    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(a);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt fOffP = 0;
    const PylithInt sOffLagrange = pylith::fekernels::_Elasticity::lagrange_sOff(sOff, numS);
    assert(sOffLagrange >= 0);

    const PylithScalar density = a[aOff[i_density]];assert(density > _Elasticity::tolerance);
    const PylithScalar* lagrange = &s[sOffLagrange];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        f0[fOffP+i] += +lagrange[i] / density;
    } // for
} // f0l_pos


// ------------------------------------------------------------------------------------------------
/** f0 for negative fault face with gravitational body force.
 *
 * Auxiliary fields: [density, gravity_field(dim)]
 */
void
pylith::fekernels::Elasticity::f0l_neg_grav(const PylithInt dim,
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
                                            PylithScalar f0[]) {
    const PylithInt _numA = 2;

    // Incoming auxiliary fields.
    const PylithInt i_gravityField = 1;

    assert(numS >= 3);
    assert(_numA <= numA);

    assert(sOff);
    assert(s);
    assert(f0);

    assert(aOff);
    assert(aOff[i_gravityField] >= 0);
    assert(a);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt fOffN = 0;

    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        f0[fOffN+i] += -gravityField[i];
    } // for
}


// ------------------------------------------------------------------------------------------------
/** f0 for positive fault face with gravitational body force.
 *
 * Auxiliary fields: [density, gravity_field(dim)]
 */
void
pylith::fekernels::Elasticity::f0l_pos_grav(const PylithInt dim,
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
                                            PylithScalar f0[]) {
    const PylithInt _numA = 2;

    // Incoming auxiliary fields.
    const PylithInt i_gravityField = 1;

    assert(numS >= 3);
    assert(_numA <= numA);

    assert(sOff);
    assert(s);
    assert(f0);

    assert(aOff);
    assert(aOff[i_gravityField] >= 0);
    assert(a);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt fOffP = 0;

    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        f0[fOffP+i] += +gravityField[i];
    } // for
}


// ------------------------------------------------------------------------------------------------
/** f0 function for negative fault face with body force.
 *
 * Auxiliary fields: [density, body_force(dim)]
 */
void
pylith::fekernels::Elasticity::f0l_neg_bodyforce(const PylithInt dim,
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
                                                 PylithScalar f0[]) {
    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_bodyForce = 1;

    assert(numS >= 3);
    assert(numA >= 2);

    assert(sOff);
    assert(s);
    assert(f0);

    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_bodyForce] >= 0);
    assert(a);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt fOffN = 0;

    const PylithScalar density = a[aOff[i_density]];assert(density > _Elasticity::tolerance);
    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        f0[fOffN+i] += -bodyForce[i] / density;
    } // for
}


// ------------------------------------------------------------------------------------------------
/** f0 function for positive fault face with body force.
 *
 * Auxiliary fields: [density, body_force(dim)]
 */
void
pylith::fekernels::Elasticity::f0l_pos_bodyforce(const PylithInt dim,
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
                                                 PylithScalar f0[]) {
    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_bodyForce = 1;

    assert(numS >= 3);
    assert(numA >= 2);

    assert(sOff);
    assert(s);
    assert(f0);

    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_bodyForce] >= 0);
    assert(a);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt fOffP = 0;

    const PylithScalar density = a[aOff[i_density]];assert(density > _Elasticity::tolerance);
    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        f0[fOffP+i] += +bodyForce[i] / density;
    } // for
}


// ------------------------------------------------------------------------------------------------
/** f0 function for negative fault face with both gravitational and body forces.
 *
 * Auxiliary fields: [density(1), body_force(dim), gravity_field(dim), ...]
 */
void
pylith::fekernels::Elasticity::f0l_neg_gravbodyforce(const PylithInt dim,
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
                                                     PylithScalar f0[]) {
    const PylithInt _numA = 3;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_bodyForce = 1;
    const PylithInt i_gravityField = 2;

    assert(_numA <= numA);
    assert(aOff);

    const PylithInt numSGrav = 0; // Number passed on to f0l_neg_grav.
    const PylithInt numAGrav = 2; // Number passed on to f0l_neg_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    f0l_neg_grav(dim, numSGrav, numAGrav, NULL, NULL, NULL, NULL, NULL, aOffGrav, NULL, a, a_t, NULL,
                 t, x, n, numConstants, constants, f0);

    const PylithInt numSBody = 0; // Number passed on to f0l_neg_bodyforce.
    const PylithInt numABody = 2; // Number passed on to f0l_neg_bodyforce.
    const PylithInt aOffBody[2] = { aOff[i_density], aOff[i_bodyForce] };
    f0l_neg_bodyforce(dim, numSBody, numABody, NULL, NULL, NULL, NULL, NULL, aOffBody, NULL, a, a_t, NULL,
                      t, x, n, numConstants, constants, f0);
}


// ------------------------------------------------------------------------------------------------
/** f0 function for positive fault face with both gravitational and body forces.
 *
 * Auxiliary fields: [density(1), body_force(dim), gravity_field(dim), ...]
 */
void
pylith::fekernels::Elasticity::f0l_pos_gravbodyforce(const PylithInt dim,
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
                                                     PylithScalar f0[]) {
    const PylithInt _numA = 3;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_bodyForce = 1;
    const PylithInt i_gravityField = 2;

    assert(_numA <= numA);
    assert(aOff);

    const PylithInt numSGrav = 0; // Number passed on to f0l_pos_grav.
    const PylithInt numAGrav = 2; // Number passed on to f0l_pos_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    f0l_pos_grav(dim, numSGrav, numAGrav, NULL, NULL, NULL, NULL, NULL, aOffGrav, NULL, a, a_t, NULL,
                 t, x, n, numConstants, constants, f0);

    const PylithInt numSBody = 0; // Number passed on to f0l_pos_bodyforce.
    const PylithInt numABody = 2; // Number passed on to f0l_pos_bodyforce.
    const PylithInt aOffBody[2] = { aOff[i_density], aOff[i_bodyForce] };
    f0l_pos_bodyforce(dim, numSBody, numABody, NULL, NULL, NULL, NULL, NULL, aOffBody, NULL, a, a_t, NULL,
                      t, x, n, numConstants, constants, f0);
}


// ------------------------------------------------------------------------------------------------
// Jf0 function for dynamic slip constraint equation for negative side of the fault.
void
pylith::fekernels::Elasticity::Jf0ll_neg(const PylithInt dim,
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
                                         PylithScalar Jf0[]) {
    const PylithInt i_density = 0;

    assert(numS >= 1);
    assert(a);
    assert(aOff);
    assert(aOff[i_density] >= 0);

    assert(numS >= 2);
    assert(Jf0);
    assert(sOff);
    assert(n);

    const PylithInt spaceDim = dim+1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithScalar density = a[aOff[i_density]];assert(density > _Elasticity::tolerance);

    const PylithInt gOffN = 0;
    const PylithInt ncols = 2*spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jf0[i*ncols+gOffN+i] += +1.0 / density;
    } // for

}


// ------------------------------------------------------------------------------------------------
// Jf0 function for dynamic slip constraint equation for positive side of the fault.
void
pylith::fekernels::Elasticity::Jf0ll_pos(const PylithInt dim,
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
                                         PylithScalar Jf0[]) {
    const PylithInt i_density = 0;

    assert(numS >= 1);
    assert(a);
    assert(aOff);
    assert(aOff[i_density] >= 0);

    assert(numS >= 2);
    assert(Jf0);
    assert(sOff);
    assert(n);

    const PylithInt spaceDim = dim+1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithScalar density = a[aOff[i_density]];assert(density > _Elasticity::tolerance);

    const PylithInt gOffP = 0;
    const PylithInt ncols = 2*spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jf0[i*ncols+gOffP+i] += +1.0 / density;
    } // for

}


// ================================================================================================
// Kernels for elasticity plane strain.
// ================================================================================================

// ------------------------------------------------------------------------------------------------
/* Calculate Cauchy strain for 2-D plane strain elasticity.
 *
 * Order of output components is xx, yy, zz, xy.
 *
 * Solution fields: [disp(dim)]
 */
void
pylith::fekernels::ElasticityPlaneStrain::cauchyStrain(const PylithInt dim,
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
                                                       PylithScalar strain[]) {
    const PylithInt _dim = 2;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(sOff_x);
    assert(s_x);
    assert(strain);

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    const PylithScalar strain_xx = disp_x[0*_dim+0];
    const PylithScalar strain_yy = disp_x[1*_dim+1];
    const PylithScalar strain_zz = 0.0;
    const PylithScalar strain_xy = 0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]);

    strain[0] = strain_xx;
    strain[1] = strain_yy;
    strain[2] = strain_zz;
    strain[3] = strain_xy;
} // cauchyStrain


// ================================================================================================
// Kernels for elasticity in 3D
// ================================================================================================

// ------------------------------------------------------------------------------------------------
/** Calculate Cauchy strain for 3-D elasticity.
 *
 * Order of output components is xx, yy, zz, xy, yz, xz.
 *
 * Solution fields: [disp(dim)]
 */
void
pylith::fekernels::Elasticity3D::cauchyStrain(const PylithInt dim,
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
                                              PylithScalar strain[]) {
    const PylithInt _dim = 3;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(sOff_x);
    assert(s_x);
    assert(strain);

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    const PylithScalar strain_xx = disp_x[0*_dim+0];
    const PylithScalar strain_yy = disp_x[1*_dim+1];
    const PylithScalar strain_zz = disp_x[2*_dim+2];
    const PylithScalar strain_xy = 0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]);
    const PylithScalar strain_yz = 0.5*(disp_x[1*_dim+2] + disp_x[2*_dim+1]);
    const PylithScalar strain_xz = 0.5*(disp_x[0*_dim+2] + disp_x[2*_dim+0]);

    strain[0] = strain_xx;
    strain[1] = strain_yy;
    strain[2] = strain_zz;
    strain[3] = strain_xy;
    strain[4] = strain_yz;
    strain[5] = strain_xz;
} // cauchyStrain


// ----------------------------------------------------------------------
// Get offset in s where Lagrange multiplier field starts.
PylithInt
pylith::fekernels::_Elasticity::lagrange_sOff(const PylithInt sOff[],
                                              const PylithInt numS) {
    PylithInt off = 0;
    const PylithInt numCount = numS - 1; // Don't include last field (Lagrange multiplier)
    for (PylithInt i = 0; i < numCount; ++i) {
        off += 2*(sOff[i+1] - sOff[i]);
    } // for
    return off;
} // lagrange_sOff


// End of file
