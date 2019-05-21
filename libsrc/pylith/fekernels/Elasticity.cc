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

#include "pylith/fekernels/Elasticity.hh"

#include <cassert> // USES assert()

// =====================================================================================================================
// Generic elasticity kernels for inertia and body forces.
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
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

    // Incoming solution fields.
    const PylithInt i_density = 0;

    // Incoming auxiliary fields.
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


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
// Jf0 function for generic elasticity terms (inertia) with implicit time stepping.
void
pylith::fekernels::Elasticity::Jf0vv_inertiaimplicit(const PylithInt dim,
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
    const PylithInt _numA = 1;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(a);
    assert(s_tshift > 0);

    const PylithScalar density = a[aOff[i_density]];

    for (PylithInt i = 0; i < dim; ++i) {
        for (PylithInt j = 0; j < dim; ++j) {
            Jf0[i*dim+j] += s_tshift * density;
        } // for
    } // for
} // Jf0vv_inertiaimplicit


// ---------------------------------------------------------------------------------------------------------------------
// Jf0 function for generic elasticity terms (inertia) with explicit time stepping.
void
pylith::fekernels::Elasticity::Jf0vv_inertiaexplicit(const PylithInt dim,
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
    const PylithInt _numA = 1;

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_vel = 1;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(aOff);
    assert(a);

    const PylithScalar density = a[aOff[i_density]];

    Jf0[i_disp*_numS+i_vel] += density;
} // Jf0vv_inertiaexplicit


// =====================================================================================================================
// Kernels for elasticity plane strain.
// =====================================================================================================================

#include <iostream>
// ---------------------------------------------------------------------------------------------------------------------
/* Calculate Cauchy strain for 2-D plane strain elasticity.
 *
 * Order of output components is xx, yy, zz, xy.
 *
 * Solution fields: [disp(dim)]
 */
void
pylith::fekernels::ElasticityPlaneStrain::cauchy_strain(const PylithInt dim,
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

    // DEBUGGING
    const PylithReal p_density = 2500.0;
    const PylithReal p_vs = 3000.0;
    const PylithReal p_vp = 5291.502622129181;
    const PylithReal p_mu = p_density*p_vs*p_vs;
    const PylithReal p_lambda = p_density*p_vp*p_vp - 2*p_mu;

    const PylithReal sxx = 1.0e+7;
    const PylithReal sxy = 0.0;
    const PylithReal syy = 0.0;
    const PylithReal szz = p_lambda/(2*p_lambda+2*p_mu)*(sxx+syy);

    const PylithReal exx = 1.0/(2*p_mu) * (sxx - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz));
    const PylithReal eyy = 1.0/(2*p_mu) * (syy - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz));
    const PylithReal ezz = 1.0/(2*p_mu) * (szz - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz));
    const PylithReal exy = 1.0/(2*p_mu) * (sxy);

    const PylithReal ux = exx*x[0] + exy*x[1];
    const PylithReal uy = exy*x[0] + eyy*x[1];

    std::cout << "x:"<<x[0]<<", y:"<<x[1]<<", disp_expected: "<<ux<<", "<<uy<<", disp_kernel_arg:"<<s[0]<<", "<<s[1]<<std::endl;

    strain[0] = strain_xx;
    strain[1] = strain_yy;
    strain[2] = strain_zz;
    strain[3] = strain_xy;
} // cauchy_strain


// =====================================================================================================================
// Kernels for elasticity in 3D
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
/** Calculate Cauchy strain for 3-D elasticity.
 *
 * Order of output components is xx, yy, zz, xy, yz, xz.
 *
 * Solution fields: [disp(dim)]
 */
void
pylith::fekernels::Elasticity3D::cauchy_strain(const PylithInt dim,
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
} // cauchy_strain


// End of file
