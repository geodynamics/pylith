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
 * Copyright (c) 2010-2022 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

/** @file libsrc/fekernels/Elasticity.hh
 *
 * Solution fields: [disp(dim), vel(dim, optional)]
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * Kernel interface.
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field.
 * @param[in] numA Number of registered subfields in auxiliary field.
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] numConstants Number of registered constants.
 * @param[in] constants Array of registered constants.
 * @param[out] f0 [dim].
 */

#if !defined(pylith_fekernels_elasticity_hh)
#define pylith_fekernels_elasticity_hh

#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::Elasticity {
public:

    typedef void (*strainfn_type) (const PylithInt,
                                   const PylithInt,
                                   const PylithInt[],
                                   const PylithInt[],
                                   const PylithScalar[],
                                   const PylithScalar[],
                                   const PylithScalar[],
                                   const PylithScalar[],
                                   PylithScalar[]);
    typedef void (*stressfn_type) (const PylithInt,
                                   const PylithInt,
                                   const PylithInt,
                                   const PylithInt[],
                                   const PylithInt[],
                                   const PylithScalar[],
                                   const PylithScalar[],
                                   const PylithScalar[],
                                   const PylithInt[],
                                   const PylithInt[],
                                   const PylithScalar[],
                                   const PylithScalar[],
                                   const PylithScalar[],
                                   const PylithReal t,
                                   const PylithScalar[],
                                   const PylithInt,
                                   const PylithScalar[],
                                   const PylithScalar[],
                                   PylithScalar[]);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    /** f0 function for elasticity equation.
     *
     * Solution fields: [disp(dim), vel(dim)]
     * Auxiliary fields: [density(1), ...]
     */
    static inline
    void f0v(const PylithInt dim,
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

    // --------------------------------------------------------------------------------------------
    /** Jf0 function for elasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), ...]
     */
    static inline
    void Jf0vv(const PylithInt dim,
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

    // --------------------------------------------------------------------------------------------
    /** g0 function for elasticity equation with gravitational body force.
     *
     * \vec{g0} = \vec{f}(t)
     *
     * Solution fields: []
     * Auxiliary fields: [density, gravity_field(dim)]
     */
    static inline
    void g0v_grav(const PylithInt dim,
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

    // --------------------------------------------------------------------------------------------
    /** g0 function for elasticity equation with body force.
     *
     * \vec{g0} = \vec{f}(t)
     *
     * Auxiliary fields: [body_force(dim)]
     */
    static inline
    void g0v_bodyforce(const PylithInt dim,
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

    // --------------------------------------------------------------------------------------------
    /** g0 function for elasticity with both gravitational and body forces.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), body_force(dim), gravity_field(dim), ...]
     */
    static inline
    void g0v_gravbodyforce(const PylithInt dim,
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

}; // Elasticity

// ------------------------------------------------------------------------------------------------
/// Kernels specific to elasticity plane strain.
class pylith::fekernels::ElasticityPlaneStrain {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    /** Calculate infinitesimal strain tensor for 2D plane strain elasticity.
     *
     * Order of output components is xx, xy, yx, yy.
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void infinitesimalStrain(const PylithInt dim,
                             const PylithInt numS,
                             const PylithInt sOff[],
                             const PylithInt sOff_x[],
                             const PylithScalar s[],
                             const PylithScalar s_t[],
                             const PylithScalar s_x[],
                             const PylithScalar x[],
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
        const PylithScalar strain_xy = 0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]);

        strain[0] = strain_xx;
        strain[1] = strain_xy;
        strain[2] = strain_xy;
        strain[3] = strain_yy;
    } // infinitesimalStrain

    // --------------------------------------------------------------------------------------------
    /** Calculate vector with infinitesimal strain for 2D plane strain elasticity.
     *
     * Order of output components is xx, yy, zz, xy.
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void infinitesimalStrain_asVector(const PylithInt dim,
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
                                      PylithScalar strainVector[]) {
        assert(strainVector);

        PylithScalar strainTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
        infinitesimalStrain(dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);

        strainVector[0] = strainTensor[0]; // xx
        strainVector[1] = strainTensor[3]; // yy
        strainVector[2] = 0.0; // zz
        strainVector[3] = strainTensor[1]; // xy
    } // infinitesimalStrain_asVector

    // --------------------------------------------------------------------------------------------
    /** Calculate vector with stress for 2D plane strain elasticity.
     *
     * Order of output components is xx, yy, zz, xy.
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void stress_asVector(const PylithInt dim,
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
                         Elasticity::strainfn_type strainFn,
                         Elasticity::stressfn_type stressFn,
                         PylithScalar stressVector[]) {
        const PylithInt _dim = 2;
        assert(stressVector);

        PylithScalar strainTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
        strainFn(dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);

        PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
        stressFn(_dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
                 t, x, numConstants, constants, strainTensor, stressTensor);

        stressVector[0] = stressTensor[0]; // xx
        stressVector[1] = stressTensor[3]; // yy
        stressVector[2] = 0.0; // zz
        stressVector[3] = stressTensor[1]; // xy
    } // cauchyStress_asVector

    // --------------------------------------------------------------------------------------------
    // f1 function for plane strain elasticity.
    static inline
    void f1v(const PylithInt dim,
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
             Elasticity::strainfn_type strainFn,
             Elasticity::stressfn_type stressFn,
             PylithScalar f1[]) {
        const PylithInt _dim = 2;

        PylithScalar strainTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
        strainFn(_dim, numS, sOff, sOff, s, s_t, s_x, x, strainTensor);

        PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
        stressFn(_dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
                 t, x, numConstants, constants, strainTensor, stressTensor);

        for (PylithInt i = 0; i < _dim*_dim; ++i) {
            f1[i] -= stressTensor[i];
        } // for
    } // f1v

}; // ElasticityPlaneStrain

// ------------------------------------------------------------------------------------------------
/// Kernels specific to elasticity in 3D.
class pylith::fekernels::Elasticity3D {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    /** Calculate infinitesimal strain for 3D elasticity.
     *
     * Order of output components is xx, xy, xz, yx, yy, yz, zx, zy, zz.
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void infinitesimalStrain(const PylithInt dim,
                             const PylithInt numS,
                             const PylithInt sOff[],
                             const PylithInt sOff_x[],
                             const PylithScalar s[],
                             const PylithScalar s_t[],
                             const PylithScalar s_x[],
                             const PylithScalar x[],
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
        strain[1] = strain_xy;
        strain[2] = strain_xz;
        strain[3] = strain_xy;
        strain[4] = strain_yy;
        strain[5] = strain_yz;
        strain[6] = strain_xz;
        strain[7] = strain_yz;
        strain[8] = strain_zz;
    } // infinitesimalStrain

    // --------------------------------------------------------------------------------------------
    /** Calculate vector with infinitesimal strain for 3D elasticity.
     *
     * Order of output components is xx, yy, zz, xy, yz, xz.
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void infinitesimalStrain_asVector(const PylithInt dim,
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
                                      PylithScalar strainVector[]) {
        assert(strainVector);

        PylithScalar strainTensor[9] = { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
        infinitesimalStrain(dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);

        strainVector[0] = strainTensor[0]; // xx
        strainVector[1] = strainTensor[4]; // yy
        strainVector[2] = strainTensor[8]; // zz
        strainVector[3] = strainTensor[1]; // xy
        strainVector[4] = strainTensor[5]; // yz
        strainVector[5] = strainTensor[2]; // xz
    } // infinitesimalStrain_asVector

    // --------------------------------------------------------------------------------------------
    /** Calculate vector with stress for 3D elasticity.
     *
     * Order of output components is xx, yy, zz, xy, yz, xz.
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void stress_asVector(const PylithInt dim,
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
                         Elasticity::strainfn_type strainFn,
                         Elasticity::stressfn_type stressFn,
                         PylithScalar stressVector[]) {
        const PylithInt _dim = 3;
        assert(stressVector);

        PylithScalar strainTensor[9] = { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
        strainFn(dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);

        PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
        stressFn(_dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
                 t, x, numConstants, constants, strainTensor, stressTensor);

        stressVector[0] = stressTensor[0]; // xx
        stressVector[1] = stressTensor[4]; // yy
        stressVector[2] = stressTensor[8]; // zz
        stressVector[3] = stressTensor[1]; // xy
        stressVector[4] = stressTensor[5]; // yz
        stressVector[5] = stressTensor[2]; // xz
    } // stress_asVector

    // --------------------------------------------------------------------------------------------
    // f1 function for 3D elasticity.
    static inline
    void f1v(const PylithInt dim,
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
             Elasticity::strainfn_type strainFn,
             Elasticity::stressfn_type stressFn,
             PylithScalar f1[]) {
        const PylithInt _dim = 3;

        PylithScalar strainTensor[9] = { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
        strainFn(dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);

        PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
        stressFn(_dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
                 t, x, numConstants, constants, strainTensor, stressTensor);

        for (PylithInt i = 0; i < _dim*_dim; ++i) {
            f1[i] -= stressTensor[i];
        } // for
    } // f1v

};

// Elasticity3D

#endif // pylith_fekernels_elasticity3d_hh

// End of file
