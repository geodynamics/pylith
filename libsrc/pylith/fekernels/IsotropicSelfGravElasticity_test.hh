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

/** @file libsrc/fekernels/IsotropicLinearIncompElasticity.hh
 *
 * Kernels for linear incompressible elasticity WITHOUT inertia.
 *
 * Solution fields: [disp(dim), potential(1)]
 *
 * Auxiliary fields:
 * - 0: density(1)
 * - 1: body_force(dim,optional)
 * - 2: gravity_field (dim, optional)
 * - 3: reference_stress(optional)
 *     2D: 4 components (stress_xx, stress_yy, stress_zz, stress_xy)
 *     3D: 6 components (stress_xx, stress_yy, stress_zz, stress_xy, stress_yz, stress_xz)
 * - 4: reference_strain(optional)
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 * - 5: shear_modulus(1)
 * - 6: bulk_modulus(1)
 *
 * The elasticity subfields come first (with required ones before optional ones) followed by the rheology subfields
 * (optional ones before required ones). The rheology fields have required fields last because we index from the back.
 *
 * :TODO: @brad Add equation here
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

#if !defined(pylith_fekernels_isotropiclinearincompelasticity_hh)
#define pylith_fekernels_isotropiclinearincompelasticity_hh

#include "fekernelsfwd.hh" // forward declarations

#include "pylith/fekernels/SelfGravElasticity.hh"        // USES SelfGravElasticity kernels
#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity kernels
#include "pylith/fekernels/Elasticity.hh"                // USES Elasticity kernels

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
class pylith::fekernels::IsotropicLinearIncompElasticity
{
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:
    struct Context
    {
        PylithReal potential;
        PylithReal shearModulus;
        PylithReal bulkModulus;
        pylith::fekernels::Tensor refStress;
        pylith::fekernels::Tensor refStrain;

        Context(void) : potential(0.0),
                        shearModulus(0.0),
                        bulkModulus(0.0) {}
    };

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:
    // --------------------------------------------------------------------------------------------
    static inline void setContext(Context *context,
                                  const PylithInt dim,
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
                                  const pylith::fekernels::TensorOps &tensorOps)
    {
        assert(context);

        const PylithInt i_potential = numS - 1;
        const PylithInt i_shearModulus = numA - 2;
        const PylithInt i_bulkModulus = numA - 1;

        assert(numS >= 2);
        assert(s);
        assert(sOff[i_potential] >= 0);
        assert(numA >= 3); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);

        context->potential = s[sOff[i_potential]];
        context->shearModulus = a[aOff[i_shearModulus]];
        assert(context->shearModulus > 0.0);
        context->bulkModulus = a[aOff[i_bulkModulus]];
        assert(context->bulkModulus > 0.0);
    } // setContext

    // --------------------------------------------------------------------------------------------
    static inline void setContext_refState(Context *context,
                                           const PylithInt dim,
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
                                           const pylith::fekernels::TensorOps &tensorOps)
    {
        assert(context);

        const PylithInt i_potential = numS - 1;
        const PylithInt i_refStress = numA - 4;
        const PylithInt i_refStrain = numA - 3;
        const PylithInt i_shearModulus = numA - 2;
        const PylithInt i_bulkModulus = numA - 1;

        assert(numS >= 2);
        assert(s);
        assert(sOff[i_potential] >= 0);
        assert(numA >= 5); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);

        context->potential = s[sOff[i_potential]];
        context->shearModulus = a[aOff[i_shearModulus]];
        assert(context->shearModulus > 0.0);
        context->bulkModulus = a[aOff[i_bulkModulus]];
        assert(context->bulkModulus > 0.0);

        tensorOps.fromVector(&a[aOff[i_refStress]], &context->refStress);
        tensorOps.fromVector(&a[aOff[i_refStrain]], &context->refStrain);
    } // createContext

    // ===========================================================================================
    // Kernels for elasticity equation
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Jf0_pp entry function for isotropic linear incompressible elasticity .
     *
     * Solution fields: [disp(dim), potential(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline void Jf3pp(const PylithInt dim,
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
                             PylithScalar Jf3[])
    {
        // Incoming auxiliary subfields
        const PylithInt i_density = 0;

        assert(numA >= 1);
        assert(aOff);
        assert(aOff[i_density] >= 0);
        assert(a);
        assert(Jf3);

        const PylithScalar density = a[aOff[i_density]];

        Jf3[0] += 1 / (4 * 3.14159 * 6.67 * pow(10, -11) * density);
        Jf3[4] += 1 / (4 * 3.14159 * 6.67 * pow(10, -11) * density);
        Jf3[8] += 1 / (4 * 3.14159 * 6.67 * pow(10, -11) * density);
        Jf3[36] += 1 / (4 * 3.14159 * 6.67 * pow(10, -11) * density);
        Jf3[40] += 1 / (4 * 3.14159 * 6.67 * pow(10, -11) * density);
        Jf3[44] += 1 / (4 * 3.14159 * 6.67 * pow(10, -11) * density);
        Jf3[72] += 1 / (4 * 3.14159 * 6.67 * pow(10, -11) * density);
        Jf3[76] += 1 / (4 * 3.14159 * 6.67 * pow(10, -11) * density);
        Jf3[80] += 1 / (4 * 3.14159 * 6.67 * pow(10, -11) * density); /// Fix

    } // Jf3pp

    // ===========================================================================================
    // Helper functions
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f0p helper function for isotropic linear incompressible elasticity WITH reference stress
     * and reference strain.
     *
     * ISA pylith::fekernels::SelfGravElasticity::incompressiblefn_type
     */
    static inline void incompressibleTerm(void *rheologyContext,
                                          const pylith::fekernels::Tensor &strain,
                                          const pylith::fekernels::TensorOps &tensorOps,
                                          PylithScalar *value)
    {
        assert(value);
        Context *context = (Context *)(rheologyContext);
        assert(context);

        const PylithScalar potential = context->potential;
        const PylithScalar bulkModulus = context->bulkModulus;

        const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        *value = strainTrace + potential / bulkModulus;
    } // incompressibleTerm

    // --------------------------------------------------------------------------------------------
    /** f0p helper function for isotropic linear incompressible elasticity WITH reference stress
     * and reference strain.
     *
     * ISA pylith::fekernels::SelfGravElasticity::incompressiblefn_type
     */
    static inline void incompressibleTerm_refState(void *rheologyContext,
                                                   const pylith::fekernels::Tensor &strain,
                                                   const pylith::fekernels::TensorOps &tensorOps,
                                                   PylithScalar *value)
    {
        assert(value);
        Context *context = (Context *)(rheologyContext);
        assert(context);

        const PylithScalar potential = context->potential;
        const PylithScalar bulkModulus = context->bulkModulus;
        const pylith::fekernels::Tensor &refStress = context->refStress;
        const pylith::fekernels::Tensor &refStrain = context->refStrain;

        const PylithReal meanRefStress = (refStress.xx + refStress.yy + refStress.zz) / 3.0;
        const PylithReal refStrainTrace = refStrain.xx + refStrain.yy + refStrain.zz;
        const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        *value = (strainTrace - refStrainTrace) + (potential + meanRefStress) / bulkModulus;
    } // incompressibleTerm_refState

    // --------------------------------------------------------------------------------------------
    /** Calculate stress for 2D plane strain isotropic linear elasticity WITHOUT a reference
     * stress and strain.
     *
     * ISA Elasticity::stressFn
     *
     * @param[in] rheologyContext IsotropicLinearElasticity context.
     * @param[in] strain Strain tensor.
     * @param[in] tensorOps Tensor operations.
     * @param[out] stress Stress tensor.
     */
    static inline void cauchyStress(void *rheologyContext,
                                    const pylith::fekernels::Tensor &strain,
                                    const pylith::fekernels::TensorOps &tensorOps,
                                    Tensor *stress)
    {
        Context *context = (Context *)(rheologyContext);
        assert(context);
        assert(stress);

        const PylithReal potential = context->potential;
        // pylith::fekernels::SelfGravElasticity::meanStress(potential, stress);

        const PylithReal shearModulus = context->shearModulus;
        pylith::fekernels::IsotropicLinearElasticity::deviatoricStress(shearModulus, strain, stress);
    } // cauchyStress

    // --------------------------------------------------------------------------------------------
    /** Calculate stress for 2-D plane strain isotropic linear incompressible elasticity WITH a
     * reference stress/strain.
     *
     * ISA Elasticity::stressFn
     *
     * @param[in] rheologyContext IsotropicLinearElasticity context.
     * @param[in] strain Strain tensor.
     * @param[in] tensorOps Tensor operations.
     * @param[out] stress Stress tensor.
     */
    static inline void cauchyStress_refState(void *rheologyContext,
                                             const pylith::fekernels::Tensor &strain,
                                             const pylith::fekernels::TensorOps &tensorOps,
                                             Tensor *stress)
    {
        Context *context = (Context *)(rheologyContext);
        assert(context);
        assert(stress);

        const pylith::fekernels::Tensor &refStress = context->refStress;
        const pylith::fekernels::Tensor &refStrain = context->refStrain;

        const PylithReal potential = context->potential;
        // pylith::fekernels::SelfGravElasticity::meanStress_refState(potential, refStress, stress);

        const PylithReal shearModulus = context->shearModulus;
        pylith::fekernels::IsotropicLinearElasticity::deviatoricStress_refState(shearModulus, refStress, refStrain, strain, stress);
    } // cauchyStress_refState

}; // IsotropicLinearIncompElasticity

// ------------------------------------------------------------------------------------------------
class pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain
{
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:
    // ===========================================================================================
    // Kernels for elasticity equation
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f0p entry function for isotropic linear incompressible elasticity with infinitesimal strain
     * WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), potential(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline void f0p(const PylithInt dim,
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
                           PylithScalar f0[])
    {
        f0[0] = -1;
    } //

    // --------------------------------------------------------------------------------------------
    /** f0p entry function for isotropic linear incompressible elasticity with infinitesimal strain
     * WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), potential(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for 2D plane strain isotropic linear incompressible elasticity with
     * infinitesimal strain WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), potential(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline void f1u_infinitesimalStrain(const PylithInt dim,
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
                                               PylithScalar f1[])
    {
        const PylithInt _dim = 2;
        assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearIncompElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearIncompElasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearIncompElasticity::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            f1);
    } // f1u_infinitesimalStrain

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for 2D plane strain isotropic linear incompressible elasticity with
     * infinitesimal strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), potential(1)]
     * Auxiliary fields: [..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
     */
    static inline void f1u_infinitesimalStrain_refState(const PylithInt dim,
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
                                                        PylithScalar f1[])
    {
        const PylithInt _dim = 2;
        assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearIncompElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearIncompElasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearIncompElasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops2D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** Jf3_uu entry function for 2-D plane strain isotropic linear incompressible elasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline void Jf3uu_infinitesimalStrain(const PylithInt dim,
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
                                                 PylithScalar Jf3[])
    {
        const PylithInt _dim = 2;
        assert(_dim == dim);

        pylith::fekernels::IsotropicLinearIncompElasticity::Context context;
        pylith::fekernels::IsotropicLinearIncompElasticity::setContext(
            &context, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        const PylithScalar shearModulus = context.shearModulus;

        const PylithReal C1111 = 4.0 / 3.0 * shearModulus;
        const PylithReal C2222 = C1111;
        const PylithReal C1122 = -2.0 / 3.0 * shearModulus;
        const PylithReal C1212 = shearModulus;

        /* j(f,g,df,dg) = C(f,df,g,dg)
         *
         * 0: j0000 = C1111
         * 1: j0001 = C1112 = 0
         * 4: j0100 = C1121, symmetry C1112 = 0
         * 5: j0101 = C1122
         *
         * 2: j0010 = C1211 = 0
         * 3: j0011 = C1212
         * 6: j0110 = C1221, symmetry C1212
         * 7: j0111 = C1222 = 0
         *
         * 8: j1000 = C2111 = 0
         * 9: j1001 = C2112, symmetry C1212
         * 12: j1100 = C2121, symmetry C1212
         * 13: j1101 = C2122, symmetry C1222 = 0
         *
         * 10: j1010 = C2211, symmetry C1122
         * 11: j1011 = C2212, symmetry C1222 = 0
         * 14: j1110 = C2221, symmetry C1222 = 0
         * 15: j1111 = C2222
         */

        Jf3[0] -= C1111;  // j0000
        Jf3[3] -= C1212;  // j0011
        Jf3[5] -= C1122;  // j0101
        Jf3[6] -= C1212;  // j0110, C1221
        Jf3[9] -= C1212;  // j1001, C2112
        Jf3[10] -= C1122; // j1010, C2211
        Jf3[12] -= C1212; // j1100, C2121
        Jf3[15] -= C2222; // j1111
    }                     // Jf3uu

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 2D plane strain isotropic linear
     * elasticity with infinitesimal strain WITHOUT a reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline void cauchyStress_infinitesimalStrain_asVector(const PylithInt dim,
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
                                                                 PylithScalar stressVector[])
    {
        const PylithInt _dim = 2;
        assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearIncompElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearIncompElasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearIncompElasticity::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            stressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 2D plane strain isotropic linear
     * elasticity with infinitesimal strain WITH a reference stress and strain.
     *
     * Used in output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline void cauchyStress_infinitesimalStrain_refState_asVector(const PylithInt dim,
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
                                                                          PylithScalar stressVector[])
    {
        const PylithInt _dim = 2;
        assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearIncompElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearIncompElasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearIncompElasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops2D,
            stressVector);
    }

}; // IsotropicLinearIncompElasticityPlaneStrain

// ------------------------------------------------------------------------------------------------
class pylith::fekernels::IsotropicLinearIncompElasticity3D
{
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:
    // ===========================================================================================
    // Kernels for elasticity equation
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f0p entry function for 3D isotropic linear incompressible elasticity with infinitesimal strain
     * WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), potential(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline void f0p(const PylithInt dim,
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
                           PylithScalar f0[])
    {
        f0[0] = -1;
    } //

    // --------------------------------------------------------------------------------------------
    /** f0p entry function for 3D isotropic linear incompressible elasticity with infinitesimal strain
     * WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), potential(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for 3D isotropic linear incompressible elasticity with infinitesimal strain
     * WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), potential(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline void f1u_infinitesimalStrain(const PylithInt dim,
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
                                               PylithScalar f1[])
    {
        const PylithInt _dim = 3;
        assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearIncompElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearIncompElasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearIncompElasticity::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            f1);
    } // f1u_infinitesimalStrain

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for 3D isotropic linear incompressible elasticity with infinitesimal
     * strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), potential(1)]
     * Auxiliary fields: [..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
     */
    static inline void f1u_infinitesimalStrain_refState(const PylithInt dim,
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
                                                        PylithScalar f1[])
    {
        const PylithInt _dim = 3;
        assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearIncompElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearIncompElasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearIncompElasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops3D,
            f1);
    } // f1u_infinitesimalStrain_refState

    void f1p_potential(const PylithInt dim,
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
                       PylithScalar f1[])
    {

        const PylithInt i_density = 0;

        assert(numA >= 1);
        assert(aOff);
        assert(aOff[i_density] >= 0);

        const PylithScalar density = a[aOff[i_density]];

        for (PetscInt d = 0; d < dim; ++d)
            f1[d] = s_x[sOff_x[1] + d] * 1 / (4 * 3.14159 * 6.67 * pow(10, -11) * density);

    } // f1p_potential

    // --------------------------------------------------------------------------------------------
    /** Jf3_uu entry function for 3D isotropic linear incompressible elasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline void Jf3uu_infinitesimalStrain(const PylithInt dim,
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
                                                 PylithScalar Jf3[])
    {
        const PylithInt _dim = 3;
        assert(_dim == dim);

        pylith::fekernels::IsotropicLinearIncompElasticity::Context context;
        pylith::fekernels::IsotropicLinearIncompElasticity::setContext(
            &context, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        const PylithScalar shearModulus = context.shearModulus;

        // All other values are either zero or equal to one of these.
        const PylithReal C1111 = 4.0 / 3.0 * shearModulus;
        const PylithReal C2222 = C1111;
        const PylithReal C3333 = C1111;
        const PylithReal C1122 = -2.0 / 3.0 * shearModulus;
        const PylithReal C1133 = C1122;
        const PylithReal C2233 = C1122;
        const PylithReal C1212 = shearModulus;
        const PylithReal C1313 = C1212;
        const PylithReal C2323 = C1212;

        /* j(f,g,df,dg) = C(f,df,g,dg)
         *
         *  0:  j0000 = C1111
         *  1:  j0001 = C1112 = 0
         *  2:  j0002 = C1113 = 0
         *  3:  j0010 = C1211 = 0
         *  4:  j0011 = C1212
         *  5:  j0012 = C1213 = 0
         *  6:  j0020 = C1311 = 0
         *  7:  j0021 = C1312 = 0
         *  8:  j0022 = C1313
         *  9:  j0100 = C1121 = 0
         * 10:  j0101 = C1122
         * 11:  j0102 = C1123 = 0
         * 12:  j0110 = C1221 = C1212
         * 13:  j0111 = C1222 = 0
         * 14:  j0112 = C1223 = 0
         * 15:  j0120 = C1321 = 0
         * 16:  j0121 = C1322 = 0
         * 17:  j0122 = C1323 = 0
         * 18:  j0200 = C1131 = 0
         * 19:  j0201 = C1132 = 0
         * 20:  j0202 = C1133
         * 21:  j0210 = C1231 = 0
         * 22:  j0211 = C1232 = 0
         * 23:  j0212 = C1233 = 0
         * 24:  j0220 = C1331 = C1313
         * 25:  j0221 = C1332 = 0
         * 26:  j0222 = C1333 = 0
         * 27:  j1000 = C2111 = 0
         * 28:  j1001 = C2112 = C1212
         * 29:  j1002 = C2113 = 0
         * 30:  j1010 = C2211 = C1122
         * 31:  j1011 = C2212 = 0
         * 32:  j1012 = C2213 = 0
         * 33:  j1020 = C2311 = 0
         * 34:  j1021 = C2312 = 0
         * 35:  j1022 = C2313 = 0
         * 36:  j1100 = C2121 = C1212
         * 37:  j1101 = C2122 = 0
         * 38:  j1102 = C2123 = 0
         * 39:  j1110 = C2221 = 0
         * 40:  j1111 = C2222
         * 41:  j1112 = C2223 = 0
         * 42:  j1120 = C2321 = 0
         * 43:  j1121 = C2322 = 0
         * 44:  j1122 = C2323
         * 45:  j1200 = C2131 = 0
         * 46:  j1201 = C2132 = 0
         * 47:  j1202 = C2133 = 0
         * 48:  j1210 = C2231 = 0
         * 49:  j1211 = C2232 = 0
         * 50:  j1212 = C2233
         * 51:  j1220 = C2331 = 0
         * 52:  j1221 = C2332 = C2323
         * 53:  j1222 = C2333 = 0
         * 54:  j2000 = C3111 = 0
         * 55:  j2001 = C3112 = 0
         * 56:  j2002 = C3113 = C1313
         * 57:  j2010 = C3211 = 0
         * 58:  j2011 = C3212 = 0
         * 59:  j2012 = C3213 = 0
         * 60:  j2020 = C3311 = C1133
         * 61:  j2021 = C3312 = 0
         * 62:  j2022 = C3313 = 0
         * 63:  j2100 = C3121 = 0
         * 64:  j2101 = C3122 = 0
         * 65:  j2102 = C3123 = 0
         * 66:  j2110 = C3221 = 0
         * 67:  j2111 = C3222 = 0
         * 68:  j2112 = C3223 = C2323
         * 69:  j2120 = C3321 = 0
         * 70:  j2121 = C3322 = C2233
         * 71:  j2122 = C3323 = 0
         * 72:  j2200 = C3131 = C1313
         * 73:  j2201 = C3132 = 0
         * 74:  j2202 = C3133 = 0
         * 75:  j2210 = C3231 = 0
         * 76:  j2211 = C3232 = C2323
         * 77:  j2212 = C3233 = 0
         * 78:  j2220 = C3331 = 0
         * 79:  j2221 = C3332 = 0
         * 80:  j2222 = C3333
         */

        // Nonzero Jacobian entries.
        Jf3[0] -= C1111;  // j0000
        Jf3[4] -= C1212;  // j0011
        Jf3[8] -= C1313;  // j0022
        Jf3[10] -= C1122; // j0101
        Jf3[12] -= C1212; // j0110
        Jf3[20] -= C1133; // j0202
        Jf3[24] -= C1313; // j0220
        Jf3[28] -= C1212; // j1001
        Jf3[30] -= C1122; // j1010
        Jf3[36] -= C1212; // j1100
        Jf3[40] -= C2222; // j1111
        Jf3[44] -= C2323; // j1122
        Jf3[50] -= C2233; // j1212
        Jf3[52] -= C2323; // j1221
        Jf3[56] -= C1313; // j2002
        Jf3[60] -= C1133; // j2020
        Jf3[68] -= C2323; // j2112
        Jf3[70] -= C2233; // j2121
        Jf3[72] -= C1313; // j2200
        Jf3[76] -= C2323; // j2211
        Jf3[80] -= C3333; // j2222
    }                     // Jf3uu

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic linear elasticity with
     * infinisteismal strain WITHOUT a reference stress and strain.
     *
     * Used in output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline void cauchyStress_infinitesimalStrain_asVector(const PylithInt dim,
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
                                                                 PylithScalar stressVector[])
    {
        const PylithInt _dim = 3;
        assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearIncompElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearIncompElasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearIncompElasticity::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            stressVector);
    } // cauchyStress_infinitesimalStrain_asVector

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D linear elasticity WITH a reference
     * stress and strain.
     *
     * Used in output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline void cauchyStress_infinitesimalStrain_refState_asVector(const PylithInt dim,
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
                                                                          PylithScalar stressVector[])
    {
        const PylithInt _dim = 3;
        assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearIncompElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearIncompElasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearIncompElasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops3D,
            stressVector);
    } // cauchyStress_infinitesimalStrain_refState_asVector

}; // IsotropicLinearIncompElasticity3D

#endif // pylith_fekernels_isotropiclinearincompelasticity_hh

// End of file
