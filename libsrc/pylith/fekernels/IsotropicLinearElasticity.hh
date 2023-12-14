/*
 * ================================================================================================
 * This code is part of PyLith, developed through the Computational Infrastructure
 * for Geodynamics (https://github.com/geodynamics/pylith).
 *
 * Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
 * All rights reserved.
 *
 * See https://mit-license.org/ and LICENSE.md and for license information.
 * =================================================================================================
 */
#pragma once

/*
 * Kernels for isotropic, linear elasticity.
 *
 * Solution fields: [disp(dim), vel(dim, optional)]
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
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * =====================================================================================================================
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

#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/FaultCohesiveKin.hh" // USES FaultCohesiveKin kernels

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linear elasticity.
class pylith::fekernels::IsotropicLinearElasticity {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    struct Context {
        PylithReal shearModulus;
        PylithReal bulkModulus;
        pylith::fekernels::Tensor refStress;
        pylith::fekernels::Tensor refStrain;

        Context(void) :
            shearModulus(0.0),
            bulkModulus(0.0) {}


    };

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    static inline
    void setContext(Context* context,
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
                    const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);

        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(numA >= 3); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);

        context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
        context->bulkModulus = a[aOff[i_bulkModulus]];assert(context->bulkModulus > 0.0);
    } // setContext

    // --------------------------------------------------------------------------------------------
    static inline
    void setContext_refState(Context* context,
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
                             const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);

        const PylithInt i_refStress = numA-4;
        const PylithInt i_refStrain = numA-3;
        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(numA >= 5); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);

        context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
        context->bulkModulus = a[aOff[i_bulkModulus]];assert(context->bulkModulus > 0.0);

        tensorOps.fromVector(&a[aOff[i_refStress]], &context->refStress);
        tensorOps.fromVector(&a[aOff[i_refStrain]], &context->refStrain);
    } // createContext

    // --------------------------------------------------------------------------------------------
    /** Helper function for calculating Cauchy stress for WITHOUT a reference stress and strain.
     *
     * ISA Elasticity::stressFn
     *
     * @param[in] rheologyContext IsotropicLinearElasticity context.
     * @param[in] strain Strain tensor.
     * @param[in] tensorOps Tensor operations.
     * @param[out] stress Stress tensor.
     */
    static inline
    void cauchyStress(void* rheologyContext,
                      const pylith::fekernels::Tensor& strain,
                      const pylith::fekernels::TensorOps& tensorOps,
                      pylith::fekernels::Tensor* stress) {
        Context* context = (Context*)(rheologyContext);
        assert(context);
        assert(stress);

        meanStress(context->bulkModulus, strain, stress);
        deviatoricStress(context->shearModulus, strain, stress);
    } // cauchyStress

    // --------------------------------------------------------------------------------------------
    /** Helper function for calculating Cauchy stress WITH reference stress/strain.
     *
     * @param[in] rheologyContext IsotropicLinearElastcity context.
     * @param[in] strain Strain tensor.
     * @param[in] tensorOps Tensor operations.
     * @param[out] stress Stress tensor.
     */
    static inline
    void cauchyStress_refState(void* rheologyContext,
                               const pylith::fekernels::Tensor& strain,
                               const pylith::fekernels::TensorOps& tensorOps,
                               pylith::fekernels::Tensor* stress) {
        Context* context = (Context*)(rheologyContext);
        assert(context);
        assert(stress);

        const pylith::fekernels::Tensor& refStress = context->refStress;
        const pylith::fekernels::Tensor& refStrain = context->refStrain;
        meanStress_refState(context->bulkModulus, refStress, refStrain, strain, stress);
        deviatoricStress_refState(context->shearModulus, refStress, refStrain, strain, stress);
    } // cauchyStress_refState

    // --------------------------------------------------------------------------------------------
    /** Calculate mean stress WITHOUT reference stress and reference strain.
     */
    static inline
    void meanStress(const PylithReal bulkModulus,
                    const pylith::fekernels::Tensor& strain,
                    pylith::fekernels::Tensor* stress) {
        assert(bulkModulus > 0.0);
        assert(stress);

        const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        const PylithReal meanStress = bulkModulus * strainTrace;

        stress->xx += meanStress;
        stress->yy += meanStress;
        stress->zz += meanStress;
    } // meanStress

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITHOUT reference stress and strain.
     */
    static inline
    void deviatoricStress(const PylithReal shearModulus,
                          const pylith::fekernels::Tensor& strain,
                          pylith::fekernels::Tensor* stress) {
        assert(shearModulus > 0.0);
        assert(stress);

        const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        const PylithReal traceTerm = -2.0/3.0*shearModulus * strainTrace;

        stress->xx += 2.0*shearModulus*strain.xx + traceTerm;
        stress->yy += 2.0*shearModulus*strain.yy + traceTerm;
        stress->zz += 2.0*shearModulus*strain.zz + traceTerm;
        stress->xy += 2.0*shearModulus*strain.xy;
        stress->yz += 2.0*shearModulus*strain.yz;
        stress->xz += 2.0*shearModulus*strain.xz;
    } // deviatoricStress

    // --------------------------------------------------------------------------------------------
    /** Calculate mean stress WITH reference stress and reference strain.
     */
    static inline
    void meanStress_refState(const PylithReal bulkModulus,
                             const pylith::fekernels::Tensor& refStress,
                             const pylith::fekernels::Tensor& refStrain,
                             const pylith::fekernels::Tensor& strain,
                             pylith::fekernels::Tensor* stress) {
        // Incoming auxiliary fields.
        assert(bulkModulus > 0.0);
        assert(stress);

        const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        const PylithReal refStrainTrace = refStrain.xx + refStrain.yy + refStrain.zz;

        const PylithReal meanRefStress = (refStress.xx + refStress.yy + refStress.zz) / 3.0;
        const PylithReal meanStress = meanRefStress + bulkModulus * (strainTrace - refStrainTrace);

        stress->xx += meanStress;
        stress->yy += meanStress;
        stress->zz += meanStress;
    } // meanStress_refState

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITH reference stress and strain.
     */
    static inline
    void deviatoricStress_refState(const PylithReal shearModulus,
                                   const pylith::fekernels::Tensor& refStress,
                                   const pylith::fekernels::Tensor& refStrain,
                                   const pylith::fekernels::Tensor& strain,
                                   pylith::fekernels::Tensor* stress) {
        assert(shearModulus > 0.0);
        assert(stress);

        const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        const PylithReal refStrainTrace = refStrain.xx + refStrain.yy + refStrain.zz;
        const PylithReal meanRefStress = (refStress.xx + refStress.yy + refStress.zz) / 3.0;
        const PylithReal traceTerm = -2.0/3.0*shearModulus * (strainTrace - refStrainTrace);

        stress->xx += refStress.xx - meanRefStress + 2.0*shearModulus*(strain.xx-refStrain.xx) + traceTerm;
        stress->yy += refStress.yy - meanRefStress + 2.0*shearModulus*(strain.yy-refStrain.yy) + traceTerm;
        stress->zz += refStress.zz - meanRefStress + 2.0*shearModulus*(strain.zz-refStrain.zz) + traceTerm;
        stress->xy += refStress.xy + 2.0*shearModulus*(strain.xy - refStrain.xy);
        stress->yz += refStress.yz + 2.0*shearModulus*(strain.yz - refStrain.yz);
        stress->xz += refStress.xz + 2.0*shearModulus*(strain.xz - refStrain.xz);
    } // deviatoricStress_refState

}; // IsotropicLinearElasticity

// ------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linearly elasticity plane strain.
class pylith::fekernels::IsotropicLinearElasticityPlaneStrain {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // ===========================================================================================
    // Kernels for elasticity equation
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for isotropic linear elasticity plane strain with infinitesimal strain
     * WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f1v_infinitesimalStrain(const PylithInt dim,
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
                                 PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for isotropic linear elasticity plane strain with infinitesimal strain WITH
     * reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f1v_infinitesimalStrain_refState(const PylithInt dim,
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
                                          PylithScalar f1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops2D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** Jf3_vu entry function for 2D plane strain isotropic linear elasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     *
     * stress_ij = C_ijkl strain_kl
     *
     * stress_11 = C1111 strain_11 + C1122 strain_22, C1111=lambda+2mu, C1122=lambda.
     *
     * stress_12 = C1212 strain_12 + C1221 strain_21. C1212 = C1221 from symmetry, so C1212 = C1221 = shearModulus.
     *
     * For reference:
     *
     * Isotropic:
     *  C_ijkl = bulkModulus * delta_ij * delta_kl
     *   + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
     */
    static inline
    void Jf3vu_infinitesimalStrain(const PylithInt dim,
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
                                   PylithScalar Jf3[]) {
        const PylithInt _dim = 2;assert(_dim == dim);
        pylith::fekernels::IsotropicLinearElasticity::Context context;
        pylith::fekernels::IsotropicLinearElasticity::setContext(
            &context, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        assert(Jf3);

        const PylithScalar shearModulus = context.shearModulus;
        const PylithScalar bulkModulus = context.bulkModulus;

        const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
        const PylithScalar lambda2mu = lambda + 2.0*shearModulus;

        const PylithReal C1111 = lambda2mu;
        const PylithReal C2222 = lambda2mu;
        const PylithReal C1122 = lambda;
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

        Jf3[ 0] -= C1111; // j0000
        Jf3[ 3] -= C1212; // j0011
        Jf3[ 5] -= C1122; // j0101
        Jf3[ 6] -= C1212; // j0110, C1221
        Jf3[ 9] -= C1212; // j1001, C2112
        Jf3[10] -= C1122; // j1010, C2211
        Jf3[12] -= C1212; // j1100, C2121
        Jf3[15] -= C2222; // j1111
    } // Jf3vu

    // ===========================================================================================
    // Kernels for fault interfaces and elasticity
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for isotropic linear elasticity plane strain with
     * infinitesimal strain WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), vel(dim), lagrange(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f0l_neg_infinitesimalStrain(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for isotropic linear elasticity plane strain with
     * infinitesimal strain WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), vel(dim), lagrange(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f0l_pos_infinitesimalStrain(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for isotropic linear elasticity plane strain with
     * infinitesimal strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), vel(dim), lagrange(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f0l_neg_infinitesimalStrain_refState(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress_refState,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for isotropic linear elasticity plane strain with
     * infinitesimal strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), vel(dim), lagrange(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f0l_pos_infinitesimalStrain_refState(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress_refState,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 2D plane strain isotropic linear
     * elasticity WITHOUT a reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void cauchyStress_infinitesimalStrain_asVector(const PylithInt dim,
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
                                                   PylithScalar stressVector[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            stressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 2D plane strain isotropic linear
     * elasticity with infinitesimal strain WITH a reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void cauchyStress_infinitesimalStrain_refState_asVector(const PylithInt dim,
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
                                                            PylithScalar stressVector[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops2D,
            stressVector);
    } // cauchyStress_infinitesimalStrain_refState_asVector

}; // IsotropicLinearElasticityPlaneStrain

// ------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linearly elasticity in 3D.
class pylith::fekernels::IsotropicLinearElasticity3D {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // ===========================================================================================
    // Kernels for elasticity equation
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for 3D isotropic linear elasticity with infinitesimal strain WITHOUT
     * reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f1v_infinitesimalStrain(const PylithInt dim,
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
                                 PylithScalar f1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for 3D isotropic linear elasticity with infinitesimal strain WITH
     * reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., refstress(6), refstrain(6), shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f1v_infinitesimalStrain_refState(const PylithInt dim,
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
                                          PylithScalar f1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops3D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** Jf3_vu entry function for 3D isotropic linear elasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void Jf3vu_infinitesimalStrain(const PylithInt dim,
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
                                   PylithScalar Jf3[]) {
        const PylithInt _dim = 3;assert(_dim == dim);
        pylith::fekernels::IsotropicLinearElasticity::Context context;
        pylith::fekernels::IsotropicLinearElasticity::setContext(
            &context, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        assert(Jf3);

        const PylithScalar shearModulus = context.shearModulus;
        const PylithScalar bulkModulus = context.bulkModulus;

        const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
        const PylithScalar lambda2mu = lambda + 2.0*shearModulus;

        // All other values are either zero or equal to one of these.
        const PylithReal C1111 = lambda2mu;
        const PylithReal C1122 = lambda;
        const PylithReal C1212 = shearModulus;
        /* j(f,g,df,dg) = C(f,df,g,dg)
         *
         *  0:  j0000 = C1111 = lambda + 2.0*shearModulus
         *  1:  j0001 = C1112 = 0
         *  2:  j0002 = C1113 = 0
         *  3:  j0010 = C1211 = 0
         *  4:  j0011 = C1212 = shearModulus
         *  5:  j0012 = C1213 = 0
         *  6:  j0020 = C1311 = 0
         *  7:  j0021 = C1312 = 0
         *  8:  j0022 = C1313 = shearModulus
         *  9:  j0100 = C1121 = 0
         * 10:  j0101 = C1122 = lambda
         * 11:  j0102 = C1123 = 0
         * 12:  j0110 = C1221 = shearModulus
         * 13:  j0111 = C1222 = 0
         * 14:  j0112 = C1223 = 0
         * 15:  j0120 = C1321 = 0
         * 16:  j0121 = C1322 = 0
         * 17:  j0122 = C1323 = 0
         * 18:  j0200 = C1131 = 0
         * 19:  j0201 = C1132 = 0
         * 20:  j0202 = C1133 = lambda
         * 21:  j0210 = C1231 = 0
         * 22:  j0211 = C1232 = 0
         * 23:  j0212 = C1233 = 0
         * 24:  j0220 = C1331 = shearModulus
         * 25:  j0221 = C1332 = 0
         * 26:  j0222 = C1333 = 0
         * 27:  j1000 = C2111 = 0
         * 28:  j1001 = C2112 = shearModulus
         * 29:  j1002 = C2113 = 0
         * 30:  j1010 = C2211 = lambda
         * 31:  j1011 = C2212 = 0
         * 32:  j1012 = C2213 = 0
         * 33:  j1020 = C2311 = 0
         * 34:  j1021 = C2312 = 0
         * 35:  j1022 = C2313 = 0
         * 36:  j1100 = C2121 = shearModulus
         * 37:  j1101 = C2122 = 0
         * 38:  j1102 = C2123 = 0
         * 39:  j1110 = C2221 = 0
         * 40:  j1111 = C2222 = lambda + 2.0*shearModulus
         * 41:  j1112 = C2223 = 0
         * 42:  j1120 = C2321 = 0
         * 43:  j1121 = C2322 = 0
         * 44:  j1122 = C2323 = shearModulus
         * 45:  j1200 = C2131 = 0
         * 46:  j1201 = C2132 = 0
         * 47:  j1202 = C2133 = 0
         * 48:  j1210 = C2231 = 0
         * 49:  j1211 = C2232 = 0
         * 50:  j1212 = C2233 = lambda
         * 51:  j1220 = C2331 = 0
         * 52:  j1221 = C2332 = shearModulus
         * 53:  j1222 = C2333 = 0
         * 54:  j2000 = C3111 = 0
         * 55:  j2001 = C3112 = 0
         * 56:  j2002 = C3113 = shearModulus
         * 57:  j2010 = C3211 = 0
         * 58:  j2011 = C3212 = 0
         * 59:  j2012 = C3213 = 0
         * 60:  j2020 = C3311 = lambda
         * 61:  j2021 = C3312 = 0
         * 62:  j2022 = C3313 = 0
         * 63:  j2100 = C3121 = 0
         * 64:  j2101 = C3122 = 0
         * 65:  j2102 = C3123 = 0
         * 66:  j2110 = C3221 = 0
         * 67:  j2111 = C3222 = 0
         * 68:  j2112 = C3223 = shearModulus
         * 69:  j2120 = C3321 = 0
         * 70:  j2121 = C3322 = lambda
         * 71:  j2122 = C3323 = 0
         * 72:  j2200 = C3131 = shearModulus
         * 73:  j2201 = C3132 = 0
         * 74:  j2202 = C3133 = 0
         * 75:  j2210 = C3231 = 0
         * 76:  j2211 = C3232 = shearModulus
         * 77:  j2212 = C3233 = 0
         * 78:  j2220 = C3331 = 0
         * 79:  j2221 = C3332 = 0
         * 80:  j2222 = C3333 = lambda + 2.0*shearModulus
         */

        // Nonzero Jacobian entries.
        Jf3[ 0] -= C1111; // j0000
        Jf3[ 4] -= C1212; // j0011
        Jf3[ 8] -= C1212; // j0022
        Jf3[10] -= C1122; // j0101
        Jf3[12] -= C1212; // j0110
        Jf3[20] -= C1122; // j0202
        Jf3[24] -= C1212; // j0220
        Jf3[28] -= C1212; // j1001
        Jf3[30] -= C1122; // j1010
        Jf3[36] -= C1212; // j1100
        Jf3[40] -= C1111; // j1111
        Jf3[44] -= C1212; // j1122
        Jf3[50] -= C1122; // j1212
        Jf3[52] -= C1212; // j1221
        Jf3[56] -= C1212; // j2002
        Jf3[60] -= C1122; // j2020
        Jf3[68] -= C1212; // j2112
        Jf3[70] -= C1122; // j2121
        Jf3[72] -= C1212; // j2200
        Jf3[76] -= C1212; // j2211
        Jf3[80] -= C1111; // j2222
    } // Jf3vu_infinitesimalStrain

    // ===========================================================================================
    // Kernels for fault interfaces and elasticity
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for 3D isotropic linear elasticity with
     * infinitesimal strain WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), vel(dim), lagrange(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f0l_neg_infinitesimalStrain(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for 3D isotropic linear elasticity with
     * infinitesimal strain WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), vel(dim), lagrange(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f0l_pos_infinitesimalStrain(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for 3D isotropic linear elasticity with
     * infinitesimal strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), vel(dim), lagrange(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f0l_neg_infinitesimalStrain_refState(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress_refState,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for 3D isotropic linear elasticity with
     * infinitesimal strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), vel(dim), lagrange(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f0l_pos_infinitesimalStrain_refState(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress_refState,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic linear elasticity with
     * infinitesimal strain WITHOUT reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void cauchyStress_infinitesimalStrain_asVector(const PylithInt dim,
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
                                                   PylithScalar stressVector[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            stressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic linear elasticity with
     * infinitesimal strain WITH a reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void cauchyStress_infinitesimalStrain_refState_asVector(const PylithInt dim,
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
                                                            PylithScalar stressVector[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearElasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearElasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearElasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops3D, stressVector);
    }

}; // IsotropicLinearElasticity3D

// End of file
