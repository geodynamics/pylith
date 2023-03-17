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

/** @file libsrc/fekernels/IsotropicLinearMaxwell.hh
 *
 * Kernels for linear Maxwell viscoelastic plane strain.
 *
 * Solution fields: [disp(dim), ...]
 *
 * Isotropic, linear Maxwell viscoelastic plane strain without reference stress/strain.
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
 * - 7: maxwell_time(1)
 * - 8: viscous_strain
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 * - 9: total_strain
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 *
 * The elasticity subfields come first (with required ones before optional ones) followed by the rheology subfields
 * (optional ones before required ones). The rheology fields have required fields last because we index from the back.
 *
 * IMPORTANT: Viscous strain must be before total strain, because viscous strain at t+dt depends on total strain at t.
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * ================================================================================================
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

#if !defined(pylith_fekernels_isotropiclinearmaxwell_hh)
#define pylith_fekernels_isotropiclinearmaxwell_hh

#include "fekernelsfwd.hh" // forward declarations
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity* kernels

#include "pylith/utils/types.hh"

// ------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linear Maxwell viscoelasticity (dimension independent).
class pylith::fekernels::IsotropicLinearMaxwell {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    struct Context {
        PylithReal shearModulus;
        PylithReal bulkModulus;
        PylithReal maxwellTime;
        PylithReal dt;
        pylith::fekernels::Tensor viscousStrain;
        pylith::fekernels::Tensor totalStrain;
        pylith::fekernels::Tensor refStress;
        pylith::fekernels::Tensor refStrain;
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

        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_bulkModulus = numA-4;
        const PylithInt i_maxwellTime = numA-3;
        const PylithInt i_viscousStrain = numA-2;
        const PylithInt i_totalStrain = numA-1;

        assert(numA >= 6); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_maxwellTime] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_totalStrain] >= 0);
        assert(1 == numConstants);
        assert(constants);

        context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
        context->bulkModulus = a[aOff[i_bulkModulus]];assert(context->bulkModulus > 0.0);
        context->maxwellTime = a[aOff[i_maxwellTime]];assert(context->maxwellTime > 0.0);
        context->dt = constants[0];

        tensorOps.fromVector(&a[aOff[i_viscousStrain]], &context->viscousStrain);
        tensorOps.fromVector(&a[aOff[i_totalStrain]], &context->totalStrain);
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

        const PylithInt i_refStress = numA-7;
        const PylithInt i_refStrain = numA-6;
        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_bulkModulus = numA-4;
        const PylithInt i_maxwellTime = numA-3;
        const PylithInt i_viscousStrain = numA-2;
        const PylithInt i_totalStrain = numA-1;

        assert(numA >= 8); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_maxwellTime] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_totalStrain] >= 0);

        context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
        context->bulkModulus = a[aOff[i_bulkModulus]];assert(context->bulkModulus > 0.0);
        context->maxwellTime = a[aOff[i_maxwellTime]];assert(context->maxwellTime > 0.0);

        tensorOps.fromVector(&a[aOff[i_viscousStrain]], &context->viscousStrain);
        tensorOps.fromVector(&a[aOff[i_totalStrain]], &context->totalStrain);
        tensorOps.fromVector(&a[aOff[i_refStress]], &context->refStress);
        tensorOps.fromVector(&a[aOff[i_refStrain]], &context->refStrain);
    } // createContext

    // --------------------------------------------------------------------------------------------
    /** Viscous strain coefficient function for Maxwell viscoelastic materials.
     *
     * @param[in] dt Time step size.
     * @param[in] maxwellTime Relaxation time for material.
     *
     * @returns Viscous strain coefficient.
     */
    static inline
    PylithReal viscousStrainCoeff(const PylithReal dt,
                                  const PylithReal maxwellTime) {
        return maxwellTime*(1.0-exp(-dt/maxwellTime))/dt;
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate viscous strain as a vector.
     *
     * Used to output of viscous strain.
     *
     */
    static inline
    void viscousStrain_asVector(const pylith::fekernels::Elasticity::StrainContext& strainContext,
                                const Context& rheologyContext,
                                pylith::fekernels::Elasticity::strainfn_type strainFn,
                                const pylith::fekernels::TensorOps& tensorOps,
                                PylithScalar viscousStrainVector[]) {
        assert(viscousStrainVector);

        Tensor strain;
        strainFn(strainContext, &strain);

        const PylithReal dt = rheologyContext.dt;
        const PylithReal maxwellTime = rheologyContext.maxwellTime;
        const pylith::fekernels::Tensor& totalStrain = rheologyContext.totalStrain;
        const pylith::fekernels::Tensor& viscousStrainPrev = rheologyContext.viscousStrain;
        pylith::fekernels::Tensor viscousStrainTensor;
        viscousStrain(maxwellTime, viscousStrainPrev, totalStrain, strain, dt, &viscousStrainTensor);

        tensorOps.toVector(viscousStrainTensor, viscousStrainVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate viscous strain.
     */
    static inline
    void viscousStrain(const PylithReal maxwellTime,
                       const pylith::fekernels::Tensor& viscousStrainPrev,
                       const pylith::fekernels::Tensor& totalStrain,
                       const pylith::fekernels::Tensor& strain,
                       const PylithReal dt,
                       pylith::fekernels::Tensor* viscousStrain) {
        assert(viscousStrain);

        pylith::fekernels::Tensor devStrain;
        pylith::fekernels::Elasticity::deviatoric(strain, &devStrain);

        pylith::fekernels::Tensor devTotalStrain;
        pylith::fekernels::Elasticity::deviatoric(totalStrain, &devTotalStrain);

        const PylithScalar dq = pylith::fekernels::IsotropicLinearMaxwell::viscousStrainCoeff(dt, maxwellTime);
        const PylithScalar expFac = exp(-dt/maxwellTime);
        viscousStrain->xx = expFac * viscousStrainPrev.xx + dq * (devStrain.xx - devTotalStrain.xx);
        viscousStrain->yy = expFac * viscousStrainPrev.yy + dq * (devStrain.yy - devTotalStrain.yy);
        viscousStrain->zz = expFac * viscousStrainPrev.zz + dq * (devStrain.zz - devTotalStrain.zz);
        viscousStrain->xy = expFac * viscousStrainPrev.xy + dq * (devStrain.xy - devTotalStrain.xy);
        viscousStrain->yz = expFac * viscousStrainPrev.yz + dq * (devStrain.yz - devTotalStrain.yz);
        viscousStrain->xz = expFac * viscousStrainPrev.xz + dq * (devStrain.xz - devTotalStrain.xz);
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate Cauchy stress WITHOUT a reference stress and strain.
     *
     * ISA pylith::fekernels::Elasticity::stressfn_type
     *
     */
    static inline
    void cauchyStress(void* rheologyContext,
                      const pylith::fekernels::Tensor& strain,
                      const pylith::fekernels::TensorOps& tensorOps,
                      pylith::fekernels::Tensor* stress) {
        assert(stress);
        Context* context = (Context*)(rheologyContext);assert(context);

        const PylithReal bulkModulus = context->bulkModulus;
        pylith::fekernels::IsotropicLinearElasticity::meanStress(bulkModulus, strain, stress);

        const PylithReal dt = context->dt;
        const PylithReal maxwellTime = context->maxwellTime;
        const pylith::fekernels::Tensor& totalStrain = context->totalStrain;
        const pylith::fekernels::Tensor& viscousStrainPrev = context->viscousStrain;
        pylith::fekernels::Tensor viscousStrain;
        pylith::fekernels::IsotropicLinearMaxwell::viscousStrain(maxwellTime, viscousStrainPrev, totalStrain, strain, dt, &viscousStrain);

        const PylithReal shearModulus = context->shearModulus;
        _deviatoricStress(shearModulus, viscousStrain, stress);

    } // cauchyStress

    // --------------------------------------------------------------------------------------------
    /** Calculate Cauchy stress WITH a reference stress/strain.
     *
     * ISA pylith::fekernels::Elasticity::stressfn_type
     *
     */
    static inline
    void cauchyStress_refState(void* rheologyContext,
                               const pylith::fekernels::Tensor& strain,
                               const pylith::fekernels::TensorOps& tensorOps,
                               pylith::fekernels::Tensor* stress) {
        assert(stress);
        Context* context = (Context*)(rheologyContext);assert(context);

        const pylith::fekernels::Tensor& refStress = context->refStress;
        const pylith::fekernels::Tensor& refStrain = context->refStrain;

        const PylithReal bulkModulus = context->bulkModulus;
        pylith::fekernels::IsotropicLinearElasticity::meanStress_refState(bulkModulus, refStress, refStrain, strain, stress);

        const PylithReal dt = context->dt;
        const PylithReal maxwellTime = context->maxwellTime;
        const pylith::fekernels::Tensor& totalStrain = context->totalStrain;
        const pylith::fekernels::Tensor& viscousStrainPrev = context->viscousStrain;
        pylith::fekernels::Tensor viscousStrain;
        pylith::fekernels::IsotropicLinearMaxwell::viscousStrain(maxwellTime, viscousStrainPrev, totalStrain, strain, dt, &viscousStrain);

        const PylithReal shearModulus = context->shearModulus;
        _deviatoricStress_refState(shearModulus, refStress, refStrain, viscousStrain, stress);
    } // cauchyStress_refState

    // --------------------------------------------------------------------------------------------
    /** Calculate Cauchy stress WITHOUT a reference stress and strain using current state variables.
     *
     * ISA pylith::fekernels::Elasticity::stressfn_type
     */
    static inline
    void cauchyStress_stateVars(void* rheologyContext,
                                const pylith::fekernels::Tensor& strain,
                                const pylith::fekernels::TensorOps& tensorOps,
                                pylith::fekernels::Tensor* stress) {
        assert(stress);
        Context* context = (Context*)(rheologyContext);assert(context);

        const PylithReal bulkModulus = context->bulkModulus;
        pylith::fekernels::IsotropicLinearElasticity::meanStress(bulkModulus, strain, stress);

        const PylithReal shearModulus = context->shearModulus;
        const pylith::fekernels::Tensor& viscousStrain = context->viscousStrain;
        _deviatoricStress(shearModulus, viscousStrain, stress);
    } // cauchyStress_stateVars

    // --------------------------------------------------------------------------------------------
    /** Calculate Cauchy stress WITH a reference stress/strain using current state variables.
     *
     * ISA pylith::fekernels::Elasticity::stressfn_type
     *
     */
    static inline
    void cauchyStress_refState_stateVars(void* rheologyContext,
                                         const pylith::fekernels::Tensor& strain,
                                         const pylith::fekernels::TensorOps& tensorOps,
                                         pylith::fekernels::Tensor* stress) {
        assert(stress);
        Context* context = (Context*)(rheologyContext);assert(context);

        const pylith::fekernels::Tensor& refStress = context->refStress;
        const pylith::fekernels::Tensor& refStrain = context->refStrain;

        const PylithReal bulkModulus = context->bulkModulus;
        pylith::fekernels::IsotropicLinearElasticity::meanStress_refState(bulkModulus, refStress, refStrain, strain, stress);

        const PylithReal shearModulus = context->shearModulus;
        const pylith::fekernels::Tensor& viscousStrain = context->viscousStrain;
        _deviatoricStress_refState(shearModulus, refStress, refStrain, viscousStrain, stress);
    } // cauchyStress_refState

private:

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITHOUT reference stress and strain.
     */
    static inline
    void _deviatoricStress(const PylithReal shearModulus,
                           const pylith::fekernels::Tensor& viscousStrain,
                           pylith::fekernels::Tensor* stress) {
        assert(shearModulus > 0.0);
        assert(stress);

        stress->xx += 2.0 * shearModulus * viscousStrain.xx;
        stress->yy += 2.0 * shearModulus * viscousStrain.yy;
        stress->zz += 2.0 * shearModulus * viscousStrain.zz;
        stress->xy += 2.0 * shearModulus * viscousStrain.xy;
        stress->yz += 2.0 * shearModulus * viscousStrain.yz;
        stress->xz += 2.0 * shearModulus * viscousStrain.xz;
    } // deviatoricStress

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITH reference stress and strain.
     */
    static inline
    void _deviatoricStress_refState(const PylithReal shearModulus,
                                    const pylith::fekernels::Tensor& refStress,
                                    const pylith::fekernels::Tensor& refStrain,
                                    const pylith::fekernels::Tensor& viscousStrain,
                                    pylith::fekernels::Tensor* stress) {
        assert(shearModulus > 0.0);
        assert(stress);

        pylith::fekernels::Tensor devRefStrain;
        pylith::fekernels::Elasticity::deviatoric(refStrain, &devRefStrain);

        pylith::fekernels::Tensor devRefStress;
        pylith::fekernels::Elasticity::deviatoric(refStress, &devRefStress);

        stress->xx += devRefStress.xx + 2.0 * shearModulus * (viscousStrain.xx - devRefStrain.xx);
        stress->yy += devRefStress.yy + 2.0 * shearModulus * (viscousStrain.yy - devRefStrain.yy);
        stress->zz += devRefStress.zz + 2.0 * shearModulus * (viscousStrain.zz - devRefStrain.zz);
        stress->xy += devRefStress.xy + 2.0 * shearModulus * (viscousStrain.xy - devRefStrain.xy);
        stress->yz += devRefStress.yz + 2.0 * shearModulus * (viscousStrain.yz - devRefStrain.yz);
        stress->xz += devRefStress.xz + 2.0 * shearModulus * (viscousStrain.xz - devRefStrain.xz);
    }

}; // IsotropicLinearMaxwell

// ------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linear Maxwell plane strain.
class pylith::fekernels::IsotropicLinearMaxwellPlaneStrain {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // ===========================================================================================
    // Kernels for elasticity equation
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for isotropic linear Maxwell plane strain with infinitesimal strain WITHOUT
     * reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for isotropic linear Maxwell plane strain with infinitesimal strain WITH
     * reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    maxwell_time(1), viscous_strain(4), total_strain(4)]
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress_refState,
            pylith::fekernels::Tensor::ops2D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** Jf3_vu entry function for 2D plane strain isotropic linear Maxwell viscoelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
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

        pylith::fekernels::IsotropicLinearMaxwell::Context context;
        pylith::fekernels::IsotropicLinearMaxwell::setContext(
            &context, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        const PylithScalar shearModulus = context.shearModulus;
        const PylithScalar bulkModulus = context.bulkModulus;
        const PylithScalar maxwellTime = context.maxwellTime;
        const PylithScalar dt = context.dt;

        const PylithScalar dq = pylith::fekernels::IsotropicLinearMaxwell::viscousStrainCoeff(dt, maxwellTime);

        // Unique components of Jacobian.
        const PylithReal C1111 = bulkModulus + 4.0/3.0 * shearModulus * dq;
        const PylithReal C1122 = bulkModulus - 2.0/3.0 * shearModulus * dq;
        const PylithReal C1212 = shearModulus * dq;

        /* j(f,g,df,dg) = C(f,df,g,dg)
         *
         * 0:  j0000 = C1111 = 1.0*bulkModulus + 1.33333333333333*delHM*shearModulus
         * 1:  j0001 = C1112 = 0
         * 2:  j0010 = C1211 = 0
         * 3:  j0011 = C1212 = 1.0*delHM*shearModulus
         * 4:  j0100 = C1121 = 0
         * 5:  j0101 = C1122 = 1.0*bulkModulus - 0.666666666666667*delHM*shearModulus
         * 6:  j0110 = C1221 = 1.0*delHM*shearModulus
         * 7:  j0111 = C1222 = 0
         * 8:  j1000 = C2111 = 0
         * 9:  j1001 = C2112 = 1.0*delHM*shearModulus
         * 10:  j1010 = C2211 = 1.0*bulkModulus - 0.666666666666667*delHM*shearModulus
         * 11:  j1011 = C2212 = 0
         * 12:  j1100 = C2121 = 1.0*delHM*shearModulus
         * 13:  j1101 = C2122 = 0
         * 14:  j1110 = C2221 = 0
         * 15:  j1111 = C2222 = 1.0*bulkModulus + 1.33333333333333*delHM*shearModulus
         */

        /* Nonzero Jacobian entries. */
        Jf3[0] -= C1111; /* j0000 */
        Jf3[3] -= C1212; /* j0011 */
        Jf3[5] -= C1122; /* j0101 */
        Jf3[6] -= C1212; /* j0110 */
        Jf3[9] -= C1212; /* j1001 */
        Jf3[10] -= C1122; /* j1010 */
        Jf3[12] -= C1212; /* j1100 */
        Jf3[15] -= C1111; /* j1111 */
    }

    // ===========================================================================================
    // Kernels for fault interfaces and elasticity
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for isotropic linear Maxwell viscoelasticity
     * plane strain with infinitesimal strain WITHOUT reference stress and reference strain.
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for isotropic linear Maxwell viscoelasticity
     * plane strain with infinitesimal strain WITHOUT reference stress and reference strain.
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for isotropic linear Maxwell viscoelasticity
     * plane strain with infinitesimal strain WITH reference stress and reference strain.
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress_refState,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for isotropic linear Maxwell viscoelasticity
     * plane strain with infinitesimal strain WITH reference stress and reference strain.
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress_refState,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // ===========================================================================================
    // Kernels for updating state variables
    // ===========================================================================================

    // Use pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector() to update total strain.

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating viscous strain as a vector for 2D plane strain isotropic
     * linear Maxwell viscoelasticity.
     *
     * Used to output viscous strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void viscousStrain_infinitesimalStrain_asVector(const PylithInt dim,
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
                                                    PylithScalar viscousStrain[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearMaxwell::viscousStrain_asVector(
            strainContext, rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::Tensor::ops2D,
            viscousStrain);
    }

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 2D plane strain isotropic linear
     * Maxwell viscoelasticity with infinitesimal strain WITHOUT a reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress_stateVars,
            pylith::fekernels::Tensor::ops2D,
            stressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 2D plane strain isotropic linear
     * Maxwell viscoelasticity with infinitesimal strain WITH a reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    maxwell_time(1), viscous_strain(4), total_strain(4)]
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress_refState_stateVars,
            pylith::fekernels::Tensor::ops2D,
            stressVector);
    }

}; // IsotropicLinearMaxwellPlaneStrain

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linear Maxwell in 3D.
class pylith::fekernels::IsotropicLinearMaxwell3D {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    // ===========================================================================================
    // Kernels for elasticity equation
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for isotropic linear Maxwell 3D with infinitesimal strain WITHOUT
     * reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), total_strain(6), viscous_strain(6)]
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for isotropic linear Maxwell 3D with infinitesimal strain WITH reference
     * stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., reference_stress(6), reference_strain(6), shear_modulus(1), bulk_modulus(1),
     *                    maxwell_time(1), total_strain(6), viscous_strain(6)]
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress_refState,
            pylith::fekernels::Tensor::ops3D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** Jf3_vu entry function for 3-D isotropic linear Maxwell viscoelasticity WITHOUT reference stress and
     * reference strain.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), total_strain(6), viscous_strain(6)]
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

        pylith::fekernels::IsotropicLinearMaxwell::Context context;
        pylith::fekernels::IsotropicLinearMaxwell::setContext(
            &context, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        const PylithScalar shearModulus = context.shearModulus;
        const PylithScalar bulkModulus = context.bulkModulus;
        const PylithScalar maxwellTime = context.maxwellTime;
        const PylithScalar dt = context.dt;

        const PylithScalar dq = pylith::fekernels::IsotropicLinearMaxwell::viscousStrainCoeff(dt, maxwellTime);

        /* Unique components of Jacobian. */
        const PylithReal C1111 = bulkModulus + 4.0*dq*shearModulus/3.0;
        const PylithReal C1122 = bulkModulus - 2.0*dq*shearModulus/3.0;
        const PylithReal C1212 = dq*shearModulus;

        /* j(f,g,df,dg) = C(f,df,g,dg)
         *
         * 0:  j0000 = C1111 = bulkModulus + 4*dq*shearModulus/3
         * 1:  j0001 = C1112 = 0
         * 2:  j0002 = C1113 = 0
         * 3:  j0010 = C1211 = 0
         * 4:  j0011 = C1212 = dq*shearModulus
         * 5:  j0012 = C1213 = 0
         * 6:  j0020 = C1311 = 0
         * 7:  j0021 = C1312 = 0
         * 8:  j0022 = C1313 = dq*shearModulus
         * 9:  j0100 = C1121 = 0
         * 10:  j0101 = C1122 = bulkModulus - 2*dq*shearModulus/3
         * 11:  j0102 = C1123 = 0
         * 12:  j0110 = C1221 = dq*shearModulus
         * 13:  j0111 = C1222 = 0
         * 14:  j0112 = C1223 = 0
         * 15:  j0120 = C1321 = 0
         * 16:  j0121 = C1322 = 0
         * 17:  j0122 = C1323 = 0
         * 18:  j0200 = C1131 = 0
         * 19:  j0201 = C1132 = 0
         * 20:  j0202 = C1133 = bulkModulus - 2*dq*shearModulus/3
         * 21:  j0210 = C1231 = 0
         * 22:  j0211 = C1232 = 0
         * 23:  j0212 = C1233 = 0
         * 24:  j0220 = C1331 = dq*shearModulus
         * 25:  j0221 = C1332 = 0
         * 26:  j0222 = C1333 = 0
         * 27:  j1000 = C2111 = 0
         * 28:  j1001 = C2112 = dq*shearModulus
         * 29:  j1002 = C2113 = 0
         * 30:  j1010 = C2211 = bulkModulus - 2*dq*shearModulus/3
         * 31:  j1011 = C2212 = 0
         * 32:  j1012 = C2213 = 0
         * 33:  j1020 = C2311 = 0
         * 34:  j1021 = C2312 = 0
         * 35:  j1022 = C2313 = 0
         * 36:  j1100 = C2121 = dq*shearModulus
         * 37:  j1101 = C2122 = 0
         * 38:  j1102 = C2123 = 0
         * 39:  j1110 = C2221 = 0
         * 40:  j1111 = C2222 = bulkModulus + 4*dq*shearModulus/3
         * 41:  j1112 = C2223 = 0
         * 42:  j1120 = C2321 = 0
         * 43:  j1121 = C2322 = 0
         * 44:  j1122 = C2323 = dq*shearModulus
         * 45:  j1200 = C2131 = 0
         * 46:  j1201 = C2132 = 0
         * 47:  j1202 = C2133 = 0
         * 48:  j1210 = C2231 = 0
         * 49:  j1211 = C2232 = 0
         * 50:  j1212 = C2233 = bulkModulus - 2*dq*shearModulus/3
         * 51:  j1220 = C2331 = 0
         * 52:  j1221 = C2332 = dq*shearModulus
         * 53:  j1222 = C2333 = 0
         * 54:  j2000 = C3111 = 0
         * 55:  j2001 = C3112 = 0
         * 56:  j2002 = C3113 = dq*shearModulus
         * 57:  j2010 = C3211 = 0
         * 58:  j2011 = C3212 = 0
         * 59:  j2012 = C3213 = 0
         * 60:  j2020 = C3311 = bulkModulus - 2*dq*shearModulus/3
         * 61:  j2021 = C3312 = 0
         * 62:  j2022 = C3313 = 0
         * 63:  j2100 = C3121 = 0
         * 64:  j2101 = C3122 = 0
         * 65:  j2102 = C3123 = 0
         * 66:  j2110 = C3221 = 0
         * 67:  j2111 = C3222 = 0
         * 68:  j2112 = C3223 = dq*shearModulus
         * 69:  j2120 = C3321 = 0
         * 70:  j2121 = C3322 = bulkModulus - 2*dq*shearModulus/3
         * 71:  j2122 = C3323 = 0
         * 72:  j2200 = C3131 = dq*shearModulus
         * 73:  j2201 = C3132 = 0
         * 74:  j2202 = C3133 = 0
         * 75:  j2210 = C3231 = 0
         * 76:  j2211 = C3232 = dq*shearModulus
         * 77:  j2212 = C3233 = 0
         * 78:  j2220 = C3331 = 0
         * 79:  j2221 = C3332 = 0
         * 80:  j2222 = C3333 = bulkModulus + 4*dq*shearModulus/3
         */

        /* Nonzero Jacobian entries. */
        Jf3[0] -= C1111; /* j0000 */
        Jf3[4] -= C1212; /* j0011 */
        Jf3[8] -= C1212; /* j0022 */
        Jf3[10] -= C1122; /* j0101 */
        Jf3[12] -= C1212; /* j0110 */
        Jf3[20] -= C1122; /* j0202 */
        Jf3[24] -= C1212; /* j0220 */
        Jf3[28] -= C1212; /* j1001 */
        Jf3[30] -= C1122; /* j1010 */
        Jf3[36] -= C1212; /* j1100 */
        Jf3[40] -= C1111; /* j1111 */
        Jf3[44] -= C1212; /* j1122 */
        Jf3[50] -= C1122; /* j1212 */
        Jf3[52] -= C1212; /* j1221 */
        Jf3[56] -= C1212; /* j2002 */
        Jf3[60] -= C1122; /* j2020 */
        Jf3[68] -= C1212; /* j2112 */
        Jf3[70] -= C1122; /* j2121 */
        Jf3[72] -= C1212; /* j2200 */
        Jf3[76] -= C1212; /* j2211 */
        Jf3[80] -= C1111; /* j2222 */
    } // Jf3vu_infinitesimalStrain

    // ===========================================================================================
    // Kernels for fault interfaces and elasticity
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for 3D isotropic linear Maxwell viscoelasticity
     * with infinitesimal strain WITHOUT reference stress and reference strain.
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for 3D isotropic linear Maxwell viscoelasticity
     * with infinitesimal strain WITHOUT reference stress and reference strain.
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for 3D isotropic linear Maxwell viscoelasticity
     * with infinitesimal strain WITH reference stress and reference strain.
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress_refState,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for 3D isotropic linear Maxwell viscoelasticity
     * with infinitesimal strain WITH reference stress and reference strain.
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress_refState,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // ===========================================================================================
    // Kernels for updating state variables
    // ===========================================================================================

    // Use pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector() to update total strain.

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating viscous strain for 3D isotropic linear Maxwell
     * viscoelasticity.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void viscousStrain_infinitesimalStrain_asVector(const PylithInt dim,
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
                                                    PylithScalar viscousStrain[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearMaxwell::viscousStrain_asVector(
            strainContext, rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::Tensor::ops3D,
            viscousStrain);
    }

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating stress for 3D isotropic linear Maxwell viscoelasticity
     * WITHOUT a reference stress and strain.
     *
     * Used in output of Cauchy stress.
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress_stateVars,
            pylith::fekernels::Tensor::ops3D,
            stressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating stress for 3D isotropic linear Maxwell WITH a reference
     * stress and strain.
     *
     * Used in output of Cauchy stress.
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

        pylith::fekernels::IsotropicLinearMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearMaxwell::cauchyStress_refState_stateVars,
            pylith::fekernels::Tensor::ops3D,
            stressVector);
    }

}; // IsotropicLinearMaxwell3D

#endif // pylith_fekernels_isotropiclinearmaxwell_hh

// End of file
