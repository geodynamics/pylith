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

/** @file libsrc/fekernels/IsotropicLinearGenMaxwell.hh
 *
 * Kernels for linear Generalized Maxwell viscoelastic 3 Maxwell elements.
 *
 * Solution fields: [disp(dim), ...]
 *
 * Isotropic, linear Generalized Maxwell viscoelastic plane strain without reference stress/strain.
 * Auxiliary fields:
 * - 0: density(1)
 * - 1: body_force(dim,optional)
 * - 2: gravity_field (dim, optional)
 * - 3: reference_stress(optional) (stress_xx, stress_yy, stress_zz, stress_xy)
 *     2D: 4 components (stress_xx, stress_yy, stress_zz, stress_xy)
 *     3D: 6 components (stress_xx, stress_yy, stress_zz, stress_xy, stress_yz, stress_xz)
 * - 4: reference_strain(optional) (strain_xx, strain_yy, strain_zz, strain_xy)
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 * - 5: shear_modulus(1)
 * - 6: bulk_modulus(1)
 * - 7: maxwell_time(3) (maxwell_time_1, maxwell_time_2, maxwell_time_3)
 * - 8: shear_modulus_ratio(3) (shear_modulus_ratio_1, shear_modulus_ratio_2, shear_modulus_ratio_3)
 * - 9: viscous_strain
 *     2D: 3*4 components (strain1_xx, strain1_yy, strain1_zz, strain1_xy, ...)
 *     3D: 3*6 components (strain1_xx, strain1_yy, strain1_zz, strain1_xy, strain1_yz, strain1_xz, ...)
 * - 10: total_strain
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 *
 * The elasticity subfields come first (with required ones before optional ones) followed by the rheology subfields
 * (optional ones before required ones). The rheology fields have required fields last because we index from the back.
 *
 * Viscous strain must be before total strain, because viscous strain at t+dt depends on total strain at t.
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

#if !defined(pylith_fekernels_isotropiclineargenmaxwell_hh)
#define pylith_fekernels_isotropiclineargenmaxwell_hh

#include "fekernelsfwd.hh" // forward declarations
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity kernels
#include "pylith/fekernels/IsotropicLinearMaxwell.hh" // USES IsotropicLinearMaxwell kernels

#include "pylith/utils/types.hh"

// ------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linear generalized Maxwell viscoelasticity (dimension independent).
class pylith::fekernels::IsotropicLinearGenMaxwell {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    static const PylithInt numParallel;

    struct Context {
        PylithReal shearModulus;
        PylithReal bulkModulus;
        const PylithReal* maxwellTime; // size is numParallel
        const PylithReal* shearModulusRatio; // size is numParallel
        PylithReal dt;
        pylith::fekernels::Tensor viscousStrain[3]; // Size must match numParallel
        pylith::fekernels::Tensor totalStrain;
        pylith::fekernels::Tensor refStress;
        pylith::fekernels::Tensor refStrain;

        Context(void) :
            shearModulus(0.0),
            bulkModulus(0.0),
            maxwellTime(nullptr),
            shearModulusRatio(nullptr),
            dt(0.0) {}


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

        const PylithInt i_shearModulus = numA-6;
        const PylithInt i_bulkModulus = numA-5;
        const PylithInt i_maxwellTime = numA-4;
        const PylithInt i_shearModulusRatio = numA-3;
        const PylithInt i_viscousStrain = numA-2;
        const PylithInt i_totalStrain = numA-1;

        assert(numA >= 7); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_maxwellTime] >= 0);
        assert(aOff[i_shearModulusRatio] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_totalStrain] >= 0);
        assert(1 == numConstants);
        assert(constants);

        context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
        context->bulkModulus = a[aOff[i_bulkModulus]];assert(context->bulkModulus > 0.0);
        context->maxwellTime = &a[aOff[i_maxwellTime]];assert(context->maxwellTime);
        context->shearModulusRatio = &a[aOff[i_shearModulusRatio]];assert(context->shearModulusRatio);
        context->dt = constants[0];

        for (PylithInt i = 0; i < numParallel; ++i) {
            const PylithInt offset = i*tensorOps.vectorSize;
            tensorOps.fromVector(&a[aOff[i_viscousStrain]+offset], &context->viscousStrain[i]);
        } // for
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

        const PylithInt i_refStress = numA-8;
        const PylithInt i_refStrain = numA-7;
        const PylithInt i_shearModulus = numA-6;
        const PylithInt i_bulkModulus = numA-5;
        const PylithInt i_maxwellTime = numA-4;
        const PylithInt i_shearModulusRatio = numA-3;
        const PylithInt i_viscousStrain = numA-2;
        const PylithInt i_totalStrain = numA-1;

        assert(numA >= 9); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_maxwellTime] >= 0);
        assert(aOff[i_shearModulusRatio] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_totalStrain] >= 0);
        assert(1 == numConstants);
        assert(constants);

        context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
        context->bulkModulus = a[aOff[i_bulkModulus]];assert(context->bulkModulus > 0.0);
        context->maxwellTime = &a[aOff[i_maxwellTime]];assert(context->maxwellTime);
        context->shearModulusRatio = &a[aOff[i_shearModulusRatio]];assert(context->shearModulusRatio);
        context->dt = constants[0];

        for (PylithInt i = 0; i < numParallel; ++i) {
            const PylithInt offset = i*tensorOps.vectorSize;
            tensorOps.fromVector(&a[aOff[i_viscousStrain]+offset], &context->viscousStrain[i]);
        } // for
        tensorOps.fromVector(&a[aOff[i_totalStrain]], &context->totalStrain);
        tensorOps.fromVector(&a[aOff[i_refStress]], &context->refStress);
        tensorOps.fromVector(&a[aOff[i_refStrain]], &context->refStrain);
    } // createContext

    // --------------------------------------------------------------------------------------------
    /** Calculate viscous strain as a vector.
     *
     * Use in output of viscous strain.
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
        const pylith::fekernels::Tensor& totalStrain = rheologyContext.totalStrain;
        for (PylithInt i = 0; i < numParallel; ++i ) {
            const PylithReal maxwellTime = rheologyContext.maxwellTime[i];
            const pylith::fekernels::Tensor& viscousStrainPrev = rheologyContext.viscousStrain[i];
            pylith::fekernels::Tensor viscousStrainTensor;

            pylith::fekernels::IsotropicLinearMaxwell::viscousStrain(maxwellTime, viscousStrainPrev, totalStrain, strain, dt, &viscousStrainTensor);

            const PylithInt offset = i*tensorOps.vectorSize;
            tensorOps.toVector(viscousStrainTensor, &viscousStrainVector[offset]);
        } // for
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

        const PylithReal shearModulus = context->shearModulus;assert(shearModulus);
        const PylithReal* shearModulusRatio = context->shearModulusRatio;
        const PylithReal dt = context->dt;
        const pylith::fekernels::Tensor& totalStrain = context->totalStrain;

        pylith::fekernels::Tensor viscousStrain[numParallel];
        for (PylithInt i = 0; i < numParallel; ++i ) {
            const PylithReal maxwellTime = context->maxwellTime[i];
            const pylith::fekernels::Tensor& viscousStrainPrev = context->viscousStrain[i];
            pylith::fekernels::IsotropicLinearMaxwell::viscousStrain(maxwellTime, viscousStrainPrev, totalStrain, strain, dt, &viscousStrain[i]);
        } // for

        pylith::fekernels::Tensor devStrain;
        pylith::fekernels::Elasticity::deviatoric(strain, &devStrain);
        _deviatoricStress(shearModulus, shearModulusRatio, devStrain, viscousStrain, stress);
    } // cauchyStress

    // --------------------------------------------------------------------------------------------
    /** Calculate Cauchy stress WITH a reference stress/strain.
     *
     * ISA pylith::fekernels::Elasticity::stressfn_type
     *
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                   maxwell_time(1), viscous_strain(4), total_strain(4)]
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

        const PylithReal shearModulus = context->shearModulus;
        const PylithReal* shearModulusRatio = context->shearModulusRatio;
        const PylithReal dt = context->dt;
        const pylith::fekernels::Tensor& totalStrain = context->totalStrain;

        pylith::fekernels::Tensor viscousStrain[numParallel];
        for (PylithInt i = 0; i < numParallel; ++i) {
            const PylithReal maxwellTime = context->maxwellTime[i];
            const pylith::fekernels::Tensor& viscousStrainPrev = context->viscousStrain[i];
            pylith::fekernels::IsotropicLinearMaxwell::viscousStrain(maxwellTime, viscousStrainPrev, totalStrain, strain, dt, &viscousStrain[i]);
        } // for

        pylith::fekernels::Tensor devStrain;
        pylith::fekernels::Elasticity::deviatoric(strain, &devStrain);
        _deviatoricStress_refState(shearModulus, shearModulusRatio, refStress, refStrain, devStrain, viscousStrain, stress);
    } // cauchyStress_refState

    // --------------------------------------------------------------------------------------------
    /** Calculate Cauchy stress WITHOUT a reference stress and strain using current state variables.
     *
     * ISA pylith::fekernels::Elasticity::stressfn_type
     *
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
        const PylithReal* shearModulusRatio = context->shearModulusRatio;
        const pylith::fekernels::Tensor* viscousStrain = context->viscousStrain;

        pylith::fekernels::Tensor devStrain;
        pylith::fekernels::Elasticity::deviatoric(strain, &devStrain);
        _deviatoricStress(shearModulus, shearModulusRatio, devStrain, viscousStrain, stress);
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
        const PylithReal* shearModulusRatio = context->shearModulusRatio;
        const pylith::fekernels::Tensor* viscousStrain = context->viscousStrain;

        pylith::fekernels::Tensor devStrain;
        pylith::fekernels::Elasticity::deviatoric(strain, &devStrain);
        _deviatoricStress_refState(shearModulus, shearModulusRatio, refStress, refStrain, devStrain, viscousStrain, stress);
    } // cauchyStress_refState

    static inline
    PylithReal shearModulusElastic(const PylithReal shearModulus,
                                   const PylithReal* shearModulusRatio) {
        PylithReal ratio = 1.0;
        for (PylithInt i = 0; i < numParallel; ++i) {
            ratio -= shearModulusRatio[i];
        } // for
        return shearModulus * ratio;
    }

private:

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITHOUT reference stress and strain.
     */
    static inline
    void _deviatoricStress(const PylithReal shearModulus,
                           const PylithReal shearModulusRatio[],
                           const pylith::fekernels::Tensor& devStrain,
                           const pylith::fekernels::Tensor viscousStrain[],
                           pylith::fekernels::Tensor* stress) {
        const size_t numParallel = pylith::fekernels::IsotropicLinearGenMaxwell::numParallel;

        assert(shearModulus > 0.0);
        assert(shearModulusRatio);
        assert(viscousStrain);
        assert(stress);

        const PylithReal shearModulusElastic = pylith::fekernels::IsotropicLinearGenMaxwell::shearModulusElastic(shearModulus, shearModulusRatio);

        stress->xx += 2.0 * shearModulusElastic * devStrain.xx;
        for (size_t i = 0; i < numParallel; ++i) {
            stress->xx += 2.0 * shearModulus * shearModulusRatio[i] * viscousStrain[i].xx;
        } // for

        stress->yy += 2.0 * shearModulusElastic * devStrain.yy;
        for (size_t i = 0; i < numParallel; ++i) {
            stress->yy += 2.0 * shearModulus * shearModulusRatio[i] * viscousStrain[i].yy;
        } // for

        stress->zz += 2.0 * shearModulusElastic * devStrain.zz;
        for (size_t i = 0; i < numParallel; ++i) {
            stress->zz += 2.0 * shearModulus * shearModulusRatio[i] * viscousStrain[i].zz;
        } // for

        stress->xy += 2.0 * shearModulusElastic * devStrain.xy;
        for (size_t i = 0; i < numParallel; ++i) {
            stress->xy += 2.0 * shearModulus * shearModulusRatio[i] * viscousStrain[i].xy;
        } // for

        stress->yz += 2.0 * shearModulusElastic * devStrain.yz;
        for (size_t i = 0; i < numParallel; ++i) {
            stress->yz += 2.0 * shearModulus * shearModulusRatio[i] * viscousStrain[i].yz;
        } // for

        stress->xz += 2.0 * shearModulusElastic * devStrain.xz;
        for (size_t i = 0; i < numParallel; ++i) {
            stress->xz += 2.0 * shearModulus * shearModulusRatio[i] * viscousStrain[i].xz;
        } // for
    } // deviatoricStress

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITH reference stress and strain.
     */
    static inline
    void _deviatoricStress_refState(const PylithReal shearModulus,
                                    const PylithReal shearModulusRatio[],
                                    const pylith::fekernels::Tensor& refStress,
                                    const pylith::fekernels::Tensor& refStrain,
                                    const pylith::fekernels::Tensor& devStrain,
                                    const pylith::fekernels::Tensor viscousStrain[],
                                    pylith::fekernels::Tensor* stress) {
        const size_t numParallel = pylith::fekernels::IsotropicLinearGenMaxwell::numParallel;

        assert(shearModulus > 0.0);
        assert(shearModulusRatio);
        assert(viscousStrain);
        assert(stress);

        pylith::fekernels::Tensor devRefStrain;
        pylith::fekernels::Elasticity::deviatoric(refStrain, &devRefStrain);

        pylith::fekernels::Tensor devRefStress;
        pylith::fekernels::Elasticity::deviatoric(refStress, &devRefStress);

        const PylithReal shearModulusElastic = pylith::fekernels::IsotropicLinearGenMaxwell::shearModulusElastic(shearModulus, shearModulusRatio);

        stress->xx += devRefStress.xx + 2.0 * (shearModulusElastic * devStrain.xx - shearModulus * devRefStrain.xx);
        for (size_t i = 0; i < numParallel; ++i) {
            stress->xx += 2.0 * shearModulus * shearModulusRatio[i] * viscousStrain[i].xx;
        } // for

        stress->yy += devRefStress.yy + 2.0 * (shearModulusElastic * devStrain.yy - shearModulus * devRefStrain.yy);
        for (size_t i = 0; i < numParallel; ++i) {
            stress->yy += 2.0 * shearModulus * shearModulusRatio[i] * viscousStrain[i].yy;
        } // for

        stress->zz += devRefStress.zz + 2.0 * (shearModulusElastic * devStrain.zz - shearModulus * devRefStrain.zz);
        for (size_t i = 0; i < numParallel; ++i) {
            stress->zz += 2.0 * shearModulus * shearModulusRatio[i] * viscousStrain[i].zz;
        } // for

        stress->xy += devRefStress.xy + 2.0 * (shearModulusElastic * devStrain.xy - shearModulus * devRefStrain.xy);
        for (size_t i = 0; i < numParallel; ++i) {
            stress->xy += 2.0 * shearModulus * shearModulusRatio[i] * viscousStrain[i].xy;
        } // for

        stress->yz += devRefStress.yz + 2.0 * (shearModulusElastic * devStrain.yz - shearModulus * devRefStrain.yz);
        for (size_t i = 0; i < numParallel; ++i) {
            stress->yz += 2.0 * shearModulus * shearModulusRatio[i] * viscousStrain[i].yz;
        } // for

        stress->xz += devRefStress.xz + 2.0 * (shearModulusElastic * devStrain.xz - shearModulus * devRefStrain.xz);
        for (size_t i = 0; i < numParallel; ++i) {
            stress->xz += 2.0 * shearModulus * shearModulusRatio[i] * viscousStrain[i].xz;
        } // for
    }

}; // IsotropicLinearGenMaxwell

// ------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linear generalized Maxwell plane strain.
class pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // ===========================================================================================
    // Kernels for elasticity equation
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for isotropic linear generalized Maxwell plane strain with infinitesimal
     * strain WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
     *                    total_strain(4), viscous_strain(12)]
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            f1);
    }

    // ------------------------------------------------------------------------------------------------
    /** f1 entry function for isotropic linear generalized Maxwell plane strain with infinitesimal
     * strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    maxwell_time(3), shear_modulus_ratio(3), total_strain(4), viscous_strain(12)]
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress_refState,
            pylith::fekernels::Tensor::ops2D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** Jf3_vu entry function for 2-D plane strain isotropic linear generalized Maxwell viscoelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
     *                    total_strain(4), viscous_strain(12)]
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
        const PylithInt numParallel = pylith::fekernels::IsotropicLinearGenMaxwell::numParallel;

        pylith::fekernels::IsotropicLinearGenMaxwell::Context context;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext(
            &context, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        const PylithScalar shearModulus = context.shearModulus;
        const PylithScalar bulkModulus = context.bulkModulus;
        const PylithReal* shearModulusRatio = context.shearModulusRatio;
        const PylithScalar dt = context.dt;

        const PylithReal shearModulusRatio_0 = 1.0 - shearModulusRatio[0] - shearModulusRatio[1] - shearModulusRatio[2];
        PylithReal shearFactor = shearModulus * shearModulusRatio_0;
        for (PylithInt i = 0; i < numParallel; ++i) {
            const PylithReal maxwellTime = context.maxwellTime[i];
            const PylithReal dq = pylith::fekernels::IsotropicLinearMaxwell::viscousStrainCoeff(dt, maxwellTime);
            shearFactor += shearModulus * dq * shearModulusRatio[i];
        } // for

        const PylithReal C1111 = bulkModulus + 4.0/3.0 * shearFactor;
        const PylithReal C1122 = bulkModulus - 2.0/3.0 * shearFactor;
        const PylithReal C1212 = shearFactor;
        /* j(f,g,df,dg) = C(f,df,g,dg)
         *
         * 0:  j0000 = C1111 = bulkModulus + 4.0/3.0*shearModulus * (dq_i*shearModulusRatio_i + shearModulusRatio_0)
         * 1:  j0001 = C1112 = 0
         * 2:  j0010 = C1211 = 0
         * 3:  j0011 = C1212 = shearModulus*(dq_i*shearModulusRatio_i + shearModulusRatio_0)
         * 4:  j0100 = C1121 = 0
         * 5:  j0101 = C1122 = bulkModulus - 2.0/3.0*shearModulus*(dq_i*shearModulusRatio_i + shearModulusRatio_0)
         * 6:  j0110 = C1221 = C1212
         * 7:  j0111 = C1222 = 0
         * 8:  j1000 = C2111 = 0
         * 9:  j1001 = C2112 = C1212
         * 10:  j1010 = C2211 = C1122
         * 11:  j1011 = C2212 = 0
         * 12:  j1100 = C2121 = C1212
         * 13:  j1101 = C2122 = 0
         * 14:  j1110 = C2221 = 0
         * 15:  j1111 = C2222 = C1111
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
    /** f0 entry function for negative fault face for isotropic linear generalized Maxwell
     * viscoelasticity plane strain with infinitesimal strain WITHOUT reference stress and
     * reference strain.
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for isotropic linear generalized Maxwell
     * viscoelasticity plane strain with infinitesimal strain WITHOUT reference stress and
     * reference strain.
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for isotropic linear generalized Maxwell
     * viscoelasticity plane strain with infinitesimal strain WITH reference stress and
     * reference strain.
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress_refState,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for isotropic linear generalized Maxwell
     * viscoelasticity plane strain with infinitesimal strain WITH reference stress and
     * reference strain.
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress_refState,
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
     * linear elasticity.
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearGenMaxwell::viscousStrain_asVector(
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
     * elasticity with infinitesimal strain WITHOUT a reference stress and strain.
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress_stateVars,
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress_refState_stateVars,
            pylith::fekernels::Tensor::ops2D,
            stressVector);
    }

}; // IsotropicLinearGenMaxwellPlaneStrain

// ------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linear generalized Maxwell viscoelastic in 3D.
class pylith::fekernels::IsotropicLinearGenMaxwell3D {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // ===========================================================================================
    // Kernels for elasticity equation
    // ===========================================================================================

    /** f1 entry function for isotropic linear generalized Maxwell 3D with infinitesimal strain WITHOUT reference stress
     * and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
     *                    total_strain(4), viscous_strain(12)]
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for isotropic linear generalized Maxwell 3D with infinitesimal strain
     * WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    maxwell_time(3), shear_modulus_ratio(3), total_strain(4), viscous_strain(12)]
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress_refState,
            pylith::fekernels::Tensor::ops3D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** Jf3_vu entry function for 3-D isotropic linear generalized Maxwell viscoelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
     *                    total_strain(4), viscous_strain(12)]
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
                                   const PylithReal utshift,
                                   const PylithScalar x[],
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar Jf3[]) {
        assert(Jf3);

        const PylithInt _dim = 3;assert(_dim == dim);
        const PylithInt numParallel = pylith::fekernels::IsotropicLinearGenMaxwell::numParallel;

        pylith::fekernels::IsotropicLinearGenMaxwell::Context context;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext(
            &context, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        const PylithReal shearModulus = context.shearModulus;
        const PylithReal bulkModulus = context.bulkModulus;
        const PylithReal* shearModulusRatio = context.shearModulusRatio;
        const PylithReal dt = context.dt;

        const PylithReal shearModulusRatio_0 = 1.0 - shearModulusRatio[0] - shearModulusRatio[1] - shearModulusRatio[2];
        PylithReal shearFactor = shearModulus * shearModulusRatio_0;
        for (PylithInt i = 0; i < numParallel; ++i) {
            const PylithReal maxwellTime = context.maxwellTime[i];
            const PylithReal dq = pylith::fekernels::IsotropicLinearMaxwell::viscousStrainCoeff(dt, maxwellTime);
            shearFactor += shearModulus * dq * shearModulusRatio[i];
        } // for

        // Unique components of Jacobian.
        const PylithReal C1111 = bulkModulus + 4.0 * shearFactor/3.0;
        const PylithReal C1122 = bulkModulus - 2.0 * shearFactor/3.0;
        const PylithReal C1212 = shearFactor;

        /* j(f,g,df,dg) = C(f,df,g,dg)
         * 0:  j0000 = C1111 = bulkModulus + 4.0/3.0*shearModulus * (dq_i*shearModulusRatio_i + shearModulusRatio_0)
         * 1:  j0001 = C1112 = 0
         * 2:  j0002 = C1113 = 0
         * 3:  j0010 = C1211 = 0
         * 4:  j0011 = C1212 = shearModulus*(dq_i*shearModulusRatio_i + shearModulusRatio_0)
         * 5:  j0012 = C1213 = 0
         * 6:  j0020 = C1311 = 0
         * 7:  j0021 = C1312 = 0
         * 8:  j0022 = C1313 = shearModulus*(dq_i*shearModulusRatio_i + shearModulusRatio_0)
         * 9:  j0100 = C1121 = 0
         * 10:  j0101 = C1122 = bulkModulus - 2.0/3.0*shearModulus*(dq_i*shearModulusRatio_i + shearModulusRatio_0)
         * 11:  j0102 = C1123 = 0
         * 12:  j0110 = C1221 = shearModulus*(dq_i*shearModulusRatio_i + shearModulusRatio_0)
         * 13:  j0111 = C1222 = 0
         * 14:  j0112 = C1223 = 0
         * 15:  j0120 = C1321 = 0
         * 16:  j0121 = C1322 = 0
         * 17:  j0122 = C1323 = 0
         * 18:  j0200 = C1131 = 0
         * 19:  j0201 = C1132 = 0
         * 20:  j0202 = C1133 = bulkModulus - 2.0/3.0*shearModulus*(dq_i*shearModulusRatio_i + shearModulusRatio_0)
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
         * 40:  j1111 = C2222 = C1111
         * 41:  j1112 = C2223 = 0
         * 42:  j1120 = C2321 = 0
         * 43:  j1121 = C2322 = 0
         * 44:  j1122 = C2323 = C1212
         * 45:  j1200 = C2131 = 0
         * 46:  j1201 = C2132 = 0
         * 47:  j1202 = C2133 = 0
         * 48:  j1210 = C2231 = 0
         * 49:  j1211 = C2232 = 0
         * 50:  j1212 = C2233 = C1122
         * 51:  j1220 = C2331 = 0
         * 52:  j1221 = C2332 = C1212
         * 53:  j1222 = C2333 = 0
         * 54:  j2000 = C3111 = 0
         * 55:  j2001 = C3112 = 0
         * 56:  j2002 = C3113 = C1212
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
         * 80:  j2222 = C3333 = C1111
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
    }

    // ===========================================================================================
    // Kernels for fault interfaces and elasticity
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for 3D isotropic linear generalized Maxwell
     * viscoelasticity with infinitesimal strain WITHOUT reference stress and
     * reference strain.
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for 3D isotropic linear generalized Maxwell
     * viscoelasticity with infinitesimal strain WITHOUT reference stress and
     * reference strain.
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for 3D isotropic linear generalized Maxwell
     * viscoelasticity with infinitesimal strain WITH reference stress and
     * reference strain.
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress_refState,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for 3D isotropic linear generalized Maxwell
     * viscoelasticity with infinitesimal strain WITH reference stress and
     * reference strain.
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress_refState,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // ===========================================================================================
    // Kernels for updating state variables
    // ===========================================================================================

    // Use pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector() to update total strain.

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating viscous strain for 3D isotropic linear elasticity.
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearGenMaxwell::viscousStrain_asVector(
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress_stateVars,
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

        pylith::fekernels::IsotropicLinearGenMaxwell::Context rheologyContext;
        pylith::fekernels::IsotropicLinearGenMaxwell::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearGenMaxwell::cauchyStress_refState_stateVars,
            pylith::fekernels::Tensor::ops3D,
            stressVector);
    }

}; // IsotropicLinearGenMaxwell3D

#endif // pylith_fekernels_isotropiclineargenmaxwell_hh

// End of file
