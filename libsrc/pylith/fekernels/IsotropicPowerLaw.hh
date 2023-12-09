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

/** @file libsrc/fekernels/IsotropicPowerLaw.hh
 *
 * Kernels for power-law viscoelastic material.
 *
 * Solution fields: [disp(dim), ...]
 *
 * Isotropic, power-law viscoelastic material.
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
 * - 7: power_law_reference_strain_rate(1)
 * - 8: power_law_reference_stress(1)
 * - 9: power_law_exponent(1)
 * -10: viscous_strain
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 * -11: deviatoric_stress
 *     2D: 4 components (devstress_xx, devstress_yy, devstress_zz, devstress_xy)
 *     3D: 6 components (devstress_xx, devstress_yy, devstress_zz, devstress_xy, devstress_yz, devstress_xz)
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

#if !defined(pylith_fekernels_isotropicpowerlaw_hh)
#define pylith_fekernels_isotropicpowerlaw_hh

#include "fekernelsfwd.hh" // forward declarations
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity* kernels

#include "pylith/utils/types.hh"

// ------------------------------------------------------------------------------------------------
/// Kernels for isotropic power-law viscoelasticity (dimension independent).
class pylith::fekernels::IsotropicPowerLaw {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    static const PylithReal powerLawAlpha;
    struct Context {
        PylithReal bulkModulus;
        PylithReal shearModulus;
        PylithReal powerLawRefStress;
        PylithReal powerLawRefStrainRate;
        PylithReal powerLawExponent;
        PylithReal dt;
        pylith::fekernels::Tensor refStress;
        pylith::fekernels::Tensor refStrain;
        pylith::fekernels::Tensor viscousStrain;
        pylith::fekernels::Tensor devStress;

        Context(void) :
            bulkModulus(0.0),
            shearModulus(0.0),
            powerLawRefStress(0.0),
            powerLawRefStrainRate(0.0),
            powerLawExponent(0.0),
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

        const PylithInt i_shearModulus = numA-7;
        const PylithInt i_bulkModulus = numA-6;
        const PylithInt i_powerLawRefStrainRate = numA-5;
        const PylithInt i_powerLawRefStress = numA-4;
        const PylithInt i_powerLawExponent = numA-3;
        const PylithInt i_viscousStrain = numA-2;
        const PylithInt i_devStress = numA-1;

        assert(numA >= 8); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_powerLawRefStrainRate] >= 0);
        assert(aOff[i_powerLawRefStress] >= 0);
        assert(aOff[i_powerLawExponent] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_devStress] >= 0);
        assert(1 == numConstants);
        assert(constants);

        context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
        context->bulkModulus = a[aOff[i_bulkModulus]];assert(context->bulkModulus > 0.0);
        context->powerLawRefStrainRate = a[aOff[i_powerLawRefStrainRate]];assert(context->powerLawRefStrainRate > 0.0);
        context->powerLawRefStress = a[aOff[i_powerLawRefStress]];assert(context->powerLawRefStress > 0.0);
        context->powerLawExponent = a[aOff[i_powerLawExponent]];assert(context->powerLawExponent > 0.0);
        context->dt = constants[0];assert(context->dt > 0.0);

        tensorOps.fromVector(&a[aOff[i_viscousStrain]], &context->viscousStrain);
        tensorOps.fromVector(&a[aOff[i_devStress]], &context->devStress);
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

        const PylithInt i_refStress = numA-9;
        const PylithInt i_refStrain = numA-8;
        const PylithInt i_shearModulus = numA-7;
        const PylithInt i_bulkModulus = numA-6;
        const PylithInt i_powerLawRefStrainRate = numA-5;
        const PylithInt i_powerLawRefStress = numA-4;
        const PylithInt i_powerLawExponent = numA-3;
        const PylithInt i_viscousStrain = numA-2;
        const PylithInt i_devStress = numA-1;

        assert(numA >= 10); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_powerLawRefStrainRate] >= 0);
        assert(aOff[i_powerLawRefStress] >= 0);
        assert(aOff[i_powerLawExponent] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_devStress] >= 0);
        assert(1 == numConstants);
        assert(constants);

        context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
        context->bulkModulus = a[aOff[i_bulkModulus]];assert(context->bulkModulus > 0.0);
        context->powerLawRefStrainRate = a[aOff[i_powerLawRefStrainRate]];assert(context->powerLawRefStrainRate > 0.0);
        context->powerLawRefStress = a[aOff[i_powerLawRefStress]];assert(context->powerLawRefStress > 0.0);
        context->powerLawExponent = a[aOff[i_powerLawExponent]];assert(context->powerLawExponent > 0.0);
        context->dt = constants[0];assert(context->dt > 0.0);

        tensorOps.fromVector(&a[aOff[i_refStress]], &context->refStress);
        tensorOps.fromVector(&a[aOff[i_refStrain]], &context->refStrain);
        tensorOps.fromVector(&a[aOff[i_viscousStrain]], &context->viscousStrain);
        tensorOps.fromVector(&a[aOff[i_devStress]], &context->devStress);
    } // createContext

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

        deviatoricStress(*context, strain, stress);
    } // cauchyStress

    // --------------------------------------------------------------------------------------------
    /** Calculate Cauchy stress WITH reference stress and strain.
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

        const PylithReal bulkModulus = context->bulkModulus;
        const pylith::fekernels::Tensor& refStress = context->refStress;
        const pylith::fekernels::Tensor& refStrain = context->refStrain;
        pylith::fekernels::IsotropicLinearElasticity::meanStress_refState(bulkModulus, refStress, refStrain, strain, stress);

        deviatoricStress_refState(*context, strain, stress);
    } // cauchyStress_refState

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITHOUT reference stress and strain.
     */
    static inline
    void deviatoricStress(const Context& context,
                          const pylith::fekernels::Tensor& strain,
                          pylith::fekernels::Tensor* stress) {
        assert(stress);

        const PylithReal shearModulus = context.shearModulus;
        const PylithReal dt = context.dt;
        const PylithReal powerLawRefStrainRate = context.powerLawRefStrainRate;
        const PylithReal powerLawRefStress = context.powerLawRefStress;
        const PylithReal powerLawExponent = context.powerLawExponent;
        const pylith::fekernels::Tensor& viscousStrain = context.viscousStrain;
        const pylith::fekernels::Tensor& devStress = context.devStress;

        const PylithReal ae = 1.0 / (2.0 * shearModulus);
        const PylithReal timeFac = dt * (1.0 - powerLawAlpha);

        const PylithScalar devStressScalarProd = pylith::fekernels::TensorOps::scalarProduct(devStress, devStress);
        const PylithScalar j2T = sqrt(0.5 * devStressScalarProd);

        pylith::fekernels::Tensor devStrain;
        pylith::fekernels::Elasticity::deviatoric(strain, &devStrain);

        pylith::fekernels::Tensor strainPP;
        strainPP.xx = devStrain.xx - viscousStrain.xx;
        strainPP.yy = devStrain.yy - viscousStrain.yy;
        strainPP.zz = devStrain.zz - viscousStrain.zz;
        strainPP.xy = devStrain.xy - viscousStrain.xy;
        strainPP.yz = devStrain.yz - viscousStrain.yz;
        strainPP.xz = devStrain.xz - viscousStrain.xz;

        const PylithReal strainPPInvar2 = 0.5 * pylith::fekernels::TensorOps::scalarProduct(strainPP, strainPP);
        const PylithReal strainStressInvar2T = pylith::fekernels::TensorOps::scalarProduct(strainPP, devStress);

        // Finish defining parameters needed for root-finding algorithm.
        const PylithReal b = strainPPInvar2;
        const PylithReal c = strainStressInvar2T * timeFac;
        const PylithReal d = timeFac * j2T;
        PylithReal j2Tpdt = 0.0;
        if ((b != 0.0) || (c != 0.0) || (d != 0.0)) {
            const PylithReal j2InitialGuess = j2T;
            const PylithReal stressScale = shearModulus;
            j2Tpdt = _effectiveStress(j2InitialGuess, stressScale, ae, b, c, d, dt, j2T,
                                      powerLawExponent, powerLawRefStrainRate,
                                      powerLawRefStress);
        } // if
        // Compute deviatoric stresses from effective stress.
        const PylithReal j2Tau = (1.0 - powerLawAlpha) * j2T + powerLawAlpha * j2Tpdt;
        const PylithReal gammaTau = powerLawRefStrainRate * pow((j2Tau / powerLawRefStress), (powerLawExponent - 1.0)) / powerLawRefStress;
        const PylithReal factor1 = 1.0 / (ae + powerLawAlpha * dt * gammaTau);
        const PylithReal factor2 = timeFac * gammaTau;

        stress->xx += factor1 * (strainPP.xx - factor2 * devStress.xx);
        stress->yy += factor1 * (strainPP.yy - factor2 * devStress.yy);
        stress->zz += factor1 * (strainPP.zz - factor2 * devStress.zz);
        stress->xy += factor1 * (strainPP.xy - factor2 * devStress.xy);
        stress->yz += factor1 * (strainPP.yz - factor2 * devStress.yz);
        stress->xz += factor1 * (strainPP.xz - factor2 * devStress.xz);
    } // deviatoricStress

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITH reference stress and strain.
     */
    static inline
    void deviatoricStress_refState(const Context& context,
                                   const pylith::fekernels::Tensor& strain,
                                   pylith::fekernels::Tensor* stress) {
        assert(stress);

        const PylithReal shearModulus = context.shearModulus;
        const PylithReal dt = context.dt;
        const PylithReal powerLawRefStrainRate = context.powerLawRefStrainRate;
        const PylithReal powerLawRefStress = context.powerLawRefStress;
        const PylithReal powerLawExponent = context.powerLawExponent;
        const pylith::fekernels::Tensor& refStress = context.refStress;
        const pylith::fekernels::Tensor& refStrain = context.refStrain;
        const pylith::fekernels::Tensor& viscousStrain = context.viscousStrain;
        const pylith::fekernels::Tensor& devStress = context.devStress;

        const PylithReal ae = 1.0 / (2.0 * shearModulus);
        const PylithReal timeFac = dt * (1.0 - powerLawAlpha);

        pylith::fekernels::Tensor devRefStress;
        pylith::fekernels::Elasticity::deviatoric(refStress, &devRefStress);
        const PylithReal j2RefStressSquared = 0.5 * pylith::fekernels::TensorOps::scalarProduct(devRefStress, devRefStress);

        const PylithScalar devStressScalarProd = pylith::fekernels::TensorOps::scalarProduct(devStress, devStress);
        const PylithScalar j2T = sqrt(0.5 * devStressScalarProd);

        pylith::fekernels::Tensor devStrain;
        pylith::fekernels::Elasticity::deviatoric(strain, &devStrain);

        pylith::fekernels::Tensor strainPP;
        strainPP.xx = devStrain.xx - viscousStrain.xx - refStrain.xx;
        strainPP.yy = devStrain.yy - viscousStrain.yy - refStrain.yy;
        strainPP.zz = devStrain.zz - viscousStrain.zz - refStrain.zz;
        strainPP.xy = devStrain.xy - viscousStrain.xy - refStrain.xy;
        strainPP.yz = devStrain.yz - viscousStrain.yz - refStrain.yz;
        strainPP.xz = devStrain.xz - viscousStrain.xz - refStrain.xz;

        const PylithReal strainPPInvar2 = 0.5 * pylith::fekernels::TensorOps::scalarProduct(strainPP, strainPP);
        const PylithReal strainStressInvar2T = pylith::fekernels::TensorOps::scalarProduct(strainPP, devStress);
        const PylithReal strainRefStressInvar2T = pylith::fekernels::TensorOps::scalarProduct(strainPP, devRefStress);
        const PylithReal stressRefStressInvar2T = pylith::fekernels::TensorOps::scalarProduct(devStress, devRefStress);

        // Finish defining parameters needed for root-finding algorithm.
        const PylithReal b = strainPPInvar2 + ae * strainRefStressInvar2T + ae * ae * j2RefStressSquared;
        const PylithReal c = strainStressInvar2T * timeFac + ae * stressRefStressInvar2T;
        const PylithReal d = timeFac * j2T;
        PylithReal j2Tpdt = 0.0;
        if ((b != 0.0) || (c != 0.0) || (d != 0.0)) {
            const PylithReal j2InitialGuess = j2T;
            const PylithReal stressScale = shearModulus;
            j2Tpdt = _effectiveStress(j2InitialGuess, stressScale, ae, b, c, d, dt, j2T,
                                      powerLawExponent, powerLawRefStrainRate,
                                      powerLawRefStress);
        } // if
        // Compute deviatoric stresses from effective stress.
        const PylithReal j2Tau = (1.0 - powerLawAlpha) * j2T + powerLawAlpha * j2Tpdt;
        const PylithReal gammaTau = powerLawRefStrainRate * pow((j2Tau / powerLawRefStress), (powerLawExponent - 1.0)) / powerLawRefStress;
        const PylithReal factor1 = 1.0 / (ae + powerLawAlpha * dt * gammaTau);
        const PylithReal factor2 = timeFac * gammaTau;

        stress->xx += factor1 * (strainPP.xx - factor2 * devStress.xx + ae * devRefStress.xx);
        stress->yy += factor1 * (strainPP.yy - factor2 * devStress.yy + ae * devRefStress.yy);
        stress->zz += factor1 * (strainPP.zz - factor2 * devStress.zz + ae * devRefStress.zz);
        stress->xy += factor1 * (strainPP.xy - factor2 * devStress.xy + ae * devRefStress.xy);
        stress->yz += factor1 * (strainPP.yz - factor2 * devStress.yz + ae * devRefStress.yz);
        stress->xz += factor1 * (strainPP.xz - factor2 * devStress.xz + ae * devRefStress.xz);
    }

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
        IsotropicLinearElasticity::meanStress(bulkModulus, strain, stress);

        const pylith::fekernels::Tensor& devStress = context->devStress;
        stress->xx += devStress.xx;
        stress->yy += devStress.yy;
        stress->zz += devStress.zz;
        stress->xy += devStress.xy;
        stress->yz += devStress.yz;
        stress->xz += devStress.xz;
    } // cauchyStress_stateVars

    // --------------------------------------------------------------------------------------------
    /** Calculate Cauchy stress WITH a reference stress/strain using current state variables.
     *
     * ISA pylith::fekernels::Elasticity::stressfn_type
     *
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                   maxwell_time(1), viscous_strain(4), total_strain(4)]
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

        const pylith::fekernels::Tensor& devStress = context->devStress;

        pylith::fekernels::Tensor devRefStress;
        pylith::fekernels::Elasticity::deviatoric(refStress, &devRefStress);

        stress->xx += devRefStress.xx + devStress.xx;
        stress->yy += devRefStress.yy + devStress.yy;
        stress->zz += devRefStress.zz + devStress.zz;
        stress->xy += devRefStress.xy + devStress.xy;
        stress->yz += devRefStress.yz + devStress.yz;
        stress->xz += devRefStress.xz + devStress.xz;
    } // cauchyStress_refState_stateVars

    // --------------------------------------------------------------------------------------------
    /** Calculate viscous strain as a vector.
     *
     * Used to update viscous strain state variable.
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

        pylith::fekernels::Tensor devStressTpdt;
        deviatoricStress(rheologyContext, strain, &devStressTpdt);

        pylith::fekernels::Tensor viscousStrainTensor;
        _viscousStrain(rheologyContext, devStressTpdt, &viscousStrainTensor);

        tensorOps.toVector(viscousStrainTensor, viscousStrainVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate viscous strain WITH reference stress and strain as a vector.
     *
     * Used to update the viscous strain state variable.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void viscousStrain_refState_asVector(const pylith::fekernels::Elasticity::StrainContext& strainContext,
                                         const Context& rheologyContext,
                                         pylith::fekernels::Elasticity::strainfn_type strainFn,
                                         const pylith::fekernels::TensorOps& tensorOps,
                                         PylithScalar viscousStrainVector[]) {
        assert(viscousStrainVector);

        Tensor strain;
        strainFn(strainContext, &strain);

        pylith::fekernels::Tensor devStressTpdt;
        deviatoricStress_refState(rheologyContext, strain, &devStressTpdt);

        pylith::fekernels::Tensor viscousStrainTensor;
        _viscousStrain(rheologyContext, devStressTpdt, &viscousStrainTensor);

        tensorOps.toVector(viscousStrainTensor, viscousStrainVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress as a vector.
     *
     * Used to update deviatoric stress state variable.
     *
     */
    static inline
    void deviatoricStress_asVector(const pylith::fekernels::Elasticity::StrainContext& strainContext,
                                   const Context& rheologyContext,
                                   pylith::fekernels::Elasticity::strainfn_type strainFn,
                                   const pylith::fekernels::TensorOps& tensorOps,
                                   PylithScalar devStressVector[]) {
        assert(devStressVector);

        Tensor strain;
        strainFn(strainContext, &strain);

        pylith::fekernels::Tensor devStressTensor;
        deviatoricStress(rheologyContext, strain, &devStressTensor);

        tensorOps.toVector(devStressTensor, devStressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITH reference stress and strain as a vector.
     *
     * Used to update the deviatoric stress state variable.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void deviatoricStress_refState_asVector(const pylith::fekernels::Elasticity::StrainContext& strainContext,
                                            const Context& rheologyContext,
                                            pylith::fekernels::Elasticity::strainfn_type strainFn,
                                            const pylith::fekernels::TensorOps& tensorOps,
                                            PylithScalar devStressVector[]) {
        assert(devStressVector);

        Tensor strain;
        strainFn(strainContext, &strain);

        pylith::fekernels::Tensor devStressTensor;
        deviatoricStress_refState(rheologyContext, strain, &devStressTensor);

        tensorOps.toVector(devStressTensor, devStressVector);
    }

private:

    // --------------------------------------------------------------------------------------------
    static inline
    void _viscousStrain(const Context& context,
                        const pylith::fekernels::Tensor& devStressTpdt,
                        pylith::fekernels::Tensor* viscousStrain) {
        assert(viscousStrain);

        const PylithReal powerLawRefStrainRate = context.powerLawRefStrainRate;
        const PylithReal powerLawRefStress = context.powerLawRefStress;
        const PylithReal powerLawExponent = context.powerLawExponent;
        const PylithReal dt = context.dt;
        const pylith::fekernels::Tensor& viscousStrainT = context.viscousStrain;
        const pylith::fekernels::Tensor& devStressT = context.devStress;

        const PylithReal devStressTpdtScalarProd = pylith::fekernels::TensorOps::scalarProduct(devStressTpdt, devStressTpdt);
        const PylithReal j2Tpdt = sqrt(0.5 * devStressTpdtScalarProd);

        // Compute stress quantities at time T.
        const PylithScalar devStressTScalarProd = pylith::fekernels::TensorOps::scalarProduct(devStressT, devStressT);
        const PylithScalar j2T = sqrt(0.5 * devStressTScalarProd);

        // Compute quantities at intermediate time.
        pylith::fekernels::Tensor devStressTau;
        devStressTau.xx = (1.0 - powerLawAlpha) * devStressT.xx + powerLawAlpha * devStressTpdt.xx;
        devStressTau.yy = (1.0 - powerLawAlpha) * devStressT.yy + powerLawAlpha * devStressTpdt.yy;
        devStressTau.zz = (1.0 - powerLawAlpha) * devStressT.zz + powerLawAlpha * devStressTpdt.zz;
        devStressTau.xy = (1.0 - powerLawAlpha) * devStressT.xy + powerLawAlpha * devStressTpdt.xy;
        devStressTau.yz = (1.0 - powerLawAlpha) * devStressT.yz + powerLawAlpha * devStressTpdt.yz;
        devStressTau.xz = (1.0 - powerLawAlpha) * devStressT.xz + powerLawAlpha * devStressTpdt.xz;

        const PylithScalar j2Tau = (1.0 - powerLawAlpha) * j2T + powerLawAlpha * j2Tpdt;
        const PylithScalar gammaTau = powerLawRefStrainRate * pow((j2Tau / powerLawRefStress), (powerLawExponent - 1.0)) / powerLawRefStress;

        viscousStrain->xx = viscousStrainT.xx + dt * gammaTau * devStressTau.xx;
        viscousStrain->yy = viscousStrainT.yy + dt * gammaTau * devStressTau.yy;
        viscousStrain->zz = viscousStrainT.zz + dt * gammaTau * devStressTau.zz;
        viscousStrain->xy = viscousStrainT.xy + dt * gammaTau * devStressTau.xy;
        viscousStrain->yz = viscousStrainT.yz + dt * gammaTau * devStressTau.yz;
        viscousStrain->xz = viscousStrainT.xz + dt * gammaTau * devStressTau.xz;
    }

    // --------------------------------------------------------------------------------------------
    /** Compute effective stress for power-law material, given an initial guess and the current parameters.
     *
     * Used to compute stress and viscous strain.
     *
     */
    static inline
    PylithReal _effectiveStress(const PylithScalar j2InitialGuess,
                                const PylithScalar stressScale,
                                const PylithScalar ae,
                                const PylithScalar b,
                                const PylithScalar c,
                                const PylithScalar d,
                                const PylithScalar dt,
                                const PylithScalar j2T,
                                const PylithScalar powerLawExponent,
                                const PylithScalar powerLawRefStrainRate,
                                const PylithScalar powerLawRefStress) {
        assert(j2InitialGuess >= 0.0);
        // If initial guess is too low, use stress scale instead.
        const PylithReal xMin = 1.0e-10;

        // Bracket the root.
        PylithReal xL = 0.0;
        PylithReal xR = 0.0;
        if (j2InitialGuess > xMin) {
            xL = 0.5 * j2InitialGuess;
            xR = 1.5 * j2InitialGuess;
        } else {
            xL = 0.5 * stressScale;
            xR = 1.5 * stressScale;
        } // else

        _bracket(&xL, &xR, ae, b, c, d, dt, j2T, powerLawExponent, powerLawRefStrainRate,
                 powerLawRefStress);

        // Find effective stress using Newton's method with bisection.
        PylithReal effStress = _search(xL, xR, ae, b, c, d, dt, j2T, powerLawExponent,
                                       powerLawRefStrainRate, powerLawRefStress);

        return effStress;
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate effective stress function for a power-law viscoelastic material.
     *
     * Used for bracketing.
     *
     */
    static inline
    PylithScalar _effectiveStressFn(const PylithReal j2Tpdt,
                                    const PylithReal ae,
                                    const PylithReal b,
                                    const PylithReal c,
                                    const PylithReal d,
                                    const PylithReal dt,
                                    const PylithReal j2T,
                                    const PylithReal powerLawExponent,
                                    const PylithReal powerLawRefStrainRate,
                                    const PylithReal powerLawRefStress) {
        const PylithReal factor1 = 1.0 - powerLawAlpha;
        const PylithReal j2Tau = factor1 * j2T + powerLawAlpha * j2Tpdt;
        const PylithReal gammaTau = powerLawRefStrainRate * pow((j2Tau / powerLawRefStress), (powerLawExponent - 1.0)) /
                                    powerLawRefStress;
        const PylithReal a = ae + powerLawAlpha * dt * gammaTau;
        const PylithReal y = a * a * j2Tpdt * j2Tpdt - b + c * gammaTau - d * d * gammaTau * gammaTau;

        return y;
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate effective stress function and its derivative for a power-law viscoelastic material.
     *
     * Used for root-finding.
     *
     */
    static inline
    void _effectiveStressFnDerivative(PylithReal *func,
                                      PylithReal *dfunc,
                                      const PylithReal j2Tpdt,
                                      const PylithReal ae,
                                      const PylithReal b,
                                      const PylithReal c,
                                      const PylithReal d,
                                      const PylithReal dt,
                                      const PylithReal j2T,
                                      const PylithReal powerLawExponent,
                                      const PylithReal powerLawRefStrainRate,
                                      const PylithReal powerLawRefStress) {
        PylithReal y = *func;
        PylithReal dy = *dfunc;

        const PylithReal factor1 = 1.0 - powerLawAlpha;
        const PylithReal j2Tau = factor1 * j2T + powerLawAlpha * j2Tpdt;
        const PylithReal gammaTau = powerLawRefStrainRate * pow((j2Tau / powerLawRefStress), (powerLawExponent - 1.0)) / powerLawRefStress;
        const PylithReal dGammaTau = powerLawRefStrainRate * powerLawAlpha * (powerLawExponent - 1.0) * pow((j2Tau / powerLawRefStress), (powerLawExponent - 2.0)) / (powerLawRefStress * powerLawRefStress);
        const PylithReal a = ae + powerLawAlpha * dt * gammaTau;
        y = a * a * j2Tpdt * j2Tpdt - b + c * gammaTau - d * d * gammaTau * gammaTau;
        dy = 2.0 * a * a * j2Tpdt + dGammaTau * (2.0 * a * powerLawAlpha * dt * j2Tpdt * j2Tpdt + c - 2.0 * d * d * gammaTau);

        *func = y;
        *dfunc = dy;
    }

    // --------------------------------------------------------------------------------------------
    /** Bracket effective stress root.
     *
     * Used to place bounds on effective stress.
     *
     */
    static inline
    void _bracket(PylithReal *px1,
                  PylithReal *px2,
                  const PylithReal ae,
                  const PylithReal b,
                  const PylithReal c,
                  const PylithReal d,
                  const PylithReal dt,
                  const PylithReal j2T,
                  const PylithReal powerLawExponent,
                  const PylithReal powerLawRefStrainRate,
                  const PylithReal powerLawRefStress) {
        const size_t maxIterations = 50;

        // Arbitrary factor by which to increase the brackets.
        const PylithReal bracketFactor = 2;
        // Minimum allowed value for effective stress.
        const PylithReal xMin = 0.0;
        PylithReal x1 = *px1;
        PylithReal x2 = *px2;

        PylithReal funcValue1 = _effectiveStressFn(x1, ae, b, c, d, dt, j2T, powerLawExponent, powerLawRefStrainRate, powerLawRefStress);
        PylithReal funcValue2 = _effectiveStressFn(x2, ae, b, c, d, dt, j2T, powerLawExponent, powerLawRefStrainRate, powerLawRefStress);

        bool bracketed = false;
        for (size_t i = 0; i < maxIterations; ++i) {
            if ((funcValue1 * funcValue2) < 0.0) {
                bracketed = true;
                break;
            } // if

            if (fabs(funcValue1) < fabs(funcValue2)) {
                x1 += bracketFactor * (x1 - x2);
                x1 = std::max(x1, xMin);
                funcValue1 = _effectiveStressFn(x1, ae, b, c, d, dt, j2T, powerLawExponent, powerLawRefStrainRate, powerLawRefStress);
            } else {
                x2 += bracketFactor * (x1 - x2);
                x2 = std::max(x2, xMin);
                funcValue2 = _effectiveStressFn(x2, ae, b, c, d, dt, j2T, powerLawExponent, powerLawRefStrainRate, powerLawRefStress);
            } // else
        } // for

        *px1 = x1;
        *px2 = x2;

        if (!bracketed) {
            throw std::runtime_error("Unable to bracket effective stress.");
        } // if
    }

    // --------------------------------------------------------------------------------------------
    /** Find zero of effective stress function using Newton's method with bisection.
     *
     * Used to find the effective stress.
     *
     */
    static inline
    PylithReal _search(PylithReal x1,
                       PylithReal x2,
                       const PylithReal ae,
                       const PylithReal b,
                       const PylithReal c,
                       const PylithReal d,
                       const PylithReal dt,
                       const PylithReal j2T,
                       const PylithReal powerLawExponent,
                       const PylithReal powerLawRefStrainRate,
                       const PylithReal powerLawRefStress) {
        const size_t maxIterations = 100;

        // Desired accuracy for root. This is a bit arbitrary for now.
        const PylithReal accuracy = 1.0e-16;

        // Organize search so that _effectiveStressFn(xLow) is less than zero.
        PylithReal funcValueLow = _effectiveStressFn(x1, ae, b, c, d, dt, j2T, powerLawExponent, powerLawRefStrainRate, powerLawRefStress);
        PylithReal funcValueHigh = _effectiveStressFn(x2, ae, b, c, d, dt, j2T, powerLawExponent, powerLawRefStrainRate, powerLawRefStress);
        assert(funcValueLow * funcValueHigh <= 0.0);

        PylithReal effStress = 0.0;
        PylithReal xLow = 0.0;
        PylithReal xHigh = 0.0;
        bool converged = false;

        if (funcValueLow < 0.0) {
            xLow = x1;
            xHigh = x2;
        } else {
            xLow = x2;
            xHigh = x1;
        } // if/else

        effStress = 0.5 * (x1 + x2);
        PylithReal dxPrevious = fabs(x2 - x1);
        PylithReal dx = dxPrevious;
        PylithReal funcValue = 0.0;
        PylithReal funcDeriv = 0.0;
        PylithReal funcXHigh = 0.0;
        PylithReal funcXLow = 0.0;
        _effectiveStressFnDerivative(&funcValue, &funcDeriv, effStress, ae, b, c, d, dt, j2T, powerLawExponent, powerLawRefStrainRate, powerLawRefStress);

        for (size_t i = 0; i < maxIterations; ++i) {
            funcXHigh = (effStress - xHigh) * funcDeriv - funcValue;
            funcXLow = (effStress - xLow) * funcDeriv - funcValue;
            if (fabs(funcValue) < accuracy) {
                converged = true;
                break;
            } // if
            // Use bisection if solution goes out of bounds.
            if (funcXHigh * funcXLow >= 0.0) {
                dx = 0.5 * (xHigh - xLow);
                effStress = xLow + dx;
            } else {
                dxPrevious = dx;
                dx = funcValue / funcDeriv;
                effStress = effStress - dx;
            } // else
            _effectiveStressFnDerivative(&funcValue, &funcDeriv, effStress, ae, b, c, d, dt, j2T, powerLawExponent, powerLawRefStrainRate, powerLawRefStress);
            if (funcValue < 0.0) {
                xLow = effStress;
            } else {
                xHigh = effStress;
            } // else
        } // for

        if (!converged) {
            throw std::runtime_error("Cannot find root of effective stress function.");
        } // if

        return effStress;
    }

}; // IsotropicPowerLaw

// ------------------------------------------------------------------------------------------------
/// Kernels for isotropic power-law plane strain.
class pylith::fekernels::IsotropicPowerLawPlaneStrain {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    struct ElasticConstants {
        PylithReal C1111;
        PylithReal C1122;
        PylithReal C1212;
        PylithReal C2211;
        PylithReal C2222;
    };

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    // ===========================================================================================
    // Kernels for elasticity equation
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for plane strain isotropic power-law with infinitesimal strain WITHOUT
     * reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     * power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), stress(4)]
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for plane strain isotropic power-law with infinitesimal strain WITH
     * reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), deviatoric_stress(4)]
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress_refState,
            pylith::fekernels::Tensor::ops2D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** Jf3_vu entry function for plane strain isotropic power-law viscoelasticity with
     * infinitesimal strain WITHOUT reference stress/strain.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     * power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), deviatoric_stress(4)]
     *
     * stress_ij = C_ijkl strain_kl
     *
     * Isotropic:
     *  C_ijkl = bulkModulus * delta_ij * delta_kl
     *   + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
     *
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

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Tensor strain;
        pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain(strainContext, &strain);

        const pylith::fekernels::Tensor& devStressT = rheologyContext.devStress;
        const PylithReal devStressScalarProd = pylith::fekernels::TensorOps::scalarProduct(devStressT, devStressT);
        const PylithReal j2T = sqrt(0.5 * devStressScalarProd);

        // Compute quantities based on stress at t = T + dt.
        pylith::fekernels::Tensor devStressTpdt;
        pylith::fekernels::IsotropicPowerLaw::deviatoricStress(rheologyContext, strain, &devStressTpdt);
        const PylithReal devStressTpdtScalarProd = pylith::fekernels::TensorOps::scalarProduct(devStressTpdt, devStressTpdt);
        const PylithReal j2Tpdt = sqrt(0.5 * devStressTpdtScalarProd);

        ElasticConstants elasticityMat;
        _elasticConstants(&elasticityMat, rheologyContext, devStressTpdt, j2T, j2Tpdt);

        /* Nonzero Jacobian entries. */
        Jf3[ 0] -= elasticityMat.C1111; /* j0000 */
        Jf3[ 3] -= elasticityMat.C1212; /* j0011 */
        Jf3[ 5] -= elasticityMat.C1122; /* j0101 */
        Jf3[ 6] -= elasticityMat.C1212; /* j0110 */
        Jf3[ 9] -= elasticityMat.C1212; /* j1001 */
        Jf3[10] -= elasticityMat.C2211; /* j1010 */
        Jf3[12] -= elasticityMat.C1212; /* j1100 */
        Jf3[15] -= elasticityMat.C2222; /* j1111 */
    }

    // --------------------------------------------------------------------------------------------
    /** Jf3_vu entry function for plane strain isotropic power-law viscoelasticity with
     * infinitesimal strain WITH reference stress/strain.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), deviatoric_stress(4)]
     */
    static inline
    void Jf3vu_infinitesimalStrain_refState(const PylithInt dim,
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

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Tensor strain;
        pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain(strainContext, &strain);

        const pylith::fekernels::Tensor& devStressT = rheologyContext.devStress;
        const PylithReal devStressScalarProd = pylith::fekernels::TensorOps::scalarProduct(devStressT, devStressT);
        const PylithReal j2T = sqrt(0.5 * devStressScalarProd);

        // Compute quantities based on stress at t = T + dt.
        pylith::fekernels::Tensor devStressTpdt;
        pylith::fekernels::IsotropicPowerLaw::deviatoricStress_refState(rheologyContext, strain, &devStressTpdt);
        const PylithReal devStressTpdtScalarProd = pylith::fekernels::TensorOps::scalarProduct(devStressTpdt, devStressTpdt);
        const PylithReal j2Tpdt = sqrt(0.5 * devStressTpdtScalarProd);

        ElasticConstants elasticityMat;
        _elasticConstants(&elasticityMat, rheologyContext, devStressTpdt, j2T, j2Tpdt);

        /* Nonzero Jacobian entries. */
        Jf3[ 0] -= elasticityMat.C1111; /* j0000 */
        Jf3[ 3] -= elasticityMat.C1212; /* j0011 */
        Jf3[ 5] -= elasticityMat.C1122; /* j0101 */
        Jf3[ 6] -= elasticityMat.C1212; /* j0110 */
        Jf3[ 9] -= elasticityMat.C1212; /* j1001 */
        Jf3[10] -= elasticityMat.C2211; /* j1010 */
        Jf3[12] -= elasticityMat.C1212; /* j1100 */
        Jf3[15] -= elasticityMat.C2222; /* j1111 */
    }

    // ===========================================================================================
    // Kernels for fault interfaces and elasticity
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for isotropic power-law viscoelasticity plane strain with
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for isotropic power-law viscoelasticity plane strain with
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for isotropic power-law viscoelasticity plane strain with
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress_refState,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for isotropic power-law viscoelasticity plane strain with
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress_refState,
            pylith::fekernels::ElasticityPlaneStrain::traction,
            pylith::fekernels::Tensor::ops2D,
            f0);
    }

    // ===========================================================================================
    // Kernels for updating state variables
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating viscous strain as a vector for plane strain isotropic
     * power-law viscoelasticity.
     *
     * Used to update the viscous strain state variable.
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicPowerLaw::viscousStrain_asVector(
            strainContext, rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::Tensor::ops2D,
            viscousStrain);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating viscous strain as a vector for plane strain isotropic
     * power-law viscoelasticity WITH reference stress and strain
     *
     * Used to update the viscous strain state variable.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void viscousStrain_infinitesimalStrain_refState_asVector(const PylithInt dim,
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicPowerLaw::viscousStrain_refState_asVector(
            strainContext, rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::Tensor::ops2D,
            viscousStrain);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating deviatoric stress as a vector for plane strain isotropic
     * power-law viscoelasticity.
     *
     * Used to update the deviatoric stress state variable.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void deviatoricStress_infinitesimalStrain_asVector(const PylithInt dim,
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
                                                       PylithScalar devStress[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicPowerLaw::deviatoricStress_asVector(
            strainContext, rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::Tensor::ops2D,
            devStress);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating deviatoric stress as a vector for plane strain isotropic
     * power-law viscoelasticity.
     *
     * Used to update the deviatoric stress state variable.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void deviatoricStress_infinitesimalStrain_refState_asVector(const PylithInt dim,
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
                                                                PylithScalar devStress[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicPowerLaw::deviatoricStress_refState_asVector(
            strainContext, rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::Tensor::ops2D,
            devStress);
    }

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for plane strain isotropic power-law
     * viscoelasticity with infinitesimal strain WITHOUT a reference stress and strain.
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress_stateVars,
            pylith::fekernels::Tensor::ops2D,
            stressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for plane strain isotropic power-law
     * viscoelasticity with infinitesimal strain WITH a reference stress and strain.
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress_refState_stateVars,
            pylith::fekernels::Tensor::ops2D,
            stressVector);
    }

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    // --------------------------------------------------------------------------------------------
    /* j(f,g,df,dg) = C(f,df,g,dg)
     *
     * 0:  j0000 = C1111 = bulkModulus + 2/(3*(alpha**2*deltaT*gammaFTau*s11**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s11*s11T*(1 - alpha)*(n -
     * 1)/(2*j2FTau*j2FTplusDt) + 1/(2*shearModulus)))
     * 1:  j0001 = C1112 = 0
     * 2:  j0010 = C1211 = 0
     * 3:  j0011 = C1212 = 1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s12*s12T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 4:  j0100 = C1121 = 0
     * 5:  j0101 = C1122 = bulkModulus - 1/(3*(alpha**2*deltaT*gammaFTau*s11**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s11*s11T*(1 - alpha)*(n -
     * 1)/(2*j2FTau*j2FTplusDt) + 1/(2*shearModulus)))
     * 6:  j0110 = C1221 = 1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s12*s12T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 7:  j0111 = C1222 = 0
     * 8:  j1000 = C2111 = 0
     * 9:  j1001 = C2112 = 1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s12*s12T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 10: j1010 = C2211 = bulkModulus - 1/(3*(alpha**2*deltaT*gammaFTau*s22**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s22*s22T*(1 - alpha)*(n -
     * 1)/(2*j2FTau*j2FTplusDt) + 1/(2*shearModulus)))
     * 11: j1011 = C2212 = 0
     * 12: j1100 = C2121 = 1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s12*s12T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 13: j1101 = C2122 = 0
     * 14: j1110 = C2221 = 0
     * 15: j1111 = C2222 = bulkModulus + 2/(3*(alpha**2*deltaT*gammaFTau*s22**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s22*s22T*(1 - alpha)*(n -
     * 1)/(2*j2FTau*j2FTplusDt) + 1/(2*shearModulus)))
     */
    static inline
    void _elasticConstants(ElasticConstants* elasticityMat,
                           pylith::fekernels::IsotropicPowerLaw::Context& context,
                           const pylith::fekernels::Tensor& devStressTpdt,
                           const PylithReal j2T,
                           const PylithReal j2Tpdt) {
        assert(elasticityMat);

        // Compute quantities at intermediate time tau.
        const PylithReal powerLawAlpha = pylith::fekernels::IsotropicPowerLaw::powerLawAlpha;
        const PylithReal powerLawRefStrainRate = context.powerLawRefStrainRate;
        const PylithReal powerLawRefStress = context.powerLawRefStress;
        const PylithReal powerLawExponent = context.powerLawExponent;
        const PylithReal dt = context.dt;
        const pylith::fekernels::Tensor& devStressT = context.devStress;

        const PylithReal j2Tau = powerLawAlpha * j2Tpdt + (1.0 - powerLawAlpha) * j2T;
        const PylithReal gammaTau = powerLawRefStrainRate * pow((j2Tau / powerLawRefStress), (powerLawExponent - 1.0)) / powerLawRefStress;

        const PylithReal bulkModulus = context.bulkModulus;
        const PylithReal shearModulus = context.shearModulus;

        if ((j2Tpdt == 0.0) && (j2Tau == 0.0)) {
            // Elastic Jacobian if effective stress is zero.
            elasticityMat->C1111 = bulkModulus + 4.0 * shearModulus / 3.0;
            elasticityMat->C1122 = bulkModulus - 2.0 * shearModulus / 3.0;
            elasticityMat->C2211 = elasticityMat->C1122;
            elasticityMat->C1212 = shearModulus;
            elasticityMat->C2222 = elasticityMat->C1111;
        } else {
            // Viscoelastic Jacobian if effective stress is nonzero.
            const PylithReal ae = 1.0 / (2.0 * shearModulus);
            const PylithReal denom = 2.0 * j2Tau * j2Tpdt;
            const PylithReal factor1 = powerLawAlpha * dt * gammaTau;
            const PylithReal factor2 = factor1 * (powerLawExponent - 1.0) / denom;
            const PylithReal factor3 = powerLawAlpha * factor2;
            const PylithReal factor4 = factor2 * (1.0 - powerLawAlpha);

            /* Unique components of Jacobian. */
            elasticityMat->C1111 = bulkModulus + 2 / (3 * (factor3 * devStressTpdt.xx * devStressTpdt.xx + factor1 +
                                                           factor4 * devStressTpdt.xx * devStressT.xx + ae));
            elasticityMat->C1122 = bulkModulus - 1 / (3 * (factor3 * devStressTpdt.xx * devStressTpdt.xx + factor1 +
                                                           factor4 * devStressTpdt.xx * devStressT.xx + ae));
            elasticityMat->C1212 = 1 / (2 * (factor3 * devStressTpdt.xy * devStressTpdt.xy + factor1 +
                                             factor4 * devStressTpdt.xy * devStressT.xy + ae));
            elasticityMat->C2211 = bulkModulus - 1 / (3 * (factor3 * devStressTpdt.yy * devStressTpdt.yy + factor1 +
                                                           factor4 * devStressTpdt.yy * devStressT.yy + ae));
            elasticityMat->C2222 = bulkModulus + 2 / (3 * (factor3 * devStressTpdt.yy * devStressTpdt.yy + factor1 +
                                                           factor4 * devStressTpdt.yy * devStressT.yy + ae));
        } // if

    }

}; // IsotropicPowerLawPlaneStrain

// ------------------------------------------------------------------------------------------------
/// Kernels for 3D isotropic power-law viscoelasticity.
class pylith::fekernels::IsotropicPowerLaw3D {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    struct ElasticConstants {
        PylithReal C1111;
        PylithReal C1122;
        PylithReal C1212;
        PylithReal C1313;
        PylithReal C2211;
        PylithReal C2222;
        PylithReal C2323;
        PylithReal C3311;
        PylithReal C3333;
    };

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    // ===========================================================================================
    // Kernels for elasticity equation
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for 3D  isotropic power-law viscoelasticity with infinitesimal strain
     * WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     * power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), stress(4)]
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** f1 entry function for 3D isotropic power-law viscoelasticity with infinitesimal strain WITH
     * reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), deviatoric_stress(4)]
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress_refState,
            pylith::fekernels::Tensor::ops3D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** Jf3_vu entry function for 3D isotropic power-law viscoelasticity with infinitesimal strain
     * WITHOUT reference stress/strain.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     * power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), deviatoric_stress(4)]
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

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Tensor strain;
        pylith::fekernels::Elasticity3D::infinitesimalStrain(strainContext, &strain);

        const pylith::fekernels::Tensor& devStressT = rheologyContext.devStress;
        const PylithReal devStressScalarProd = pylith::fekernels::TensorOps::scalarProduct(devStressT, devStressT);
        const PylithReal j2T = sqrt(0.5 * devStressScalarProd);

        // Compute quantities based on stress at t = T + dt.
        pylith::fekernels::Tensor devStressTpdt;
        pylith::fekernels::IsotropicPowerLaw::deviatoricStress(rheologyContext, strain, &devStressTpdt);
        const PylithReal devStressTpdtScalarProd = pylith::fekernels::TensorOps::scalarProduct(devStressTpdt, devStressTpdt);
        const PylithReal j2Tpdt = sqrt(0.5 * devStressTpdtScalarProd);

        ElasticConstants elasticityMat;
        _elasticConstants(&elasticityMat, rheologyContext, devStressTpdt, j2T, j2Tpdt);

        /* Nonzero Jacobian entries. */
        Jf3[ 0] -= elasticityMat.C1111;  /* j0000 */
        Jf3[ 4] -= elasticityMat.C1212;  /* j0011 */
        Jf3[ 8] -= elasticityMat.C1313;  /* j0022 */
        Jf3[10] -= elasticityMat.C1122; /* j0101 */
        Jf3[12] -= elasticityMat.C1212; /* j0110 */
        Jf3[20] -= elasticityMat.C1122; /* j0202 */
        Jf3[24] -= elasticityMat.C1313; /* j0220 */
        Jf3[28] -= elasticityMat.C1212; /* j1001 */
        Jf3[30] -= elasticityMat.C2211; /* j1010 */
        Jf3[36] -= elasticityMat.C1212; /* j1100 */
        Jf3[40] -= elasticityMat.C2222; /* j1111 */
        Jf3[44] -= elasticityMat.C2323; /* j1122 */
        Jf3[50] -= elasticityMat.C2211; /* j1212 */
        Jf3[52] -= elasticityMat.C2323; /* j1221 */
        Jf3[56] -= elasticityMat.C1313; /* j2002 */
        Jf3[60] -= elasticityMat.C3311; /* j2020 */
        Jf3[68] -= elasticityMat.C2323; /* j2112 */
        Jf3[70] -= elasticityMat.C3311; /* j2121 */
        Jf3[72] -= elasticityMat.C1313; /* j2200 */
        Jf3[76] -= elasticityMat.C2323; /* j2211 */
        Jf3[80] -= elasticityMat.C3333; /* j2222 */
    }

    // --------------------------------------------------------------------------------------------
    /** Jf3_vu entry function for 3D isotropic power-law viscoelasticity with infinitesimal strain
     * WITH reference stress/strain.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), deviatoric_stress(4)]
     */
    static inline
    void Jf3vu_infinitesimalStrain_refState(const PylithInt dim,
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

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Tensor strain;
        pylith::fekernels::Elasticity3D::infinitesimalStrain(strainContext, &strain);

        const pylith::fekernels::Tensor& devStressT = rheologyContext.devStress;
        const PylithReal devStressScalarProd = pylith::fekernels::TensorOps::scalarProduct(devStressT, devStressT);
        const PylithReal j2T = sqrt(0.5 * devStressScalarProd);

        // Compute quantities based on stress at t = T + dt.
        pylith::fekernels::Tensor devStressTpdt;
        pylith::fekernels::IsotropicPowerLaw::deviatoricStress_refState(rheologyContext, strain, &devStressTpdt);
        const PylithReal devStressTpdtScalarProd = pylith::fekernels::TensorOps::scalarProduct(devStressTpdt, devStressTpdt);
        const PylithReal j2Tpdt = sqrt(0.5 * devStressTpdtScalarProd);

        ElasticConstants elasticityMat;
        _elasticConstants(&elasticityMat, rheologyContext, devStressTpdt, j2T, j2Tpdt);

        /* Nonzero Jacobian entries. */
        Jf3[ 0] -= elasticityMat.C1111;  /* j0000 */
        Jf3[ 4] -= elasticityMat.C1212;  /* j0011 */
        Jf3[ 8] -= elasticityMat.C1313;  /* j0022 */
        Jf3[10] -= elasticityMat.C1122; /* j0101 */
        Jf3[12] -= elasticityMat.C1212; /* j0110 */
        Jf3[20] -= elasticityMat.C1122; /* j0202 */
        Jf3[24] -= elasticityMat.C1313; /* j0220 */
        Jf3[28] -= elasticityMat.C1212; /* j1001 */
        Jf3[30] -= elasticityMat.C2211; /* j1010 */
        Jf3[36] -= elasticityMat.C1212; /* j1100 */
        Jf3[40] -= elasticityMat.C2222; /* j1111 */
        Jf3[44] -= elasticityMat.C2323; /* j1122 */
        Jf3[50] -= elasticityMat.C2211; /* j1212 */
        Jf3[52] -= elasticityMat.C2323; /* j1221 */
        Jf3[56] -= elasticityMat.C1313; /* j2002 */
        Jf3[60] -= elasticityMat.C3311; /* j2020 */
        Jf3[68] -= elasticityMat.C2323; /* j2112 */
        Jf3[70] -= elasticityMat.C3311; /* j2121 */
        Jf3[72] -= elasticityMat.C1313; /* j2200 */
        Jf3[76] -= elasticityMat.C2323; /* j2211 */
        Jf3[80] -= elasticityMat.C3333; /* j2222 */
    }

    // ===========================================================================================
    // Kernels for fault interfaces and elasticity
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for 3D isotropic power-law viscoelasticity with
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for 3D isotropic power-law viscoelasticity with
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for negative fault face for 3D isotropic power-law viscoelasticity with
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_neg(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress_refState,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // --------------------------------------------------------------------------------------------
    /** f0 entry function for positive fault face for 3D isotropic power-law viscoelasticity with
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::FaultCohesiveKin::f0l_pos(
            dim, numS, sOff, s, n,
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress_refState,
            pylith::fekernels::Elasticity3D::traction,
            pylith::fekernels::Tensor::ops3D,
            f0);
    }

    // ===========================================================================================
    // Kernels for updating state variables
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating viscous strain as a vector for 3D isotropic power-law
     * viscoelasticity.
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
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicPowerLaw::viscousStrain_asVector(
            strainContext, rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::Tensor::ops3D,
            viscousStrain);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating viscous strain as a vector for 3D isotropic power-law
     * viscoelasticity WITH reference stress and strain.
     *
     * Used to update the viscous strain state variable.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void viscousStrain_infinitesimalStrain_refState_asVector(const PylithInt dim,
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

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicPowerLaw::viscousStrain_refState_asVector(
            strainContext, rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::Tensor::ops3D,
            viscousStrain);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating deviatoric stress as a vector for 3D isotropic
     * power-law viscoelasticity.
     *
     * Used to update the deviatoric stress state variable.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void deviatoricStress_infinitesimalStrain_asVector(const PylithInt dim,
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
                                                       PylithScalar devStress[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicPowerLaw::deviatoricStress_asVector(
            strainContext, rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::Tensor::ops3D,
            devStress);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating deviatoric stress as a vector for 3D isotropic
     * power-law viscoelasticity.
     *
     * Used to update the deviatoric stress state variable.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void deviatoricStress_infinitesimalStrain_refState_asVector(const PylithInt dim,
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
                                                                PylithScalar devStress[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicPowerLaw::deviatoricStress_refState_asVector(
            strainContext, rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::Tensor::ops3D,
            devStress);
    }

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic power-law viscoelasticity
     * with infinitesimal strain WITHOUT a reference stress and strain.
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
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress_stateVars,
            pylith::fekernels::Tensor::ops3D,
            stressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic power-law viscoelasticity
     * with infinitesimal strain WITH a reference stress and strain.
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
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicPowerLaw::Context rheologyContext;
        pylith::fekernels::IsotropicPowerLaw::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicPowerLaw::cauchyStress_refState_stateVars,
            pylith::fekernels::Tensor::ops3D,
            stressVector);
    }

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    /* j(f,g,df,dg) = C(f,df,g,dg)
     *
     * 0:  j0000 = C1111 = bulkModulus + 2/(3*(alpha**2*deltaT*gammaFTau*s11**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s11*s11T*(1 - alpha)*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 1:  j0001 = C1112 = 0
     * 2:  j0002 = C1113 = 0
     * 3:  j0010 = C1211 = 0
     * 4:  j0011 = C1212 = 1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s12*s12T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 5:  j0012 = C1213 = 0
     * 6:  j0020 = C1311 = 0
     * 7:  j0021 = C1312 = 0
     * 8:  j0022 = C1313 = 1/(2*(alpha**2*deltaT*gammaFTau*s13**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s13*s13T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 9:  j0100 = C1121 = 0
     * 10: j0101 = C1122 = bulkModulus - 1/(3*(alpha**2*deltaT*gammaFTau*s11**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s11*s11T*(1 - alpha)*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 11: j0102 = C1123 = 0
     * 12: j0110 = C1221 = 1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s12*s12T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 13: j0111 = C1222 = 0
     * 14: j0112 = C1223 = 0
     * 15: j0120 = C1321 = 0
     * 16: j0121 = C1322 = 0
     * 17: j0122 = C1323 = 0
     * 18: j0200 = C1131 = 0
     * 19: j0201 = C1132 = 0
     * 20: j0202 = C1133 = bulkModulus - 1/(3*(alpha**2*deltaT*gammaFTau*s11**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s11*s11T*(1 - alpha)*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 21: j0210 = C1231 = 0
     * 22: j0211 = C1232 = 0
     * 23: j0212 = C1233 = 0
     * 24: j0220 = C1331 = 1/(2*(alpha**2*deltaT*gammaFTau*s13**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s13*s13T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 25: j0221 = C1332 = 0
     * 26: j0222 = C1333 = 0
     * 27: j1000 = C2111 = 0
     * 28: j1001 = C2112 = 1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s12*s12T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 29: j1002 = C2113 = 0
     * 30: j1010 = C2211 = bulkModulus - 1/(3*(alpha**2*deltaT*gammaFTau*s22**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s22*s22T*(1 - alpha)*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 31: j1011 = C2212 = 0
     * 32: j1012 = C2213 = 0
     * 33: j1020 = C2311 = 0
     * 34: j1021 = C2312 = 0
     * 35: j1022 = C2313 = 0
     * 36: j1100 = C2121 = 1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s12*s12T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 37: j1101 = C2122 = 0
     * 38: j1102 = C2123 = 0
     * 39: j1110 = C2221 = 0
     * 40: j1111 = C2222 = bulkModulus + 2/(3*(alpha**2*deltaT*gammaFTau*s22**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s22*s22T*(1 - alpha)*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 41: j1112 = C2223 = 0
     * 42: j1120 = C2321 = 0
     * 43: j1121 = C2322 = 0
     * 44: j1122 = C2323 = 1/(2*(alpha**2*deltaT*gammaFTau*s23**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s23*s23T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 45: j1200 = C2131 = 0
     * 46: j1201 = C2132 = 0
     * 47: j1202 = C2133 = 0
     * 48: j1210 = C2231 = 0
     * 49: j1211 = C2232 = 0
     * 50: j1212 = C2233 = bulkModulus - 1/(3*(alpha**2*deltaT*gammaFTau*s22**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s22*s22T*(1 - alpha)*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 51: j1220 = C2331 = 0
     * 52: j1221 = C2332 = 1/(2*(alpha**2*deltaT*gammaFTau*s23**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s23*s23T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 53: j1222 = C2333 = 0
     * 54: j2000 = C3111 = 0
     * 55: j2001 = C3112 = 0
     * 56: j2002 = C3113 = 1/(2*(alpha**2*deltaT*gammaFTau*s13**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s13*s13T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 57: j2010 = C3211 = 0
     * 58: j2011 = C3212 = 0
     * 59: j2012 = C3213 = 0
     * 60: j2020 = C3311 = bulkModulus - 1/(3*(alpha**2*deltaT*gammaFTau*s33**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s33*s33T*(1 - alpha)*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 61: j2021 = C3312 = 0
     * 62: j2022 = C3313 = 0
     * 63: j2100 = C3121 = 0
     * 64: j2101 = C3122 = 0
     * 65: j2102 = C3123 = 0
     * 66: j2110 = C3221 = 0
     * 67: j2111 = C3222 = 0
     * 68: j2112 = C3223 = 1/(2*(alpha**2*deltaT*gammaFTau*s23**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s23*s23T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 69: j2120 = C3321 = 0
     * 70: j2121 = C3322 = bulkModulus - 1/(3*(alpha**2*deltaT*gammaFTau*s33**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s33*s33T*(1 - alpha)*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 71: j2122 = C3323 = 0
     * 72: j2200 = C3131 = 1/(2*(alpha**2*deltaT*gammaFTau*s13**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s13*s13T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 73: j2201 = C3132 = 0
     * 74: j2202 = C3133 = 0
     * 75: j2210 = C3231 = 0
     * 76: j2211 = C3232 = 1/(2*(alpha**2*deltaT*gammaFTau*s23**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau +
     *                        alpha*deltaT*gammaFTau*s23*s23T*(1 - alpha)*(n - 1)/(j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     * 77: j2212 = C3233 = 0
     * 78: j2220 = C3331 = 0
     * 79: j2221 = C3332 = 0
     * 80: j2222 = C3333 = bulkModulus + 2/(3*(alpha**2*deltaT*gammaFTau*s33**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau +
     *                                      alpha*deltaT*gammaFTau*s33*s33T*(1 - alpha)*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * 1/(2*shearModulus)))
     */
    // --------------------------------------------------------------------------------------------
    static inline
    void _elasticConstants(ElasticConstants* elasticityMat,
                           pylith::fekernels::IsotropicPowerLaw::Context& context,
                           const pylith::fekernels::Tensor& devStressTpdt,
                           const PylithReal j2T,
                           const PylithReal j2Tpdt) {
        assert(elasticityMat);

        // Compute quantities at intermediate time tau.
        const PylithReal powerLawAlpha = pylith::fekernels::IsotropicPowerLaw::powerLawAlpha;
        const PylithReal powerLawRefStrainRate = context.powerLawRefStrainRate;
        const PylithReal powerLawRefStress = context.powerLawRefStress;
        const PylithReal powerLawExponent = context.powerLawExponent;
        const PylithReal dt = context.dt;
        const pylith::fekernels::Tensor& devStressT = context.devStress;

        const PylithReal j2Tau = powerLawAlpha * j2Tpdt + (1.0 - powerLawAlpha) * j2T;
        const PylithReal gammaTau = powerLawRefStrainRate * pow((j2Tau / powerLawRefStress), (powerLawExponent - 1.0)) / powerLawRefStress;

        const PylithReal bulkModulus = context.bulkModulus;
        const PylithReal shearModulus = context.shearModulus;

        if ((j2Tpdt == 0.0) && (j2Tau == 0.0)) {
            // Elastic Jacobian if effective stress is zero.
            elasticityMat->C1111 = bulkModulus + 4.0 * shearModulus / 3.0;
            elasticityMat->C1122 = bulkModulus - 2.0 * shearModulus / 3.0;
            elasticityMat->C1212 = shearModulus;
            elasticityMat->C1313 = shearModulus;
            elasticityMat->C2211 = elasticityMat->C1122;
            elasticityMat->C2222 = elasticityMat->C1111;
            elasticityMat->C2323 = elasticityMat->C1212;
            elasticityMat->C3311 = elasticityMat->C1122;
            elasticityMat->C3333 = elasticityMat->C1111;
        } else {
            // Viscoelastic Jacobian if effective stress is nonzero.
            const PylithReal ae = 1.0 / (2.0 * shearModulus);
            const PylithReal denom = 2.0 * j2Tau * j2Tpdt;
            const PylithReal factor1 = powerLawAlpha * dt * gammaTau;
            const PylithReal factor2 = factor1 * (powerLawExponent - 1.0) / denom;
            const PylithReal factor3 = powerLawAlpha * factor2;
            const PylithReal factor4 = factor2 * (1.0 - powerLawAlpha);

            /* Unique components of Jacobian. */
            elasticityMat->C1111 = bulkModulus + 2 / (3 * (factor3 * devStressTpdt.xx * devStressTpdt.xx + factor1 +
                                                           factor4 * devStressTpdt.xx * devStressT.xx + ae));
            elasticityMat->C1122 = bulkModulus - 1 / (3 * (factor3 * devStressTpdt.xx * devStressTpdt.xx + factor1 +
                                                           factor4 * devStressTpdt.xx * devStressT.xx + ae));
            elasticityMat->C1212 = 1 / (2 * (factor3 * devStressTpdt.xy * devStressTpdt.xy + factor1 +
                                             factor4 * devStressTpdt.xy * devStressT.xy + ae));
            elasticityMat->C1313 = 1 / (2 * (factor3 * devStressTpdt.xz * devStressTpdt.xz + factor1 +
                                             factor4 * devStressTpdt.xz * devStressT.xz + ae));
            elasticityMat->C2211 = bulkModulus - 1 / (3 * (factor3 * devStressTpdt.yy * devStressTpdt.yy + factor1 +
                                                           factor4 * devStressTpdt.yy * devStressT.yy + ae));
            elasticityMat->C2222 = bulkModulus + 2 / (3 * (factor3 * devStressTpdt.yy * devStressTpdt.yy + factor1 +
                                                           factor4 * devStressTpdt.yy * devStressT.yy + ae));
            elasticityMat->C2323 = 1 / (2 * (factor3 * devStressTpdt.yz * devStressTpdt.yz + factor1 +
                                             factor4 * devStressTpdt.yz * devStressT.yz + ae));
            elasticityMat->C3311 = bulkModulus - 1 / (3 * (factor3 * devStressTpdt.zz * devStressTpdt.zz + factor1 +
                                                           factor4 * devStressTpdt.zz * devStressT.zz + ae));
            elasticityMat->C3333 = bulkModulus + 2 / (3 * (factor3 * devStressTpdt.zz * devStressTpdt.zz + factor1 +
                                                           factor4 * devStressTpdt.zz * devStressT.zz + ae));
        } // if
    }

}; // IsotropicPowerLaw3D

#endif // pylith_fekernels_isotropicpowerlaw_hh

// End of file
