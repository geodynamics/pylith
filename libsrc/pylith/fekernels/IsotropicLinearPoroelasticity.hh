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

/** @file libsrc/fekernels/IsotropicLinearPoroelasticity.hh
 *
 * Kernels for linear poroelasticity plane strain.
 *
 * Solution fields: [disp(dim), pressure(1),trace_strain(1) ] (QS)
 * OR
 * Solution fields: [disp(dim), pressure(1),velocity(dim) ] (DYN)
 *
 * Auxiliary fields:
 * -- numA : number of auxiliary fields
 ***** Required fields(govening equations) + option fields + required fields (rheology)
 * - 0: solid_density(1)
 * - 1: fluid_density(1)
 * - 2: fluid_viscosity(1)
 * - 3: porosity(1)
 *
 ** Optional fields
 * - +1: gravity_field (dim, optional) (4)
 * - +1: body_force(dim,optional) (4,5)
 * - +1: source_density(1,optional) (4,5,6)
 * - +1: reference_stress(4,optional) (stress_xx, stress_yy, stress_xy, stress_zz)  // numA - 7
 *     2D: 4 components (stress_xx, stress_yy, stress_zz, stress_xy)
 *     3D: 6 components (stress_xx, stress_yy, stress_zz, stress_xy, stress_yz, stress_xz)
 * - +1: reference_strain(4,optional) (strain_xx, strain_yy, strain_xy, strain_zz) // numA - 6
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 *
 ** Rheological fields
 * - numA - 5: addShearModulus(1)
 * - numA - 4: addDrainedBulkModulus(1)
 * - numA - 3: addBiotCoefficient(1)
 *      Isotropic: numA - 2: addIsotropicPermeability(1)
 *      Tensor: numA - 2: addTensorPermeability(4,optional) (permeability_xx, permeability_yy, permeability_zz,
 * permeability_xy)
 *          2D: 4 components (permeability_xx, permeability_yy, permeability_zz, permeability_xy)
 *          3D: 6 components (permeability_xx, permeability_yy, permeability_zz, permeability_xy, permeability_yz,
 * permeability_xz)
 * - numA - 1: addFluidBulkModulus(1)
 *
 * The poroelasticity subfields come first (with required ones before optional ones) followed by the rheology subfields
 * (optional ones before required ones). The rheology fields have required fields last because we index from the back.
 *
 * :TODO: @robert Add equation here
 *
 *
 * ======================================================================
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

#if !defined(pylith_fekernels_isotropiclinearporoelasticity_hh)
#define pylith_fekernels_isotropiclinearporoelasticity_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations
#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linear poroelasticity.
class pylith::fekernels::IsotropicLinearPoroelasticity {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    struct Context {
        PylithReal pressure;
        PylithReal trace_strain;
        PylithReal trace_strain_t;
        PylithReal fluidViscosity; // Poroelastic Auxiliaries
        PylithReal shearModulus; // Rheologic Auxiliaries
        PylithReal drainedBulkModulus;
        PylithReal biotCoefficient;
        PylithReal biotModulus;
        pylith::fekernels::Tensor permeability;
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

        // Incoming solution subfields
        const PylithInt i_pressure = 1;

        // Incoming poroelastic auxiliary subfields
        const PylithInt i_fluidViscosity = 2;

        // Incoming rheology auxiliary subfields.
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_biotModulus = numA - 2;

        assert(numA >= 3); // also have density
        assert(s);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(a);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_drainedBulkModulus] >= 0);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(aOff[i_fluidViscosity] >= 0);

        // Solution Variables
        context->pressure = s[sOff[i_pressure]];

        // Poroelastic Auxiliary Variables
        context->fluidViscosity = a[aOff[i_fluidViscosity]];assert(context->fluidViscosity > 0.0);

        // Rheology Specific Auxiliary Variables
        context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
        context->drainedBulkModulus = a[aOff[i_drainedBulkModulus]];assert(context->drainedBulkModulus > 0.0);
        context->biotCoefficient = a[aOff[i_biotCoefficient]];assert(context->biotCoefficient > 0.0);
        context->biotModulus = a[aOff[i_biotModulus]];assert(context->biotModulus > 0.0);

    } // setContext

    // --------------------------------------------------------------------------------------------
    static inline
    void setContextIsotropicPerm(Context* context,
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

        // Incoming auxiliary fields.
        const PylithInt i_isotropicPermeability = numA - 1;

        // Using isotropic permeability
        tensorOps.fromScalar(a[aOff[i_isotropicPermeability]], &context->permeability);

    } // setContextTensorPerm

    // --------------------------------------------------------------------------------------------
    static inline
    void setContextTensorPerm(Context* context,
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

        // Incoming auxiliary fields.
        const PylithInt i_tensorPermeability = numA - 1;

        // Using tensor permeability
        tensorOps.fromVector(&a[aOff[i_tensorPermeability]], &context->permeability);

    } // setContextTensorPerm

    // --------------------------------------------------------------------------------------------
    static inline
    void setContextRefState(Context* context,
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
        // Incoming auxiliary fields.
        const PylithInt i_refStress = numA - 7;
        const PylithInt i_refStrain = numA - 6;

        // Reference stress and strain
        tensorOps.fromVector(&a[aOff[i_refStress]], &context->refStress);
        tensorOps.fromVector(&a[aOff[i_refStrain]], &context->refStrain);

    } // setContextTensorPerm

    // --------------------------------------------------------------------------------------------
    static inline
    void setContextQuasistatic(Context* context,
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
        // Incoming solution fields.
        const PylithInt i_trace_strain = 2;
        assert(sOff[i_trace_strain] >= 0);

        // Variables &c
        context->trace_strain = s[sOff[i_trace_strain]];

    } // setContextQuasistatic

    // --------------------------------------------------------------------------------------------
    static inline
    void setContextDynamic(Context* context,
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
        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        // const PylithInt i_pressure = 1;
        const PylithInt i_velocity = 2;

        assert(sOff[i_displacement] >= 0);
        // assert(sOff[i_pressure] >= 0);
        assert(sOff[i_velocity] >= 0);
        // Improvised values

        const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar *velocity_x = &s_x[sOff_x[i_velocity]];

        PylithScalar trace_strain = 0.0;
        for (PylithInt d = 0; d < dim; ++d) {
            trace_strain += displacement_x[d * dim + d];
        }

        PylithScalar trace_strain_t = 0.0;
        for (PylithInt d = 0; d < dim; ++d) {
            trace_strain_t += velocity_x[d * dim + d];
        }

        // Variables
        context->trace_strain = trace_strain;
        context->trace_strain_t = trace_strain_t;

    } // setContextDynamic

    // --------------------------------------------------------------------------------------------
    /** Helper function for calculating Cauchy stress for WITHOUT a reference stress and strain.
     *
     * ISA Poroelasticity::stressFn
     *
     * @param[in] rheologyContext IsotropicLinearPoroelasticity context.
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

        meanStress(context->pressure, context->trace_strain, context->drainedBulkModulus, context->biotCoefficient, strain, stress);
        deviatoricStress(context->trace_strain, context->shearModulus, strain, stress);
    } // cauchyStress

    // --------------------------------------------------------------------------------------------
    /** Helper function for calculating Cauchy stress WITH reference stress/strain.
     *
     * @param[in] rheologyContext IsotropicLinearPoroelastcity context.
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
        meanStress_refState(context->pressure, context->trace_strain, context->drainedBulkModulus, context->biotCoefficient, refStress, refStrain, strain, stress);
        deviatoricStress_refState(context->trace_strain, context->shearModulus, refStress, refStrain, strain, stress);
    } // cauchyStress_refState

    // --------------------------------------------------------------------------------------------
    /** Calculate mean (volumetric) stress WITHOUT reference stress and reference strain.
     */
    static inline
    void meanStress(const PylithReal pressure,
                    const PylithReal strainTrace,
                    const PylithReal drainedBulkModulus,
                    const PylithReal biotCoefficient,
                    const pylith::fekernels::Tensor& strain,
                    pylith::fekernels::Tensor* stress) {
        assert(drainedBulkModulus > 0.0);
        assert(stress);

        // const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        const PylithReal meanStress = drainedBulkModulus * strainTrace - biotCoefficient*pressure;

        stress->xx += meanStress;
        stress->yy += meanStress;
        stress->zz += meanStress;
    } // meanStress

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITHOUT reference stress and strain.
     */
    static inline
    void deviatoricStress(const PylithReal strainTrace,
                          const PylithReal shearModulus,
                          const pylith::fekernels::Tensor& strain,
                          pylith::fekernels::Tensor* stress) {
        assert(shearModulus > 0.0);
        assert(stress);

        // const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
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
    void meanStress_refState(const PylithReal pressure,
                             const PylithReal strainTrace,
                             const PylithReal drainedBulkModulus,
                             const PylithReal biotCoefficient,
                             const pylith::fekernels::Tensor& refStress,
                             const pylith::fekernels::Tensor& refStrain,
                             const pylith::fekernels::Tensor& strain,
                             pylith::fekernels::Tensor* stress) {
        // Incoming auxiliary fields.
        assert(drainedBulkModulus > 0.0);
        assert(stress);

        // const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        const PylithReal refStrainTrace = refStrain.xx + refStrain.yy + refStrain.zz;

        const PylithReal meanRefStress = (refStress.xx + refStress.yy + refStress.zz) / 3.0;
        const PylithReal meanStress = meanRefStress + drainedBulkModulus * (strainTrace - refStrainTrace) - biotCoefficient*pressure;

        stress->xx += meanStress;
        stress->yy += meanStress;
        stress->zz += meanStress;
    } // meanStress_refState

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITH reference stress and strain.
     */
    static inline
    void deviatoricStress_refState(const PylithReal strainTrace,
                                   const PylithReal shearModulus,
                                   const pylith::fekernels::Tensor& refStress,
                                   const pylith::fekernels::Tensor& refStrain,
                                   const pylith::fekernels::Tensor& strain,
                                   pylith::fekernels::Tensor* stress) {
        assert(shearModulus > 0.0);
        assert(stress);

        // const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
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

    // ================================ Kernels ====================================
    // --------------------------------------------------------------------------------------------
    /** Helper function for calculating Cauchy stress for WITHOUT a reference stress and strain.
     *
     * ISA Poroelasticity::fluxrateFn
     *
     * @param[in] rheologyContext IsotropicLinearElasticity context.
     * @param[in] strain Strain tensor.
     * @param[in] tensorOps Tensor operations.
     * @param[out] stress Stress tensor.
     */
    static inline
    void darcyFluxRate(const pylith::fekernels::Poroelasticity::Context& poroelasticContext,
                       void* rheologyContext,
                       const pylith::fekernels::TensorOps& tensorOps,
                       pylith::fekernels::Tensor* fluxRate) {
        Context* context = (Context*)(rheologyContext);
        assert(context);
        assert(fluxRate);

        const PylithInt dim = poroelasticContext.dim;
        // Solution Variables
        const PylithScalar *pressure_x = poroelasticContext.pressure_x;

        // Poroelastic Auxiliaries
        const PylithScalar fluidDensity = poroelasticContext.fluidDensity;
        const PylithScalar fluidViscosity = poroelasticContext.fluidViscosity;
        const PylithScalar *bodyForce = poroelasticContext.bodyForce;
        const PylithScalar *gravityField = poroelasticContext.gravityField;

        // Rheological Auxiliaries

        PylithScalar tensorPermeability[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        PylithScalar fluxRateVector[3] = {0.0, 0.0, 0.0};
        tensorOps.toTensor(context->permeability, tensorPermeability);

        for (PylithInt i = 0; i < dim; ++i) {
            for (PylithInt j = 0; j < dim; j++) {
                fluxRateVector[i] += (tensorPermeability[i * dim + j] / fluidViscosity) * (pressure_x[j] - bodyForce[j] - fluidDensity * gravityField[j]);
            } // for
        } // for

        fluxRate->xx = fluxRateVector[0];
        fluxRate->yy = fluxRateVector[1];
        fluxRate->zz = fluxRateVector[2];
        fluxRate->xy = 0.0;
        fluxRate->yz = 0.0;
        fluxRate->xz = 0.0;
    } // darcyFluxRate

}; // IsotropicLinearPoroelasticity

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linear poroelasticity in Plane Strain.
class pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain {
    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    // ================================= LHS =======================================

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_explicit(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Solution Variables
        const PylithScalar pressure_t = poroelasticContext.pressure_t;

        // Rheological Auxiliaries
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        f0[0] += pressure_t / biotModulus;
    } // f0p_explicit

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms.
    static inline
    void f0p_implicit(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Solution Variables
        PylithScalar pressure_t = (poroelasticContext.pressure_t ? poroelasticContext.pressure_t : 0.0);
        PylithScalar trace_strain_t = (poroelasticContext.trace_strain_t ? poroelasticContext.trace_strain_t : 0.0);

        // Rheological Auxiliaries
        const PylithReal biotCoefficient = rheologyContext.biotCoefficient;
        const PylithReal biotModulus = rheologyContext.biotModulus;

        f0[0] += s_t ? (biotCoefficient * trace_strain_t) : 0.0;
        f0[0] += s_t ? (pressure_t / biotModulus) : 0.0;

    } // f0p_implicit

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Solution Variables
        const PylithReal pressure_t = poroelasticContext.pressure_t;
        const PylithReal trace_strain_t = poroelasticContext.trace_strain_t;

        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static
    void f0p_implicit_source_body(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextBodyForceSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Solution Variables
        const PylithReal pressure_t = poroelasticContext.pressure_t;
        const PylithReal trace_strain_t = poroelasticContext.trace_strain_t;

        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source_body

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_grav(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextGravitySourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Solution Variables
        const PylithReal pressure_t = poroelasticContext.pressure_t;
        const PylithReal trace_strain_t = poroelasticContext.trace_strain_t;

        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source_grav

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_grav_body(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextGravityBodyForceSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Solution Variables
        const PylithReal pressure_t = poroelasticContext.pressure_t;
        const PylithReal trace_strain_t = poroelasticContext.trace_strain_t;

        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source_grav_body

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Quasi - Static Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static inline
    void f1u(const PylithInt dim,
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

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextQuasistatic(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1u

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static inline
    void f1u_refstate(const PylithInt dim,
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

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextQuasistatic(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextRefState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1u_refstate

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity
     *
     */
    static inline
    void f1p(const PylithInt dim,
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
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_tensor_permeability(const PylithInt dim,
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
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body(const PylithInt dim,
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
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_body

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_tensor_permeability(const PylithInt dim,
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
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_body_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_gravity(const PylithInt dim,
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
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_gravity

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_gravity_tensor_permeability(const PylithInt dim,
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

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_gravity_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_gravity(const PylithInt dim,
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

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravityBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_body_gravity

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_gravity_tensor_permeability(const PylithInt dim,
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

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravityBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            f1);

    } // f1p_body_gravity_tensor_permeability

    // =========================== LHS Jacobian ============================

    // ----------------------------------------------------------------------
    /* Jf3_uu entry function for isotropic linear poroelasticity.
     */
    static inline
    void Jf3uu(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Rheological Auxiliaries
        const PylithScalar shearModulus = rheologyContext.shearModulus;

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; ++j) {
                Jf3[((i * _dim + i) * _dim + j) * _dim + j] -= shearModulus;
                Jf3[((i * _dim + j) * _dim + j) * _dim + i] -= shearModulus;
            }
        }
    } // Jf3uu

    // ----------------------------------------------------------------------
    /** Jf2_up entry function for isotropic linear poroelasticity.
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf2up(const PylithInt dim,
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
               PylithScalar Jf2[]) {
        const PylithInt _dim = 2;

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Rheology Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf2[d * _dim + d] += biotCoefficient;
        } // for
    } // Jf2up

    // -----------------------------------------------------------------------------
    // Jf2ue function for isotropic linear poroelasticity.
    static inline
    void Jf2ue(const PylithInt dim,
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
               PylithScalar Jf2[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Rheological Auxiliaries
        const PylithScalar shearModulus = rheologyContext.shearModulus;
        const PylithScalar drainedBulkModulus = rheologyContext.drainedBulkModulus;

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf2[d * _dim + d] -= drainedBulkModulus - (2.0 * shearModulus) / 3.0;
        } // for
    } // Jf2ue

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf3pp(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Poroelastic Auxiliaries
        const PylithScalar fluidViscosity = rheologyContext.fluidViscosity;

        // Rheological Auxiliaries
        // const PylithScalar *tensorPermeablity = rheologyContext.permeability;
        PylithScalar tensorPermeability[4] = {0.0, 0.0, 0.0, 0.0};
        pylith::fekernels::Tensor::ops2D.toTensor(rheologyContext.permeability, tensorPermeability);

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                Jf3[i * _dim + j] += tensorPermeability[i * _dim + j] / fluidViscosity;
            } // for
        } // for

    } // Jf3pp

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf3pp_tensor_permeability(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Poroelastic Auxiliaries
        const PylithScalar fluidViscosity = rheologyContext.fluidViscosity;

        // Rheological Auxiliaries
        PylithScalar tensorPermeability[4] = {0.0, 0.0, 0.0, 0.0};
        pylith::fekernels::Tensor::ops2D.toTensor(rheologyContext.permeability, tensorPermeability);

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                Jf3[i * _dim + j] += tensorPermeability[i * _dim + j] / fluidViscosity;
            } // for
        } // for
    } // Jf3pp_tensor_permeability

    // ----------------------------------------------------------------------
    /** Jf0_pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0pp(const PylithInt dim,
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
               PylithScalar Jf0[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Rheological Auxiliaries
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        Jf0[0] += utshift / biotModulus;

    } // Jf0pp

    // ----------------------------------------------------------------------
    /** Jf0_pe entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0pe(const PylithInt dim,
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
               PylithScalar Jf0[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        Jf0[0] += utshift * biotCoefficient;
    } // Jf0pe

    // ----------------------------------------------------------------------
    /** Jf0_ppdot entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0ppdot(const PylithInt dim,
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
                  PylithScalar Jf0[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Rheological Auxiliaries
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        Jf0[0] += 1.0 / biotModulus;
    } // Jf0ppdot

    // ----------------------------------------------------------------------
    /** Jf0_pedot entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0pedot(const PylithInt dim,
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
                  PylithScalar Jf0[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        Jf0[0] += biotCoefficient;
    } // Jf0pedot

    // ============================== RHS Residual =================================

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Solution Variables
        const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;

        // Rheology Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_implicit

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Solution Variables
        const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;

        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;

        // Rheologic Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source_body(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextBodyForceSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Solution Variables
        const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;

        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;

        // Rheologic Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source_body

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source_grav(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravitySourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Solution Variables
        const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;

        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;

        // Rheologic Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source_grav

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source_grav_body(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravityBodyForceSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Solution Variables
        const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;

        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;

        // Rheologic Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source_grav_body

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p(const PylithInt dim,
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
             PylithScalar g1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            g1);

    } // g1p

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p_tensor_permeability(const PylithInt dim,
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
                                 PylithScalar g1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            g1);

    } // g1p_tensor_permeability

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p_gravity(const PylithInt dim,
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
                     PylithScalar g1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            g1);

    } // g1p_gravity

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p_gravity_tensor_permeability(const PylithInt dim,
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
                                         PylithScalar g1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops2D,
            g1);

    } // g1p_gravity_tensor_permeability

    // -----------------------------------------------------------------------------
    /** g1v function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Dynamic Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static inline
    void g1v(const PylithInt dim,
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
             PylithScalar g1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextDynamic(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            g1);

    } // g1v

    // -----------------------------------------------------------------------------
    /** g1v function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static inline
    void g1v_refstate(const PylithInt dim,
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
                      PylithScalar g1[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextRefState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextDynamic(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops2D,
            g1);

    } // g1v_refstate

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic linear poroelasticity with
     * infinitesimal strain WITHOUT reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., biot_coefficient(1), shear_modulus(1), drained_bulk_modulus(1)]
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

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        // Using dynamic formulation for trace strain, assuming that it will be equal to the variable
        // for QS
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextDynamic(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            stressVector);

    } // cauchyStress_infinitesimalStrain_asVector

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic linear poroelasticity with
     * infinitesimal strain WITH a reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., biot_coefficient(1), shear_modulus(1), drained_bulk_modulus(1)]
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

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextRefState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);
        // Using dynamic formulation for trace strain, assuming that it will be equal to the variable
        // for QS
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextDynamic(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops2D, stressVector);
    } // cauchyStress_infinitesimalStrain_refState_asVector

    // ========================== Update Kernels ===================================

    // ---------------------------------------------------------------------------------------------------------------------
    /* Update porosity for a linear poroelastic material, implicit.
     */
    static inline
    void updatePorosityImplicit(const PylithInt dim,
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
                                PylithScalar porosity[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Incoming solution fields.
        const PylithInt i_pressure_t = 4;
        const PylithInt i_trace_strain_t = 5;

        // Incoming re-packed auxiliary field.

        // Poroelasticity
        const PylithInt i_porosityPrev = 3;

        // IsotropicLinearPoroelasticity
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        // Constants
        const PylithScalar dt = constants[0];

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 6);
        assert(numA >= 6);
        assert(aOff);
        assert(aOff[i_porosityPrev] >= 0);
        assert(porosity);

        // Do stuff
        const PylithScalar pressure_t = s ? s[sOff[i_pressure_t]] : 0.0;
        const PylithScalar trace_strain_t = s ? s[sOff[i_trace_strain_t]] : 0.0;

        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar porosityPrev = a[aOff[i_porosityPrev]];

        // Update porosity
        porosity[0] = porosityPrev + dt * ((biotCoefficient - porosityPrev) * trace_strain_t +
                                           ((1.0 - biotCoefficient) * (biotCoefficient - porosityPrev)) /
                                           drainedBulkModulus * pressure_t);
        porosity[0] = std::max(0.0, std::min(1.0, porosity[0]));

        // setContextQuasistatic sets pressure_t - this throws a segfault / memory out of range error.

        // // Poroelastic Context
        // pylith::fekernels::Poroelasticity::Context poroelasticContext;
        // pylith::fekernels::Poroelasticity::setContextQuasistatic(
        //     &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        // pylith::fekernels::Poroelasticity::setContextQS_sixField(
        //     &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // // Rheology Context
        // pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        // pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
        //     &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
        //     t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // // Solution Variables
        // const PylithScalar pressure_t = poroelasticContext.pressure_t;
        // const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;

        // // Poroelastic Auxiliaries
        // const PylithScalar porosityPrev = poroelasticContext.porosity;

        // // Rheologic Auxiliaries
        // const PylithScalar drainedBulkModulus = rheologyContext.drainedBulkModulus;
        // const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        // // Constants
        // const PylithScalar dt = constants[0];

        // // Update porosity
        // porosity[0] = porosityPrev + dt * ((biotCoefficient - porosityPrev) * trace_strain_t +
        //                                    ((1.0 - biotCoefficient) * (biotCoefficient - porosityPrev)) /
        //                                    drainedBulkModulus * pressure_t);

    } // updatePorosityImplicit

    // ---------------------------------------------------------------------------------------------------------------------
    /* Update porosity for a linear poroelastic material, explicit.
     */
    static inline
    void updatePorosityExplicit(const PylithInt dim,
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
                                PylithScalar porosity[]) {
        const PylithInt _dim = 2;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // Solution Variables
        const PylithScalar pressure_t = poroelasticContext.pressure_t;
        const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;

        // Poroelastic Auxiliaries
        const PylithScalar porosityPrev = poroelasticContext.porosity;

        // Rheologic Auxiliaries
        const PylithScalar drainedBulkModulus = rheologyContext.drainedBulkModulus;
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        // Constants
        const PylithScalar dt = constants[0];

        // Update porosity
        porosity[0] = porosityPrev + dt * ((biotCoefficient - porosityPrev) * trace_strain_t +
                                           ((1.0 - biotCoefficient) * (biotCoefficient - porosityPrev)) /
                                           drainedBulkModulus * pressure_t);
        porosity[0] = std::max(0.0, std::min(1.0, porosity[0]));
    } // updatePorosityExplicit

}; // IsotropicLinearPoroelasticityPlaneStrain

// ------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linear poroelasticity in 3D.
class pylith::fekernels::IsotropicLinearPoroelasticity3D {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // ================================= LHS =======================================

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_explicit(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Solution Variables
        const PylithScalar pressure_t = poroelasticContext.pressure_t;

        // Rheology Auxiliaries
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        f0[0] += pressure_t / biotModulus;
    } // f0p_explicit

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms.
    static inline
    void f0p_implicit(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Solution Variables
        PylithScalar pressure_t = (poroelasticContext.pressure_t ? poroelasticContext.pressure_t : 0.0);
        PylithScalar trace_strain_t = (poroelasticContext.trace_strain_t ? poroelasticContext.trace_strain_t : 0.0);

        // Rheological Auxiliaries
        const PylithReal biotCoefficient = rheologyContext.biotCoefficient;
        const PylithReal biotModulus = rheologyContext.biotModulus;

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;

    } // f0p_implicit

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Solution Variables
        PylithScalar pressure_t = (poroelasticContext.pressure_t ? poroelasticContext.pressure_t : 0.0);
        PylithScalar trace_strain_t = (poroelasticContext.trace_strain_t ? poroelasticContext.trace_strain_t : 0.0);

        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_body(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextBodyForceSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Solution Variables
        PylithScalar pressure_t = (poroelasticContext.pressure_t ? poroelasticContext.pressure_t : 0.0);
        PylithScalar trace_strain_t = (poroelasticContext.trace_strain_t ? poroelasticContext.trace_strain_t : 0.0);

        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source_body

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_grav(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravitySourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Solution Variables
        PylithScalar pressure_t = (poroelasticContext.pressure_t ? poroelasticContext.pressure_t : 0.0);
        PylithScalar trace_strain_t = (poroelasticContext.trace_strain_t ? poroelasticContext.trace_strain_t : 0.0);

        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source_grav

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_grav_body(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravityBodyForceSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Solution Variables
        PylithScalar pressure_t = (poroelasticContext.pressure_t ? poroelasticContext.pressure_t : 0.0);
        PylithScalar trace_strain_t = (poroelasticContext.trace_strain_t ? poroelasticContext.trace_strain_t : 0.0);

        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source_grav_body

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Quasi - Static Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static inline
    void f1u(const PylithInt dim,
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

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextQuasistatic(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1u

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */

    static inline
    void f1u_refstate(const PylithInt dim,
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

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextRefState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextQuasistatic(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1u_refstate

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p(const PylithInt dim,
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

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_tensor_permeability(const PylithInt dim,
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

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body(const PylithInt dim,
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

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_gravity

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_tensor_permeability(const PylithInt dim,
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

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextBodyForce(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_body_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_gravity(const PylithInt dim,
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

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_gravity

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_gravity_tensor_permeability(const PylithInt dim,
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

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_gravity_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_gravity(const PylithInt dim,
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

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);
    } // f1p_body_gravity

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline void f1p_body_gravity_tensor_permeability(const PylithInt dim,
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

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            f1);

    } // f1p_body_gravity_tensor_permeability

    // =========================== LHS Jacobian ============================

    // ----------------------------------------------------------------------
    /* Jf3_uu entry function for isotropic linear poroelasticity.
     */
    static inline
    void Jf3uu(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Rheological Auxiliaries
        const PylithScalar shearModulus = rheologyContext.shearModulus;

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; ++j) {
                Jf3[((i * _dim + i) * _dim + j) * _dim + j] -= shearModulus;
                Jf3[((i * _dim + j) * _dim + j) * _dim + i] -= shearModulus;
            } // for
        } // for

    } // Jf3uu

    // ----------------------------------------------------------------------
    /** Jf2_up entry function for isotropic linear poroelasticity.
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf2up(const PylithInt dim,
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
               PylithScalar Jf2[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf2[d * _dim + d] += biotCoefficient;
        } // for

    } // Jf2up

    // -----------------------------------------------------------------------------
    // Jf2ue function for isotropic linear poroelasticity.
    static inline
    void Jf2ue(const PylithInt dim,
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
               PylithScalar Jf2[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Rheological Auxiliaries
        const PylithScalar shearModulus = rheologyContext.shearModulus;
        const PylithScalar drainedBulkModulus = rheologyContext.drainedBulkModulus;

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf2[d * _dim + d] -= drainedBulkModulus - (2.0 * shearModulus) / 3.0;
        } // for

    } // Jf2ue

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf3pp(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Poroelastic Auxiliaries
        const PylithScalar fluidViscosity = rheologyContext.fluidViscosity;

        // Rheological Auxiliaries
        PylithScalar tensorPermeability[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        pylith::fekernels::Tensor::ops3D.toTensor(rheologyContext.permeability, tensorPermeability);

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                Jf3[i * _dim + j] += tensorPermeability[i * _dim + j] / fluidViscosity;
            } // for
        } // for

    } // Jf3pp

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf3pp_tensor_permeability(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Poroelastic Auxiliaries
        const PylithScalar fluidViscosity = rheologyContext.fluidViscosity;

        // Rheological Auxiliaries
        PylithScalar tensorPermeability[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        pylith::fekernels::Tensor::ops3D.toTensor(rheologyContext.permeability, tensorPermeability);

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                Jf3[i * _dim + j] += tensorPermeability[i * _dim + j] / fluidViscosity;
            } // for
        } // for
    } // Jf3pp_tensor_permeability

    // ----------------------------------------------------------------------
    /** Jf0_pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0pp(const PylithInt dim,
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
               PylithScalar Jf0[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Rheological Auxiliaries
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        Jf0[0] = utshift / biotModulus;

    } // Jf0pp

    // ----------------------------------------------------------------------
    /** Jf0_pe entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0pe(const PylithInt dim,
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
               PylithScalar Jf0[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        Jf0[0] += utshift * biotCoefficient;

    } // Jf0pe

    // ----------------------------------------------------------------------
    /** Jf0_ppdot entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0ppdot(const PylithInt dim,
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
                  PylithScalar Jf0[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Rheological Auxiliaries
        const PylithScalar biotModulus = rheologyContext.biotModulus;

        Jf0[0] += 1.0 / biotModulus;
    } // Jf0ppdot

    // ----------------------------------------------------------------------
    /** Jf0_pedot entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0pedot(const PylithInt dim,
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
                  PylithScalar Jf0[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Rheological Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        Jf0[0] += biotCoefficient;
    } // Jf0pedot

    // ============================== RHS Residual =================================

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p(const PylithInt dim,
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
        const PylithInt _dim = 3;

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Solution Variables
        const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;

        // Rheology Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_implicit

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source(const PylithInt dim,
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
        const PylithInt _dim = 3;

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Solution Variables
        const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;
        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;
        // Rheologic Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source_body(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextBodyForceSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Solution Variables
        const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;
        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;
        // Rheologic Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source_body

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source_grav(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravitySourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Solution Variables
        const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;
        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;
        // Rheologic Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source_grav

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source_grav_body(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravitySourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Solution Variables
        const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;
        // Poroelastic Auxiliaries
        const PylithScalar source = poroelasticContext.sourceDensity;
        // Rheologic Auxiliaries
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source_grav_body

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / without gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p(const PylithInt dim,
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
             PylithScalar g1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            g1);

    } // g1p

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p_tensor_permeability(const PylithInt dim,
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
                                 PylithScalar g1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            g1);

    } // g1p_tensor_permeability

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p_gravity(const PylithInt dim,
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
                     PylithScalar g1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextIsotropicPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            g1);

    } // g1p_gravity

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p_gravity_tensor_permeability(const PylithInt dim,
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
                                         PylithScalar g1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextTensorPerm(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Use f1p / fluxrate / darcy function
        pylith::fekernels::Poroelasticity::f1p(
            poroelasticContext, &rheologyContext,
            pylith::fekernels::IsotropicLinearPoroelasticity::darcyFluxRate,
            pylith::fekernels::Tensor::ops3D,
            g1);

    } // g1p_gravity_tensor_permeability

    // -----------------------------------------------------------------------------
    /** g1v function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Dynamic Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static inline
    void g1v(const PylithInt dim,
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
             PylithScalar g1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextDynamic(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            g1);

    } // g1v

    // -----------------------------------------------------------------------------
    /** g1v function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static inline
    void g1v_refstate(const PylithInt dim,
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
                      PylithScalar g1[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextRefState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextDynamic(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops3D,
            g1);

    } // g1v_refstate

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic linear poroelasticity with
     * infinitesimal strain WITHOUT reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., biot_coefficient(1), shear_modulus(1), drained_bulk_modulus(1)]
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

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheological Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        // Using dynamic formulation for trace strain, assuming that it will be equal to the variable
        // for QS
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextDynamic(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            stressVector);

    } // cauchyStress_infinitesimalStrain_asVector

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic linear poroelasticity with
     * infinitesimal strain WITH a reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., biot_coefficient(1), shear_modulus(1), drained_bulk_modulus(1)]
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

        // Strain Context
        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextRefState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);
        // Using dynamic formulation for trace strain, assuming that it will be equal to the variable
        // for QS
        pylith::fekernels::IsotropicLinearPoroelasticity::setContextDynamic(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops3D, stressVector);
    } // cauchyStress_infinitesimalStrain_refState_asVector

    // ========================== Update Kernels ===================================

    // ---------------------------------------------------------------------------------------------------------------------
    /* Update porosity for a linear poroelastic material, implicit.
     */
    static inline
    void updatePorosityImplicit(const PylithInt dim,
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
                                PylithScalar porosity[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Incoming solution fields.
        const PylithInt i_pressure_t = 4;
        const PylithInt i_trace_strain_t = 5;

        // Incoming re-packed auxiliary field.

        // Poroelasticity
        const PylithInt i_porosityPrev = 3;

        // IsotropicLinearPoroelasticity
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        // Constants
        const PylithScalar dt = constants[0];

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 6);
        assert(numA >= 6);
        assert(aOff);
        assert(aOff[i_porosityPrev] >= 0);
        assert(porosity);

        // Do stuff
        const PylithScalar pressure_t = s ? s[sOff[i_pressure_t]] : 0.0;
        const PylithScalar trace_strain_t = s ? s[sOff[i_trace_strain_t]] : 0.0;

        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar porosityPrev = a[aOff[i_porosityPrev]];

        // Update porosity
        porosity[0] = porosityPrev + dt * ((biotCoefficient - porosityPrev) * trace_strain_t +
                                           ((1.0 - biotCoefficient) * (biotCoefficient - porosityPrev)) /
                                           drainedBulkModulus * pressure_t);
        porosity[0] = std::max(0.0, std::min(1.0, porosity[0]));

        // setContextQuasistatic sets pressure_t - this throws a segfault / memory out of range error.

        // // Poroelastic Context
        // pylith::fekernels::Poroelasticity::Context poroelasticContext;
        // pylith::fekernels::Poroelasticity::setContextQuasistatic(
        //     &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        // pylith::fekernels::Poroelasticity::setContextQS_sixField(
        //     &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // // Rheology Context
        // pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        // pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
        //     &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
        //     t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        // // Solution Variables
        // const PylithScalar pressure_t = poroelasticContext.pressure_t;
        // const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;

        // // Poroelastic Auxiliaries
        // const PylithScalar porosityPrev = poroelasticContext.porosity;

        // // Rheologic Auxiliaries
        // const PylithScalar drainedBulkModulus = rheologyContext.drainedBulkModulus;
        // const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        // // Constants
        // const PylithScalar dt = constants[0];

        // // Update porosity
        // porosity[0] = porosityPrev + dt * ((biotCoefficient - porosityPrev) * trace_strain_t +
        //                                    ((1.0 - biotCoefficient) * (biotCoefficient - porosityPrev)) /
        //                                    drainedBulkModulus * pressure_t);
    } // updatePorosityImplicit

    // ---------------------------------------------------------------------------------------------------------------------
    /* Update porosity for a linear poroelastic material, explicit.
     */
    static inline
    void updatePorosityExplicit(const PylithInt dim,
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
                                PylithScalar porosity[]) {
        const PylithInt _dim = 3;assert(_dim == dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroelasticContext;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        pylith::fekernels::Poroelasticity::setContextSourceDensity(
            &poroelasticContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Rheology Context
        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        // Solution Variables
        const PylithScalar pressure_t = poroelasticContext.pressure_t;
        const PylithScalar trace_strain_t = poroelasticContext.trace_strain_t;

        // Poroelastic Auxiliaries
        const PylithScalar porosityPrev = poroelasticContext.porosity;

        // Rheologic Auxiliaries
        const PylithScalar drainedBulkModulus = rheologyContext.drainedBulkModulus;
        const PylithScalar biotCoefficient = rheologyContext.biotCoefficient;

        // Constants
        const PylithScalar dt = constants[0];

        // Update porosity
        porosity[0] = porosityPrev + dt * ((biotCoefficient - porosityPrev) * trace_strain_t +
                                           ((1.0 - biotCoefficient) * (biotCoefficient - porosityPrev)) /
                                           drainedBulkModulus * pressure_t);
        porosity[0] = std::max(0.0, std::min(1.0, porosity[0]));

    } // updatePorosityExplicit

}; // IsotropicLinearPoroelasticity3D

#endif // pylith_fekernels_isotropiclinearporoelasticity_hh

// End of file
