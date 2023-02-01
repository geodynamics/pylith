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
 * - 1: body_force(2,optional)
 * - 2: gravity_field (2, optional)
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
 * Viscous strain must be before total strain, because viscous strain at t+dt depends on total strain at t.
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
#include "pylith/fekernels/Viscoelasticity.hh" // USES Viscoelasticity kernels

#include "pylith/utils/types.hh"

// ------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linear Maxwell plane strain.
class pylith::fekernels::IsotropicLinearMaxwellPlaneStrain {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    /** f1 function for isotropic linear Maxwell plane strain WITHOUT reference stress and reference strain.
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
        pylith::fekernels::ElasticityPlaneStrain::f1v(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain, _cauchyStress,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** f1 function for isotropic linear Maxwell plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void f1v_infinitesimalStrain_refstate(const PylithInt dim,
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
        pylith::fekernels::ElasticityPlaneStrain::f1v(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain, _cauchyStress_refstate,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate stress for 2D plane strain isotropic linear elasticity WITHOUT a reference stress and strain.
     *
     * Used in output of Cauchy stress.
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
        assert(stressVector);

        PylithScalar strainTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
        ElasticityPlaneStrain::infinitesimalStrain(dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);
        _cauchyStress_stateVars_asVector(dim, numA, aOff, a, strainTensor, stressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate stress for 2D plane strain isotropic linear elasticity WITH a reference stress and strain.
     *
     * Used in output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void cauchyStress_infinitesimalStrain_refstate_asVector(const PylithInt dim,
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
        assert(stressVector);

        PylithScalar strainTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
        ElasticityPlaneStrain::infinitesimalStrain(dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);
        _cauchyStress_refstate_stateVars_asVector(dim, numA, aOff, a, strainTensor, stressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate viscous strain for 2D plane strain isotropic linear elasticity.
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
        assert(viscousStrain);

        PylithScalar strainTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
        ElasticityPlaneStrain::infinitesimalStrain(dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);

        // Incoming auxiliary fields.
        const PylithInt i_maxwellTime = numA-3;
        const PylithInt i_viscousStrain = numA-2;
        const PylithInt i_totalStrain = numA-1;

        assert(numA >= 6);
        assert(aOff);
        assert(aOff[i_maxwellTime] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_totalStrain] >= 0);
        assert(a);
        assert(1 == numConstants);
        assert(constants);

        const PylithReal maxwellTime = a[aOff[i_maxwellTime]];
        const PylithReal* viscousStrainPrev = &a[aOff[i_viscousStrain]];
        const PylithReal* totalStrain = &a[aOff[i_totalStrain]];
        const PylithReal dt = constants[0];
        _viscousStrain_asVector(dim, maxwellTime, viscousStrainPrev, totalStrain, dt, strainTensor, viscousStrain);
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
        const PylithInt _dim = 2;

        // Incoming auxiliary fields.
        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_bulkModulus = numA-4;
        const PylithInt i_maxwellTime = numA-3;

        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 5);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_maxwellTime] >= 0);
        assert(a);
        assert(Jf3);
        assert(numConstants == 1);
        assert(constants);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];
        const PylithScalar dt = constants[0];

        const PylithScalar dq = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime);

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

private:

    // --------------------------------------------------------------------------------------------
    /** Calculate stress tensor for 2D plane strain isotropic linear elasticity
     * WITHOUT a reference stress and strain.
     *
     * ISA Elasticity::stressfn_type
     *
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1),
     * viscous_strain(4), total_strain(4)]
     */
    static inline
    void _cauchyStress(const PylithInt dim,
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
                       const PylithScalar strain[],
                       PylithScalar stress[]) {
        const PylithInt _dim = 2;

        // Auxiliary fields used.
        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_bulkModulus = numA-4;
        const PylithInt i_maxwellTime = numA-3;
        const PylithInt i_viscousStrain = numA-2;
        const PylithInt i_totalStrain = numA-1;

        assert(_dim == dim);
        assert(numA >= 6);
        assert(a);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_maxwellTime] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_totalStrain] >= 0);
        assert(1 == numConstants);
        assert(constants);
        assert(strain);
        assert(stress);

        const PylithReal bulkModulus = a[aOff[i_bulkModulus]];
        IsotropicLinearElasticityPlaneStrain::meanStress(_dim, bulkModulus, strain, stress);

        const PylithReal shearModulus = a[aOff[i_shearModulus]];
        const PylithReal maxwellTime = a[aOff[i_maxwellTime]];
        const PylithReal* viscousStrainPrev = &a[aOff[i_viscousStrain]];
        const PylithReal* totalStrain = &a[aOff[i_totalStrain]];
        const PylithReal dt = constants[0];
        _deviatoricStress(_dim, shearModulus, maxwellTime, viscousStrainPrev, totalStrain, strain, dt, stress);
    } // cauchyStress

    // --------------------------------------------------------------------------------------------
    /** Calculate stress vector for 2D plane strain isotropic linear elasticity
     * WITHOUT a reference stress and strain while making use of state variables.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void _cauchyStress_stateVars_asVector(const PylithInt dim,
                                          const PylithInt numA,
                                          const PylithInt aOff[],
                                          const PylithScalar a[],
                                          const PylithScalar strain[],
                                          PylithScalar stressVector[]) {
        const PylithInt _dim = 2;

        // Incoming auxiliary fields.
        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_bulkModulus = numA-4;
        const PylithInt i_viscousStrain = numA-2;

        assert(numA >= 3);
        assert(aOff);

        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_viscousStrain] >= 0);

        PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 }; // Full stress tensor in vector form.
        const PylithReal bulkModulus = a[aOff[i_bulkModulus]];
        IsotropicLinearElasticityPlaneStrain::meanStress(_dim, bulkModulus, strain, stressTensor);
        const PylithScalar meanStress = stressTensor[0];

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar* viscousStrain = &a[aOff[i_viscousStrain]];
        stressVector[0] = meanStress + 2.0 * shearModulus * viscousStrain[0]; // stress_xx
        stressVector[1] = meanStress + 2.0 * shearModulus * viscousStrain[1]; // stress_yy
        stressVector[2] = meanStress + 2.0 * shearModulus * viscousStrain[2]; // stress_zz
        stressVector[3] = 2.0 * shearModulus * viscousStrain[3]; // stress_xy
    } // cauchyStress_stateVars_asVector

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress for 2D plane strain isotropic linear
     * Maxwell viscoelasticity WITHOUT reference stress and strain.
     */
    static inline
    void _deviatoricStress(const PylithInt dim,
                           const PylithReal shearModulus,
                           const PylithReal maxwellTime,
                           const PylithReal viscousStrainPrev[],
                           const PylithReal strainPrev[],
                           const PylithScalar strain[],
                           const PylithReal dt,
                           PylithScalar stress[]) {
        const PylithInt _dim = 2;

        assert(_dim == dim);
        assert(shearModulus > 0.0);
        assert(maxwellTime > 0.0);
        assert(viscousStrainPrev);
        assert(strainPrev);
        assert(strain);
        assert(dt > 0.0);
        assert(stress);

        PylithScalar viscousStrain[4] = {0.0, 0.0, 0.0, 0.0}; // Viscous strain vector.
        _viscousStrain_asVector(_dim, maxwellTime, viscousStrainPrev, strainPrev, dt, strain, viscousStrain);

        stress[0] += 2.0 * shearModulus * viscousStrain[0]; // stress_xx
        stress[1] += 2.0 * shearModulus * viscousStrain[3]; // stress_xy
        stress[2] += 2.0 * shearModulus * viscousStrain[3]; // stress_yx
        stress[3] += 2.0 * shearModulus * viscousStrain[1]; // stress_yy
    } // deviatoricStress

    // --------------------------------------------------------------------------------------------
    /** Calculate stress tensor for 2D plane strain isotropic linear
     * Maxwell WITH a reference stress/strain.
     *
     * ISA Elasticity::stressfn_type
     *
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                   maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void _cauchyStress_refstate(const PylithInt dim,
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
                                const PylithScalar strain[],
                                PylithScalar stress[]) {
        const PylithInt _dim = 2;

        // Auxiliary fields used.
        const PylithInt i_refStress = numA-7;
        const PylithInt i_refStrain = numA-6;
        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_bulkModulus = numA-4;
        const PylithInt i_maxwellTime = numA-3;
        const PylithInt i_viscousStrain = numA-2;
        const PylithInt i_totalStrain = numA-1;

        assert(_dim == dim);
        assert(numA >= 8);
        assert(a);
        assert(aOff);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_maxwellTime] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_totalStrain] >= 0);
        assert(1 == numConstants);
        assert(constants);
        assert(strain);
        assert(stress);

        const PylithReal* refStress = &a[aOff[i_refStress]];
        const PylithReal* refStrain = &a[aOff[i_refStrain]];

        const PylithReal bulkModulus = a[aOff[i_bulkModulus]];
        IsotropicLinearElasticityPlaneStrain::meanStress_refstate(_dim, bulkModulus, refStress, refStrain, strain, stress);

        const PylithReal shearModulus = a[aOff[i_shearModulus]];
        const PylithReal maxwellTime = a[aOff[i_maxwellTime]];
        const PylithReal* viscousStrainPrev = &a[aOff[i_viscousStrain]];
        const PylithReal* totalStrain = &a[aOff[i_totalStrain]];
        const PylithReal dt = constants[0];
        _deviatoricStress_refstate(_dim, shearModulus, maxwellTime, refStress, refStrain, viscousStrainPrev, totalStrain, strain, dt, stress);
    } // cauchyStress_refstate

    // --------------------------------------------------------------------------------------------
    /** Calculate stress vector for 2D plane strain isotropic linear
     * Maxwell WITH a reference stress/strain using state variables.
     *
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                   maxwell_time(1) viscous_strain(4), total_strain(4)]
     */
    static inline
    void _cauchyStress_refstate_stateVars_asVector(const PylithInt dim,
                                                   const PylithInt numA,
                                                   const PylithInt aOff[],
                                                   const PylithScalar a[],
                                                   const PylithScalar strain[],
                                                   PylithScalar stressVector[]) {
        const PylithInt _dim = 2;

        // Incoming auxiliary fields.
        const PylithInt i_refStress = numA-7;
        const PylithInt i_refStrain = numA-6;
        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_bulkModulus = numA-4;
        const PylithInt i_viscousStrain = numA-2;

        assert(numA >= 8);
        assert(aOff);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(stressVector);

        const PylithScalar* refStress = &a[aOff[i_refStress]]; // xx, yy, zz, xy
        const PylithScalar* refStrain = &a[aOff[i_refStrain]]; // xx, yy, zz, xy

        PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 }; // Full stress tensor in vector form.
        const PylithReal bulkModulus = a[aOff[i_bulkModulus]];assert(bulkModulus > 0.0);
        IsotropicLinearElasticityPlaneStrain::meanStress_refstate(_dim, bulkModulus, refStress, refStrain, strain, stressTensor);
        const PylithScalar meanStress = stressTensor[0];

        const PylithReal meanRefStrain = (refStrain[0] + refStrain[1] + refStrain[2]) / 3.0;
        const PylithScalar devRefStrain[4] = {refStrain[0] - meanRefStrain,
                                              refStrain[1] - meanRefStrain,
                                              refStrain[2] - meanRefStrain,
                                              refStrain[3]};

        const PylithReal meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
        const PylithScalar devRefStress[4] = {refStress[0] - meanRefStress,
                                              refStress[1] - meanRefStress,
                                              refStress[2] - meanRefStress,
                                              refStress[3]};

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];assert(shearModulus > 0.0);
        const PylithScalar* viscousStrain = &a[aOff[i_viscousStrain]];
        stressVector[0] = meanStress + devRefStress[0] + 2.0 * shearModulus * (viscousStrain[0] - devRefStrain[0]); // stress_xx
        stressVector[1] = meanStress + devRefStress[1] + 2.0 * shearModulus * (viscousStrain[1] - devRefStrain[1]); // stress_yy
        stressVector[2] = meanStress + devRefStress[2] + 2.0 * shearModulus * (viscousStrain[2] - devRefStrain[2]); // stress_zz
        stressVector[3] = devRefStress[3] + 2.0 * shearModulus * (viscousStrain[3] - devRefStrain[3]); // stress_xy
    } // cauchyStress_refstate_stateVars_asVector

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress tensor for 2D plane strain isotropic linear
     * Maxwell viscoelasticity WITH reference stress and strain.
     */
    static inline
    void _deviatoricStress_refstate(const PylithInt dim,
                                    const PylithReal shearModulus,
                                    const PylithReal maxwellTime,
                                    const PylithReal refStress[],
                                    const PylithReal refStrain[],
                                    const PylithReal viscousStrainPrev[],
                                    const PylithReal strainPrev[],
                                    const PylithScalar strain[],
                                    const PylithReal dt,
                                    PylithScalar stress[]) {
        const PylithInt _dim = 2;

        assert(_dim == dim);
        assert(shearModulus > 0.0);
        assert(maxwellTime > 0.0);
        assert(refStress);
        assert(refStrain);
        assert(viscousStrainPrev);
        assert(strainPrev);
        assert(strain);
        assert(dt > 0.0);
        assert(stress);

        PylithScalar viscousStrain[4] = {0.0, 0.0, 0.0, 0.0}; // Viscous strain vector.
        _viscousStrain_asVector(_dim, maxwellTime, viscousStrainPrev, strainPrev, dt, strain, viscousStrain);

        // Compute reference deviatoric values.
        const PylithReal meanRefStrain = (refStrain[0] + refStrain[1] + refStrain[2]) / 3.0;
        const PylithScalar devRefStrain[4] = {
            refStrain[0] - meanRefStrain,
            refStrain[1] - meanRefStrain,
            refStrain[2] - meanRefStrain,
            refStrain[3],
        };

        const PylithReal meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
        const PylithScalar devRefStress[4] = {
            refStress[0] - meanRefStress,
            refStress[1] - meanRefStress,
            refStress[2] - meanRefStress,
            refStress[3],
        };

        stress[0] += devRefStress[0] + 2.0 * shearModulus * (viscousStrain[0] - devRefStrain[0]); // stress_xx
        stress[1] += devRefStress[3] + 2.0 * shearModulus * (viscousStrain[3] - devRefStrain[3]); // stress_xy
        stress[2] += devRefStress[3] + 2.0 * shearModulus * (viscousStrain[3] - devRefStrain[3]); // stress_yx
        stress[3] += devRefStress[1] + 2.0 * shearModulus * (viscousStrain[1] - devRefStrain[1]); // stress_yy
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate viscous strain vector at t+dt for 2D plane strain isotropic linear
     * Maxwell viscoelasticity.
     *
     * Note: viscousStrainPrev, totalStrain, and viscousStrain are all stored as vectors.
     */
    static inline
    void _viscousStrain_asVector(const PylithInt dim,
                                 const PylithReal maxwellTime,
                                 const PylithReal viscousStrainPrev[],
                                 const PylithReal totalStrain[],
                                 const PylithReal dt,
                                 const PylithScalar strainTensor[],
                                 PylithScalar viscousStrain[]) {
        const PylithInt _dim = 2;

        assert(_dim == dim);
        assert(maxwellTime > 0.0);
        assert(viscousStrainPrev);
        assert(totalStrain);
        assert(viscousStrain);

        const PylithReal meanStrain = (strainTensor[0] + strainTensor[3]) / 3.0;
        const PylithScalar devStrain[4] = {
            strainTensor[0] - meanStrain,
            strainTensor[3] - meanStrain,
            0.0 - meanStrain,
            strainTensor[1]
        };

        const PylithReal meanTotalStrain = (totalStrain[0] + totalStrain[1]) / 3.0;
        const PylithScalar devTotalStrain[4] = {
            totalStrain[0] - meanTotalStrain,
            totalStrain[1] - meanTotalStrain,
            totalStrain[2] - meanTotalStrain,
            totalStrain[3]
        };

        const PylithScalar dq = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime);
        const PylithScalar expFac = exp(-dt/maxwellTime);
        for (int iComp = 0; iComp < 4; ++iComp) {
            viscousStrain[iComp] = expFac * viscousStrainPrev[iComp] + dq * (devStrain[iComp] - devTotalStrain[iComp]);
        } // for
    }

}; // IsotropicLinearMaxwellPlaneStrain

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linear Maxwell in 3D.
class pylith::fekernels::IsotropicLinearMaxwell3D {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    /** f1 function for isotropic linear Maxwell 3D WITHOUT reference stress and reference strain.
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
        pylith::fekernels::Elasticity3D::f1v(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            pylith::fekernels::Elasticity3D::infinitesimalStrain, _cauchyStress,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** f1 function for isotropic linear Maxwell 3D WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., reference_stress(6), reference_strain(6), shear_modulus(1), bulk_modulus(1),
     *                    maxwell_time(1), total_strain(6), viscous_strain(6)]
     */
    static inline
    void f1v_infinitesimalStrain_refstate(const PylithInt dim,
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
        pylith::fekernels::Elasticity3D::f1v(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            pylith::fekernels::Elasticity3D::infinitesimalStrain, _cauchyStress_refstate,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate stress for 3D isotropic linear Maxwell viscoelasticity WITHOUT a reference stress and strain.
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
        assert(stressVector);

        PylithScalar strainTensor[9] = { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
        Elasticity3D::infinitesimalStrain(dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);
        _cauchyStress_stateVars_asVector(dim, numA, aOff, a, strainTensor, stressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate stress for 3D isotropic linear Maxwell WITH a reference stress and strain.
     *
     * Used in output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void cauchyStress_infinitesimalStrain_refstate_asVector(const PylithInt dim,
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
        assert(stressVector);

        PylithScalar strainTensor[9] = { 0.0, 0.0, 0.0,  0.0, 0, 0.0,  0.0, 0.0, 0.0 };
        Elasticity3D::infinitesimalStrain(dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);
        _cauchyStress_refstate_stateVars_asVector(dim, numA, aOff, a, strainTensor, stressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate viscous strain for 3D isotropic linear elasticity.
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
        assert(viscousStrain);

        PylithScalar strainTensor[9] = { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
        Elasticity3D::infinitesimalStrain(dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);

        // Incoming auxiliary fields.
        const PylithInt i_maxwellTime = numA-3;
        const PylithInt i_viscousStrain = numA-2;
        const PylithInt i_totalStrain = numA-1;

        assert(numA >= 6);
        assert(aOff);
        assert(aOff[i_maxwellTime] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_totalStrain] >= 0);
        assert(a);
        assert(1 == numConstants);
        assert(constants);

        const PylithReal maxwellTime = a[aOff[i_maxwellTime]];
        const PylithReal* viscousStrainPrev = &a[aOff[i_viscousStrain]];
        const PylithReal* totalStrain = &a[aOff[i_totalStrain]];
        const PylithReal dt = constants[0];
        _viscousStrain_asVector(dim, maxwellTime, viscousStrainPrev, totalStrain, dt, strainTensor, viscousStrain);
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
        const PylithInt _dim = 3;

        // Incoming auxiliary fields.
        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_bulkModulus = numA-4;
        const PylithInt i_maxwellTime = numA-3;

        assert(_dim == dim);
        assert(numA >= 5);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_maxwellTime] >= 0);
        assert(a);
        assert(Jf3);
        assert(numConstants == 1);
        assert(constants);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];
        const PylithScalar dt = constants[0];

        const PylithScalar dq = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime);

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

private:

    // --------------------------------------------------------------------------------------------
    /** Calculate stress tensor for 3D isotropic linear elasticity
     * WITHOUT a reference stress and strain.
     *
     * ISA Elasticity::stressfn_type
     *
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(1),
     * viscous_strain(4), total_strain(4)]
     */
    static inline
    void _cauchyStress(const PylithInt dim,
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
                       const PylithScalar strain[],
                       PylithScalar stress[]) {
        const PylithInt _dim = 3;

        // Auxiliary fields used.
        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_bulkModulus = numA-4;
        const PylithInt i_maxwellTime = numA-3;
        const PylithInt i_viscousStrain = numA-2;
        const PylithInt i_totalStrain = numA-1;

        assert(_dim == dim);
        assert(numA >= 6);
        assert(a);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_maxwellTime] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_totalStrain] >= 0);
        assert(1 == numConstants);
        assert(constants);
        assert(strain);
        assert(stress);

        const PylithReal bulkModulus = a[aOff[i_bulkModulus]];
        IsotropicLinearElasticity3D::meanStress(_dim, bulkModulus, strain, stress);

        const PylithReal shearModulus = a[aOff[i_shearModulus]];
        const PylithReal maxwellTime = a[aOff[i_maxwellTime]];
        const PylithReal* viscousStrainPrev = &a[aOff[i_viscousStrain]];
        const PylithReal* totalStrain = &a[aOff[i_totalStrain]];
        const PylithReal dt = constants[0];
        _deviatoricStress(_dim, shearModulus, maxwellTime, viscousStrainPrev, totalStrain, strain, dt, stress);
    } // cauchyStress

    // --------------------------------------------------------------------------------------------
    /** Calculate stress vector for 3D isotropic linear elasticity
     * WITHOUT a reference stress and strain while making use of state variables.
     *
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void _cauchyStress_stateVars_asVector(const PylithInt dim,
                                          const PylithInt numA,
                                          const PylithInt aOff[],
                                          const PylithScalar a[],
                                          const PylithScalar strain[],
                                          PylithScalar stressVector[]) {
        const PylithInt _dim = 3;

        // Incoming auxiliary fields.
        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_bulkModulus = numA-4;
        const PylithInt i_viscousStrain = numA-2;

        assert(numA >= 6);
        assert(aOff);

        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_viscousStrain] >= 0);

        PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 }; // Stress tensor as vector.
        const PylithReal bulkModulus = a[aOff[i_bulkModulus]];
        IsotropicLinearElasticity3D::meanStress(_dim, bulkModulus, strain, stressTensor);
        const PylithScalar meanStress = stressTensor[0];

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar* viscousStrain = &a[aOff[i_viscousStrain]];
        stressVector[0] = meanStress + 2.0 * shearModulus * viscousStrain[0]; // stress_xx
        stressVector[1] = meanStress + 2.0 * shearModulus * viscousStrain[1]; // stress_yy
        stressVector[2] = meanStress + 2.0 * shearModulus * viscousStrain[2]; // stress_zz
        stressVector[3] = 2.0 * shearModulus * viscousStrain[3]; // stress_xy
        stressVector[4] = 2.0 * shearModulus * viscousStrain[4]; // stress_yz
        stressVector[5] = 2.0 * shearModulus * viscousStrain[5]; // stress_xz
    } // cauchyStress_stateVars_asVector

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress tensor for 3D isotropic linear
     * Maxwell viscoelasticity WITHOUT reference stress and strain.
     */
    static inline
    void _deviatoricStress(const PylithInt dim,
                           const PylithReal shearModulus,
                           const PylithReal maxwellTime,
                           const PylithReal viscousStrainPrev[],
                           const PylithReal strainPrev[],
                           const PylithScalar strain[],
                           const PylithReal dt,
                           PylithScalar stress[]) {
        const PylithInt _dim = 3;

        assert(_dim == dim);
        assert(shearModulus > 0.0);
        assert(maxwellTime > 0.0);
        assert(viscousStrainPrev);
        assert(strainPrev);
        assert(strain);
        assert(dt > 0.0);
        assert(stress);

        PylithScalar viscousStrain[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Viscous strain vector.
        _viscousStrain_asVector(_dim, maxwellTime, viscousStrainPrev, strainPrev, dt, strain, viscousStrain);

        stress[0] += 2.0 * shearModulus * viscousStrain[0]; // xx
        stress[1] += 2.0 * shearModulus * viscousStrain[3]; // xy
        stress[2] += 2.0 * shearModulus * viscousStrain[5]; // xz
        stress[3] += 2.0 * shearModulus * viscousStrain[3]; // yx
        stress[4] += 2.0 * shearModulus * viscousStrain[1]; // yy
        stress[5] += 2.0 * shearModulus * viscousStrain[4]; // yz
        stress[6] += 2.0 * shearModulus * viscousStrain[5]; // zx
        stress[7] += 2.0 * shearModulus * viscousStrain[4]; // zy
        stress[8] += 2.0 * shearModulus * viscousStrain[2]; // zx
    } // deviatoricStress

    // --------------------------------------------------------------------------------------------
    /** Calculate stress tensor for 3D isotropic linear
     * Maxwell WITH a reference stress/strain.
     *
     * ISA Elasticity::stressfn_type
     *
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                   maxwell_time(1), viscous_strain(4), total_strain(4)]
     */
    static inline
    void _cauchyStress_refstate(const PylithInt dim,
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
                                const PylithScalar strain[],
                                PylithScalar stress[]) {
        const PylithInt _dim = 3;

        // Auxiliary fields used.
        const PylithInt i_refStress = numA-7;
        const PylithInt i_refStrain = numA-6;
        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_bulkModulus = numA-4;
        const PylithInt i_maxwellTime = numA-3;
        const PylithInt i_viscousStrain = numA-2;
        const PylithInt i_totalStrain = numA-1;

        assert(_dim == dim);
        assert(numA >= 8);
        assert(a);
        assert(aOff);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_maxwellTime] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(aOff[i_totalStrain] >= 0);
        assert(1 == numConstants);
        assert(constants);
        assert(strain);
        assert(stress);

        const PylithReal* refStress = &a[aOff[i_refStress]];
        const PylithReal* refStrain = &a[aOff[i_refStrain]];

        const PylithReal bulkModulus = a[aOff[i_bulkModulus]];
        IsotropicLinearElasticity3D::meanStress_refstate(_dim, bulkModulus, refStress, refStrain, strain, stress);

        const PylithReal shearModulus = a[aOff[i_shearModulus]];
        const PylithReal maxwellTime = a[aOff[i_maxwellTime]];
        const PylithReal* viscousStrainPrev = &a[aOff[i_viscousStrain]];
        const PylithReal* totalStrain = &a[aOff[i_totalStrain]];
        const PylithReal dt = constants[0];
        _deviatoricStress_refstate(_dim, shearModulus, maxwellTime, refStress, refStrain, viscousStrainPrev, totalStrain, strain, dt, stress);
    } // cauchyStress_refstate

    // --------------------------------------------------------------------------------------------
    /** Calculate stress vector for 3D isotropic linear
     * Maxwell WITH a reference stress/strain using state variables.
     *
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                   maxwell_time(1) viscous_strain(4), total_strain(4)]
     */
    static inline
    void _cauchyStress_refstate_stateVars_asVector(const PylithInt dim,
                                                   const PylithInt numA,
                                                   const PylithInt aOff[],
                                                   const PylithScalar a[],
                                                   const PylithScalar strain[],
                                                   PylithScalar stressVector[]) {
        const PylithInt _dim = 3;

        // Incoming auxiliary fields.
        const PylithInt i_refStress = numA-7;
        const PylithInt i_refStrain = numA-6;
        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_bulkModulus = numA-4;
        const PylithInt i_viscousStrain = numA-2;

        assert(numA >= 8);
        assert(aOff);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_viscousStrain] >= 0);
        assert(stressVector);

        const PylithScalar* refStress = &a[aOff[i_refStress]]; // xx, yy, zz, xy, yz, xz
        const PylithScalar* refStrain = &a[aOff[i_refStrain]]; // xx, yy, zz, xy, yz, xz

        PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
        const PylithReal bulkModulus = a[aOff[i_bulkModulus]];assert(bulkModulus > 0.0);
        IsotropicLinearElasticity3D::meanStress_refstate(_dim, bulkModulus, refStress, refStrain, strain, stressTensor);
        const PylithScalar meanStress = stressTensor[0];

        const PylithReal meanRefStrain = (refStrain[0] + refStrain[1] + refStrain[2]) / 3.0;
        const PylithScalar devRefStrain[6] = {
            refStrain[0] - meanRefStrain,
            refStrain[1] - meanRefStrain,
            refStrain[2] - meanRefStrain,
            refStrain[3],
            refStrain[4],
            refStrain[5],
        };

        const PylithReal meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
        const PylithScalar devRefStress[6] = {
            refStress[0] - meanRefStress,
            refStress[1] - meanRefStress,
            refStress[2] - meanRefStress,
            refStress[3],
            refStress[4],
            refStress[5],
        };

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar* viscousStrain = &a[aOff[i_viscousStrain]];
        stressVector[0] = meanStress + devRefStress[0] + 2.0 * shearModulus * (viscousStrain[0] - devRefStrain[0]); // xx
        stressVector[1] = meanStress + devRefStress[1] + 2.0 * shearModulus * (viscousStrain[1] - devRefStrain[1]); // yy
        stressVector[2] = meanStress + devRefStress[2] + 2.0 * shearModulus * (viscousStrain[2] - devRefStrain[2]); // zz
        stressVector[3] = devRefStress[3] + 2.0 * shearModulus * (viscousStrain[3] - devRefStrain[3]); // xy
        stressVector[4] = devRefStress[4] + 2.0 * shearModulus * (viscousStrain[4] - devRefStrain[4]); // yz
        stressVector[5] = devRefStress[5] + 2.0 * shearModulus * (viscousStrain[4] - devRefStrain[5]); // xz
    } // cauchyStress_refstate_stateVars_asVector

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress for 3D isotropic linear
     * Maxwell viscoelasticity WITH reference stress and strain.
     */
    static inline
    void _deviatoricStress_refstate(const PylithInt dim,
                                    const PylithReal shearModulus,
                                    const PylithReal maxwellTime,
                                    const PylithReal refStress[],
                                    const PylithReal refStrain[],
                                    const PylithReal viscousStrainPrev[],
                                    const PylithReal strainPrev[],
                                    const PylithScalar strain[],
                                    const PylithReal dt,
                                    PylithScalar stress[]) {
        const PylithInt _dim = 3;

        assert(_dim == dim);
        assert(shearModulus > 0.0);
        assert(maxwellTime > 0.0);
        assert(refStress);
        assert(refStrain);
        assert(viscousStrainPrev);
        assert(strainPrev);
        assert(strain);
        assert(dt > 0.0);
        assert(stress);

        PylithScalar viscousStrain[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Viscous strain vector.
        _viscousStrain_asVector(_dim, maxwellTime, viscousStrainPrev, strainPrev, dt, strain, viscousStrain);

        // Compute reference deviatoric values.
        const PylithReal meanRefStrain = (refStrain[0] + refStrain[1] + refStrain[2]) / 3.0;
        const PylithScalar devRefStrain[6] = {
            refStrain[0] - meanRefStrain,
            refStrain[1] - meanRefStrain,
            refStrain[2] - meanRefStrain,
            refStrain[3],
            refStrain[4],
            refStrain[5],
        };

        const PylithReal meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
        const PylithScalar devRefStress[6] = {
            refStress[0] - meanRefStress,
            refStress[1] - meanRefStress,
            refStress[2] - meanRefStress,
            refStress[3],
            refStress[4],
            refStress[5],
        };

        stress[0] += devRefStress[0] + 2.0 * shearModulus * (viscousStrain[0] - devRefStrain[0]); // stress_xx
        stress[1] += devRefStress[3] + 2.0 * shearModulus * (viscousStrain[3] - devRefStrain[3]); // stress_xy
        stress[2] += devRefStress[5] + 2.0 * shearModulus * (viscousStrain[5] - devRefStrain[5]); // stress_xz
        stress[3] += devRefStress[3] + 2.0 * shearModulus * (viscousStrain[3] - devRefStrain[3]); // stress_yx
        stress[4] += devRefStress[1] + 2.0 * shearModulus * (viscousStrain[1] - devRefStrain[1]); // stress_yy
        stress[5] += devRefStress[4] + 2.0 * shearModulus * (viscousStrain[4] - devRefStrain[4]); // stress_yz
        stress[6] += devRefStress[5] + 2.0 * shearModulus * (viscousStrain[5] - devRefStrain[5]); // stress_zx
        stress[7] += devRefStress[4] + 2.0 * shearModulus * (viscousStrain[4] - devRefStrain[4]); // stress_zy
        stress[8] += devRefStress[2] + 2.0 * shearModulus * (viscousStrain[2] - devRefStrain[2]); // stress_yzz
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate viscous strain at t+dt for 3D isotropic linear
     * Maxwell viscoelasticity.
     */
    static inline
    void _viscousStrain_asVector(const PylithInt dim,
                                 const PylithReal maxwellTime,
                                 const PylithReal viscousStrainPrev[],
                                 const PylithReal totalStrain[],
                                 const PylithReal dt,
                                 const PylithScalar strainTensor[],
                                 PylithScalar viscousStrain[]) {
        const PylithInt _dim = 3;

        assert(_dim == dim);
        assert(maxwellTime > 0.0);
        assert(viscousStrainPrev);
        assert(totalStrain);
        assert(viscousStrain);

        const PylithReal meanStrain = (strainTensor[0] + strainTensor[4] + strainTensor[8]) / 3.0;
        const PylithScalar devStrain[6] = {
            strainTensor[0] - meanStrain,
            strainTensor[4] - meanStrain,
            strainTensor[8] - meanStrain,
            strainTensor[1],
            strainTensor[5],
            strainTensor[2],
        };

        const PylithReal meanTotalStrain = (totalStrain[0] + totalStrain[1] + totalStrain[2]) / 3.0;
        const PylithScalar devTotalStrain[6] = {
            totalStrain[0] - meanTotalStrain,
            totalStrain[1] - meanTotalStrain,
            totalStrain[2] - meanTotalStrain,
            totalStrain[3],
            totalStrain[4],
            totalStrain[5],
        };

        const PylithScalar dq = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime);
        const PylithScalar expFac = exp(-dt/maxwellTime);
        for (int iComp = 0; iComp < 6; ++iComp) {
            viscousStrain[iComp] = expFac * viscousStrainPrev[iComp] + dq * (devStrain[iComp] - devTotalStrain[iComp]);
        } // for
    }

}; // IsotropicLinearMaxwell3D

#endif // pylith_fekernels_isotropiclinearmaxwell_hh

// End of file
