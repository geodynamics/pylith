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

/** @file libsrc/fekernels/IsotropicLinearIncompElasticityPlaneStrain.hh
 *
 * Kernels for linear incompressible elasticity plane strain WITHOUT inertia.
 *
 * Solution fields: [disp(dim), pressure(1)]
 *
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

#include "pylith/fekernels/IncompressibleElasticity.hh" // USES IncompressibleElasticity kernels
#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
class pylith::fekernels::IsotropicLinearIncompElasticity {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Jf0_pp function for isotropic linear incompressible elasticity .
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf0[]) {
        // Incoming auxiliary subfields
        const PylithInt i_bulkModulus = numA-1;

        assert(numA >= 1);
        assert(aOff);
        assert(aOff[i_bulkModulus] >= 0);
        assert(a);
        assert(Jf0);

        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

        Jf0[0] += 1.0 / bulkModulus;
    } // Jf0pp

    /** f0p function for isotropic linear incompressible elasticity WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void incompressibleTerm(const PylithInt dim,
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
                            PylithScalar* value) {
        // Incoming solution subfields
        const PylithInt i_pres = 1;

        // Incoming auxiliary subfields
        const PylithInt i_bulkModulus = numA-1;

        assert(numS >= 2);
        assert(sOff);
        assert(sOff[i_pres] >= 0);
        assert(s);
        assert(numA >= 1);
        assert(aOff);
        assert(aOff[i_bulkModulus] >= 0);
        assert(a);
        assert(strain);
        assert(value);

        const PylithScalar pressure = s[sOff[i_pres]];

        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

        PylithScalar strainTrace = 0;
        for (PylithInt i = 0; i < dim; ++i) {
            strainTrace += strain[i*dim+i];
        } // for
        *value = strainTrace + pressure / bulkModulus;
    } // incompressibleTerm

    /** f0p function for isotropic linear incompressible elasticity WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void incompressibleTerm_refstate(const PylithInt dim,
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
                                     PylithScalar* value) {
        // Incoming solution subfields
        const PylithInt i_pres = 1;

        // Incoming auxiliary subfields
        const PylithInt i_rstress = numA-4;
        const PylithInt i_rstrain = numA-3;
        const PylithInt i_bulkModulus = numA-1;

        assert(numS >= 2);
        assert(sOff);
        assert(sOff[i_pres] >= 0);
        assert(s);
        assert(s_x);
        assert(numA >= 4);
        assert(aOff);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_rstress] >= 0);
        assert(aOff[i_rstrain] >= 0);
        assert(a);
        assert(strain);
        assert(value);

        const PylithScalar pressure = s[sOff[i_pres]];

        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, ...
        const PylithScalar* refStrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, ...

        const PylithScalar meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
        const PylithScalar refStrainTrace = (refStrain[0] + refStrain[1] + refStrain[2]);

        PylithScalar strainTrace = 0;
        for (PylithInt i = 0; i < dim; ++i) {
            strainTrace += strain[i*dim+i];
        } // for
        *value = (strainTrace - refStrainTrace) + (pressure + meanRefStress) / bulkModulus;
    } // incompressibleTerm_refstate

}; // IsotropicLinearIncompElasticity

// ---------------------------------------------------------------------------------------------------------------------
class pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** f0p function for isotropic linear incompressible elasticity WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f0p_infinitesimalStrain(const PylithInt dim,
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
        pylith::fekernels::IncompressibleElasticityPlaneStrain::f0p(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            ElasticityPlaneStrain::infinitesimalStrain,
            IsotropicLinearIncompElasticity::incompressibleTerm,
            f0);
    } // f0p

    /** f0p function for isotropic linear incompressible elasticity WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f0p_infinitesimalStrain_refstate(const PylithInt dim,
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
        pylith::fekernels::IncompressibleElasticityPlaneStrain::f0p(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            ElasticityPlaneStrain::infinitesimalStrain,
            IsotropicLinearIncompElasticity::incompressibleTerm_refstate,
            f0);
    } // f0p

    /** f1 function for 2-D plane strain isotropic linear incompressible elasticity WITHOUT reference stress and
     * reference strain.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f1u_infinitesimalStrain(const PylithInt dim,
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
            ElasticityPlaneStrain::infinitesimalStrain, _cauchyStress, f1);
    } // f1u_infinitesimalStrain

    /** f1 function for 2-D plane strain isotropic linear incompressible elasticity WITH reference stress and reference
     * strain.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f1u_infinitesimalStrain_refstate(const PylithInt dim,
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
            ElasticityPlaneStrain::infinitesimalStrain, _cauchyStress_refstate, f1);
    }

    /** Jf3_uu entry function for 2-D plane strain isotropic linear incompressible elasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void Jf3uu_infinitesimalStrain(const PylithInt dim,
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
        const PylithInt i_shearModulus = numA-2;

        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 2);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(a);
        assert(Jf3);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];

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

        Jf3[ 0] -= C1111; // j0000
        Jf3[ 3] -= C1212; // j0011
        Jf3[ 5] -= C1122; // j0101
        Jf3[ 6] -= C1212; // j0110, C1221
        Jf3[ 9] -= C1212; // j1001, C2112
        Jf3[10] -= C1122; // j1010, C2211
        Jf3[12] -= C1212; // j1100, C2121
        Jf3[15] -= C2222; // j1111
    } // Jf3uu

    /** Calculate stress for 2D plane strain isotropic linear elasticity WITHOUT a reference stress and strain.
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
        pylith::fekernels::ElasticityPlaneStrain::stress_asVector(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain, _cauchyStress, stressVector);

        // Incoming auxiliary fields.
        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);

        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
        const PylithScalar stress_zz = 0.5*lambda/(lambda+shearModulus) * (stressVector[0] + stressVector[1]);
        stressVector[2] = stress_zz;
    }

    /** Calculate stress for 2D plane strain isotropic linear elasticity WITH a reference stress and strain.
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
        // Incoming auxiliary fields.
        const PylithInt i_rstress = numA-4;
        const PylithInt i_rstrain = numA-3;
        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(numS >= 1);
        assert(numA >= 5);
        assert(sOff);
        assert(sOff_x);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_rstress] >= 0);
        assert(aOff[i_rstrain] >= 0);
        assert(stressVector);

        PylithScalar strainTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
        ElasticityPlaneStrain::infinitesimalStrain(dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);

        PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
        _cauchyStress_refstate(dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
                               t, x, numConstants, constants, strainTensor, stressTensor);

        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
        const PylithScalar* rstress = &a[aOff[i_rstress]];
        const PylithScalar stress_zz = rstress[2] +
                                       0.5*lambda/(lambda+shearModulus) *
                                       (stressTensor[0]-rstress[0] + stressTensor[3]-rstress[1]);

        stressVector[0] = stressTensor[0]; // stress_xx
        stressVector[1] = stressTensor[3]; // stress_yy
        stressVector[2] = stress_zz;
        stressVector[3] = stressTensor[1]; // stress_xy
    } // cauchyStress_infinitesimalStrain_refstate_asVector

    /** Calculate stress for 2D plane strain isotropic linear elasticity WITHOUT a reference stress and strain.
     *
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
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

        // Incoming solution fields.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary fields.
        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(_dim == dim);
        assert(numA >= 3);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(s);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(a);
        assert(strain);
        assert(stress);

        const PylithReal pressure = s[sOff[i_pressure]];
        IncompressibleElasticity::meanStress(_dim, pressure, stress);

        const PylithReal shearModulus = a[aOff[i_shearModulus]];
        IsotropicLinearElasticityPlaneStrain::deviatoricStress(_dim, shearModulus, strain, stress);
    } // cauchyStress

    /** Calculate stress for 2-D plane strain isotropic linear incompressible elasticity WITH a reference stress/strain.
     *
     * Used to output the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
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

        // Incoming solution fields.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary fields.
        const PylithInt i_refStress = numA-4;
        const PylithInt i_refStrain = numA-3;
        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(_dim == dim);
        assert(numA >= 4);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(stress);

        const PylithReal pressure = s[sOff[i_pressure]];
        const PylithReal* refStress = &a[aOff[i_refStress]];
        IncompressibleElasticity::meanStress_refstate(_dim, pressure, refStress, stress);

        const PylithReal shearModulus = a[sOff[i_shearModulus]];
        const PylithReal* refStrain = &a[aOff[i_refStrain]];
        IsotropicLinearElasticityPlaneStrain::deviatoricStress_refstate(_dim, shearModulus, refStress, refStrain, strain, stress);
    } // cauchyStress_refstate

}; // IsotropicLinearIncompElasticityPlaneStrain

// ------------------------------------------------------------------------------------------------
class pylith::fekernels::IsotropicLinearIncompElasticity3D {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    /** f0p function for isotropic linear incompressible elasticity WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f0p_infinitesimalStrain(const PylithInt dim,
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
        pylith::fekernels::IncompressibleElasticity3D::f0p(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            Elasticity3D::infinitesimalStrain,
            IsotropicLinearIncompElasticity::incompressibleTerm,
            f0);
    } // f0p_infinitesimalStrain

    /** f0p function for isotropic linear incompressible elasticity WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f0p_infinitesimalStrain_refstate(const PylithInt dim,
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
        pylith::fekernels::IncompressibleElasticity3D::f0p(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants,
            Elasticity3D::infinitesimalStrain,
            IsotropicLinearIncompElasticity::incompressibleTerm_refstate,
            f0);
    } // f0p_infinitesimalStrain_refstate

    /** f1 function for 3-D isotropic linear incompressible elasticity WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f1u_infinitesimalStrain(const PylithInt dim,
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
            Elasticity3D::infinitesimalStrain, _cauchyStress, f1);
    } // f1u_infinitesimalStrain

    /** f1 function for 3D isotropic linear incompressible elasticity WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void f1u_infinitesimalStrain_refstate(const PylithInt dim,
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
            Elasticity3D::infinitesimalStrain, _cauchyStress_refstate, f1);
    } // f1u_infinitesimalStrain_refstate

    /** Jf3_uu entry function for 3-D isotropic linear incompressible elasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static inline
    void Jf3uu_infinitesimalStrain(const PylithInt dim,
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
        const PylithInt i_shearModulus = numA-2;

        assert(_dim == dim);
        assert(numA >= 2);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(a);
        assert(Jf3);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];

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
        Jf3[ 0] -= C1111; // j0000
        Jf3[ 4] -= C1212; // j0011
        Jf3[ 8] -= C1313; // j0022
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
    } // Jf3uu

    /** Calculate stress for 3D isotropic linear elasticity WITHOUT a reference stress and strain.
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
        pylith::fekernels::Elasticity3D::stress_asVector(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            pylith::fekernels::Elasticity3D::infinitesimalStrain, _cauchyStress, stressVector);
    } // cauchyStress_infinitesimalStrain_asVector

    /** Calculate stress for 3D linear elasticity WITH a reference stress and strain.
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
        pylith::fekernels::Elasticity3D::stress_asVector(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            pylith::fekernels::Elasticity3D::infinitesimalStrain, _cauchyStress_refstate, stressVector);
    } // cauchyStress_infinitesimalStrain_refstate_asVector

private:

    /** Calculate stress for 3D isotropic linear elasticity WITHOUT a reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
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

        // Incoming solution fields.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary fields.
        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(_dim == dim);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(strain);
        assert(stress);

        const PylithReal pressure = s[sOff[i_pressure]];
        IncompressibleElasticity::meanStress(_dim, pressure, stress);

        const PylithReal shearModulus = a[aOff[i_shearModulus]];
        IsotropicLinearElasticity3D::deviatoricStress(_dim, shearModulus, strain, stress);
    } // cauchyStress

    /** Calculate stress for 3D isotropic linear incompressible elasticity WITH a reference stress/strain.
     *
     * Used to output the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
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

        // Incoming solution fields.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary fields.
        const PylithInt i_refStress = numA-4;
        const PylithInt i_refStrain = numA-3;
        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(_dim == dim);
        assert(numA >= 4);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(stress);

        const PylithReal pressure = s[sOff[i_pressure]];
        const PylithReal* refStress = &a[aOff[i_refStress]];
        IncompressibleElasticity::meanStress_refstate(_dim, pressure, refStress, stress);

        const PylithReal shearModulus = a[sOff[i_shearModulus]];
        const PylithReal* refStrain = &a[aOff[i_refStrain]];
        IsotropicLinearElasticity3D::deviatoricStress_refstate(_dim, shearModulus, refStress, refStrain, strain, stress);
    } // cauchyStress_refstate

}; // IsotropicLinearIncompElasticity3D

#endif // pylith_fekernels_isotropiclinearincompelasticity_hh

// End of file
