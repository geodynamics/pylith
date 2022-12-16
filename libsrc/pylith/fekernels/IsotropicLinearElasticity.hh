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

/** @file libsrc/fekernels/IsotropicLinearElasticity.hh
 *
 * Kernels for isotropic, linear elasticity.
 *
 * Solution fields: [disp(dim), vel(dim, optional)]
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

#if !defined(pylith_fekernels_isotropiclinearelasticity_hh)
#define pylith_fekernels_isotropiclinearelasticity_hh

#include "fekernelsfwd.hh" // forward declarations
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linearly elasticity plane strain.
class pylith::fekernels::IsotropicLinearElasticityPlaneStrain {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    /** f1 function for isotropic linear elasticity plane strain WITHOUT reference stress and reference strain.
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
        pylith::fekernels::ElasticityPlaneStrain::f1v(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain, _cauchyStress,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** f1 function for isotropic linear elasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
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

    // --------------------------------------------------------------------------------------------
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
        pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain(
            dim, numS, sOff, sOff_x, s, s_t, s_x, x, strainTensor);

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

    // --------------------------------------------------------------------------------------------
    /** Jf3_vu entry function for 2-D plane strain isotropic linear elasticity.
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
        const PylithInt _dim = 2;

        // Incoming auxiliary fields.
        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 2);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(a);
        assert(Jf3);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

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

    // --------------------------------------------------------------------------------------------
    /** Calculate mean stress for 2-D plane strain isotropic linear
     * elasticity WITHOUT reference stress and reference strain.
     */
    static inline
    void meanStress(const PylithInt dim,
                    const PylithReal bulkModulus,
                    const PylithScalar strain[],
                    PylithScalar stress[]) {
        const PylithInt _dim = 2;

        assert(_dim == dim);
        assert(bulkModulus > 0.0);
        assert(strain);
        assert(stress);

        const PylithReal strainTrace = strain[0*_dim+0] + strain[1*_dim+1];
        const PylithReal meanStress = bulkModulus * strainTrace;

        for (PylithInt i = 0; i < _dim; ++i) {
            stress[i*_dim+i] += meanStress;
        } // for
    } // meanStress

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress for 2-D plane strain isotropic linear
     * elasticity WITHOUT reference stress and strain.
     */
    static inline
    void deviatoricStress(const PylithInt dim,
                          const PylithReal shearModulus,
                          const PylithScalar strain[],
                          PylithScalar stress[]) {
        const PylithInt _dim = 2;

        assert(_dim == dim);
        assert(shearModulus > 0.0);
        assert(stress);

        const PylithReal strainTrace = strain[0*_dim+0] + strain[1*_dim+1];
        const PylithReal traceTerm = -2.0/3.0*shearModulus * strainTrace;

        const PylithScalar stress_xx = 2.0*shearModulus*strain[0*_dim+0] + traceTerm;
        const PylithScalar stress_yy = 2.0*shearModulus*strain[1*_dim+1] + traceTerm;
        const PylithScalar stress_xy = 2.0*shearModulus*strain[0*_dim+1];

        stress[0*_dim+0] += stress_xx;
        stress[1*_dim+1] += stress_yy;
        stress[0*_dim+1] += stress_xy;
        stress[1*_dim+0] += stress_xy;
    } // deviatoricStress

    // --------------------------------------------------------------------------------------------
    /** Calculate mean stress for 2-D plane strain isotropic linear
     * elasticity WITH reference stress and reference strain.
     */
    static inline
    void meanStress_refstate(const PylithInt dim,
                             const PylithReal bulkModulus,
                             const PylithReal refStress[],
                             const PylithReal refStrain[],
                             const PylithScalar strain[],
                             PylithScalar stress[]) {
        const PylithInt _dim = 2;

        assert(_dim == dim);
        assert(bulkModulus > 0.0);
        assert(refStress);
        assert(refStrain);
        assert(strain);
        assert(stress);

        const PylithReal strainTrace = strain[0*_dim+0] + strain[1*_dim+1];
        const PylithReal refStrainTrace = refStrain[0] + refStrain[1] + refStrain[2];

        const PylithReal meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
        const PylithReal meanStress = meanRefStress + bulkModulus * (strainTrace - refStrainTrace);

        for (PylithInt i = 0; i < _dim; ++i) {
            stress[i*_dim+i] += meanStress;
        } // for
    } // meanStress_refstate

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress for 2-D plane strain isotropic linear
     * elasticity WITH reference stress and strain.
     *
     * devStress_ij = stress_ij - meanStress*delta_ij
     *
     * i==j
     * devStress_ii = refstress_ii - meanRefstress
     *  + 2*shearModulus*(strain_ii - refStrain_ii) - 2/3*shearModulus*(strain_kk - refstrain_kk)
     *
     * i!=j
     * devStress_ij = refstress_ij + 2*shearModulus*(strain_ij - refstrain_ij)
     */
    static inline
    void deviatoricStress_refstate(const PylithInt dim,
                                   const PylithReal shearModulus,
                                   const PylithReal refStress[],
                                   const PylithReal refStrain[],
                                   const PylithScalar strain[],
                                   PylithScalar stress[]) {
        const PylithInt _dim = 2;

        assert(_dim == dim);
        assert(shearModulus > 0.0);
        assert(refStress);
        assert(refStrain);
        assert(strain);
        assert(stress);

        const PylithReal strainTrace = strain[0*_dim+0] + strain[1*_dim+1];
        const PylithReal refStrainTrace = refStrain[0] + refStrain[1] + refStrain[2];
        const PylithReal meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
        const PylithReal traceTerm = -2.0/3.0*shearModulus * (strainTrace - refStrainTrace);

        const PylithScalar stress_xx = refStress[0] - meanRefStress + 2.0*shearModulus*(strain[0*_dim+0]-refStrain[0]) + traceTerm;
        const PylithScalar stress_yy = refStress[1] - meanRefStress + 2.0*shearModulus*(strain[1*_dim+1]-refStrain[1]) + traceTerm;
        const PylithScalar stress_xy = refStress[3] + 2.0*shearModulus * (strain[0*_dim+1] - refStrain[3]);

        stress[0*_dim+0] += stress_xx;
        stress[1*_dim+1] += stress_yy;
        stress[0*_dim+1] += stress_xy;
        stress[1*_dim+0] += stress_xy;
    } // deviatoricStress_refstate

private:

    // --------------------------------------------------------------------------------------------
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

        const PylithReal bulkModulus = a[aOff[i_bulkModulus]];
        meanStress(_dim, bulkModulus, strain, stress);

        const PylithReal shearModulus = a[aOff[i_shearModulus]];
        deviatoricStress(_dim, shearModulus, strain, stress);
    } // cauchyStress

    // --------------------------------------------------------------------------------------------
    /** Calculate stress for 2-D plane strain isotropic linear
     * elasticity WITH a reference stress/strain.
     *
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

        // Incoming auxiliary fields.
        const PylithInt i_rstress = numA-4;
        const PylithInt i_rstrain = numA-3;
        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(_dim == dim);
        assert(numA >= 5);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_rstress] >= 0);
        assert(aOff[i_rstrain] >= 0);
        assert(strain);
        assert(stress);

        const PylithScalar* refStress = &a[aOff[i_rstress]];
        const PylithScalar* refStrain = &a[aOff[i_rstrain]];

        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        meanStress_refstate(_dim, bulkModulus, refStress, refStrain, strain, stress);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        deviatoricStress_refstate(_dim, shearModulus, refStress, refStrain, strain, stress);
    } // cauchyStress_refstate

}; // IsotropicLinearElasticityPlaneStrain

// ------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linearly elasticity in 3D.
class pylith::fekernels::IsotropicLinearElasticity3D {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    /** f1 function for 3D isotropic linear elasticity WITHOUT reference stress and reference strain.
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
        pylith::fekernels::Elasticity3D::f1v(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            pylith::fekernels::Elasticity3D::infinitesimalStrain, _cauchyStress,
            f1);
    }

    // --------------------------------------------------------------------------------------------
    /** f1 function for 3D isotropic linear elasticity WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., refstress(6), refstrain(6), shear_modulus(1), bulk_modulus(1)]
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
        return pylith::fekernels::Elasticity3D::stress_asVector(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
            numConstants, constants,
            pylith::fekernels::Elasticity3D::infinitesimalStrain, _cauchyStress,
            stressVector);
    }

    // --------------------------------------------------------------------------------------------
    /** Calculate stress for 3D isotropic linear elasticity WITH a reference stress and strain.
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
        return pylith::fekernels::Elasticity3D::stress_asVector(dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
                                                                numConstants, constants,
                                                                Elasticity3D::infinitesimalStrain, _cauchyStress_refstate, stressVector);
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
        const PylithInt _dim = 3;

        // Incoming auxiliary fields.
        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(_dim == dim);
        assert(numA >= 2);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(a);
        assert(Jf3);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

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

    // --------------------------------------------------------------------------------------------
    /** Calculate mean stress for 3D isotropic linear
     * elasticity WITHOUT reference stress and reference strain.
     */
    static inline
    void meanStress(const PylithInt dim,
                    const PylithReal bulkModulus,
                    const PylithScalar strain[],
                    PylithScalar stress[]) {
        const PylithInt _dim = 3;

        assert(_dim == dim);
        assert(bulkModulus > 0.0);
        assert(stress);

        const PylithReal strainTrace = strain[0*_dim+0] + strain[1*_dim+1] + strain[2*_dim+2];
        const PylithReal meanStress = bulkModulus * strainTrace;

        for (PylithInt i = 0; i < _dim; ++i) {
            stress[i*_dim+i] += meanStress;
        } // for
    } // meanStress

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress for 3D isotropic linear
     * elasticity WITHOUT reference stress and strain.
     */
    static inline
    void deviatoricStress(const PylithInt dim,
                          const PylithReal shearModulus,
                          const PylithScalar strain[],
                          PylithScalar stress[]) {
        const PylithInt _dim = 3;

            << << <<< HEAD
            // Incoming solution field.
            const PylithInt i_disp = 0;

        // Incoming auxiliary field.
        const PylithInt i_shearModulus = 0;

        assert(_dim == dim);
        assert(1 == numS);
        assert(1 == numA);
        assert(sOff_x);
        assert(sOff_x[i_disp] >= 0);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        == == == =
            assert(_dim == dim);
        assert(shearModulus > 0.0);
        >> >> >>> 4d 311233e (Update linear elasticity.)
        assert(stress);

        const PylithReal strainTrace = strain[0*_dim+0] + strain[1*_dim+1] + strain[2*_dim+2];
        const PylithReal traceTerm = -2.0/3.0*shearModulus * strainTrace;

        const PylithScalar stress_xx = 2.0*shearModulus*strain[0*_dim+0] + traceTerm;
        const PylithScalar stress_yy = 2.0*shearModulus*strain[1*_dim+1] + traceTerm;
        const PylithScalar stress_zz = 2.0*shearModulus*strain[2*_dim+2] + traceTerm;
        const PylithScalar stress_xy = 2.0*shearModulus*strain[0*_dim+1];
        const PylithScalar stress_yz = 2.0*shearModulus*strain[1*_dim+2];
        const PylithScalar stress_xz = 2.0*shearModulus*strain[0*_dim+2];

        stress[0*_dim+0] += stress_xx;
        stress[0*_dim+1] += stress_xy;
        stress[0*_dim+2] += stress_xz;
        stress[1*_dim+0] += stress_xy;
        stress[1*_dim+1] += stress_yy;
        stress[1*_dim+2] += stress_yz;
        stress[2*_dim+0] += stress_xz;
        stress[2*_dim+1] += stress_yz;
        stress[2*_dim+2] += stress_zz;
    } // deviatoricStress

    // --------------------------------------------------------------------------------------------
    /** Calculate mean stress for 3D isotropic linear
     * elasticity WITH reference stress and reference strain.
     */
    static inline
    void meanStress_refstate(const PylithInt dim,
                             const PylithReal bulkModulus,
                             const PylithReal refStress[],
                             const PylithReal refStrain[],
                             const PylithScalar strain[],
                             PylithScalar stress[]) {
        const PylithInt _dim = 3;

        // Incoming auxiliary fields.
        assert(_dim == dim);
        assert(bulkModulus > 0.0);
        assert(refStress);
        assert(refStrain);
        assert(stress);

        const PylithReal strainTrace = strain[0*_dim+0] + strain[1*_dim+1] + strain[2*_dim+2];
        const PylithReal refStrainTrace = refStrain[0] + refStrain[1] + refStrain[2];

        const PylithReal meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
        const PylithReal meanStress = meanRefStress + bulkModulus * (strainTrace - refStrainTrace);

        for (PylithInt i = 0; i < _dim; ++i) {
            stress[i*_dim+i] += meanStress;
        } // for
    } // meanStress_refstate

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress for 3D isotropic linear
     * elasticity WITH reference stress and strain.
     */
    static inline
    void deviatoricStress_refstate(const PylithInt dim,
                                   const PylithReal shearModulus,
                                   const PylithReal refStress[],
                                   const PylithReal refStrain[],
                                   const PylithScalar strain[],
                                   PylithScalar stress[]) {
        const PylithInt _dim = 3;

        assert(_dim == dim);
        assert(shearModulus > 0.0);
        assert(refStress);
        assert(refStrain);
        assert(strain);
        assert(stress);

        const PylithReal strainTrace = strain[0*_dim+0] + strain[1*_dim+1] + strain[2*_dim+2];
        const PylithReal refStrainTrace = refStrain[0] + refStrain[1] + refStrain[2];
        const PylithReal meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
        const PylithReal traceTerm = -2.0/3.0*shearModulus * (strainTrace - refStrainTrace);

        const PylithScalar stress_xx = refStress[0] - meanRefStress + 2.0*shearModulus*(strain[0*_dim+0]-refStrain[0]) + traceTerm;
        const PylithScalar stress_yy = refStress[1] - meanRefStress + 2.0*shearModulus*(strain[1*_dim+1]-refStrain[1]) + traceTerm;
        const PylithScalar stress_zz = refStress[2] - meanRefStress + 2.0*shearModulus*(strain[2*_dim+2]-refStrain[2]) + traceTerm;
        const PylithScalar stress_xy = refStress[3] + 2.0*shearModulus*(strain[0*_dim+1] - refStrain[3]);
        const PylithScalar stress_yz = refStress[4] + 2.0*shearModulus*(strain[1*_dim+2] - refStrain[4]);
        const PylithScalar stress_xz = refStress[5] + 2.0*shearModulus*(strain[0*_dim+2] - refStrain[5]);

        stress[0*_dim+0] += stress_xx;
        stress[1*_dim+1] += stress_yy;
        stress[2*_dim+2] += stress_zz;
        stress[0*_dim+1] += stress_xy;
        stress[1*_dim+0] += stress_xy;
        stress[1*_dim+2] += stress_yz;
        stress[2*_dim+1] += stress_yz;
        stress[0*_dim+2] += stress_xz;
        stress[2*_dim+0] += stress_xz;
    } // deviatoricStress_refstate

private:

    // --------------------------------------------------------------------------------------------
    /** Calculate stress for 3D isotropic linear elasticity WITHOUT a reference stress and strain.
     *
     * Used to output the stress field.
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

        // Incoming auxiliary fields.
        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(_dim == dim);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);

        const PylithReal bulkModulus = a[aOff[i_bulkModulus]];
        meanStress(_dim, bulkModulus, strain, stress);

        const PylithReal shearModulus = a[aOff[i_shearModulus]];
        deviatoricStress(_dim, shearModulus, strain, stress);
    } // cauchyStress

    // --------------------------------------------------------------------------------------------
    /** Calculate stress for 3D isotropic linear elasticity WITH a reference stress/strain.
     *
     * Used to output the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., refstress(6), refstrain(6), shear_modulus(1), bulk_modulus(1)]
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

        // Incoming auxiliary fields.
        const PylithInt i_rstress = numA-4;
        const PylithInt i_rstrain = numA-3;
        const PylithInt i_shearModulus = numA-2;
        const PylithInt i_bulkModulus = numA-1;

        assert(_dim == dim);
        assert(numA >= 5);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_rstress] >= 0);
        assert(aOff[i_rstrain] >= 0);

        const PylithScalar* refStress = &a[aOff[i_rstress]];
        const PylithScalar* refStrain = &a[aOff[i_rstrain]];

        const PylithReal bulkModulus = a[aOff[i_bulkModulus]];
        meanStress_refstate(_dim, bulkModulus, refStress, refStrain, strain, stress);

        const PylithReal shearModulus = a[aOff[i_shearModulus]];
        deviatoricStress_refstate(_dim, shearModulus, refStress, refStrain, strain, stress);
    } // cauchyStress_refstate

}; // IsotropicLinearElasticity3D

#endif // pylith_fekernels_isotropiclinearelasticity_hh

// End of file
