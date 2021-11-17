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
 * Copyright (c) 2010-2021 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "IsotropicLinearElasticity.hh" // Implementation of object methods

#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include <cassert> // USES assert()

// ================================================================================================
// Kernels for isotropic, linear elasticity plane strain.
// ================================================================================================

// ------------------------------------------------------------------------------------------------
// f1 function for isotropic linear elasticity plane strain WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1v(const PylithInt dim,
                                                             const PylithInt numS,
                                                             const PylithInt numA,
                                                             const PylithInt sOff[],
                                                             const PylithInt sOff_x[],
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
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 2);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(f1);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v


// ----------------------------------------------------------------------
// g1 function for isotropic linear elasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1v_refstate(const PylithInt dim,
                                                                      const PylithInt numS,
                                                                      const PylithInt numA,
                                                                      const PylithInt sOff[],
                                                                      const PylithInt sOff_x[],
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
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f1);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v_refstate


// ------------------------------------------------------------------------------------------------
/* Jf3_vu entry function for 2-D plane strain isotropic linear elasticity WITHOUT reference stress and reference strain.
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
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jf3vu(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
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
} // Jg3vu


// ------------------------------------------------------------------------------------------------
// f0 function for negative fault face.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_neg(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff_x);
    assert(aOff_x[i_density] >= 0);
    assert(f0);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);
    const PylithScalar* density_x = &a_x[aOff_x[i_density]];

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);
    for (PylithInt j = 0; j < _dim; ++j) {
        for (PylithInt i = 0; i < _dim; ++i) {
            f0[j] += density_x[i] / (density*density) * stressTensor[i*_dim+j];
        } // for
    } // for
} // f0l_neg


// ------------------------------------------------------------------------------------------------
// f0 function for positive fault face.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_pos(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff_x);
    assert(aOff_x[i_density] >= 0);
    assert(f0);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);
    const PylithScalar* density_x = &a_x[aOff_x[i_density]];

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);
    for (PylithInt j = 0; j < _dim; ++j) {
        for (PylithInt i = 0; i < _dim; ++i) {
            f0[j] -= density_x[i] / (density*density) * stressTensor[i*_dim+j];
        } // for
    } // for
} // f0l_pos


// ------------------------------------------------------------------------------------------------
// f0 function for negative fault face.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_refstate_neg(const PylithInt dim,
                                                                          const PylithInt numS,
                                                                          const PylithInt numA,
                                                                          const PylithInt sOff[],
                                                                          const PylithInt sOff_x[],
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
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f0);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);
    const PylithScalar* density_x = &a_x[aOff_x[i_density]];

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt j = 0; j < _dim; ++j) {
        for (PylithInt i = 0; i < _dim; ++i) {
            f0[j] += density_x[i] / (density*density) * stressTensor[i*_dim+j];
        } // for
    } // for
} // f0l_refstate_neg


// ------------------------------------------------------------------------------------------------
// f0 function for positive fault face.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_refstate_pos(const PylithInt dim,
                                                                          const PylithInt numS,
                                                                          const PylithInt numA,
                                                                          const PylithInt sOff[],
                                                                          const PylithInt sOff_x[],
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
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f0);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);
    const PylithScalar* density_x = &a_x[aOff_x[i_density]];

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt j = 0; j < _dim; ++j) {
        for (PylithInt i = 0; i < _dim; ++i) {
            f0[j] -= density_x[i] / (density*density) * stressTensor[i*_dim+j];
        } // for
    } // for
} // f0l_refstate_pos


// ------------------------------------------------------------------------------------------------
// f1 function for negative fault face.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1l_neg(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
                                                                 PylithScalar f1[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(f1);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] += stressTensor[i] / (density*density);
    } // for
} // f1l_neg


// ------------------------------------------------------------------------------------------------
// f1 function for positive fault face.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1l_pos(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
                                                                 PylithScalar f1[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(f1);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i] / (density*density);
    } // for
} // f1l_pos


// ------------------------------------------------------------------------------------------------
// f1 function with reference state for negative fault face.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1l_refstate_neg(const PylithInt dim,
                                                                          const PylithInt numS,
                                                                          const PylithInt numA,
                                                                          const PylithInt sOff[],
                                                                          const PylithInt sOff_x[],
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
                                                                          PylithScalar f1[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f1);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] += stressTensor[i] / (density*density);
    } // for
} // f1l_refstate_neg


// ------------------------------------------------------------------------------------------------
// f1 function with reference state for positive fault face.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1l_refstate_pos(const PylithInt dim,
                                                                          const PylithInt numS,
                                                                          const PylithInt numA,
                                                                          const PylithInt sOff[],
                                                                          const PylithInt sOff_x[],
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
                                                                          PylithScalar f1[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f1);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i] / (density*density);
    } // for
} // f1l_refstate_pos


// ------------------------------------------------------------------------------------------------
// Jf1lu function for dynamic slip constraint equation for negative fault face.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jf1lu_neg(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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
                                                                   const PylithReal n[],
                                                                   const PylithInt numConstants,
                                                                   const PylithScalar constants[],
                                                                   PylithScalar Jf1[]) {
    const PylithInt _dim = 2;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 2);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff_x);
    assert(aOff_x[i_density] >= 0);
    assert(a);
    assert(Jf1);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);
    const PylithScalar* density_x = &a_x[aOff_x[i_density]];

    const size_t size = 16;
    PylithScalar elasticConstants[size];
    for (size_t i = 0; i < size; ++i) {
        elasticConstants[i] = 0.0;
    } // for
    Jf3vu(_dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, s_tshift, x,
          numConstants, constants, elasticConstants);

    // j(f,g,dg) = rho,df * C(f,df,g,dg)
    const PylithInt dim3 = _dim*_dim*_dim;
    const PylithInt dim2 = _dim*_dim;
    const PylithInt dim1 = _dim;
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            for (PylithInt l = 0; l < _dim; ++l) {
                for (PylithInt k = 0; k < _dim; ++k) {
                    Jf1[i*dim3+j*dim1+l] += density_x[k] / (density*density) * elasticConstants[i*dim3+k*dim2+j*dim1+l];
                } // for`
            } // for`
        } // for`
    } // for`
} // Jf1lu_neg


// ------------------------------------------------------------------------------------------------
// Jf1lu function for dynamic slip constraint equation for positive fault face.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jf1lu_pos(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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
                                                                   const PylithReal n[],
                                                                   const PylithInt numConstants,
                                                                   const PylithScalar constants[],
                                                                   PylithScalar Jf1[]) {
    const PylithInt _dim = 2;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 2);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff_x);
    assert(aOff_x[i_density] >= 0);
    assert(a);
    assert(Jf1);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);
    const PylithScalar* density_x = &a_x[aOff_x[i_density]];

    const size_t size = 16;
    PylithScalar elasticConstants[size];
    for (size_t i = 0; i < size; ++i) {
        elasticConstants[i] = 0.0;
    } // for
    Jf3vu(_dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, s_tshift, x,
          numConstants, constants, elasticConstants);

    // j(f,g,dg) = rho,df * C(f,df,g,dg)
    const PylithInt dim3 = _dim*_dim*_dim;
    const PylithInt dim2 = _dim*_dim;
    const PylithInt dim1 = _dim;
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            for (PylithInt l = 0; l < _dim; ++l) {
                for (PylithInt k = 0; k < _dim; ++k) {
                    Jf1[i*dim3+j*dim1+l] -= density_x[k] / (density*density) * elasticConstants[i*dim3+k*dim2+j*dim1+l];
                } // for`
            } // for`
        } // for`
    } // for`
} // Jf1lu_pos


// ------------------------------------------------------------------------------------------------
// Jf3lu function for dynamic slip constraint equation for negative fault face.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jf3lu_neg(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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
                                                                   const PylithReal n[],
                                                                   const PylithInt numConstants,
                                                                   const PylithScalar constants[],
                                                                   PylithScalar Jf3[]) {
    const PylithInt _dim = 2;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 2);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(a);
    assert(Jf3);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);

    const size_t size = 16;
    PylithScalar elasticConstants[size];
    for (size_t i = 0; i < size; ++i) {
        elasticConstants[i] = 0.0;
    } // for
    Jf3vu(_dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, s_tshift, x,
          numConstants, constants, elasticConstants);
    for (size_t i = 0; i < size; ++i) {
        Jf3[i] += elasticConstants[i] / density;
    } // for
} // Jf3lu_neg


// ------------------------------------------------------------------------------------------------
// Jf3lu function for dynamic slip constraint equation for positive fault face.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jf3lu_pos(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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
                                                                   const PylithReal n[],
                                                                   const PylithInt numConstants,
                                                                   const PylithScalar constants[],
                                                                   PylithScalar Jf3[]) {
    const PylithInt _dim = 2;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 2);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(a);
    assert(Jf3);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);

    const size_t size = 16;
    PylithScalar elasticConstants[size];
    for (size_t i = 0; i < size; ++i) {
        elasticConstants[i] = 0.0;
    } // for
    Jf3vu(_dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, s_tshift, x,
          numConstants, constants, elasticConstants);
    for (size_t i = 0; i < size; ++i) {
        Jf3[i] -= elasticConstants[i] / density;
    } // for
} // Jf3lu_pos


// ------------------------------------------------------------------------------------------------
// Calculate stress vector for 2-D plane strain isotropic linear
// elasticity WITHOUT a reference stress and strain.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::cauchyStress(const PylithInt dim,
                                                                      const PylithInt numS,
                                                                      const PylithInt numA,
                                                                      const PylithInt sOff[],
                                                                      const PylithInt sOff_x[],
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
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(stressVector);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
    meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);

    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar stress_zz = 0.5*lambda/(lambda+shearModulus) * (stressTensor[0*_dim+0] + stressTensor[1*_dim+1]);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stress_zz;
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
} // cauchyStress


// ------------------------------------------------------------------------------------------------
/* Calculate stress vector for 2-D plane strain isotropic linear
 * elasticity WITH a reference stress/strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim), ...]
 * Auxiliary fields: [density(1), ..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
 */
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::cauchyStress_refstate(const PylithInt dim,
                                                                               const PylithInt numS,
                                                                               const PylithInt numA,
                                                                               const PylithInt sOff[],
                                                                               const PylithInt sOff_x[],
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
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 5);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(stressVector);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
    meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);

    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar* rstress = &a[aOff[i_rstress]];
    const PylithScalar stress_zz = rstress[2] +
                                   0.5*lambda/(lambda+shearModulus) *
                                   (stressTensor[0*_dim+0]-rstress[0] + stressTensor[1*_dim+1]-rstress[1]);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stress_zz;
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
} // cauchyStress_refstate


// ------------------------------------------------------------------------------------------------
/* Calculate mean stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT reference stress and strain.
 *
 * meanStress = bulkModulus * strain_kk
 *
 * stress += meanStress * delta_ij
 */
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::meanStress(const PylithInt dim,
                                                                    const PylithInt numS,
                                                                    const PylithInt numA,
                                                                    const PylithInt sOff[],
                                                                    const PylithInt sOff_x[],
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
                                                                    PylithScalar stress[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_bulkModulus = 0;

    assert(_dim == dim);
    assert(1 == numS);
    assert(1 == numA);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_bulkModulus] >= 0);
    assert(stress);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal meanStress = bulkModulus * strainTrace;

    for (PylithInt i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += meanStress;
    } // for
} // meanStress


// ------------------------------------------------------------------------------------------------
/* Calculate mean stress for 2-D plane strain isotropic linear
 * elasticity WITH reference stress and reference strain.
 *
 * We compute the stress relative to a reference stress/strain state.
 *
 * meanStress = meanRefStress + bulkModulus * (strain_kk - refstrain_kk)
 *
 * stress += meanStress * delta_ij
 */
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::meanStress_refstate(const PylithInt dim,
                                                                             const PylithInt numS,
                                                                             const PylithInt numA,
                                                                             const PylithInt sOff[],
                                                                             const PylithInt sOff_x[],
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
                                                                             PylithScalar stress[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = 0;
    const PylithInt i_rstrain = 1;
    const PylithInt i_bulkModulus = 2;

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(stress);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[2];

    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithReal meanStress = meanrstress + bulkModulus * (strainTrace - refstrainTrace);

    for (PylithInt i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += meanStress;
    } // for
} // meanStress_refstate


// ------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * i==j
 * devStress_ii = 2*shearModulus*strain_ii - 2/3*shearModulus*strain_kk
 *
 * i!=j
 * devStress_ij = 2*shearModulus*strain_ij
 */
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::deviatoricStress(const PylithInt dim,
                                                                          const PylithInt numS,
                                                                          const PylithInt numA,
                                                                          const PylithInt sOff[],
                                                                          const PylithInt sOff_x[],
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
                                                                          PylithScalar stress[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;

    assert(_dim == dim);
    assert(1 == numS);
    assert(1 == numA);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(stress);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal traceTerm = -2.0/3.0*shearModulus * strainTrace;

    const PylithScalar stress_xx = 2.0*shearModulus*disp_x[0*_dim+0] + traceTerm;
    const PylithScalar stress_yy = 2.0*shearModulus*disp_x[1*_dim+1] + traceTerm;
    const PylithScalar stress_xy = shearModulus * (disp_x[0*_dim+1] + disp_x[1*_dim+0]);

    stress[0*_dim+0] += stress_xx;
    stress[1*_dim+1] += stress_yy;
    stress[0*_dim+1] += stress_xy;
    stress[1*_dim+0] += stress_xy;
} // deviatoricStress


// ------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITH reference stress and reference strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * i==j
 * devStress_ii = refstress_ii - meanRefstress
 *  + 2*shearModulus*(strain_ii - refstrain_ii) - 2/3*shearModulus*(strain_kk - refstrain_kk)
 *
 * i!=j
 * devStress_ij = refstress_ij + 2*shearModulus*(strain_ij - refstrain_ij)
 */
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::deviatoricStress_refstate(const PylithInt dim,
                                                                                   const PylithInt numS,
                                                                                   const PylithInt numA,
                                                                                   const PylithInt sOff[],
                                                                                   const PylithInt sOff_x[],
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
                                                                                   PylithScalar stress[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = 0;
    const PylithInt i_rstrain = 1;
    const PylithInt i_shearModulus = 2;

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(stress);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[2];
    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithReal traceTerm = -2.0/3.0*shearModulus * (strainTrace - refstrainTrace);

    const PylithScalar stress_xx = refstress[0] - meanrstress + 2.0*shearModulus*(disp_x[0*_dim+0]-refstrain[0]) + traceTerm;
    const PylithScalar stress_yy = refstress[1] - meanrstress + 2.0*shearModulus*(disp_x[1*_dim+1]-refstrain[1]) + traceTerm;
    const PylithScalar stress_xy = refstress[3] + 2.0*shearModulus * (0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]) - refstrain[3]);

    stress[0*_dim+0] += stress_xx;
    stress[1*_dim+1] += stress_yy;
    stress[0*_dim+1] += stress_xy;
    stress[1*_dim+0] += stress_xy;
} // deviatoricStress_refstate


// ================================================================================================
// Kernels for isotropic, linear elatsicity in 3D.
// ================================================================================================
// ------------------------------------------------------------------------------------------------
// f1 function for isotropic linear elasticity 3D WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearElasticity3D::f1v(const PylithInt dim,
                                                    const PylithInt numS,
                                                    const PylithInt numA,
                                                    const PylithInt sOff[],
                                                    const PylithInt sOff_x[],
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
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 2);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(f1);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v


// ------------------------------------------------------------------------------------------------
// g1 function for isotropic linear elasticity 3D with reference stress and strain.
void
pylith::fekernels::IsotropicLinearElasticity3D::f1v_refstate(const PylithInt dim,
                                                             const PylithInt numS,
                                                             const PylithInt numA,
                                                             const PylithInt sOff[],
                                                             const PylithInt sOff_x[],
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
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f1);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v_refstate


// ------------------------------------------------------------------------------------------------
/* Jf3_vu entry function for 3-D isotropic linear elasticity WITHOUT reference stress and reference strain.
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
void
pylith::fekernels::IsotropicLinearElasticity3D::Jf3vu(const PylithInt dim,
                                                      const PylithInt numS,
                                                      const PylithInt numA,
                                                      const PylithInt sOff[],
                                                      const PylithInt sOff_x[],
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
} // Jg3vu


// ------------------------------------------------------------------------------------------------
// f0 function for negative fault face.
void
pylith::fekernels::IsotropicLinearElasticity3D::f0l_neg(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
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
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff_x);
    assert(aOff_x[i_density] >= 0);
    assert(f0);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);
    const PylithScalar* density_x = &a_x[aOff_x[i_density]];

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);
    for (PylithInt j = 0; j < _dim; ++j) {
        for (PylithInt i = 0; i < _dim; ++i) {
            f0[j] += density_x[i] / (density*density) * stressTensor[i*_dim+j];
        } // for
    } // for
} // f0l_neg


// ------------------------------------------------------------------------------------------------
// f0 function for positive fault face.
void
pylith::fekernels::IsotropicLinearElasticity3D::f0l_pos(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
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
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff_x);
    assert(aOff_x[i_density] >= 0);
    assert(f0);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);
    const PylithScalar* density_x = &a_x[aOff_x[i_density]];

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);
    for (PylithInt j = 0; j < _dim; ++j) {
        for (PylithInt i = 0; i < _dim; ++i) {
            f0[j] -= density_x[i] / (density*density) * stressTensor[i*_dim+j];
        } // for
    } // for
} // f0l_pos


// ------------------------------------------------------------------------------------------------
// f0 function for negative fault face.
void
pylith::fekernels::IsotropicLinearElasticity3D::f0l_refstate_neg(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f0);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);
    const PylithScalar* density_x = &a_x[aOff_x[i_density]];

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt j = 0; j < _dim; ++j) {
        for (PylithInt i = 0; i < _dim; ++i) {
            f0[j] += density_x[i] / (density*density) * stressTensor[i*_dim+j];
        } // for
    } // for
} // f0l_refstate_neg


// ------------------------------------------------------------------------------------------------
// f0 function for positive fault face.
void
pylith::fekernels::IsotropicLinearElasticity3D::f0l_refstate_pos(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f0);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);
    const PylithScalar* density_x = &a_x[aOff_x[i_density]];

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt j = 0; j < _dim; ++j) {
        for (PylithInt i = 0; i < _dim; ++i) {
            f0[j] -= density_x[i] / (density*density) * stressTensor[i*_dim+j];
        } // for
    } // for
} // f0l_refstate_pos


// ------------------------------------------------------------------------------------------------
// f1 function for negative fault face.
void
pylith::fekernels::IsotropicLinearElasticity3D::f1l_neg(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
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
                                                        PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(f1);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] += stressTensor[i] / density;
    } // for
} // f1l_neg


// ------------------------------------------------------------------------------------------------
// f1 function for positive fault face.
void
pylith::fekernels::IsotropicLinearElasticity3D::f1l_pos(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
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
                                                        PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(f1);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i] / density;
    } // for
} // f1l_pos


// ------------------------------------------------------------------------------------------------
// f1 function with reference state for negative fault face.
void
pylith::fekernels::IsotropicLinearElasticity3D::f1l_refstate_neg(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
                                                                 PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f1);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] += stressTensor[i] / density;
    } // for
} // f1l_refstate_neg


// ------------------------------------------------------------------------------------------------
// f1 function with reference state for positive fault face.
void
pylith::fekernels::IsotropicLinearElasticity3D::f1l_refstate_pos(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
                                                                 PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f1);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i] / density;
    } // for
} // f1l_refstate_pos


// ------------------------------------------------------------------------------------------------
// Jf1lu function for dynamic slip constraint equation for negative fault face.
void
pylith::fekernels::IsotropicLinearElasticity3D::Jf1lu_neg(const PylithInt dim,
                                                          const PylithInt numS,
                                                          const PylithInt numA,
                                                          const PylithInt sOff[],
                                                          const PylithInt sOff_x[],
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
                                                          const PylithReal n[],
                                                          const PylithInt numConstants,
                                                          const PylithScalar constants[],
                                                          PylithScalar Jf1[]) {
    const PylithInt _dim = 3;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 2);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff_x);
    assert(aOff_x[i_density] >= 0);
    assert(a);
    assert(Jf1);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);
    const PylithScalar* density_x = &a_x[aOff_x[i_density]];

    const size_t size = 32;
    PylithScalar elasticConstants[size];
    for (size_t i = 0; i < size; ++i) {
        elasticConstants[i] = 0.0;
    } // for
    Jf3vu(_dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, s_tshift, x,
          numConstants, constants, elasticConstants);

    // j(f,g,dg) = rho,df * C(f,df,g,dg)
    const PylithInt dim3 = _dim*_dim*_dim;
    const PylithInt dim2 = _dim*_dim;
    const PylithInt dim1 = _dim;
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            for (PylithInt l = 0; l < _dim; ++l) {
                for (PylithInt k = 0; k < _dim; ++k) {
                    Jf1[i*dim3+j*dim1+l] += density_x[k] / (density*density) * elasticConstants[i*dim3+k*dim2+j*dim1+l];
                } // for`
            } // for`
        } // for`
    } // for`
} // Jf1lu_neg


// ------------------------------------------------------------------------------------------------
// Jf1lu function for dynamic slip constraint equation for positive fault face.
void
pylith::fekernels::IsotropicLinearElasticity3D::Jf1lu_pos(const PylithInt dim,
                                                          const PylithInt numS,
                                                          const PylithInt numA,
                                                          const PylithInt sOff[],
                                                          const PylithInt sOff_x[],
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
                                                          const PylithReal n[],
                                                          const PylithInt numConstants,
                                                          const PylithScalar constants[],
                                                          PylithScalar Jf1[]) {
    const PylithInt _dim = 3;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 2);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff_x);
    assert(aOff_x[i_density] >= 0);
    assert(a);
    assert(Jf1);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);
    const PylithScalar* density_x = &a_x[aOff_x[i_density]];

    const size_t size = 32;
    PylithScalar elasticConstants[size];
    for (size_t i = 0; i < size; ++i) {
        elasticConstants[i] = 0.0;
    } // for
    Jf3vu(_dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, s_tshift, x,
          numConstants, constants, elasticConstants);

    // j(f,g,dg) = rho,df * C(f,df,g,dg)
    const PylithInt dim3 = _dim*_dim*_dim;
    const PylithInt dim2 = _dim*_dim;
    const PylithInt dim1 = _dim;
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            for (PylithInt l = 0; l < _dim; ++l) {
                for (PylithInt k = 0; k < _dim; ++k) {
                    Jf1[i*dim3+j*dim1+l] -= density_x[k] / (density*density) * elasticConstants[i*dim3+k*dim2+j*dim1+l];
                } // for`
            } // for`
        } // for`
    } // for`
} // Jf1lu_pos


// ------------------------------------------------------------------------------------------------
// Jf3lu function for dynamic slip constraint equation for negative fault face.
void
pylith::fekernels::IsotropicLinearElasticity3D::Jf3lu_neg(const PylithInt dim,
                                                          const PylithInt numS,
                                                          const PylithInt numA,
                                                          const PylithInt sOff[],
                                                          const PylithInt sOff_x[],
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
                                                          const PylithReal n[],
                                                          const PylithInt numConstants,
                                                          const PylithScalar constants[],
                                                          PylithScalar Jf3[]) {
    const PylithInt _dim = 3;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 2);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(a);
    assert(Jf3);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);

    const size_t size = 32;
    PylithScalar elasticConstants[size];
    for (size_t i = 0; i < size; ++i) {
        elasticConstants[i] = 0.0;
    } // for
    Jf3vu(_dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, s_tshift, x,
          numConstants, constants, elasticConstants);
    for (size_t i = 0; i < size; ++i) {
        Jf3[i] += elasticConstants[i] / density;
    } // for
} // Jf3lu_neg


// ------------------------------------------------------------------------------------------------
// Jf3lu function for dynamic slip constraint equation for positive fault face.
void
pylith::fekernels::IsotropicLinearElasticity3D::Jf3lu_pos(const PylithInt dim,
                                                          const PylithInt numS,
                                                          const PylithInt numA,
                                                          const PylithInt sOff[],
                                                          const PylithInt sOff_x[],
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
                                                          const PylithReal n[],
                                                          const PylithInt numConstants,
                                                          const PylithScalar constants[],
                                                          PylithScalar Jf3[]) {
    const PylithInt _dim = 3;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_dim == dim + 1);
    assert(numS >= 1);
    assert(numA >= 2);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(a);
    assert(Jf3);

    const PylithScalar density = a[aOff[i_density]];assert(density > 0.0);

    const size_t size = 32;
    PylithScalar elasticConstants[size];
    for (size_t i = 0; i < size; ++i) {
        elasticConstants[i] = 0.0;
    } // for
    Jf3vu(_dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, s_tshift, x,
          numConstants, constants, elasticConstants);
    for (size_t i = 0; i < size; ++i) {
        Jf3[i] -= elasticConstants[i] / density;
    } // for
} // Jf3lu_pos


// ------------------------------------------------------------------------------------------------
/* Calculate mean stress for 3-D isotropic linear
 * elasticity WITHOUT reference stress and strain.
 *
 * meanStress = bulkModulus * strain_kk
 *
 * stress += meanStress * delta_ij
 */
void
pylith::fekernels::IsotropicLinearElasticity3D::meanStress(const PylithInt dim,
                                                           const PylithInt numS,
                                                           const PylithInt numA,
                                                           const PylithInt sOff[],
                                                           const PylithInt sOff_x[],
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
                                                           PylithScalar stress[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_disp = 0;

    // Incoming auxiliary field.
    const PylithInt i_bulkModulus = 0;

    assert(_dim == dim);
    assert(1 == numS);
    assert(1 == numA);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_bulkModulus] >= 0);
    assert(stress);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1] + disp_x[2*_dim+2];
    const PylithReal meanStress = bulkModulus * strainTrace;

    for (PylithInt i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += meanStress;
    } // for
} // meanStress


// ------------------------------------------------------------------------------------------------
/* Calculate mean stress for 3-D isotropic linear
 * elasticity WITH reference stress and reference strain.
 *
 * We compute the stress relative to a reference stress/strain state.
 *
 * meanStress = meanRefStress + bulkModulus * (strain_kk - refstrain_kk)
 *
 * stress += meanStress * delta_ij
 */
void
pylith::fekernels::IsotropicLinearElasticity3D::meanStress_refstate(const PylithInt dim,
                                                                    const PylithInt numS,
                                                                    const PylithInt numA,
                                                                    const PylithInt sOff[],
                                                                    const PylithInt sOff_x[],
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
                                                                    PylithScalar stress[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = 0;
    const PylithInt i_rstrain = 1;
    const PylithInt i_bulkModulus = 2;

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(stress);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy, stress_yz,
                                                         // stress_xz
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy, strain_yz,
                                                         // strain_xz
    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1] + disp_x[2*_dim+2];
    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[2];

    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithReal meanStress = meanrstress + bulkModulus * (strainTrace - refstrainTrace);

    for (PylithInt i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += meanStress;
    } // for
} // meanStress_refstate


// ------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 3-D isotropic linear
 * elasticity WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * i==j
 * devStress_ii = 2*shearModulus*strain_ii - 2/3*shearModulus*strain_kk
 *
 * i!=j
 * devStress_ij = 2*shearModulus*strain_ij
 */
void
pylith::fekernels::IsotropicLinearElasticity3D::deviatoricStress(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
                                                                 PylithScalar stress[]) {
    const PylithInt _dim = 3;

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
    assert(stress);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1] + disp_x[2*_dim+2];
    const PylithReal traceTerm = -2.0/3.0*shearModulus * strainTrace;
    const PylithReal twomu = 2.0*shearModulus;

    const PylithScalar stress_xx = twomu*disp_x[0*_dim+0] + traceTerm;
    const PylithScalar stress_yy = twomu*disp_x[1*_dim+1] + traceTerm;
    const PylithScalar stress_zz = twomu*disp_x[2*_dim+2] + traceTerm;
    const PylithScalar stress_xy = shearModulus * (disp_x[0*_dim+1] + disp_x[1*_dim+0]);
    const PylithScalar stress_yz = shearModulus * (disp_x[1*_dim+2] + disp_x[2*_dim+1]);
    const PylithScalar stress_xz = shearModulus * (disp_x[0*_dim+2] + disp_x[2*_dim+0]);

    stress[0*_dim+0] += stress_xx;
    stress[1*_dim+1] += stress_yy;
    stress[2*_dim+2] += stress_zz;
    stress[0*_dim+1] += stress_xy;
    stress[1*_dim+0] += stress_xy;
    stress[2*_dim+1] += stress_yz;
    stress[1*_dim+2] += stress_yz;
    stress[2*_dim+0] += stress_xz;
    stress[0*_dim+2] += stress_xz;
} // deviatoricStress


// ------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 3-D isotropic linear
 * elasticity WITH reference stress and reference strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * i==j
 * devStress_ii = refstress_ii - meanRefstress
 *  + 2*shearModulus*(strain_ii - refstrain_ii) - 2/3*shearModulus*(strain_kk - refstrain_kk)
 *
 * i!=j
 * devStress_ij = refstress_ij + 2*shearModulus*(strain_ij - refstrain_ij)
 */
void
pylith::fekernels::IsotropicLinearElasticity3D::deviatoricStress_refstate(const PylithInt dim,
                                                                          const PylithInt numS,
                                                                          const PylithInt numA,
                                                                          const PylithInt sOff[],
                                                                          const PylithInt sOff_x[],
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
                                                                          PylithScalar stress[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = 0;
    const PylithInt i_rstrain = 1;
    const PylithInt i_shearModulus = 2;

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(stress);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy, stress_yz,
                                                         // stress_xz
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy, strain_yz,
                                                         // strain_xz

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1] + disp_x[2*_dim+2];
    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[2];
    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithReal traceTerm = -2.0/3.0*shearModulus * (strainTrace - refstrainTrace);
    const PylithReal twomu = 2.0*shearModulus;

    const PylithScalar stress_xx = refstress[0] - meanrstress + twomu*(disp_x[0*_dim+0]-refstrain[0]) + traceTerm;
    const PylithScalar stress_yy = refstress[1] - meanrstress + twomu*(disp_x[1*_dim+1]-refstrain[1]) + traceTerm;
    const PylithScalar stress_zz = refstress[2] - meanrstress + twomu*(disp_x[2*_dim+2]-refstrain[2]) + traceTerm;
    const PylithScalar stress_xy = refstress[3] + twomu * (0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]) - refstrain[3]);
    const PylithScalar stress_yz = refstress[4] + twomu * (0.5*(disp_x[1*_dim+2] + disp_x[2*_dim+1]) - refstrain[4]);
    const PylithScalar stress_xz = refstress[5] + twomu * (0.5*(disp_x[0*_dim+2] + disp_x[2*_dim+0]) - refstrain[5]);

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


// ------------------------------------------------------------------------------------------------
/* Calculate stress vector for 3-D isotropic linear
 * elasticity WITHOUT a reference stress and strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim), ...]
 * Auxiliary fields: [density(1), ... shear_modulus(1), bulk_modulus(1)]
 */
void
pylith::fekernels::IsotropicLinearElasticity3D::cauchyStress(const PylithInt dim,
                                                             const PylithInt numS,
                                                             const PylithInt numA,
                                                             const PylithInt sOff[],
                                                             const PylithInt sOff_x[],
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
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(stressVector);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stressTensor[2*_dim+2]; // stress_zz
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
    stressVector[4] = stressTensor[1*_dim+2]; // stress_yz
    stressVector[5] = stressTensor[0*_dim+2]; // stress_xz
} // stress


// ------------------------------------------------------------------------------------------------
/* Calculate stress vector for 3-D isotropic linear
 * elasticity WITH a reference stress/strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim), ...]
 * Auxiliary fields: [density(1), ..., refstress(6), refstrain(6), shear_modulus(1), bulk_modulus(1)]
 */
void
pylith::fekernels::IsotropicLinearElasticity3D::cauchyStress_refstate(const PylithInt dim,
                                                                      const PylithInt numS,
                                                                      const PylithInt numA,
                                                                      const PylithInt sOff[],
                                                                      const PylithInt sOff_x[],
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
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 5);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(stressVector);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stressTensor[2*_dim+2]; // stress_zz
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
    stressVector[4] = stressTensor[1*_dim+2]; // stress_yz
    stressVector[5] = stressTensor[0*_dim+2]; // stress_xz
} // stress_refstate


// End of file
