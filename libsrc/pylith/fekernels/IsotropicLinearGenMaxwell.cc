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

#include "pylith/fekernels/IsotropicLinearGenMaxwell.hh" // Implementation of object methods.

#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/Viscoelasticity.hh" // USES Viscoelasticity kernels

#include <cassert> // USES assert()
#include <cmath> // USES exp()

// =====================================================================================================================
// Kernels for isotropic, linear generalized Maxwell viscoelastic plane strain.
// ====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
// f1 function for isotropic linear generalized Maxwell plane strain material WITHOUT reference
// stress and strain.
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::f1v(const PylithInt dim,
                                                             const PylithInt numS,
                                                             const PylithInt numA,
                                                             const PylithInt sOff[],
                                                             const PylithInt sOff_x[],
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
    const PylithInt i_shearModulus = numA-6;
    const PylithInt i_bulkModulus = numA-5;
    const PylithInt i_maxwellTime = numA-4;
    const PylithInt i_shearModulusRatio = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 6);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(f1);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 5; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[5] = {
        aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_shearModulusRatio], aOff[i_viscousStrain], aOff[i_totalStrain]
    };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    IsotropicLinearElasticityPlaneStrain::meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                     aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v


// ---------------------------------------------------------------------------------------------------------------------
// f1 function for isotropic linear generalized Maxwell viscoelastic plane strain material with
// reference stress and strain.
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::f1v_refstate(const PylithInt dim,
                                                                      const PylithInt numS,
                                                                      const PylithInt numA,
                                                                      const PylithInt sOff[],
                                                                      const PylithInt sOff_x[],
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
    const PylithInt i_rstress = numA-8;
    const PylithInt i_rstrain = numA-7;
    const PylithInt i_shearModulus = numA-6;
    const PylithInt i_bulkModulus = numA-5;
    const PylithInt i_maxwellTime = numA-4;
    const PylithInt i_shearModulusRatio = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 8);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f1);

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 14);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Pass bulk modulus, reference stress, and reference strain.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 7; // Pass shear modulus, Maxwell times, shear modulus ratios,
                                 // total strain, viscous strains,
                                 // reference stress, and reference strain.
    const PylithInt aOffDev[7] = {
        aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_shearModulusRatio],
        aOff[i_viscousStrain], aOff[i_totalStrain],
    };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    IsotropicLinearElasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                              aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants,
                                                              stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Jf3_vu entry function for 2-D plane strain isotropic linear generalized Maxwell
 * viscoelastic material.
 *
 * stress_ij = C_ijkl strain_kl
 *
 * dq = maxwellTime * (1.0 - exp(-dt/maxwellTime))/dt
 *
 * stress_11 = C1111 strain_11 + C1122 strain_22,
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
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::Jf3vu(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
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
    const PylithInt i_shearModulus = numA-6;
    const PylithInt i_bulkModulus = numA-5;
    const PylithInt i_maxwellTime = numA-4;
    const PylithInt i_shearModulusRatio = numA-3;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 6);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(a);
    assert(numConstants == 1);
    assert(constants);
    assert(Jf3);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar maxwellTime_1 = a[aOff[i_maxwellTime]+0];
    const PylithScalar maxwellTime_2 = a[aOff[i_maxwellTime]+1];
    const PylithScalar maxwellTime_3 = a[aOff[i_maxwellTime]+2];
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio] + 1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio] + 2];
    const PylithScalar dt = constants[0];

    const PylithScalar dq_1 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_1);
    const PylithScalar dq_2 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_2);
    const PylithScalar dq_3 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_3);

    // Unique components of Jacobian.
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2 - shearModulusRatio_3;
    const PylithScalar shearFactor = shearModulus *
                                     (dq_1 * shearModulusRatio_1 + dq_2 * shearModulusRatio_2 + dq_3 * shearModulusRatio_3 +
                                      shearModulusRatio_0);

    const PylithReal C1111 = bulkModulus + 4.0 * shearFactor/3.0;
    const PylithReal C1122 = bulkModulus - 2.0 * shearFactor/3.0;
    const PylithReal C1212 = shearFactor;
    /* j(f,g,df,dg) = C(f,df,g,dg)
     *
     * 0:  j0000 = C1111 = 1.0*bulkModulus + 2.0*shearModulus*(0.666666666666667*dq_1*shearModulusRatio_1 +
     * 0.666666666666667*dq_2*shearModulusRatio_2 + 0.666666666666667*dq_3*shearModulusRatio_3 +
     * 0.666666666666667*shearModulusRatio_0)
     * 1:  j0001 = C1112 = 0
     * 2:  j0010 = C1211 = 0
     * 3:  j0011 = C1212 = 2.0*shearModulus*(0.5*dq_1*shearModulusRatio_1 + 0.5*dq_2*shearModulusRatio_2 +
     * 0.5*dq_3*shearModulusRatio_3 + 0.5*shearModulusRatio_0)
     * 4:  j0100 = C1121 = 0
     * 5:  j0101 = C1122 = 1.0*bulkModulus + 2.0*shearModulus*(-0.333333333333333*dq_1*shearModulusRatio_1 -
     * 0.333333333333333*dq_2*shearModulusRatio_2 - 0.333333333333333*dq_3*shearModulusRatio_3 -
     * 0.333333333333333*shearModulusRatio_0)
     * 6:  j0110 = C1221 = 2.0*shearModulus*(0.5*dq_1*shearModulusRatio_1 + 0.5*dq_2*shearModulusRatio_2 +
     * 0.5*dq_3*shearModulusRatio_3 + 0.5*shearModulusRatio_0)
     * 7:  j0111 = C1222 = 0
     * 8:  j1000 = C2111 = 0
     * 9:  j1001 = C2112 = 2.0*shearModulus*(0.5*dq_1*shearModulusRatio_1 + 0.5*dq_2*shearModulusRatio_2 +
     * 0.5*dq_3*shearModulusRatio_3 + 0.5*shearModulusRatio_0)
     * 10:  j1010 = C2211 = 1.0*bulkModulus + 2.0*shearModulus*(-0.333333333333333*dq_1*shearModulusRatio_1 -
     * 0.333333333333333*dq_2*shearModulusRatio_2 - 0.333333333333333*dq_3*shearModulusRatio_3 -
     * 0.333333333333333*shearModulusRatio_0)
     * 11:  j1011 = C2212 = 0
     * 12:  j1100 = C2121 = 2.0*shearModulus*(0.5*dq_1*shearModulusRatio_1 + 0.5*dq_2*shearModulusRatio_2 +
     * 0.5*dq_3*shearModulusRatio_3 + 0.5*shearModulusRatio_0)
     * 13:  j1101 = C2122 = 0
     * 14:  j1110 = C2221 = 0
     * 15:  j1111 = C2222 = 1.0*bulkModulus + 2.0*shearModulus*(0.666666666666667*dq_1*shearModulusRatio_1 +
     * 0.666666666666667*dq_2*shearModulusRatio_2 + 0.666666666666667*dq_3*shearModulusRatio_3 +
     * 0.666666666666667*shearModulusRatio_0)
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
} // Jf3vu


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * generalized Maxwell viscoelastic material WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = 2*shearModulus*visStrain_ij
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::deviatoricStress(const PylithInt dim,
                                                                          const PylithInt numS,
                                                                          const PylithInt numA,
                                                                          const PylithInt sOff[],
                                                                          const PylithInt sOff_x[],
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
    const PylithInt _strainSize = 4;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_maxwellTime = 1;
    const PylithInt i_shearModulusRatio = 2;
    const PylithInt i_viscousStrain = 3;
    const PylithInt i_totalStrain = 4;

    assert(_dim == dim);
    assert(1 == numS);
    assert(5 == numA);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s_x);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithInt _numS = 1; // Number passed on to visStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    // Viscous strain.
    const PylithInt numAVis = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain] };

    PylithScalar visStrainTpdt[12] = {0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis, NULL, a, a_t, NULL,
                         t, x, numConstants, constants, visStrainTpdt);

    // Shear modulus ratio factors.
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]+0];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio]+1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio]+2];
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2 - shearModulusRatio_3;

    // Deviatoric strains.
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar strainTpdt[4] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        0.0,
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0])
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0;

    const PylithScalar devStrainTpdt[4] = {
        strainTpdt[0] - meanStrainTpdt,
        strainTpdt[1] - meanStrainTpdt,
        strainTpdt[2] - meanStrainTpdt,
        strainTpdt[3]
    };

    // Stresses.
    stress[0] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[0] +
                                       shearModulusRatio_1 * visStrainTpdt[0] +
                                       shearModulusRatio_2 * visStrainTpdt[_strainSize] +
                                       shearModulusRatio_3 * visStrainTpdt[2*_strainSize]); // stress_xx
    stress[1] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[3] +
                                       shearModulusRatio_1 * visStrainTpdt[3] +
                                       shearModulusRatio_2 * visStrainTpdt[3 + _strainSize] +
                                       shearModulusRatio_3 * visStrainTpdt[3 + 2*_strainSize]); // stress_xy
    stress[2] += stress[1]; // stress_yx
    stress[3] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[1] +
                                       shearModulusRatio_1 * visStrainTpdt[1] +
                                       shearModulusRatio_2 * visStrainTpdt[1 + _strainSize] +
                                       shearModulusRatio_3 * visStrainTpdt[1 + 2*_strainSize]); // stress_yy
} // deviatoricStress


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * generalized Maxwell viscoelastic material WITH reference stress and reference strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = devrefstress_ij + 2*shearModulus*(visstrain_ij - devrefstrain_ij)
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::deviatoricStress_refstate(const PylithInt dim,
                                                                                   const PylithInt numS,
                                                                                   const PylithInt numA,
                                                                                   const PylithInt sOff[],
                                                                                   const PylithInt sOff_x[],
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
    const PylithInt _strainSize = 4;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = 0;
    const PylithInt i_rstrain = 1;
    const PylithInt i_shearModulus = 2;
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_shearModulusRatio = 4;
    const PylithInt i_viscousStrain = 5;
    const PylithInt i_totalStrain = 6;

    assert(_dim == dim);
    assert(1 == numS);
    assert(7 == numA);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s_x);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithInt _numS = 1; // Number passed on to visStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    // Viscous strain.
    const PylithInt numAVis = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain] };

    PylithScalar visStrainTpdt[12] = {0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis, NULL, a, a_t, NULL,
                         t, x, numConstants, constants, visStrainTpdt);

    // Compute reference deviatoric values.
    const PylithReal meanRefStrain = (refstrain[0] + refstrain[1] + refstrain[2])/3.0;
    const PylithScalar devRefStrain[4] = {refstrain[0] - meanRefStrain,
                                          refstrain[1] - meanRefStrain,
                                          refstrain[2] - meanRefStrain,
                                          refstrain[3]};
    const PylithReal meanRefStress = (refstress[0] + refstress[1] + refstress[2])/3.0;
    const PylithScalar devRefStress[4] = {refstress[0] - meanRefStress,
                                          refstress[1] - meanRefStress,
                                          refstress[2] - meanRefStress,
                                          refstress[3]};

    // Shear modulus ratio factors.
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]+0];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio]+1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio]+2];
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2 - shearModulusRatio_3;

    // Deviatoric strains.
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar strainTpdt[4] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        0.0,
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0])
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0;

    const PylithScalar devStrainTpdt[4] = {
        strainTpdt[0] - meanStrainTpdt,
        strainTpdt[1] - meanStrainTpdt,
        strainTpdt[2] - meanStrainTpdt,
        strainTpdt[3]
    };

    // Compute stress components -- note that we are including reference deviatoric
    // stress for now. This may need to be removed after testing.
    const PylithReal twomu = 2.0*shearModulus;
    stress[0] += devRefStress[0] + twomu * (shearModulusRatio_0 * devStrainTpdt[0] +
                                            shearModulusRatio_1 * visStrainTpdt[0] +
                                            shearModulusRatio_2 * visStrainTpdt[_strainSize] +
                                            shearModulusRatio_3 * visStrainTpdt[2*_strainSize] -
                                            devRefStrain[0]); // stress_xx
    stress[1] += devRefStress[3] + twomu * (shearModulusRatio_0 * devStrainTpdt[3] +
                                            shearModulusRatio_1 * visStrainTpdt[3] +
                                            shearModulusRatio_2 * visStrainTpdt[3 + _strainSize] +
                                            shearModulusRatio_3 * visStrainTpdt[3 + 2*_strainSize] -
                                            devRefStrain[3]); // stress_xy
    stress[2] += stress[1]; // stress_yx
    stress[3] += devRefStress[1] + twomu * (shearModulusRatio_0 * devStrainTpdt[1] +
                                            shearModulusRatio_1 * visStrainTpdt[1] +
                                            shearModulusRatio_2 * visStrainTpdt[1 + _strainSize] +
                                            shearModulusRatio_3 * visStrainTpdt[1 + 2*_strainSize] -
                                            devRefStrain[3]); // stress_yy
} // deviatoricStress_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate viscous strain for a 2D plane strain generalized Maxwell viscoelastic material.
 *
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::computeViscousStrain(const PylithInt dim,
                                                                              const PylithInt numS,
                                                                              const PylithInt numA,
                                                                              const PylithInt sOff[],
                                                                              const PylithInt sOff_x[],
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
                                                                              PylithScalar visStrainTpdt[]) {
    const PylithInt _dim = 2;
    const PylithInt _strainSize = 4;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = 0;
    const PylithInt i_viscousStrain = 1;
    const PylithInt i_totalStrain = 2;

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s_x);
    assert(aOff);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(a);
    assert(numConstants == 1);
    assert(constants);
    assert(visStrainTpdt);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    const PylithScalar maxwellTime_1 = a[aOff[i_maxwellTime]+0];
    const PylithScalar maxwellTime_2 = a[aOff[i_maxwellTime]+1];
    const PylithScalar maxwellTime_3 = a[aOff[i_maxwellTime]+2];
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];
    const PylithScalar* viscousStrain_1 = &a[aOff[i_viscousStrain]];
    const PylithScalar* viscousStrain_2 = &a[aOff[i_viscousStrain] + _strainSize];
    const PylithScalar* viscousStrain_3 = &a[aOff[i_viscousStrain] + 2*_strainSize];

    const PylithScalar dt = constants[0];

    const PylithScalar dq_1 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_1);
    const PylithScalar expFac_1 = exp(-dt/maxwellTime_1);

    const PylithScalar dq_2 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_2);
    const PylithScalar expFac_2 = exp(-dt/maxwellTime_2);

    const PylithScalar dq_3 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_3);
    const PylithScalar expFac_3 = exp(-dt/maxwellTime_3);

    const PylithScalar strainTpdt[4] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        0.0,
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0])
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0;
    const PylithReal meanStrainT = (totalStrain[0] + totalStrain[1])/3.0;

    const PylithScalar devStrainTpdt[4] = {
        strainTpdt[0] - meanStrainTpdt,
        strainTpdt[1] - meanStrainTpdt,
        strainTpdt[2] - meanStrainTpdt,
        strainTpdt[3]
    };

    const PylithScalar devStrainT[4] = {
        totalStrain[0] - meanStrainT,
        totalStrain[1] - meanStrainT,
        totalStrain[2] - meanStrainT,
        totalStrain[3]
    };

    PylithScalar strainDiff = 0.0;
    for (int iComp = 0; iComp < 4; ++iComp) {
        strainDiff = devStrainTpdt[iComp] - devStrainT[iComp];
        visStrainTpdt[iComp] = expFac_1*viscousStrain_1[iComp] + dq_1*strainDiff;
        visStrainTpdt[iComp + _strainSize] = expFac_2*viscousStrain_2[iComp] + dq_2*strainDiff;
        visStrainTpdt[iComp + 2*_strainSize] = expFac_3*viscousStrain_3[iComp] + dq_3*strainDiff;
    } // for
} // computeViscousStrain


// ---------------------------------------------------------------------------------------------------------------------
/* Update total strain for a 2D plane strain Maxwell viscoelastic material.
 *
 * IMPORTANT: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::updateTotalStrain(const PylithInt dim,
                                                                           const PylithInt numS,
                                                                           const PylithInt numA,
                                                                           const PylithInt sOff[],
                                                                           const PylithInt sOff_x[],
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
                                                                           PylithScalar totalStrain[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

#if 0 // :DEBUG:
    std::cout << "dim:  " << dim << std::endl;
    std::cout << "numS:  " << numS << std::endl;
    std::cout << "numA:  " << numA << std::endl;
    std::cout << "sOff[0]:  " << sOff[0] << std::endl;
    std::cout << "sOff_x[0]:  " << sOff_x[0] << std::endl;
    std::cout << "s[0]:  " << s[0] << std::endl;
    std::cout << "aOff[0]:  " << aOff[0] << std::endl;
    std::cout << "a[0]:  " << a[0] << std::endl;
    std::cout << "t:  " << t << std::endl;
    std::cout << "x[0]:  " << x[0] << std::endl;
    std::cout << "numConstants:  " << numConstants << std::endl;
    std::cout << "totalStrain[0]:  " << totalStrain[0] << std::endl;
#endif

    assert(_dim == dim);
    assert(numS >= 1);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s_x);
    assert(totalStrain);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    totalStrain[0] = disp_x[0*_dim+0];
    totalStrain[1] = disp_x[1*_dim+1];
    totalStrain[2] = 0.0;
    totalStrain[3] = 0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]);
} // updateTotalStrain


// ---------------------------------------------------------------------------------------------------------------------
/* Update viscous strain for 2D plane strain generalized Maxwell.
 *
 * IMPORTANT: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::updateViscousStrain(const PylithInt dim,
                                                                             const PylithInt numS,
                                                                             const PylithInt numA,
                                                                             const PylithInt sOff[],
                                                                             const PylithInt sOff_x[],
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
                                                                             PylithScalar visStrain[]) {
    const PylithInt _dim = 2;
    const PylithInt _strainSize = 4;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = numA-4;
    const PylithInt i_shearModulusRatio = numA-3;
    const PylithInt i_viscousStrainPrev = numA-2;
    const PylithInt i_totalStrainPrev = numA-1;

#if 0 // :DEBUG:
    std::cout << "dim:  " << dim << std::endl;
    std::cout << "numS:  " << numS << std::endl;
    std::cout << "numA:  " << numA << std::endl;
    std::cout << "sOff[0]:  " << sOff[0] << std::endl;
    std::cout << "sOff_x[0]:  " << sOff_x[0] << std::endl;
    std::cout << "s[0]:  " << s[0] << std::endl;
    std::cout << "aOff[0]:  " << aOff[0] << std::endl;
    std::cout << "a[0]:  " << a[0] << std::endl;
    std::cout << "t:  " << t << std::endl;
    std::cout << "x[0]:  " << x[0] << std::endl;
    std::cout << "numConstants:  " << numConstants << std::endl;
    std::cout << "visStrain[0]:  " << visStrain[0] << std::endl;
#endif

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s_x);
    assert(aOff);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrainPrev] >= 0);
    assert(a);
    assert(visStrain);
    assert(constants);
    assert(1 == numConstants);

    // Compute strain, deviatoric strain, etc.
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    const PylithScalar maxwellTime_1 = a[aOff[i_maxwellTime]+0];
    const PylithScalar maxwellTime_2 = a[aOff[i_maxwellTime]+1];
    const PylithScalar maxwellTime_3 = a[aOff[i_maxwellTime]+2];
    const PylithScalar* viscousStrain_1 = &a[aOff[i_viscousStrainPrev]];
    const PylithScalar* viscousStrain_2 = &a[aOff[i_viscousStrainPrev] + _strainSize];
    const PylithScalar* viscousStrain_3 = &a[aOff[i_viscousStrainPrev] + 2*_strainSize];
    const PylithScalar* totalStrainPrev = &a[aOff[i_totalStrainPrev]];

    const PylithScalar dt = constants[0];

    const PylithScalar dq_1 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_1);
    const PylithScalar expFac_1 = exp(-dt/maxwellTime_1);

    const PylithScalar dq_2 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_2);
    const PylithScalar expFac_2 = exp(-dt/maxwellTime_2);

    const PylithScalar dq_3 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_3);
    const PylithScalar expFac_3 = exp(-dt/maxwellTime_3);

    const PylithScalar strain[4] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        0.0,
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0])
    };
    const PylithReal meanStrain = (strain[0] + strain[1])/3.0;

    const PylithScalar devStrain[4] = {
        strain[0] - meanStrain,
        strain[1] - meanStrain,
        strain[2] - meanStrain,
        strain[3]
    };

    const PylithReal meanStrainPrev = (totalStrainPrev[0] + totalStrainPrev[1])/3.0;

    const PylithScalar devStrainPrev[4] = {
        totalStrainPrev[0] - meanStrainPrev,
        totalStrainPrev[1] - meanStrainPrev,
        totalStrainPrev[2] - meanStrainPrev,
        totalStrainPrev[3]
    };

    PylithScalar strainDiff = 0.0;
    for (int iComp = 0; iComp < 4; ++iComp) {
        strainDiff = devStrain[iComp] - devStrainPrev[iComp];
        visStrain[iComp] = expFac_1*viscousStrain_1[iComp] + dq_1*strainDiff;
        visStrain[iComp + _strainSize] = expFac_2*viscousStrain_2[iComp] + dq_2*strainDiff;
        visStrain[iComp + 2*_strainSize] = expFac_3*viscousStrain_3[iComp] + dq_3*strainDiff;
    } // for
} // updateViscousStrain


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate stress for 2-D plane strain isotropic linear generalized
 * Maxwell material WITHOUT a reference stress and strain.
 *
 * Used in outputing the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
 * total_strain(4), viscous_strain(4)]
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::cauchyStress(const PylithInt dim,
                                                                      const PylithInt numS,
                                                                      const PylithInt numA,
                                                                      const PylithInt sOff[],
                                                                      const PylithInt sOff_x[],
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
    const PylithInt _strainSize = 4;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-6;
    const PylithInt i_bulkModulus = numA-5;
    const PylithInt i_maxwellTime = numA-4;
    const PylithInt i_shearModulusRatio = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 6);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    // Shear modulus and (updated) viscous strains.
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* visStrainTpdt_1 = &a[aOff[i_viscousStrain]];
    const PylithScalar* visStrainTpdt_2 = &a[aOff[i_viscousStrain] + _strainSize];
    const PylithScalar* visStrainTpdt_3 = &a[aOff[i_viscousStrain] + 2*_strainSize];
    const PylithScalar* totalStrainTpdt = &a[aOff[i_totalStrain]];

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    IsotropicLinearElasticityPlaneStrain::meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                     aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);

    const PylithScalar meanStress = stressTensor[0];

    const PylithReal meanStrainTpdt = (totalStrainTpdt[0] + totalStrainTpdt[1] + totalStrainTpdt[2])/3.0;

    const PylithScalar devStrainTpdt[4] = {
        totalStrainTpdt[0] - meanStrainTpdt,
        totalStrainTpdt[1] - meanStrainTpdt,
        totalStrainTpdt[2] - meanStrainTpdt,
        totalStrainTpdt[3]
    };

    // Shear modulus ratio factors.
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio] + 1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio] + 2];
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2 - shearModulusRatio_3;

    // Stresses.
    stressVector[0] = meanStress + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[0] +
                                                     shearModulusRatio_1*visStrainTpdt_1[0] +
                                                     shearModulusRatio_2*visStrainTpdt_2[0] +
                                                     shearModulusRatio_3*visStrainTpdt_3[0]); // stress_xx
    stressVector[1] = meanStress + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[1] +
                                                     shearModulusRatio_1*visStrainTpdt_1[1] +
                                                     shearModulusRatio_2*visStrainTpdt_2[1] +
                                                     shearModulusRatio_3*visStrainTpdt_3[1]); // stress_yy
    stressVector[2] = meanStress + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[2] +
                                                     shearModulusRatio_1*visStrainTpdt_1[2] +
                                                     shearModulusRatio_2*visStrainTpdt_2[2] +
                                                     shearModulusRatio_3*visStrainTpdt_3[2]); // stress_zz
    stressVector[3] = 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[3] +
                                        shearModulusRatio_1*visStrainTpdt_1[3] +
                                        shearModulusRatio_2*visStrainTpdt_2[3] +
                                        shearModulusRatio_3*visStrainTpdt_3[3]); // stress_xy

} // cauchyStress


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate stress for 2-D plane strain isotropic linear generalized
 * Maxwell material WITH a reference stress/strain.
 *
 * Used in outputing the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
 *                   maxwell_time(3), shear_modulus_ratio(3), total_strain(4), viscous_strain(4)]
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::cauchyStress_refstate(const PylithInt dim,
                                                                               const PylithInt numS,
                                                                               const PylithInt numA,
                                                                               const PylithInt sOff[],
                                                                               const PylithInt sOff_x[],
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
    const PylithInt _strainSize = 4;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-8;
    const PylithInt i_rstrain = numA-7;
    const PylithInt i_shearModulus = numA-6;
    const PylithInt i_bulkModulus = numA-5;
    const PylithInt i_maxwellTime = numA-4;
    const PylithInt i_shearModulusRatio = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 8);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    // Shear modulus and (updated) viscous strains.
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* visStrainTpdt_1 = &a[aOff[i_viscousStrain]];
    const PylithScalar* visStrainTpdt_2 = &a[aOff[i_viscousStrain] + _strainSize];
    const PylithScalar* visStrainTpdt_3 = &a[aOff[i_viscousStrain] + 2*_strainSize];
    const PylithScalar* totalStrainTpdt = &a[aOff[i_totalStrain]];

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    IsotropicLinearElasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                              aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);

    const PylithScalar meanStress = stressTensor[0];

    // Reference stress and strain.
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    // Compute reference deviatoric values.
    const PylithReal meanRefStrain = (refstrain[0] + refstrain[1] + refstrain[2])/3.0;
    const PylithScalar devRefStrain[4] = {refstrain[0] - meanRefStrain,
                                          refstrain[1] - meanRefStrain,
                                          refstrain[2] - meanRefStrain,
                                          refstrain[3]};
    const PylithReal meanRefStress = (refstress[0] + refstress[1] + refstress[2])/3.0;
    const PylithScalar devRefStress[4] = {refstress[0] - meanRefStress,
                                          refstress[1] - meanRefStress,
                                          refstress[2] - meanRefStress,
                                          refstress[3]};

    const PylithReal meanStrainTpdt = (totalStrainTpdt[0] + totalStrainTpdt[1] + totalStrainTpdt[2])/3.0;

    const PylithScalar devStrainTpdt[4] = {
        totalStrainTpdt[0] - meanStrainTpdt,
        totalStrainTpdt[1] - meanStrainTpdt,
        totalStrainTpdt[2] - meanStrainTpdt,
        totalStrainTpdt[3]
    };

    // Shear modulus ratio factors.
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio] + 1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio] + 2];
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2 - shearModulusRatio_3;

    // Stresses.
    stressVector[0] = meanStress + devRefStress[0] + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[0] +
                                                                       shearModulusRatio_1*visStrainTpdt_1[0] +
                                                                       shearModulusRatio_2*visStrainTpdt_2[0] +
                                                                       shearModulusRatio_3*visStrainTpdt_3[0] -
                                                                       devRefStrain[0]); // stress_xx
    stressVector[1] = meanStress + devRefStress[1] + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[1] +
                                                                       shearModulusRatio_1*visStrainTpdt_1[1] +
                                                                       shearModulusRatio_2*visStrainTpdt_2[1] +
                                                                       shearModulusRatio_3*visStrainTpdt_3[1] -
                                                                       devRefStrain[1]); // stress_yy
    stressVector[2] = meanStress + devRefStress[2] + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[2] +
                                                                       shearModulusRatio_1*visStrainTpdt_1[2] +
                                                                       shearModulusRatio_2*visStrainTpdt_2[2] +
                                                                       shearModulusRatio_3*visStrainTpdt_3[2] -
                                                                       devRefStrain[2]); // stress_zz
    stressVector[3] = devRefStress[3] + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[3] +
                                                          shearModulusRatio_1*visStrainTpdt_1[3] +
                                                          shearModulusRatio_2*visStrainTpdt_2[3] +
                                                          shearModulusRatio_3*visStrainTpdt_3[3] -
                                                          devRefStrain[3]); // stress_xy

} // cauchyStress_refstate


// =====================================================================================================================
// Kernels for isotropic, linear generalized Maxwell viscoelastic 3D material.
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
/* f1 function for isotropic linear generalized Maxwell 3D WITHOUT reference stress and reference strain.
 *
 * Solution fields: [disp(dim), ...]
 * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
 *                    total_strain(4), viscous_strain(12)]
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::f1v(const PylithInt dim,
                                                    const PylithInt numS,
                                                    const PylithInt numA,
                                                    const PylithInt sOff[],
                                                    const PylithInt sOff_x[],
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
    const PylithInt i_shearModulus = numA-6;
    const PylithInt i_bulkModulus = numA-5;
    const PylithInt i_maxwellTime = numA-4;
    const PylithInt i_shearModulusRatio = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 6);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 5; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[5] = {
        aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_shearModulusRatio], aOff[i_viscousStrain], aOff[i_totalStrain]
    };

    PylithScalar stressTensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    IsotropicLinearElasticity3D::meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                            aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v


// ---------------------------------------------------------------------------------------------------------------------
/* f1 function for isotropic linear generalized Maxwell 3D WITH reference stress and reference strain.
 *
 * Solution fields: [disp(dim), ...]
 * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
 *                    maxwell_time(3), shear_modulus_ratio(3), total_strain(4), viscous_strain(12)]
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::f1v_refstate(const PylithInt dim,
                                                             const PylithInt numS,
                                                             const PylithInt numA,
                                                             const PylithInt sOff[],
                                                             const PylithInt sOff_x[],
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
    const PylithInt i_rstress = numA-8;
    const PylithInt i_rstrain = numA-7;
    const PylithInt i_shearModulus = numA-6;
    const PylithInt i_bulkModulus = numA-5;
    const PylithInt i_maxwellTime = numA-4;
    const PylithInt i_shearModulusRatio = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 8);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 14);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Pass bulk modulus, reference stress, and reference strain.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 7; // Pass shear modulus, Maxwell times, shear modulus ratios,
                                 // total strain, viscous strains,
                                 // reference stress, and reference strain.
    const PylithInt aOffDev[7] = {
        aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_shearModulusRatio],
        aOff[i_viscousStrain], aOff[i_totalStrain],
    };

    PylithScalar stressTensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    IsotropicLinearElasticity3D::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                     aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants,
                                                     stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for

} // f1v_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Jf3_vu entry function for 3-D isotropic linear generalized Maxwell viscoelasticity.
 *
 * Solution fields: [...]
 * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
 *                    total_strain(4), viscous_strain(12)]
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::Jf3vu(const PylithInt dim,
                                                      const PylithInt numS,
                                                      const PylithInt numA,
                                                      const PylithInt sOff[],
                                                      const PylithInt sOff_x[],
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
    const PylithInt _dim = 3;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-6;
    const PylithInt i_bulkModulus = numA-5;
    const PylithInt i_maxwellTime = numA-4;
    const PylithInt i_shearModulusRatio = numA-3;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 6);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(a);
    assert(numConstants == 1);
    assert(constants);
    assert(Jf3);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar maxwellTime_1 = a[aOff[i_maxwellTime]];
    const PylithScalar maxwellTime_2 = a[aOff[i_maxwellTime] + 1];
    const PylithScalar maxwellTime_3 = a[aOff[i_maxwellTime] + 2];
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]+0];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio]+1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio]+2];
    const PylithScalar dt = constants[0];

    const PylithScalar dq_1 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_1);
    const PylithScalar dq_2 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_2);
    const PylithScalar dq_3 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_3);

    // Unique components of Jacobian.
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2 - shearModulusRatio_3;
    const PylithScalar shearFactor = shearModulus *
                                     (dq_1 * shearModulusRatio_1 + dq_2 * shearModulusRatio_2 + dq_3 * shearModulusRatio_3 +
                                      shearModulusRatio_0);

    const PylithReal C1111 = bulkModulus + 4.0 * shearFactor/3.0;
    const PylithReal C1122 = bulkModulus - 2.0 * shearFactor/3.0;
    const PylithReal C1212 = shearFactor;

    /* j(f,g,df,dg) = C(f,df,g,dg)
     * 0:  j0000 = C1111 = bulkModulus + 2*shearModulus*(2*dq_1*shearModulusRatio_1/3 + 2*dq_2*shearModulusRatio_2/3 +
     * 2*dq_3*shearModulusRatio_3/3 + 2*shearModulusRatio_0/3)
     * 1:  j0001 = C1112 = 0
     * 2:  j0002 = C1113 = 0
     * 3:  j0010 = C1211 = 0
     * 4:  j0011 = C1212 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 +
     * dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
     * 5:  j0012 = C1213 = 0
     * 6:  j0020 = C1311 = 0
     * 7:  j0021 = C1312 = 0
     * 8:  j0022 = C1313 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 +
     * dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
     * 9:  j0100 = C1121 = 0
     * 10:  j0101 = C1122 = bulkModulus + 2*shearModulus*(-dq_1*shearModulusRatio_1/3 - dq_2*shearModulusRatio_2/3 -
     * dq_3*shearModulusRatio_3/3 - shearModulusRatio_0/3)
     * 11:  j0102 = C1123 = 0
     * 12:  j0110 = C1221 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 +
     * dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
     * 13:  j0111 = C1222 = 0
     * 14:  j0112 = C1223 = 0
     * 15:  j0120 = C1321 = 0
     * 16:  j0121 = C1322 = 0
     * 17:  j0122 = C1323 = 0
     * 18:  j0200 = C1131 = 0
     * 19:  j0201 = C1132 = 0
     * 20:  j0202 = C1133 = bulkModulus + 2*shearModulus*(-dq_1*shearModulusRatio_1/3 - dq_2*shearModulusRatio_2/3 -
     * dq_3*shearModulusRatio_3/3 - shearModulusRatio_0/3)
     * 21:  j0210 = C1231 = 0
     * 22:  j0211 = C1232 = 0
     * 23:  j0212 = C1233 = 0
     * 24:  j0220 = C1331 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 +
     * dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
     * 25:  j0221 = C1332 = 0
     * 26:  j0222 = C1333 = 0
     * 27:  j1000 = C2111 = 0
     * 28:  j1001 = C2112 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 +
     * dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
     * 29:  j1002 = C2113 = 0
     * 30:  j1010 = C2211 = bulkModulus + 2*shearModulus*(-dq_1*shearModulusRatio_1/3 - dq_2*shearModulusRatio_2/3 -
     * dq_3*shearModulusRatio_3/3 - shearModulusRatio_0/3)
     * 31:  j1011 = C2212 = 0
     * 32:  j1012 = C2213 = 0
     * 33:  j1020 = C2311 = 0
     * 34:  j1021 = C2312 = 0
     * 35:  j1022 = C2313 = 0
     * 36:  j1100 = C2121 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 +
     * dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
     * 37:  j1101 = C2122 = 0
     * 38:  j1102 = C2123 = 0
     * 39:  j1110 = C2221 = 0
     * 40:  j1111 = C2222 = bulkModulus + 2*shearModulus*(2*dq_1*shearModulusRatio_1/3 + 2*dq_2*shearModulusRatio_2/3 +
     * 2*dq_3*shearModulusRatio_3/3 + 2*shearModulusRatio_0/3)
     * 41:  j1112 = C2223 = 0
     * 42:  j1120 = C2321 = 0
     * 43:  j1121 = C2322 = 0
     * 44:  j1122 = C2323 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 +
     * dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
     * 45:  j1200 = C2131 = 0
     * 46:  j1201 = C2132 = 0
     * 47:  j1202 = C2133 = 0
     * 48:  j1210 = C2231 = 0
     * 49:  j1211 = C2232 = 0
     * 50:  j1212 = C2233 = bulkModulus + 2*shearModulus*(-dq_1*shearModulusRatio_1/3 - dq_2*shearModulusRatio_2/3 -
     * dq_3*shearModulusRatio_3/3 - shearModulusRatio_0/3)
     * 51:  j1220 = C2331 = 0
     * 52:  j1221 = C2332 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 +
     * dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
     * 53:  j1222 = C2333 = 0
     * 54:  j2000 = C3111 = 0
     * 55:  j2001 = C3112 = 0
     * 56:  j2002 = C3113 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 +
     * dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
     * 57:  j2010 = C3211 = 0
     * 58:  j2011 = C3212 = 0
     * 59:  j2012 = C3213 = 0
     * 60:  j2020 = C3311 = bulkModulus + 2*shearModulus*(-dq_1*shearModulusRatio_1/3 - dq_2*shearModulusRatio_2/3 -
     * dq_3*shearModulusRatio_3/3 - shearModulusRatio_0/3)
     * 61:  j2021 = C3312 = 0
     * 62:  j2022 = C3313 = 0
     * 63:  j2100 = C3121 = 0
     * 64:  j2101 = C3122 = 0
     * 65:  j2102 = C3123 = 0
     * 66:  j2110 = C3221 = 0
     * 67:  j2111 = C3222 = 0
     * 68:  j2112 = C3223 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 +
     * dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
     * 69:  j2120 = C3321 = 0
     * 70:  j2121 = C3322 = bulkModulus + 2*shearModulus*(-dq_1*shearModulusRatio_1/3 - dq_2*shearModulusRatio_2/3 -
     * dq_3*shearModulusRatio_3/3 - shearModulusRatio_0/3)
     * 71:  j2122 = C3323 = 0
     * 72:  j2200 = C3131 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 +
     * dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
     * 73:  j2201 = C3132 = 0
     * 74:  j2202 = C3133 = 0
     * 75:  j2210 = C3231 = 0
     * 76:  j2211 = C3232 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 +
     * dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
     * 77:  j2212 = C3233 = 0
     * 78:  j2220 = C3331 = 0
     * 79:  j2221 = C3332 = 0
     * 80:  j2222 = C3333 = bulkModulus + 2*shearModulus*(2*dq_1*shearModulusRatio_1/3 + 2*dq_2*shearModulusRatio_2/3 +
     * 2*dq_3*shearModulusRatio_3/3 + 2*shearModulusRatio_0/3)
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
} // Jf3vu


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 3-D isotropic linear
 * generalized Maxwell viscoelasticity WITHOUT reference stress and strain.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [..., shear_modulus(1), maxwell_time(3), shear_modulus_ratio(3), total_strain(4),
 *                    viscous_strain(12)]
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::deviatoricStress(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
    const PylithInt _strainSize = 6;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_maxwellTime = 1;
    const PylithInt i_shearModulusRatio = 2;
    const PylithInt i_viscousStrain = 3;
    const PylithInt i_totalStrain = 4;

    assert(_dim == dim);
    assert(1 == numS);
    assert(5 == numA);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s_x);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithInt _numS = 1; // Number passed on to visStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    // Viscous strain.
    const PylithInt numAVis = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain],
                                   aOff[i_totalStrain] };

    PylithScalar visStrainTpdt[18] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis, NULL, a, a_t, NULL,
                         t, x, numConstants, constants, visStrainTpdt);

    // Shear modulus ratio factors.
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio] + 1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio] + 2];
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2 - shearModulusRatio_3;

    // Deviatoric strains.
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar strainTpdt[6] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        disp_x[2*_dim+2],
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]),
        0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]),
        0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0])
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1] + strainTpdt[2])/3.0;

    const PylithScalar devStrainTpdt[6] = {
        strainTpdt[0] - meanStrainTpdt,
        strainTpdt[1] - meanStrainTpdt,
        strainTpdt[2] - meanStrainTpdt,
        strainTpdt[3],
        strainTpdt[4],
        strainTpdt[5]
    };

    // Stresses.
    stress[0*_dim+0] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[0] +
                                              shearModulusRatio_1 * visStrainTpdt[0] +
                                              shearModulusRatio_2 * visStrainTpdt[_strainSize] +
                                              shearModulusRatio_3 * visStrainTpdt[2 * _strainSize]); // stress_xx
    stress[0*_dim+1] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[3] +
                                              shearModulusRatio_1 * visStrainTpdt[3] +
                                              shearModulusRatio_2 * visStrainTpdt[3 + _strainSize] +
                                              shearModulusRatio_3 * visStrainTpdt[3 + 2 * _strainSize]); // stress_xy
    stress[0*_dim+2] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[5] +
                                              shearModulusRatio_1 * visStrainTpdt[5] +
                                              shearModulusRatio_2 * visStrainTpdt[5 + _strainSize] +
                                              shearModulusRatio_3 * visStrainTpdt[5 + 2 * _strainSize]); // stress_xz
    stress[1*_dim+0] += stress[1]; // stress_yx
    stress[1*_dim+1] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[1] +
                                              shearModulusRatio_1 * visStrainTpdt[1] +
                                              shearModulusRatio_2 * visStrainTpdt[1 + _strainSize] +
                                              shearModulusRatio_3 * visStrainTpdt[1 + 2 * _strainSize]); // stress_yy
    stress[1*_dim+2] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[4] +
                                              shearModulusRatio_1 * visStrainTpdt[4] +
                                              shearModulusRatio_2 * visStrainTpdt[4 + _strainSize] +
                                              shearModulusRatio_3 * visStrainTpdt[4 + 2 * _strainSize]); // stress_yz
    stress[2*_dim+0] += stress[2]; // stress_zx
    stress[2*_dim+1] += stress[5]; // stress_zy
    stress[2*_dim+2] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[2] +
                                              shearModulusRatio_1 * visStrainTpdt[2] +
                                              shearModulusRatio_2 * visStrainTpdt[2 + _strainSize] +
                                              shearModulusRatio_3 * visStrainTpdt[2 + 2 * _strainSize]); // stress_zz
} // deviatoricStress


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 3-D isotropic linear
 * generalized Maxwell viscoelasticity WITH reference stress and strain.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), maxwell_time(3),
 *                    shear_modulus_ratio(3), total_strain(4), viscous_strain(12)]
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::deviatoricStress_refstate(const PylithInt dim,
                                                                          const PylithInt numS,
                                                                          const PylithInt numA,
                                                                          const PylithInt sOff[],
                                                                          const PylithInt sOff_x[],
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
    const PylithInt _strainSize = 6;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = 0;
    const PylithInt i_rstrain = 1;
    const PylithInt i_shearModulus = 2;
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_shearModulusRatio = 4;
    const PylithInt i_viscousStrain = 5;
    const PylithInt i_totalStrain = 6;

    assert(_dim == dim);
    assert(1 == numS);
    assert(6 == numA);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s_x);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy, stress_yz,
                                                         // stress_xz
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy, strain_yz,
                                                         // strain_xz

    const PylithInt _numS = 1; // Number passed on to visStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    // Viscous strain.
    const PylithInt numAVis = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain] };

    PylithScalar visStrainTpdt[18] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis, NULL, a, a_t, NULL,
                         t, x, numConstants, constants, visStrainTpdt);

    // Compute reference deviatoric values.
    const PylithReal meanRefStrain = (refstrain[0] + refstrain[1] + refstrain[2])/3.0;
    const PylithScalar devRefStrain[6] = {refstrain[0] - meanRefStrain,
                                          refstrain[1] - meanRefStrain,
                                          refstrain[2] - meanRefStrain,
                                          refstrain[3],
                                          refstrain[4],
                                          refstrain[5]};
    const PylithReal meanRefStress = (refstress[0] + refstress[1] + refstress[2])/3.0;
    const PylithScalar devRefStress[6] = {refstress[0] - meanRefStress,
                                          refstress[1] - meanRefStress,
                                          refstress[2] - meanRefStress,
                                          refstress[3],
                                          refstress[4],
                                          refstress[5]};

    // Shear modulus ratio factors.
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio] + 1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio] + 2];
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2 - shearModulusRatio_3;

    // Deviatoric strains.
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar strainTpdt[6] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        disp_x[2*_dim+2],
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]),
        0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]),
        0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0])
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1] + strainTpdt[2])/3.0;

    const PylithScalar devStrainTpdt[6] = {
        strainTpdt[0] - meanStrainTpdt,
        strainTpdt[1] - meanStrainTpdt,
        strainTpdt[2] - meanStrainTpdt,
        strainTpdt[3],
        strainTpdt[4],
        strainTpdt[5]
    };

    // Compute stress components -- note that we are including reference deviatoric
    // stress for now. This may need to be removed after testing.
    const PylithReal twomu = 2.0*shearModulus;
    stress[0*_dim+0] += devRefStress[0] + twomu * (shearModulusRatio_0 * devStrainTpdt[0] +
                                                   shearModulusRatio_1 * visStrainTpdt[0] +
                                                   shearModulusRatio_2 * visStrainTpdt[_strainSize] +
                                                   shearModulusRatio_3 * visStrainTpdt[2 * _strainSize] -
                                                   devRefStrain[0]); // stress_xx
    stress[0*_dim+1] += devRefStress[3] + twomu * (shearModulusRatio_0 * devStrainTpdt[3] +
                                                   shearModulusRatio_1 * visStrainTpdt[3] +
                                                   shearModulusRatio_2 * visStrainTpdt[3 + _strainSize] +
                                                   shearModulusRatio_3 * visStrainTpdt[3 + 2 * _strainSize] -
                                                   devRefStrain[3]); // stress_xy
    stress[0*_dim+2] += devRefStress[5] + twomu * (shearModulusRatio_0 * devStrainTpdt[5] +
                                                   shearModulusRatio_1 * visStrainTpdt[5] +
                                                   shearModulusRatio_2 * visStrainTpdt[5 + _strainSize] +
                                                   shearModulusRatio_3 * visStrainTpdt[5 + 2 * _strainSize] -
                                                   devRefStrain[5]); // stress_xz
    stress[1*_dim+0] += stress[1]; // stress_yx
    stress[1*_dim+1] += devRefStress[1] + twomu * (shearModulusRatio_0 * devStrainTpdt[1] +
                                                   shearModulusRatio_1 * visStrainTpdt[1] +
                                                   shearModulusRatio_2 * visStrainTpdt[1 + _strainSize] +
                                                   shearModulusRatio_3 * visStrainTpdt[1 + 2 * _strainSize] -
                                                   devRefStrain[1]); // stress_yy
    stress[1*_dim+2] += devRefStress[4] + twomu * (shearModulusRatio_0 * devStrainTpdt[4] +
                                                   shearModulusRatio_1 * visStrainTpdt[4] +
                                                   shearModulusRatio_2 * visStrainTpdt[4 + _strainSize] +
                                                   shearModulusRatio_3 * visStrainTpdt[4 + 2 * _strainSize] -
                                                   devRefStrain[4]); // stress_yz
    stress[2*_dim+0] += stress[2]; // stress_zx
    stress[2*_dim+1] += stress[5]; // stress_zy
    stress[2*_dim+2] += devRefStress[2] + twomu * (shearModulusRatio_0 * devStrainTpdt[2] +
                                                   shearModulusRatio_1 * visStrainTpdt[2] +
                                                   shearModulusRatio_2 * visStrainTpdt[2 + _strainSize] +
                                                   shearModulusRatio_3 * visStrainTpdt[2 + 2 * _strainSize] -
                                                   devRefStrain[2]); // stress_zz
} // deviatoricStress_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate viscous strain at t+dt for 3-D isotropic linear
 * generalized Maxwell viscoelasticity.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [maxwell_time(3), viscous_strain(18), total_strain(6)]
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::computeViscousStrain(const PylithInt dim,
                                                                     const PylithInt numS,
                                                                     const PylithInt numA,
                                                                     const PylithInt sOff[],
                                                                     const PylithInt sOff_x[],
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
                                                                     PylithScalar visStrainTpdt[]) {
    const PylithInt _dim = 3;
    const PylithInt _strainSize = 6;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = 0;
    const PylithInt i_viscousStrain = 1;
    const PylithInt i_totalStrain = 2;

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s_x);
    assert(aOff);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(a);
    assert(numConstants == 1);
    assert(constants);
    assert(visStrainTpdt);

    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    const PylithScalar maxwellTime_1 = a[aOff[i_maxwellTime]];
    const PylithScalar maxwellTime_2 = a[aOff[i_maxwellTime] + 1];
    const PylithScalar maxwellTime_3 = a[aOff[i_maxwellTime] + 2];
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];
    const PylithScalar* viscousStrain_1 = &a[aOff[i_viscousStrain]];
    const PylithScalar* viscousStrain_2 = &a[aOff[i_viscousStrain] + _strainSize];
    const PylithScalar* viscousStrain_3 = &a[aOff[i_viscousStrain] + 2 * _strainSize];

    const PylithScalar dt = constants[0];

    const PylithScalar dq_1 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_1);
    const PylithScalar expFac_1 = exp(-dt/maxwellTime_1);

    const PylithScalar dq_2 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_2);
    const PylithScalar expFac_2 = exp(-dt/maxwellTime_2);

    const PylithScalar dq_3 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_3);
    const PylithScalar expFac_3 = exp(-dt/maxwellTime_3);

    const PylithScalar strainTpdt[6] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        disp_x[2*_dim+2],
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]),
        0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]),
        0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0])
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1] + strainTpdt[2])/3.0;
    const PylithReal meanStrainT = (totalStrain[0] + totalStrain[1] + totalStrain[2])/3.0;

    const PylithScalar devStrainTpdt[6] = {
        strainTpdt[0] - meanStrainTpdt,
        strainTpdt[1] - meanStrainTpdt,
        strainTpdt[2] - meanStrainTpdt,
        strainTpdt[3],
        strainTpdt[4],
        strainTpdt[5]
    };

    const PylithScalar devStrainT[6] = {
        totalStrain[0] - meanStrainT,
        totalStrain[1] - meanStrainT,
        totalStrain[2] - meanStrainT,
        totalStrain[3],
        totalStrain[4],
        totalStrain[5]
    };

    PylithScalar strainDiff = 0.0;
    for (int iComp = 0; iComp < _strainSize; ++iComp) {
        strainDiff = devStrainTpdt[iComp] - devStrainT[iComp];
        visStrainTpdt[iComp] = expFac_1 * viscousStrain_1[iComp] + dq_1 * strainDiff;
        visStrainTpdt[iComp + _strainSize] = expFac_2 * viscousStrain_2[iComp] + dq_2 * strainDiff;
        visStrainTpdt[iComp + 2 * _strainSize] = expFac_3 * viscousStrain_3[iComp] + dq_3 * strainDiff;
    } // for
} // computeViscousStrain


// ---------------------------------------------------------------------------------------------------------------------
/* Update total strain for 3-D isotropic linear generalized Maxwell viscoelasticity.
 *
 * IMPORTANT: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::updateTotalStrain(const PylithInt dim,
                                                                  const PylithInt numS,
                                                                  const PylithInt numA,
                                                                  const PylithInt sOff[],
                                                                  const PylithInt sOff_x[],
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
                                                                  PylithScalar totalStrain[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 2;

#if 0 // :DEBUG:
    std::cout << "dim:  " << dim << std::endl;
    std::cout << "numS:  " << numS << std::endl;
    std::cout << "numA:  " << numA << std::endl;
    std::cout << "sOff[0]:  " << sOff[0] << std::endl;
    std::cout << "sOff_x[0]:  " << sOff_x[0] << std::endl;
    std::cout << "s[0]:  " << s[0] << std::endl;
    std::cout << "aOff[0]:  " << aOff[0] << std::endl;
    std::cout << "a[0]:  " << a[0] << std::endl;
    std::cout << "t:  " << t << std::endl;
    std::cout << "x[0]:  " << x[0] << std::endl;
    std::cout << "numConstants:  " << numConstants << std::endl;
    std::cout << "totalStrain[0]:  " << totalStrain[0] << std::endl;
#endif

    assert(_dim == dim);
    assert(numS >= 1);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s_x);
    assert(totalStrain);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    totalStrain[0] = disp_x[0*_dim+0];
    totalStrain[1] = disp_x[1*_dim+1];
    totalStrain[2] = disp_x[2*_dim+2];
    totalStrain[3] = 0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]);
    totalStrain[4] = 0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]);
    totalStrain[5] = 0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0]);
} // updateTotalStrain


// ---------------------------------------------------------------------------------------------------------------------
/* Update viscous strain for 3-D isotropic linear generalized Maxwell viscoelasticity.
 *
 * IMPORTANT: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::updateViscousStrain(const PylithInt dim,
                                                                    const PylithInt numS,
                                                                    const PylithInt numA,
                                                                    const PylithInt sOff[],
                                                                    const PylithInt sOff_x[],
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
                                                                    PylithScalar visStrain[]) {
    const PylithInt _dim = 3;
    const PylithInt _strainSize = 6;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = numA-4;
    const PylithInt i_shearModulusRatio = numA-3;
    const PylithInt i_viscousStrainPrev = numA-2;
    const PylithInt i_totalStrainPrev = numA-1;

#if 0 // :DEBUG:
    std::cout << "dim:  " << dim << std::endl;
    std::cout << "numS:  " << numS << std::endl;
    std::cout << "numA:  " << numA << std::endl;
    std::cout << "sOff[0]:  " << sOff[0] << std::endl;
    std::cout << "sOff_x[0]:  " << sOff_x[0] << std::endl;
    std::cout << "s[0]:  " << s[0] << std::endl;
    std::cout << "aOff[0]:  " << aOff[0] << std::endl;
    std::cout << "a[0]:  " << a[0] << std::endl;
    std::cout << "t:  " << t << std::endl;
    std::cout << "x[0]:  " << x[0] << std::endl;
    std::cout << "numConstants:  " << numConstants << std::endl;
    std::cout << "visStrain[0]:  " << visStrain[0] << std::endl;
#endif

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s_x);
    assert(aOff);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrainPrev] >= 0);
    assert(aOff[i_totalStrainPrev] >= 0);
    assert(a);
    assert(visStrain);
    assert(constants);
    assert(1 == numConstants);

    // Compute strain, deviatoric strain, etc.
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    const PylithScalar maxwellTime_1 = a[aOff[i_maxwellTime]+0];
    const PylithScalar maxwellTime_2 = a[aOff[i_maxwellTime]+1];
    const PylithScalar maxwellTime_3 = a[aOff[i_maxwellTime]+2];
    const PylithScalar* viscousStrain_1 = &a[aOff[i_viscousStrainPrev]];
    const PylithScalar* viscousStrain_2 = &a[aOff[i_viscousStrainPrev] + _strainSize];
    const PylithScalar* viscousStrain_3 = &a[aOff[i_viscousStrainPrev] + 2*_strainSize];
    const PylithScalar* totalStrainPrev = &a[aOff[i_totalStrainPrev]];

    const PylithScalar dt = constants[0];

    const PylithScalar dq_1 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_1);
    const PylithScalar expFac_1 = exp(-dt/maxwellTime_1);

    const PylithScalar dq_2 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_2);
    const PylithScalar expFac_2 = exp(-dt/maxwellTime_2);

    const PylithScalar dq_3 = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime_3);
    const PylithScalar expFac_3 = exp(-dt/maxwellTime_3);

    const PylithScalar strain[6] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        disp_x[2*_dim+2],
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]),
        0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]),
        0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0])
    };
    const PylithReal meanStrain = (strain[0] + strain[1] + strain[2]) / 3.0;

    const PylithScalar devStrain[6] = {
        strain[0] - meanStrain,
        strain[1] - meanStrain,
        strain[2] - meanStrain,
        strain[3],
        strain[4],
        strain[5]
    };

    const PylithReal meanStrainPrev = (totalStrainPrev[0] + totalStrainPrev[1] + totalStrainPrev[2])/3.0;

    const PylithScalar devStrainPrev[6] = {
        totalStrainPrev[0] - meanStrainPrev,
        totalStrainPrev[1] - meanStrainPrev,
        totalStrainPrev[2] - meanStrainPrev,
        totalStrainPrev[3],
        totalStrainPrev[4],
        totalStrainPrev[5]
    };

    // Compute viscous strain.
    PylithScalar strainDiff = 0.0;
    for (int iComp = 0; iComp < 6; ++iComp) {
        strainDiff = devStrain[iComp] - devStrainPrev[iComp];
        visStrain[iComp] = expFac_1*viscousStrain_1[iComp] + dq_1*strainDiff;
        visStrain[iComp + _strainSize] = expFac_2*viscousStrain_2[iComp] + dq_2*strainDiff;
        visStrain[iComp + 2*_strainSize] = expFac_3*viscousStrain_3[iComp] + dq_3*strainDiff;
    } // for
} // updateViscousStrain


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate stress for 3-D  isotropic linear generalized
 * Maxwell WITHOUT a reference stress and strain.
 *
 * Used in outputing the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), maxwell_time(3), shear_modulus_ratio(3),
 * total_strain(4), viscous_strain(4)]
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::cauchyStress(const PylithInt dim,
                                                             const PylithInt numS,
                                                             const PylithInt numA,
                                                             const PylithInt sOff[],
                                                             const PylithInt sOff_x[],
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
    const PylithInt _strainSize = 6;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-6;
    const PylithInt i_bulkModulus = numA-5;
    const PylithInt i_maxwellTime = numA-4;
    const PylithInt i_shearModulusRatio = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 6);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    // Shear modulus and (updated) viscous strains.
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* visStrainTpdt_1 = &a[aOff[i_viscousStrain]];
    const PylithScalar* visStrainTpdt_2 = &a[aOff[i_viscousStrain] + _strainSize];
    const PylithScalar* visStrainTpdt_3 = &a[aOff[i_viscousStrain] + 2*_strainSize];
    const PylithScalar* totalStrainTpdt = &a[aOff[i_totalStrain]];

    PylithScalar stressTensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    IsotropicLinearElasticity3D::meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                            aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);

    const PylithScalar meanStress = stressTensor[0];

    const PylithReal meanStrainTpdt = (totalStrainTpdt[0] + totalStrainTpdt[1] + totalStrainTpdt[2])/3.0;

    const PylithScalar devStrainTpdt[6] = {
        totalStrainTpdt[0] - meanStrainTpdt,
        totalStrainTpdt[1] - meanStrainTpdt,
        totalStrainTpdt[2] - meanStrainTpdt,
        totalStrainTpdt[3],
        totalStrainTpdt[4],
        totalStrainTpdt[5]
    };

    // Shear modulus ratio factors.
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio] + 1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio] + 2];
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2 - shearModulusRatio_3;

    // Stresses.
    stressVector[0] = meanStress + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[0] +
                                                     shearModulusRatio_1*visStrainTpdt_1[0] +
                                                     shearModulusRatio_2*visStrainTpdt_2[0] +
                                                     shearModulusRatio_3*visStrainTpdt_3[0]); // stress_xx
    stressVector[1] = meanStress + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[1] +
                                                     shearModulusRatio_1*visStrainTpdt_1[1] +
                                                     shearModulusRatio_2*visStrainTpdt_2[1] +
                                                     shearModulusRatio_3*visStrainTpdt_3[1]); // stress_yy
    stressVector[2] = meanStress + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[2] +
                                                     shearModulusRatio_1*visStrainTpdt_1[2] +
                                                     shearModulusRatio_2*visStrainTpdt_2[2] +
                                                     shearModulusRatio_3*visStrainTpdt_3[2]); // stress_zz
    stressVector[3] = 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[3] +
                                        shearModulusRatio_1*visStrainTpdt_1[3] +
                                        shearModulusRatio_2*visStrainTpdt_2[3] +
                                        shearModulusRatio_3*visStrainTpdt_3[3]); // stress_xy
    stressVector[4] = 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[4] +
                                        shearModulusRatio_1*visStrainTpdt_1[4] +
                                        shearModulusRatio_2*visStrainTpdt_2[4] +
                                        shearModulusRatio_3*visStrainTpdt_3[4]); // stress_yz
    stressVector[5] = 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[5] +
                                        shearModulusRatio_1*visStrainTpdt_1[5] +
                                        shearModulusRatio_2*visStrainTpdt_2[5] +
                                        shearModulusRatio_3*visStrainTpdt_3[5]); // stress_xz

} // cauchyStress


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate stress for 3-D isotropic linear generalized
 * Maxwell WITH a reference stress/strain.
 *
 * Used in outputing the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
 *                   maxwell_time(3), shear_modulus_ratio(3), total_strain(4), viscous_strain(4)]
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::cauchyStress_refstate(const PylithInt dim,
                                                                      const PylithInt numS,
                                                                      const PylithInt numA,
                                                                      const PylithInt sOff[],
                                                                      const PylithInt sOff_x[],
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
    const PylithInt _strainSize = 6;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-8;
    const PylithInt i_rstrain = numA-7;
    const PylithInt i_shearModulus = numA-6;
    const PylithInt i_bulkModulus = numA-5;
    const PylithInt i_maxwellTime = numA-4;
    const PylithInt i_shearModulusRatio = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 8);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_shearModulusRatio] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Pass bulk modulus, reference stress, and reference strain.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    // Shear modulus and (updated) viscous strains.
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* visStrainTpdt_1 = &a[aOff[i_viscousStrain]];
    const PylithScalar* visStrainTpdt_2 = &a[aOff[i_viscousStrain] + _strainSize];
    const PylithScalar* visStrainTpdt_3 = &a[aOff[i_viscousStrain] + 2*_strainSize];
    const PylithScalar* totalStrainTpdt = &a[aOff[i_totalStrain]];

    PylithScalar stressTensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    IsotropicLinearElasticity3D::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                     aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants,
                                                     stressTensor);

    const PylithScalar meanStress = stressTensor[0];

    // Reference stress and strain.
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy, stress_yz,
                                                         // stress_xz
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy, strain_yz,
                                                         // strain_xz

    // Compute reference deviatoric values.
    const PylithReal meanRefStrain = (refstrain[0] + refstrain[1] + refstrain[2])/3.0;
    const PylithScalar devRefStrain[6] = {refstrain[0] - meanRefStrain,
                                          refstrain[1] - meanRefStrain,
                                          refstrain[2] - meanRefStrain,
                                          refstrain[3],
                                          refstrain[4],
                                          refstrain[5]};
    const PylithReal meanRefStress = (refstress[0] + refstress[1] + refstress[2])/3.0;
    const PylithScalar devRefStress[6] = {refstress[0] - meanRefStress,
                                          refstress[1] - meanRefStress,
                                          refstress[2] - meanRefStress,
                                          refstress[3],
                                          refstress[4],
                                          refstress[5]};

    const PylithReal meanStrainTpdt = (totalStrainTpdt[0] + totalStrainTpdt[1] + totalStrainTpdt[2])/3.0;

    const PylithScalar devStrainTpdt[6] = {
        totalStrainTpdt[0] - meanStrainTpdt,
        totalStrainTpdt[1] - meanStrainTpdt,
        totalStrainTpdt[2] - meanStrainTpdt,
        totalStrainTpdt[3],
        totalStrainTpdt[4],
        totalStrainTpdt[5]
    };

    // Shear modulus ratio factors.
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio] + 1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio] + 2];
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2 - shearModulusRatio_3;

    // Stresses.
    stressVector[0] = meanStress + devRefStress[0] + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[0] +
                                                                       shearModulusRatio_1*visStrainTpdt_1[0] +
                                                                       shearModulusRatio_2*visStrainTpdt_2[0] +
                                                                       shearModulusRatio_3*visStrainTpdt_3[0] -
                                                                       devRefStrain[0]); // stress_xx
    stressVector[1] = meanStress + devRefStress[1] + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[1] +
                                                                       shearModulusRatio_1*visStrainTpdt_1[1] +
                                                                       shearModulusRatio_2*visStrainTpdt_2[1] +
                                                                       shearModulusRatio_3*visStrainTpdt_3[1] -
                                                                       devRefStrain[1]); // stress_yy
    stressVector[2] = meanStress + devRefStress[2] + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[2] +
                                                                       shearModulusRatio_1*visStrainTpdt_1[2] +
                                                                       shearModulusRatio_2*visStrainTpdt_2[2] +
                                                                       shearModulusRatio_3*visStrainTpdt_3[2] -
                                                                       devRefStrain[2]); // stress_zz
    stressVector[3] = devRefStress[3] + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[3] +
                                                          shearModulusRatio_1*visStrainTpdt_1[3] +
                                                          shearModulusRatio_2*visStrainTpdt_2[3] +
                                                          shearModulusRatio_3*visStrainTpdt_3[3] -
                                                          devRefStrain[3]); // stress_xy
    stressVector[4] = devRefStress[4] + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[4] +
                                                          shearModulusRatio_1*visStrainTpdt_1[4] +
                                                          shearModulusRatio_2*visStrainTpdt_2[4] +
                                                          shearModulusRatio_3*visStrainTpdt_3[4] -
                                                          devRefStrain[4]); // stress_yz
    stressVector[5] = devRefStress[5] + 2.0*shearModulus*(shearModulusRatio_0*devStrainTpdt[5] +
                                                          shearModulusRatio_1*visStrainTpdt_1[5] +
                                                          shearModulusRatio_2*visStrainTpdt_2[5] +
                                                          shearModulusRatio_3*visStrainTpdt_3[5] -
                                                          devRefStrain[5]); // stress_xz

} // cauchyStress_refstate


// End of file
