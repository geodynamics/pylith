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

#include "IsotropicLinearMaxwell.hh" // Implementation of object methods.

#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity* kernels
#include "pylith/fekernels/Viscoelasticity.hh" // USES Viscoelasticity kernels

#include <cassert> // USES assert()
#include <cmath> // USES exp()

// =====================================================================================================================
// Kernels for isotropic, linear Maxwell viscoelastic plane strain.
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
// f1 function for isotropic linear Maxwell plane strain material WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f1v(const PylithInt dim,
                                                          const PylithInt numS,
                                                          const PylithInt numA,
                                                          const PylithInt sOff[],
                                                          const PylithInt sOff_x[],
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
    const PylithInt i_shearModulus = numA-5;
    const PylithInt i_bulkModulus = numA-4;
    const PylithInt i_maxwellTime = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

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
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(f1);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 4; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[4] = {
        aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain]
    };

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
    IsotropicLinearElasticityPlaneStrain::meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                     aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v


// ---------------------------------------------------------------------------------------------------------------------
// f1 function for isotropic linear Maxwell viscoelastic plane strain material with reference
// stress and strain.
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f1v_refstate(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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
    const PylithInt i_rstress = numA-7;
    const PylithInt i_rstrain = numA-6;
    const PylithInt i_shearModulus = numA-5;
    const PylithInt i_bulkModulus = numA-4;
    const PylithInt i_maxwellTime = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 7);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f1);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Pass bulk modulus, reference stress, and reference strain.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 6; // Pass shear modulus, Maxwell time, viscous strain, total strain,
                                 // reference stress, and reference strain.
    const PylithInt aOffDev[6] = {aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_maxwellTime],
                                  aOff[i_viscousStrain], aOff[i_totalStrain] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    IsotropicLinearElasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                              aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Jf3_vu entry function for 2-D plane strain isotropic linear Maxwell viscoelastic material WITHOUT reference
 * stress/strain.
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
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::Jf3vu(const PylithInt dim,
                                                            const PylithInt numS,
                                                            const PylithInt numA,
                                                            const PylithInt sOff[],
                                                            const PylithInt sOff_x[],
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

} // Jf3vu


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * Maxwell viscoelastic material WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = 2*shearModulus*visStrain_ij
 */
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::deviatoricStress(const PylithInt dim,
                                                                       const PylithInt numS,
                                                                       const PylithInt numA,
                                                                       const PylithInt sOff[],
                                                                       const PylithInt sOff_x[],
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
    const PylithInt i_maxwellTime = 1;
    const PylithInt i_viscousStrain = 2;
    const PylithInt i_totalStrain = 3;

    assert(_dim == dim);
    assert(1 == numS);
    assert(4 == numA);
    assert(sOff);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithInt _numS = 1; // Number passed on to computeViscousStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; // Number passed on to computeViscousStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain] };

    PylithScalar visStrainTpdt[4] = {0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.
    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis, NULL, a, a_t, NULL,
                         t, x, numConstants, constants, visStrainTpdt);
    stress[0] += 2.0 * shearModulus * visStrainTpdt[0]; // stress_xx
    stress[1] += 2.0 * shearModulus * visStrainTpdt[3]; // stress_xy
    stress[2] += 2.0 * shearModulus * visStrainTpdt[3]; // stress_yx
    stress[3] += 2.0 * shearModulus * visStrainTpdt[1]; // stress_yy
} // deviatoricStress


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * Maxwell viscoelastic material WITH reference stress and reference strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = devrefstress_ij + 2*shearModulus*(visstrain_ij - devrefstrain_ij)
 */
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::deviatoricStress_refstate(const PylithInt dim,
                                                                                const PylithInt numS,
                                                                                const PylithInt numA,
                                                                                const PylithInt sOff[],
                                                                                const PylithInt sOff_x[],
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
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_viscousStrain = 4;
    const PylithInt i_totalStrain = 5;

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
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithInt _numS = 1; // Number passed on to computeViscousStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; // Number passed on to computeViscousStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain] };

    PylithScalar visStrainTpdt[4] = {0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    // Compute viscous strain for current time step.
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

    // Compute stress components -- note that we are including reference deviatoric
    // stress for now. This may need to be removed after testing.
    const PylithReal twomu = 2.0*shearModulus;

    const PylithScalar stress_xx = devRefStress[0] + twomu * (visStrainTpdt[0] - devRefStrain[0]);
    const PylithScalar stress_yy = devRefStress[1] + twomu * (visStrainTpdt[1] - devRefStrain[1]);
    const PylithScalar stress_xy = devRefStress[3] + twomu * (visStrainTpdt[3] - devRefStrain[3]);

    stress[0*_dim+0] += stress_xx;
    stress[0*_dim+1] += stress_xy;
    stress[1*_dim+0] += stress_xy;
    stress[1*_dim+1] += stress_yy;
} // deviatoricStress_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate viscous strain for a Maxwell plane strain viscoelastic material.
 *
 */
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::computeViscousStrain(const PylithInt dim,
                                                                           const PylithInt numS,
                                                                           const PylithInt numA,
                                                                           const PylithInt sOff[],
                                                                           const PylithInt sOff_x[],
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
    assert(visStrainTpdt);
    assert(1 == numConstants);
    assert(constants);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];
    const PylithScalar* viscousStrain = &a[aOff[i_viscousStrain]];
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];

    const PylithScalar dt = constants[0];

    const PylithScalar dq = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime);
    const PylithScalar expFac = exp(-dt/maxwellTime);

    const PylithScalar strainTpdt[4] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        0.0,
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0])
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0;
    const PylithReal meanStrainT = (totalStrain[0] + totalStrain[1])/3.0;
#if 0 // :DEBUG:
    std::cout << "totalStrain[0]:  " << totalStrain[0] << std::endl;
    std::cout << "strainTpdt[0]:  " << strainTpdt[0] << std::endl;
#endif

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

    for (int iComp = 0; iComp < 4; ++iComp) {
        visStrainTpdt[iComp] = expFac * viscousStrain[iComp] + dq * (devStrainTpdt[iComp] - devStrainT[iComp]);
    } // for
} // computeViscousStrain


// ---------------------------------------------------------------------------------------------------------------------
/* Update total strain for a Maxwell plane strain viscoelastic material.
 *
 * IMPORTANT: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 *
 */
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::updateTotalStrain(const PylithInt dim,
                                                                        const PylithInt numS,
                                                                        const PylithInt numA,
                                                                        const PylithInt sOff[],
                                                                        const PylithInt sOff_x[],
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
/* Update viscous strain for a Maxwell plane strain viscoelastic material.
 *
 * :ATTENSION: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 */
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::updateViscousStrain(const PylithInt dim,
                                                                          const PylithInt numS,
                                                                          const PylithInt numA,
                                                                          const PylithInt sOff[],
                                                                          const PylithInt sOff_x[],
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = numA-3;
    const PylithInt i_viscousStrainPrev = numA-2;
    const PylithInt i_totalStrainPrev = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrainPrev] >= 0);
    assert(aOff[i_totalStrainPrev] >= 0);
    assert(s_x);
    assert(a);
    assert(visStrain);
    assert(constants);
    assert(1 == numConstants);

    // Compute strain, deviatoric strain, etc.
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];
    const PylithScalar* viscousStrainPrev = &a[aOff[i_viscousStrainPrev]];
    const PylithScalar* totalStrainPrev = &a[aOff[i_totalStrainPrev]];

    const PylithScalar dt = constants[0];

    const PylithScalar dq = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime);
    const PylithScalar expFac = exp(-dt/maxwellTime);

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

    for (int iComp = 0; iComp < 4; ++iComp) {
        visStrain[iComp] = expFac * viscousStrainPrev[iComp] + dq * (devStrain[iComp] - devStrainPrev[iComp]);
    } // for

} // updateViscousStrain


// ---------------------------------------------------------------------------------------------------------------------
// Calculate stress for 2-D plane strain isotropic linear Maxwell material WITHOUT a reference stress and strain.
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::cauchyStress(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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
    const PylithInt i_shearModulus = numA-5;
    const PylithInt i_bulkModulus = numA-4;
    const PylithInt i_maxwellTime = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

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
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    // Shear modulus and (updated) viscous strain.
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* visStrainTpdt = &a[aOff[i_viscousStrain]];

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 }; // Full stress tensor in vector form.
    IsotropicLinearElasticityPlaneStrain::meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                     aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    const PylithScalar meanStress = stressTensor[0];

    stressVector[0] = meanStress + 2.0 * shearModulus * visStrainTpdt[0]; // stress_xx
    stressVector[1] = meanStress + 2.0 * shearModulus * visStrainTpdt[1]; // stress_yy
    stressVector[2] = meanStress + 2.0 * shearModulus * visStrainTpdt[2]; // stress_zz
    stressVector[3] = 2.0 * shearModulus * visStrainTpdt[3]; // stress_xy

} // cauchyStress


// ---------------------------------------------------------------------------------------------------------------------
// Calculate stress for 2-D plane strain isotropic linear Maxwell material WITH a reference stress and strain.
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::cauchyStress_refstate(const PylithInt dim,
                                                                            const PylithInt numS,
                                                                            const PylithInt numA,
                                                                            const PylithInt sOff[],
                                                                            const PylithInt sOff_x[],
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
    const PylithInt i_rstress = numA-7;
    const PylithInt i_rstrain = numA-6;
    const PylithInt i_shearModulus = numA-5;
    const PylithInt i_bulkModulus = numA-4;
    const PylithInt i_maxwellTime = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 7);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    // Shear modulus and (updated) viscous strain.
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* visStrainTpdt = &a[aOff[i_viscousStrain]];

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 }; // Full stress tensor in vector form.
    IsotropicLinearElasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                              aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    const PylithScalar meanStress = stressTensor[0];

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

    stressVector[0] = meanStress + devRefStress[0] + 2.0 * shearModulus * (visStrainTpdt[0] - devRefStrain[0]); // stress_xx
    stressVector[1] = meanStress + devRefStress[1] + 2.0 * shearModulus * (visStrainTpdt[1] - devRefStrain[1]); // stress_yy
    stressVector[2] = meanStress + devRefStress[2] + 2.0 * shearModulus * (visStrainTpdt[2] - devRefStrain[2]); // stress_zz
    stressVector[3] = devRefStress[3] + 2.0 * shearModulus * (visStrainTpdt[3] - devRefStrain[3]); // stress_xy

} // cauchyStress_refstate


// =====================================================================================================================
// Kernels for isotropic, linear Maxwell viscoelastic 3D material.
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
// f1 function for isotropic linear Maxwell 3D material WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearMaxwell3D::f1v(const PylithInt dim,
                                                 const PylithInt numS,
                                                 const PylithInt numA,
                                                 const PylithInt sOff[],
                                                 const PylithInt sOff_x[],
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
    const PylithInt i_shearModulus = numA-5;
    const PylithInt i_bulkModulus = numA-4;
    const PylithInt i_maxwellTime = numA-3;
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
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 4; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[4] = {
        aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain]
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
// f1 function for isotropic linear Maxwell viscoelastic 3D material with reference
// stress and strain.
void
pylith::fekernels::IsotropicLinearMaxwell3D::f1v_refstate(const PylithInt dim,
                                                          const PylithInt numS,
                                                          const PylithInt numA,
                                                          const PylithInt sOff[],
                                                          const PylithInt sOff_x[],
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
    const PylithInt i_rstress = numA-7;
    const PylithInt i_rstrain = numA-6;
    const PylithInt i_shearModulus = numA-5;
    const PylithInt i_bulkModulus = numA-4;
    const PylithInt i_maxwellTime = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 7);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Pass bulk modulus, reference stress, and reference strain.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 6; // Pass shear modulus, Maxwell time, viscous strain, total strain,
                                 // reference stress, and reference strain.
    const PylithInt aOffDev[6] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_maxwellTime],
                                   aOff[i_viscousStrain], aOff[i_totalStrain] };

    PylithScalar stressTensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    IsotropicLinearElasticity3D::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                     aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Jf3_vu entry function for 3-D isotropic linear Maxwell viscoelastic material.
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
pylith::fekernels::IsotropicLinearMaxwell3D::Jf3vu(const PylithInt dim,
                                                   const PylithInt numS,
                                                   const PylithInt numA,
                                                   const PylithInt sOff[],
                                                   const PylithInt sOff_x[],
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
} // Jf3vu


// ----------------------------------------------------------------------
/* Calculate deviatoric stress for 3-D isotropic linear
 * Maxwell viscoelastic material WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = 2*shearModulus*visStrain_ij
 */
void
pylith::fekernels::IsotropicLinearMaxwell3D::deviatoricStress(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_maxwellTime = 1;
    const PylithInt i_viscousStrain = 2;
    const PylithInt i_totalStrain = 3;

    assert(_dim == dim);
    assert(1 == numS);
    assert(4 == numA);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s_x);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithInt _numS = 1; // Number passed on to visStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain] };

    PylithScalar visStrainTpdt[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor (vector).

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis, NULL, a, a_t, NULL,
                         t, x, numConstants, constants, visStrainTpdt);

    stress[0] += 2.0 * shearModulus * visStrainTpdt[0]; // stress_xx
    stress[1] += 2.0 * shearModulus * visStrainTpdt[3]; // stress_xy
    stress[2] += 2.0 * shearModulus * visStrainTpdt[5]; // stress_xz
    stress[3] += stress[1]; // stress_yx
    stress[4] += 2.0 * shearModulus * visStrainTpdt[1]; // stress_yy
    stress[5] += 2.0 * shearModulus * visStrainTpdt[4]; // stress_yz
    stress[6] += stress[2]; // stress_zx
    stress[7] += stress[5]; // stress_zy
    stress[8] += 2.0 * shearModulus * visStrainTpdt[2]; // stress_zz
} // deviatoricStress


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 3-D isotropic linear
 * Maxwell viscoelastic material WITH reference stress and reference strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = devrefstress_ij + 2*shearModulus*(visstrain_ij - devrefstrain_ij)
 */
void
pylith::fekernels::IsotropicLinearMaxwell3D::deviatoricStress_refstate(const PylithInt dim,
                                                                       const PylithInt numS,
                                                                       const PylithInt numA,
                                                                       const PylithInt sOff[],
                                                                       const PylithInt sOff_x[],
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = 0;
    const PylithInt i_rstrain = 1;
    const PylithInt i_shearModulus = 2;
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_viscousStrain = 4;
    const PylithInt i_totalStrain = 5;

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
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
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

    const PylithInt numAVis = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain] };

    PylithScalar visStrainTpdt[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor (vector).

    // Compute viscous strain for current time step.
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

    // Compute stress components -- note that we are including reference deviatoric
    // stress for now. This may need to be removed after testing.
    const PylithReal twomu = 2.0*shearModulus;

    const PylithScalar stress_xx = devRefStress[0] + twomu * (visStrainTpdt[0] - devRefStrain[0]);
    const PylithScalar stress_yy = devRefStress[1] + twomu * (visStrainTpdt[1] - devRefStrain[1]);
    const PylithScalar stress_zz = devRefStress[2] + twomu * (visStrainTpdt[2] - devRefStrain[2]);
    const PylithScalar stress_xy = devRefStress[3] + twomu * (visStrainTpdt[3] - devRefStrain[3]);
    const PylithScalar stress_yz = devRefStress[4] + twomu * (visStrainTpdt[4] - devRefStrain[4]);
    const PylithScalar stress_xz = devRefStress[5] + twomu * (visStrainTpdt[5] - devRefStrain[5]);

    stress[0*_dim+0] += stress_xx;
    stress[0*_dim+1] += stress_xy;
    stress[0*_dim+2] += stress_xz;
    stress[1*_dim+0] += stress_xy;
    stress[1*_dim+1] += stress_yy;
    stress[1*_dim+2] += stress_yz;
    stress[2*_dim+0] += stress_xz;
    stress[2*_dim+1] += stress_yz;
    stress[2*_dim+2] += stress_zz;
} // deviatoricStress_refstate


// ---------------------------------------------------------------------------------------------------------------------
// Calculate viscous strain for a 3D Maxwell viscoelastic material.
void
pylith::fekernels::IsotropicLinearMaxwell3D::computeViscousStrain(const PylithInt dim,
                                                                  const PylithInt numS,
                                                                  const PylithInt numA,
                                                                  const PylithInt sOff[],
                                                                  const PylithInt sOff_x[],
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
    assert(1 == numConstants);
    assert(constants);
    assert(visStrainTpdt);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];
    const PylithScalar* viscousStrain = &a[aOff[i_viscousStrain]];
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];

    const PylithScalar dt = constants[0];

    const PylithScalar dq = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime);
    const PylithScalar expFac = exp(-dt/maxwellTime);

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
#if 0 // :DEBUG:
    std::cout << "totalStrain[0]:  " << totalStrain[0] << std::endl;
    std::cout << "strainTpdt[0]:  " << strainTpdt[0] << std::endl;
#endif

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

    for (int iComp = 0; iComp < 6; ++iComp) {
        visStrainTpdt[iComp] = expFac * viscousStrain[iComp] + dq * (devStrainTpdt[iComp] - devStrainT[iComp]);
    } // for
} // computeViscousStrain


// ---------------------------------------------------------------------------------------------------------------------
/* Update total strain for a 3D Maxwell viscoelastic material.
 *
 * :ATTENSION: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 */
void
pylith::fekernels::IsotropicLinearMaxwell3D::updateTotalStrain(const PylithInt dim,
                                                               const PylithInt numS,
                                                               const PylithInt numA,
                                                               const PylithInt sOff[],
                                                               const PylithInt sOff_x[],
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
    const PylithInt i_disp = 0;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(sOff);
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
/* Update viscous strain for a Maxwell viscoelastic material.
 *
 * :ATTENSION: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 */
void
pylith::fekernels::IsotropicLinearMaxwell3D::updateViscousStrain(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = numA-3;
    const PylithInt i_viscousStrainPrev = numA-2;
    const PylithInt i_totalStrainPrev = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrainPrev] >= 0);
    assert(aOff[i_totalStrainPrev] >= 0);
    assert(s_x);
    assert(a);
    assert(visStrain);
    assert(constants);
    assert(1 == numConstants);

    // Compute strain, deviatoric strain, etc.
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];
    const PylithScalar* viscousStrainPrev = &a[aOff[i_viscousStrainPrev]];
    const PylithScalar* totalStrainPrev = &a[aOff[i_totalStrainPrev]];

    const PylithScalar dt = constants[0];

    const PylithScalar dq = pylith::fekernels::Viscoelasticity::maxwellViscousStrainCoeff(dt, maxwellTime);
    const PylithScalar expFac = exp(-dt/maxwellTime);

    const PylithScalar strain[6] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        disp_x[2*_dim+2],
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]),
        0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]),
        0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0])
    };
    const PylithReal meanStrain = (strain[0] + strain[1] + strain[2])/3.0;

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
    for (int iComp = 0; iComp < 6; ++iComp) {
        visStrain[iComp] = expFac * viscousStrainPrev[iComp] + dq * (devStrain[iComp] - devStrainPrev[iComp]);
    } // for

} // updateViscousStrain


// ---------------------------------------------------------------------------------------------------------------------
// Calculate stress for 3-D isotropic linear Maxwell material WITHOUT a reference stress and strain.
void
pylith::fekernels::IsotropicLinearMaxwell3D::cauchyStress(const PylithInt dim,
                                                          const PylithInt numS,
                                                          const PylithInt numA,
                                                          const PylithInt sOff[],
                                                          const PylithInt sOff_x[],
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
    const PylithInt i_shearModulus = numA-5;
    const PylithInt i_bulkModulus = numA-4;
    const PylithInt i_maxwellTime = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

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
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    // Shear modulus and (updated) viscous strain.
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* visStrainTpdt = &a[aOff[i_viscousStrain]];

    PylithScalar stressTensor[9] = { // Full stress tensor
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
    };
    IsotropicLinearElasticity3D::meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                            aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);

    const PylithScalar meanStress = stressTensor[0];

    stressVector[0] = meanStress + 2.0 * shearModulus * visStrainTpdt[0]; // stress_xx
    stressVector[1] = meanStress + 2.0 * shearModulus * visStrainTpdt[1]; // stress_yy
    stressVector[2] = meanStress + 2.0 * shearModulus * visStrainTpdt[2]; // stress_zz
    stressVector[3] = 2.0 * shearModulus * visStrainTpdt[3]; // stress_xy
    stressVector[4] = 2.0 * shearModulus * visStrainTpdt[4]; // stress_yz
    stressVector[5] = 2.0 * shearModulus * visStrainTpdt[5]; // stress_xz

} // stress


// ---------------------------------------------------------------------------------------------------------------------
// Calculate stress for 2-D plane strain isotropic linear Maxwell material WITH a reference stress/strain.
void
pylith::fekernels::IsotropicLinearMaxwell3D::cauchyStress_refstate(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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
    const PylithInt i_rstress = numA-7;
    const PylithInt i_rstrain = numA-6;
    const PylithInt i_shearModulus = numA-5;
    const PylithInt i_bulkModulus = numA-4;
    const PylithInt i_maxwellTime = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_totalStrain = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 7);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Pass bulk modulus, reference stress, and reference strain.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    // Shear modulus and (updated) viscous strain.
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* visStrainTpdt = &a[aOff[i_viscousStrain]];

    PylithScalar stressTensor[9] = { // Full stress tensor
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
    };
    IsotropicLinearElasticity3D::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                     aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);

    const PylithScalar meanStress = stressTensor[0];

    // Reference stress and strain values.
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

    // Current stress values.
    stressVector[0] = meanStress + devRefStress[0] + 2.0*shearModulus*(visStrainTpdt[0] - devRefStrain[0]); // stress_xx
    stressVector[1] = meanStress + devRefStress[1] + 2.0*shearModulus*(visStrainTpdt[1] - devRefStrain[1]); // stress_yy
    stressVector[2] = meanStress + devRefStress[2] + 2.0*shearModulus*(visStrainTpdt[2] - devRefStrain[2]); // stress_zz
    stressVector[3] = devRefStress[3] + 2.0*shearModulus*(visStrainTpdt[3] - devRefStrain[3]); // stress_xy
    stressVector[4] = devRefStress[4] + 2.0*shearModulus*(visStrainTpdt[4] - devRefStrain[4]); // stress_yz
    stressVector[5] = devRefStress[5] + 2.0*shearModulus*(visStrainTpdt[5] - devRefStrain[5]); // stress_xz
} // stress


// End of file
