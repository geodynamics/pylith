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

#include "IsotropicPowerLaw.hh" // Implementation of object methods.

#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity* kernels
#include "pylith/fekernels/Viscoelasticity.hh" // USES Viscoelasticity kernels

#include <cassert> // USES assert()
#include <cmath> // USES exp()
#include <stdexcept> // USES runtime_error

// =====================================================================================================================
// Kernels for isotropic power-law viscoelastic material.
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
// f1 function for isotropic power-law plane strain WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::f1v(const PylithInt dim,
                                                     const PylithInt numS,
                                                     const PylithInt numA,
                                                     const PylithInt sOff[],
                                                     const PylithInt sOff_x[],
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
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

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
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(f1);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 6; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[6] = {aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate], aOff[i_powerLawReferenceStress],
                                  aOff[i_powerLawExponent], aOff[i_viscousStrain], aOff[i_stress]};

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
// f1 function for isotropic power-law viscoelastic plane strain with reference
// stress and strain.
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::f1v_refstate(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
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
    const PylithInt i_rstress = numA-9;
    const PylithInt i_rstrain = numA-8;
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f1);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 8; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[8] = {aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate],
                                  aOff[i_powerLawReferenceStress], aOff[i_powerLawExponent], aOff[i_viscousStrain],
                                  aOff[i_stress]};

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    IsotropicLinearElasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean,
                                                              NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Jf3_vu entry function for 2-D plane strain isotropic power-law viscoelastic WITHOUT reference stress/strain.
 *
 * stress_ij = C_ijkl strain_kl
 *
 * For reference:
 *
 * Isotropic:
 *  C_ijkl = bulkModulus * delta_ij * delta_kl
 *   + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
 */
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::Jf3vu(const PylithInt dim,
                                                       const PylithInt numS,
                                                       const PylithInt numA,
                                                       const PylithInt sOff[],
                                                       const PylithInt sOff_x[],
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

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
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(constants);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 6; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[6] = {aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate], aOff[i_powerLawReferenceStress],
                                  aOff[i_powerLawExponent], aOff[i_viscousStrain], aOff[i_stress]};

    // Compute deviatoric stress (4 components).
    PylithScalar devStressTpdt[4] = { 0.0, 0.0, 0.0, 0.0 };
    deviatoricStress4(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                      t, x, numConstants, constants, devStressTpdt);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar powerLawReferenceStrainRate = a[aOff[i_powerLawReferenceStrainRate]];
    const PylithScalar powerLawReferenceStress = a[aOff[i_powerLawReferenceStress]];
    const PylithScalar powerLawExponent = a[aOff[i_powerLawExponent]];
    // ****** I think powerLawAlpha should actually be a kernel constant.
    const PylithScalar powerLawAlpha = 0.5;
    // const PylithScalar powerLawAlpha = a[aOff[i_powerLawAlpha]];
    const PylithScalar dt = constants[0];

    // Compute quantities based on stress at t = T.
    const PylithScalar* stressT = &a[aOff[i_stress]]; // stress_xx, stress_yy, stress_zz, stress_xy at t = T.
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[4] = {stressT[0] - meanStressT,
                                        stressT[1] - meanStressT,
                                        stressT[2] - meanStressT,
                                        stressT[3]};
    const PylithScalar devStressProdT = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(devStressT, devStressT);
    const PylithScalar j2T = sqrt(0.5*devStressProdT);

    // Compute quantities based on stress at t = T + dt.
    const PylithScalar devStressProdTpdt = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(devStressTpdt, devStressTpdt);
    const PylithScalar j2Tpdt = sqrt(0.5*devStressProdTpdt);

    // Compute quantities at intermediate time tau.
    const PylithScalar j2Tau = powerLawAlpha*j2Tpdt + (1.0 - powerLawAlpha)*j2T;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress),
                                                                  (powerLawExponent - 1.0))/powerLawReferenceStress;

    /* Unique components of Jacobian. */
    const PylithReal C1111 = bulkModulus -
                             2.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[0]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressT[0]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C1122 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[1]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressT[1]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C1212 =
        -1.0/(2.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[3]*devStressTpdt[3]*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) +
                   powerLawAlpha*dt*gammaTau +
                   powerLawAlpha*dt*gammaTau*devStressTpdt[3]*devStressT[3]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C2211 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[1]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressT[0]*devStressTpdt[1]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C2222 = bulkModulus -
                             2.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressTpdt[1]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressT[1]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    /* j(f,g,df,dg) = C(f,df,g,dg)
     * 0:  j0000 = C1111 = bulkModulus - 2/(3*(alpha**2*deltaT*gammaFTau*s11**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau + alpha*deltaT*gammaFTau*s11*s11T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 1:  j0001 = C1112 = 0
     * 2:  j0010 = C1211 = 0
     * 3:  j0011 = C1212 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 4:  j0100 = C1121 = 0
     * 5:  j0101 = C1122 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s11*s22*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s11*s22T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 6:  j0110 = C1221 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 7:  j0111 = C1222 = 0
     * 8:  j1000 = C2111 = 0
     * 9:  j1001 = C2112 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 10:  j1010 = C2211 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s11*s22*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s11T*s22*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 11:  j1011 = C2212 = 0
     * 12:  j1100 = C2121 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 13:  j1101 = C2122 = 0
     * 14:  j1110 = C2221 = 0
     * 15:  j1111 = C2222 = bulkModulus - 2/(3*(alpha**2*deltaT*gammaFTau*s22**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau + alpha*deltaT*gammaFTau*s22*s22T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt) + 1/(2*mu)))
     */

    /* Nonzero Jacobian entries. */
    Jf3[0] -= C1111; /* j0000 */
    Jf3[3] -= C1212; /* j0011 */
    Jf3[5] -= C1122; /* j0101 */
    Jf3[6] -= C1212; /* j0110 */
    Jf3[9] -= C1212; /* j1001 */
    Jf3[10] -= C2211; /* j1010 */
    Jf3[12] -= C1212; /* j1100 */
    Jf3[15] -= C2222; /* j1111 */

} // Jf3vu


// ---------------------------------------------------------------------------------------------------------------------
/* Jf3_vu entry function for 2-D plane strain isotropic power-law viscoelastic WITH reference stress/strain.
 *
 * stress_ij = C_ijkl strain_kl
 *
 * For reference:
 *
 * Isotropic:
 *  C_ijkl = bulkModulus * delta_ij * delta_kl
 *   + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
 */
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::Jf3vu_refstate(const PylithInt dim,
                                                                const PylithInt numS,
                                                                const PylithInt numA,
                                                                const PylithInt sOff[],
                                                                const PylithInt sOff_x[],
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-9;
    const PylithInt i_rstrain = numA-8;
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 8; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[8] = {aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate],
                                  aOff[i_powerLawReferenceStress], aOff[i_powerLawExponent], aOff[i_viscousStrain],
                                  aOff[i_stress]};

    // Compute deviatoric stress vector (4 components).
    PylithScalar devStressTpdt[4] = { 0.0, 0.0, 0.0, 0.0 };
    deviatoricStress4_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                               t, x, numConstants, constants, devStressTpdt);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar powerLawReferenceStrainRate = a[aOff[i_powerLawReferenceStrainRate]];
    const PylithScalar powerLawReferenceStress = a[aOff[i_powerLawReferenceStress]];
    const PylithScalar powerLawExponent = a[aOff[i_powerLawExponent]];
    // ****** I think powerLawAlpha should actually be a kernel constant.
    const PylithScalar powerLawAlpha = 0.5;
    // const PylithScalar powerLawAlpha = a[aOff[i_powerLawAlpha]];
    const PylithScalar dt = constants[0];

    // Compute quantities based on stress at t = T.
    const PylithScalar* stressT = &a[aOff[i_stress]]; // stress_xx, stress_yy, stress_zz, stress_xy at t = T.
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[4] = {stressT[0] - meanStressT,
                                        stressT[1] - meanStressT,
                                        stressT[2] - meanStressT,
                                        stressT[3]};
    const PylithScalar devStressProdT = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(devStressT, devStressT);
    const PylithScalar j2T = sqrt(0.5*devStressProdT);

    // Compute quantities based on stress at t = T + dt.
    const PylithScalar devStressProdTpdt = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(devStressTpdt, devStressTpdt);
    const PylithScalar j2Tpdt = sqrt(0.5*devStressProdTpdt);

    // Compute quantities at intermediate time tau.
    const PylithScalar j2Tau = powerLawAlpha*j2Tpdt + (1.0 - powerLawAlpha)*j2T;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress),
                                                                  (powerLawExponent - 1.0))/powerLawReferenceStress;

    /* Unique components of Jacobian. */
    const PylithReal C1111 = bulkModulus -
                             2.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[0]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressT[0]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C1122 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[1]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressT[1]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C1212 =
        -1.0/(2.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[3]*devStressTpdt[3]*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) +
                   powerLawAlpha*dt*gammaTau +
                   powerLawAlpha*dt*gammaTau*devStressTpdt[3]*devStressT[3]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C2211 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[1]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressT[0]*devStressTpdt[1]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C2222 = bulkModulus -
                             2.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressTpdt[1]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressT[1]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    /* j(f,g,df,dg) = C(f,df,g,dg)
     * 0:  j0000 = C1111 = bulkModulus - 2/(3*(alpha**2*deltaT*gammaFTau*s11**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau + alpha*deltaT*gammaFTau*s11*s11T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 1:  j0001 = C1112 = 0
     * 2:  j0010 = C1211 = 0
     * 3:  j0011 = C1212 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 4:  j0100 = C1121 = 0
     * 5:  j0101 = C1122 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s11*s22*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s11*s22T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 6:  j0110 = C1221 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 7:  j0111 = C1222 = 0
     * 8:  j1000 = C2111 = 0
     * 9:  j1001 = C2112 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 10:  j1010 = C2211 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s11*s22*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s11T*s22*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 11:  j1011 = C2212 = 0
     * 12:  j1100 = C2121 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 13:  j1101 = C2122 = 0
     * 14:  j1110 = C2221 = 0
     * 15:  j1111 = C2222 = bulkModulus - 2/(3*(alpha**2*deltaT*gammaFTau*s22**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau + alpha*deltaT*gammaFTau*s22*s22T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt) + 1/(2*mu)))
     */

    /* Nonzero Jacobian entries. */
    Jf3[0] -= C1111; /* j0000 */
    Jf3[3] -= C1212; /* j0011 */
    Jf3[5] -= C1122; /* j0101 */
    Jf3[6] -= C1212; /* j0110 */
    Jf3[9] -= C1212; /* j1001 */
    Jf3[10] -= C2211; /* j1010 */
    Jf3[12] -= C1212; /* j1100 */
    Jf3[15] -= C2222; /* j1111 */

} // Jf3vu_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic
 * power-law viscoelastic WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = 2*shearModulus*visStrain_ij
 */
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::deviatoricStress(const PylithInt dim,
                                                                  const PylithInt numS,
                                                                  const PylithInt numA,
                                                                  const PylithInt sOff[],
                                                                  const PylithInt sOff_x[],
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
    const PylithInt i_powerLawReferenceStrainRate = 1;
    const PylithInt i_powerLawReferenceStress = 2;
    const PylithInt i_powerLawExponent = 3;
    const PylithInt i_viscousStrain = 4;
    const PylithInt i_stress = 5;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 6);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 6; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[6] = {aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate], aOff[i_powerLawReferenceStress],
                                  aOff[i_powerLawExponent], aOff[i_viscousStrain], aOff[i_stress]};

    // Compute deviatoric stress vector (4 components).
    PylithScalar devStressTpdt[4] = { 0.0, 0.0, 0.0, 0.0 };
    deviatoricStress4(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                      t, x, numConstants, constants, devStressTpdt);

    stress[0] += devStressTpdt[0]; // stress_xx
    stress[1] += devStressTpdt[3]; // stress_xy
    stress[2] += devStressTpdt[3]; // stress_yx
    stress[3] += devStressTpdt[1]; // stress_yy
} // deviatoricStress


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic
 * power-law viscoelastic WITH reference stress and reference strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = devrefstress_ij + 2*shearModulus*(visstrain_ij - devrefstrain_ij)
 */
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::deviatoricStress_refstate(const PylithInt dim,
                                                                           const PylithInt numS,
                                                                           const PylithInt numA,
                                                                           const PylithInt sOff[],
                                                                           const PylithInt sOff_x[],
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
    const PylithInt i_powerLawReferenceStrainRate = 3;
    const PylithInt i_powerLawReferenceStress = 4;
    const PylithInt i_powerLawExponent = 5;
    const PylithInt i_viscousStrain = 6;
    const PylithInt i_stress = 7;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 8);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 8; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[8] = {aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate],
                                  aOff[i_powerLawReferenceStress], aOff[i_powerLawExponent], aOff[i_viscousStrain],
                                  aOff[i_stress]};

    // Compute deviatoric stress vector (4 components).
    PylithScalar devStressTpdt[4] = { 0.0, 0.0, 0.0, 0.0 };
    deviatoricStress4_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                               t, x, numConstants, constants, devStressTpdt);

    stress[0] += devStressTpdt[0]; // stress_xx
    stress[1] += devStressTpdt[3]; // stress_xy
    stress[2] += devStressTpdt[3]; // stress_yx
    stress[3] += devStressTpdt[1]; // stress_yy
} // deviatoricStress_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress including stress_zz for 2-D plane strain isotropic power-law
 * viscoelastic WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = 2*shearModulus*visStrain_ij
 */
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::deviatoricStress4(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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
                                                                   PylithScalar devStressTpdt[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_powerLawReferenceStrainRate = 1;
    const PylithInt i_powerLawReferenceStress = 2;
    const PylithInt i_powerLawExponent = 3;
    const PylithInt i_viscousStrain = 4;
    const PylithInt i_stress = 5;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 6);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(constants);

    // Constants.
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar powerLawReferenceStrainRate = a[aOff[i_powerLawReferenceStrainRate]];
    const PylithScalar powerLawReferenceStress = a[aOff[i_powerLawReferenceStress]];
    const PylithScalar powerLawExponent = a[aOff[i_powerLawExponent]];
    // ****** I think powerLawAlpha should actually be a kernel constant.
    const PylithScalar powerLawAlpha = 0.5;
    // const PylithScalar powerLawAlpha = a[aOff[i_powerLawAlpha]];
    const PylithScalar dt = constants[0];
    const PylithScalar ae = 1.0/(2.0*shearModulus);
    const PylithScalar timeFac = dt*(1.0 - powerLawAlpha);

    // Compute quantities based on stress at t = T.
    const PylithScalar* stressT = &a[aOff[i_stress]]; // stress_xx, stress_yy, stress_zz, stress_xy at t = T.
    const PylithScalar* visStrainT = &a[aOff[i_viscousStrain]]; // visStrain_xx, visStrain_yy, visStrain_zz,
                                                                // visStrain_xy at t = T.
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[4] = {
        stressT[0] - meanStressT,
        stressT[1] - meanStressT,
        stressT[2] - meanStressT,
        stressT[3],
    };
    const PylithScalar devStressProdT = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(devStressT, devStressT);
    const PylithScalar j2T = sqrt(0.5*devStressProdT);

    // Compute quantities based on strain at t = T + dt.
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar strainTpdt[4] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        0.0,
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]),
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0;
    const PylithScalar devStrainTpdt[4] = {
        strainTpdt[0] - meanStrainTpdt,
        strainTpdt[1] - meanStrainTpdt,
        -meanStrainTpdt,
        strainTpdt[2],
    };
    const PylithScalar strainPPTpdt[4] = {
        devStrainTpdt[0] - visStrainT[0],
        devStrainTpdt[1] - visStrainT[1],
        devStrainTpdt[2] - visStrainT[2],
        devStrainTpdt[3] - visStrainT[3],
    };
    const PylithScalar strainPPInvar2Tpdt = 0.5*pylith::fekernels::Viscoelasticity::scalarProduct2DPS(strainPPTpdt, strainPPTpdt);
    const PylithScalar strainStressInvar2T = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(strainPPTpdt, devStressT);

    // Finish defining parameters needed for root-finding algorithm.
    const PylithScalar b = strainPPInvar2Tpdt;
    const PylithScalar c = strainStressInvar2T*timeFac;
    const PylithScalar d = timeFac*j2T;
    PylithScalar j2Tpdt = 0.0;
    if ((b != 0.0) || (c != 0.0) || (d != 0.0)) {
        const PylithScalar j2InitialGuess = j2T;
        const PylithScalar stressScale = shearModulus;
        j2Tpdt = IsotropicPowerLawEffectiveStress::computeEffectiveStress(j2InitialGuess, stressScale, ae, b, c, d, powerLawAlpha,
                                                                          dt, j2T, powerLawExponent, powerLawReferenceStrainRate,
                                                                          powerLawReferenceStress);
    } // if
    // Compute deviatoric stresses from effective stress.
    const PylithScalar j2Tau = (1.0 - powerLawAlpha) * j2T + powerLawAlpha * j2Tpdt;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress),
                                                                  (powerLawExponent - 1.0))/powerLawReferenceStress;
    const PylithScalar factor1 = 1.0/(ae + powerLawAlpha*dt*gammaTau);
    const PylithScalar factor2 = timeFac*gammaTau;
    devStressTpdt[0] += factor1*(strainPPTpdt[0] - factor2*devStressT[0]);
    devStressTpdt[1] += factor1*(strainPPTpdt[1] - factor2*devStressT[1]);
    devStressTpdt[2] += factor1*(strainPPTpdt[2] - factor2*devStressT[2]);
    devStressTpdt[3] += factor1*(strainPPTpdt[3] - factor2*devStressT[3]);

} // deviatoricStress4


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress including stress_zz for 2-D plane strain isotropic power-law
 * viscoelastic WITH reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = 2*shearModulus*visStrain_ij
 */
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::deviatoricStress4_refstate(const PylithInt dim,
                                                                            const PylithInt numS,
                                                                            const PylithInt numA,
                                                                            const PylithInt sOff[],
                                                                            const PylithInt sOff_x[],
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
                                                                            PylithScalar devStressTpdt[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = 0;
    const PylithInt i_rstrain = 1;
    const PylithInt i_shearModulus = 2;
    const PylithInt i_powerLawReferenceStrainRate = 3;
    const PylithInt i_powerLawReferenceStress = 4;
    const PylithInt i_powerLawExponent = 5;
    const PylithInt i_viscousStrain = 6;
    const PylithInt i_stress = 7;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 8);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(constants);

    // Constants.
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar powerLawReferenceStrainRate = a[aOff[i_powerLawReferenceStrainRate]];
    const PylithScalar powerLawReferenceStress = a[aOff[i_powerLawReferenceStress]];
    const PylithScalar powerLawExponent = a[aOff[i_powerLawExponent]];
    // ****** I think powerLawAlpha should actually be a kernel constant.
    const PylithScalar powerLawAlpha = 0.5;
    // const PylithScalar powerLawAlpha = a[aOff[i_powerLawAlpha]];
    const PylithScalar dt = constants[0];
    const PylithScalar ae = 1.0/(2.0*shearModulus);
    const PylithScalar timeFac = dt*(1.0 - powerLawAlpha);

    // Compute quantities based on reference stress.
    const PylithScalar* stressR = &a[aOff[i_rstress]];
    const PylithScalar meanStressR = (stressR[0] + stressR[1] + stressR[2])/3.0;
    const PylithScalar devStressR[4] = {stressR[0] - meanStressR,
                                        stressR[1] - meanStressR,
                                        stressR[2] - meanStressR,
                                        stressR[3]};
    const PylithScalar j2RSquared = 0.5*pylith::fekernels::Viscoelasticity::scalarProduct2DPS(devStressR, devStressR);

    // Compute quantities based on reference strain.
    const PylithScalar* strainR = &a[aOff[i_rstrain]];
    const PylithScalar meanStrainR = (strainR[0] + strainR[1])/3.0;

    // Compute quantities based on stress at t = T.
    const PylithScalar* stressT = &a[aOff[i_stress]]; // stress_xx, stress_yy, stress_zz, stress_xy at t = T.
    const PylithScalar* visStrainT = &a[aOff[i_viscousStrain]]; // visStrain_xx, visStrain_yy, visStrain_zz,
                                                                // visStrain_xy at t = T.
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[4] = {stressT[0] - meanStressT,
                                        stressT[1] - meanStressT,
                                        stressT[2] - meanStressT,
                                        stressT[3]};
    const PylithScalar devStressProdT = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(devStressT, devStressT);
    const PylithScalar j2T = sqrt(0.5*devStressProdT);

    // Compute quantities based on strain at t = T + dt.
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar strainTpdt[4] = {disp_x[0*_dim+0],
                                        disp_x[1*_dim+1],
                                        0.0,
                                        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0])};
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0 - meanStrainR;
    const PylithScalar strainPPTpdt[4] = {strainTpdt[0] - meanStrainTpdt - visStrainT[0] - strainR[0],
                                          strainTpdt[1] - meanStrainTpdt - visStrainT[1] - strainR[1],
                                          -meanStrainTpdt - visStrainT[2],
                                          strainTpdt[3] - visStrainT[3] - strainR[2]};
    const PylithScalar strainPPInvar2Tpdt = 0.5*pylith::fekernels::Viscoelasticity::scalarProduct2DPS(strainPPTpdt, strainPPTpdt);
    const PylithScalar strainStressInvar2T = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(strainPPTpdt, devStressT);
    const PylithScalar strainStressInvar2R = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(strainPPTpdt, devStressR);
    const PylithScalar stressInvar2RT = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(devStressT, devStressR);

    // Finish defining parameters needed for root-finding algorithm.
    const PylithScalar b = strainPPInvar2Tpdt + ae*strainStressInvar2R + ae*ae*j2RSquared;
    const PylithScalar c = strainStressInvar2T*timeFac + ae*stressInvar2RT;
    const PylithScalar d = timeFac*j2T;
    PylithScalar j2Tpdt = 0.0;
    if ((b != 0.0) || (c != 0.0) || (d != 0.0)) {
        const PylithScalar j2InitialGuess = j2T;
        const PylithScalar stressScale = shearModulus;
        j2Tpdt = IsotropicPowerLawEffectiveStress::computeEffectiveStress(j2InitialGuess, stressScale, ae, b, c, d, powerLawAlpha,
                                                                          dt, j2T, powerLawExponent, powerLawReferenceStrainRate,
                                                                          powerLawReferenceStress);
    } // if
    // Compute deviatoric stresses from effective stress.
    const PylithScalar j2Tau = (1.0 - powerLawAlpha) * j2T + powerLawAlpha * j2Tpdt;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress),
                                                                  (powerLawExponent - 1.0))/powerLawReferenceStress;
    const PylithScalar factor1 = 1.0/(ae + powerLawAlpha*dt*gammaTau);
    const PylithScalar factor2 = timeFac*gammaTau;
    devStressTpdt[0] += factor1*(strainPPTpdt[0] - factor2*devStressT[0] + ae*devStressR[0]);
    devStressTpdt[1] += factor1*(strainPPTpdt[1] - factor2*devStressT[1] + ae*devStressR[1]);
    devStressTpdt[2] += factor1*(strainPPTpdt[2] - factor2*devStressT[2] + ae*devStressR[2]);
    devStressTpdt[3] += factor1*(strainPPTpdt[3] - factor2*devStressT[3] + ae*devStressR[3]);

} // deviatoricStress4_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Update stress for a plane strain power-law viscoelastic material.
 *
 * IMPORTANT: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 *
 */
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::updateStress(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
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
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

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
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 6; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[6] = {aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate], aOff[i_powerLawReferenceStress],
                                  aOff[i_powerLawExponent], aOff[i_viscousStrain], aOff[i_stress]};

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
    IsotropicLinearElasticityPlaneStrain::meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL,
                                                     a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    stress[0] += stressTensor[0];
    stress[1] += stressTensor[0];
    stress[2] += stressTensor[0];

    // Compute deviatoric stress (4 components).
    deviatoricStress4(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                      t, x, numConstants, constants, stress);

} // updateStress


// ---------------------------------------------------------------------------------------------------------------------
/* Update stress for a plane strain power-law viscoelastic material WITH reference stress and strain.
 *
 * IMPORTANT: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 *
 */
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::updateStress_refstate(const PylithInt dim,
                                                                       const PylithInt numS,
                                                                       const PylithInt numA,
                                                                       const PylithInt sOff[],
                                                                       const PylithInt sOff_x[],
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
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-9;
    const PylithInt i_rstrain = numA-8;
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 8; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[8] = {aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate],
                                  aOff[i_powerLawReferenceStress], aOff[i_powerLawExponent], aOff[i_viscousStrain],
                                  aOff[i_stress]};

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
    IsotropicLinearElasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean,
                                                              NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    stress[0] += stressTensor[0];
    stress[1] += stressTensor[0];
    stress[2] += stressTensor[0];

    // Compute deviatoric stress vector (4 components).
    deviatoricStress4_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                               t, x, numConstants, constants, stress);

} // updateStress_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Update viscous strain for a power-law viscoelastic material.
 *
 * :ATTENTION: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 */
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::updateViscousStrain(const PylithInt dim,
                                                                     const PylithInt numS,
                                                                     const PylithInt numA,
                                                                     const PylithInt sOff[],
                                                                     const PylithInt sOff_x[],
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
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

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
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(constants);

    // Constants.
    const PylithScalar powerLawReferenceStrainRate = a[aOff[i_powerLawReferenceStrainRate]];
    const PylithScalar powerLawReferenceStress = a[aOff[i_powerLawReferenceStress]];
    const PylithScalar powerLawExponent = a[aOff[i_powerLawExponent]];
    // ****** I think powerLawAlpha should actually be a kernel constant.
    const PylithScalar powerLawAlpha = 0.5;
    // const PylithScalar powerLawAlpha = a[aOff[i_powerLawAlpha]];
    const PylithScalar dt = constants[0];

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 6; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[6] = {aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate], aOff[i_powerLawReferenceStress],
                                  aOff[i_powerLawExponent], aOff[i_viscousStrain], aOff[i_stress]};

    // Compute current deviatoric stress.
    PylithScalar devStressTpdt[4] = { 0.0, 0.0, 0.0, 0.0 };
    deviatoricStress4(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                      t, x, numConstants, constants, devStressTpdt);
    const PylithScalar devStressProdTpdt = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(devStressTpdt, devStressTpdt);
    const PylithScalar j2Tpdt = sqrt(0.5*devStressProdTpdt);

    // Compute stress quantities at time T.
    const PylithScalar* stressT = &a[aOff[i_stress]]; // stress_xx, stress_yy, stress_zz, stress_xy at t = T.
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[4] = {stressT[0] - meanStressT,
                                        stressT[1] - meanStressT,
                                        stressT[2] - meanStressT,
                                        stressT[3]};
    const PylithScalar devStressProdT = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(devStressT, devStressT);
    const PylithScalar j2T = sqrt(0.5*devStressProdT);

    // Compute quantities at intermediate time.
    const PylithScalar devStressTau[4] = {(1.0 - powerLawAlpha)*devStressT[0] + powerLawAlpha*devStressTpdt[0],
                                          (1.0 - powerLawAlpha)*devStressT[1] + powerLawAlpha*devStressTpdt[1],
                                          (1.0 - powerLawAlpha)*devStressT[2] + powerLawAlpha*devStressTpdt[2],
                                          (1.0 - powerLawAlpha)*devStressT[3] + powerLawAlpha*devStressTpdt[3]};
    const PylithScalar j2Tau = (1.0 - powerLawAlpha)*j2T + powerLawAlpha*j2Tpdt;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress),
                                                                  (powerLawExponent - 1.0))/powerLawReferenceStress;

    // Update viscous strain.
    visStrain[0] += dt*gammaTau*devStressTau[0];
    visStrain[1] += dt*gammaTau*devStressTau[1];
    visStrain[2] += dt*gammaTau*devStressTau[2];
    visStrain[3] += dt*gammaTau*devStressTau[3];

} // updateViscousStrain


// ---------------------------------------------------------------------------------------------------------------------
/* Update viscous strain for a power-law viscoelastic material WITH reference stress and strain.
 *
 * :ATTENTION: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 */
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::updateViscousStrain_refstate(const PylithInt dim,
                                                                              const PylithInt numS,
                                                                              const PylithInt numA,
                                                                              const PylithInt sOff[],
                                                                              const PylithInt sOff_x[],
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
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-9;
    const PylithInt i_rstrain = numA-8;
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(constants);

    // Constants.
    const PylithScalar powerLawReferenceStrainRate = a[aOff[i_powerLawReferenceStrainRate]];
    const PylithScalar powerLawReferenceStress = a[aOff[i_powerLawReferenceStress]];
    const PylithScalar powerLawExponent = a[aOff[i_powerLawExponent]];
    // ****** I think powerLawAlpha should actually be a kernel constant.
    const PylithScalar powerLawAlpha = 0.5;
    // const PylithScalar powerLawAlpha = a[aOff[i_powerLawAlpha]];
    const PylithScalar dt = constants[0];

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 8; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[8] = {aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate],
                                  aOff[i_powerLawReferenceStress], aOff[i_powerLawExponent], aOff[i_viscousStrain],
                                  aOff[i_stress]};

    // Compute current deviatoric stress.
    PylithScalar devStressTpdt[4] = { 0.0, 0.0, 0.0, 0.0 };
    deviatoricStress4_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                               t, x, numConstants, constants, devStressTpdt);
    const PylithScalar devStressProdTpdt = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(devStressTpdt, devStressTpdt);
    const PylithScalar j2Tpdt = sqrt(0.5*devStressProdTpdt);

    // Compute stress quantities at time T.
    const PylithScalar* stressT = &a[aOff[i_stress]]; // stress_xx, stress_yy, stress_zz, stress_xy at t = T.
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[4] = {stressT[0] - meanStressT,
                                        stressT[1] - meanStressT,
                                        stressT[2] - meanStressT,
                                        stressT[3]};
    const PylithScalar devStressProdT = pylith::fekernels::Viscoelasticity::scalarProduct2DPS(devStressT, devStressT);
    const PylithScalar j2T = sqrt(0.5*devStressProdT);

    // Compute quantities at intermediate time.
    const PylithScalar devStressTau[4] = {(1.0 - powerLawAlpha)*devStressT[0] + powerLawAlpha*devStressTpdt[0],
                                          (1.0 - powerLawAlpha)*devStressT[1] + powerLawAlpha*devStressTpdt[1],
                                          (1.0 - powerLawAlpha)*devStressT[2] + powerLawAlpha*devStressTpdt[2],
                                          (1.0 - powerLawAlpha)*devStressT[3] + powerLawAlpha*devStressTpdt[3]};
    const PylithScalar j2Tau = (1.0 - powerLawAlpha)*j2T + powerLawAlpha*j2Tpdt;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress),
                                                                  (powerLawExponent - 1.0))/powerLawReferenceStress;

    // Update viscous strain.
    visStrain[0] += dt*gammaTau*devStressTau[0];
    visStrain[1] += dt*gammaTau*devStressTau[1];
    visStrain[2] += dt*gammaTau*devStressTau[2];
    visStrain[3] += dt*gammaTau*devStressTau[3];

} // updateViscousStrain_refstate


// ---------------------------------------------------------------------------------------------------------------------
// Calculate stress for 2-D plane strain isotropic power-law viscoelastic material WITHOUT a reference stress and
// strain.
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::cauchyStress(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
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
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

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
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 6; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[6] = {aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate], aOff[i_powerLawReferenceStress],
                                  aOff[i_powerLawExponent], aOff[i_viscousStrain], aOff[i_stress]};

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
    IsotropicLinearElasticityPlaneStrain::meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL,
                                                     a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    stressVector[0] = stressTensor[0];
    stressVector[1] = stressTensor[0];
    stressVector[2] = stressTensor[0];

    // Compute deviatoric stress vector (4 components).
    deviatoricStress4(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                      t, x, numConstants, constants, stressVector);

} // stress


// ---------------------------------------------------------------------------------------------------------------------
// Calculate stress for 2-D plane strain isotropic linear power-law viscoelastic material WITH a reference
// stress/strain.
void
pylith::fekernels::IsotropicPowerLawPlaneStrain::cauchyStress_refstate(const PylithInt dim,
                                                                       const PylithInt numS,
                                                                       const PylithInt numA,
                                                                       const PylithInt sOff[],
                                                                       const PylithInt sOff_x[],
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
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-9;
    const PylithInt i_rstrain = numA-8;
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Pass bulk modulus, reference stress, and reference strain.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 8; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[8] = {aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate],
                                  aOff[i_powerLawReferenceStress], aOff[i_powerLawExponent], aOff[i_viscousStrain],
                                  aOff[i_stress]};

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
    IsotropicLinearElasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean,
                                                              NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    stressVector[0] = stressTensor[0];
    stressVector[1] = stressTensor[0];
    stressVector[2] = stressTensor[0];

    // Compute deviatoric stress vector (4 components).
    deviatoricStress4_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                               t, x, numConstants, constants, stressVector);

} // stress_refstate


// =====================================================================================================================
// Kernels for isotropic, power-law viscoelastic 3D material.
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
// f1 function for isotropic power-law 3D viscoelastic material WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicPowerLaw3D::f1v(const PylithInt dim,
                                            const PylithInt numS,
                                            const PylithInt numA,
                                            const PylithInt sOff[],
                                            const PylithInt sOff_x[],
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
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

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
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 6; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[6] = {aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate], aOff[i_powerLawReferenceStress],
                                  aOff[i_powerLawExponent], aOff[i_viscousStrain], aOff[i_stress]};

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    IsotropicLinearElasticity3D::meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                            aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v


// ---------------------------------------------------------------------------------------------------------------------
// f1 function for isotropic power-law viscoelastic 3D material with reference
// stress and strain.
void
pylith::fekernels::IsotropicPowerLaw3D::f1v_refstate(const PylithInt dim,
                                                     const PylithInt numS,
                                                     const PylithInt numA,
                                                     const PylithInt sOff[],
                                                     const PylithInt sOff_x[],
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
    const PylithInt i_rstress = numA-9;
    const PylithInt i_rstrain = numA-8;
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 8; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[8] = {aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate],
                                  aOff[i_powerLawReferenceStress], aOff[i_powerLawExponent], aOff[i_viscousStrain],
                                  aOff[i_stress]};

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    IsotropicLinearElasticity3D::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                     aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1v_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Jf3_vu entry function for 3-D isotropic power-law viscoelastic material WITHOUT reference stress/strain.
 *
 * stress_ij = C_ijkl strain_kl
 *
 * For reference:
 *
 * Isotropic:
 *  C_ijkl = bulkModulus * delta_ij * delta_kl
 *   + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
 */
void
pylith::fekernels::IsotropicPowerLaw3D::Jf3vu(const PylithInt dim,
                                              const PylithInt numS,
                                              const PylithInt numA,
                                              const PylithInt sOff[],
                                              const PylithInt sOff_x[],
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

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
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(constants);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 6; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[6] = {aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate], aOff[i_powerLawReferenceStress],
                                  aOff[i_powerLawExponent], aOff[i_viscousStrain], aOff[i_stress]};

    // Compute deviatoric stress tensor (9 components).
    PylithScalar devStressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, devStressTensor);
    const PylithScalar devStressTpdt[6] = {
        devStressTensor[0],
        devStressTensor[4],
        devStressTensor[8],
        devStressTensor[1],
        devStressTensor[5],
        devStressTensor[2],
    };

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar powerLawReferenceStrainRate = a[aOff[i_powerLawReferenceStrainRate]];
    const PylithScalar powerLawReferenceStress = a[aOff[i_powerLawReferenceStress]];
    const PylithScalar powerLawExponent = a[aOff[i_powerLawExponent]];
    // ****** I think powerLawAlpha should actually be a kernel constant.
    const PylithScalar powerLawAlpha = 0.5;
    // const PylithScalar powerLawAlpha = a[aOff[i_powerLawAlpha]];
    const PylithScalar dt = constants[0];

    // Compute quantities based on stress at t = T.
    const PylithScalar* stressT = &a[aOff[i_stress]]; // stress_xx, stress_yy, stress_zz, stress_xy, stress_yz,
                                                      // stress_xz at t = T.
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[6] = {
        stressT[0] - meanStressT,
        stressT[1] - meanStressT,
        stressT[2] - meanStressT,
        stressT[3],
        stressT[4],
        stressT[5],
    };
    const PylithScalar devStressProdT = pylith::fekernels::Viscoelasticity::scalarProduct3D(devStressT, devStressT);
    const PylithScalar j2T = sqrt(0.5*devStressProdT);

    // Compute quantities based on stress at t = T + dt.
    const PylithScalar devStressProdTpdt = pylith::fekernels::Viscoelasticity::scalarProduct3D(devStressTpdt, devStressTpdt);
    const PylithScalar j2Tpdt = sqrt(0.5*devStressProdTpdt);

    // Compute quantities at intermediate time tau.
    const PylithScalar j2Tau = powerLawAlpha*j2Tpdt + (1.0 - powerLawAlpha)*j2T;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress),
                                                                  (powerLawExponent - 1.0))/powerLawReferenceStress;

    /* Unique components of Jacobian. */
    const PylithReal C1111 = bulkModulus -
                             2.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[0]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressT[0]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C1122 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[1]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressT[1]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C1133 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[2]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressT[2]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C1212 =
        -1.0/(2.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[3]*devStressTpdt[3]*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) +
                   powerLawAlpha*dt*gammaTau +
                   powerLawAlpha*dt*gammaTau*devStressTpdt[3]*devStressT[3]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C1313 =
        -1.0/(2.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[5]*devStressTpdt[5]*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) +
                   powerLawAlpha*dt*gammaTau +
                   powerLawAlpha*dt*gammaTau*devStressTpdt[5]*devStressT[5]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C2211 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[1]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressT[0]*devStressTpdt[1]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C2222 = bulkModulus -
                             2.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressTpdt[1]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressT[1]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C2233 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressTpdt[2]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressT[2]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C2323 =
        -1.0/(2.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[4]*devStressTpdt[4]*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) +
                   powerLawAlpha*dt*gammaTau +
                   powerLawAlpha*dt*gammaTau*devStressTpdt[4]*devStressT[4]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C3311 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[2]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressT[0]*devStressTpdt[2]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C3322 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressTpdt[2]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressT[1]*devStressTpdt[2]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C3333 = bulkModulus -
                             2.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[2]*devStressTpdt[2]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[2]*devStressT[2]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    /* j(f,g,df,dg) = C(f,df,g,dg)
     *
     * 0:  j0000 = C1111 = bulkModulus - 2/(3*(alpha**2*deltaT*gammaFTau*s11**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau + alpha*deltaT*gammaFTau*s11*s11T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 1:  j0001 = C1112 = 0
     * 2:  j0002 = C1113 = 0
     * 3:  j0010 = C1211 = 0
     * 4:  j0011 = C1212 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 5:  j0012 = C1213 = 0
     * 6:  j0020 = C1311 = 0
     * 7:  j0021 = C1312 = 0
     * 8:  j0022 = C1313 = -1/(2*(alpha**2*deltaT*gammaFTau*s13**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s13*s13T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 9:  j0100 = C1121 = 0
     * 10:  j0101 = C1122 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s11*s22*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s11*s22T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 11:  j0102 = C1123 = 0
     * 12:  j0110 = C1221 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 13:  j0111 = C1222 = 0
     * 14:  j0112 = C1223 = 0
     * 15:  j0120 = C1321 = 0
     * 16:  j0121 = C1322 = 0
     * 17:  j0122 = C1323 = 0
     * 18:  j0200 = C1131 = 0
     * 19:  j0201 = C1132 = 0
     * 20:  j0202 = C1133 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s11*s33*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s11*s33T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 21:  j0210 = C1231 = 0
     * 22:  j0211 = C1232 = 0
     * 23:  j0212 = C1233 = 0
     * 24:  j0220 = C1331 = -1/(2*(alpha**2*deltaT*gammaFTau*s13**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s13*s13T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 25:  j0221 = C1332 = 0
     * 26:  j0222 = C1333 = 0
     * 27:  j1000 = C2111 = 0
     * 28:  j1001 = C2112 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 29:  j1002 = C2113 = 0
     * 30:  j1010 = C2211 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s11*s22*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s11T*s22*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 31:  j1011 = C2212 = 0
     * 32:  j1012 = C2213 = 0
     * 33:  j1020 = C2311 = 0
     * 34:  j1021 = C2312 = 0
     * 35:  j1022 = C2313 = 0
     * 36:  j1100 = C2121 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 37:  j1101 = C2122 = 0
     * 38:  j1102 = C2123 = 0
     * 39:  j1110 = C2221 = 0
     * 40:  j1111 = C2222 = bulkModulus - 2/(3*(alpha**2*deltaT*gammaFTau*s22**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau + alpha*deltaT*gammaFTau*s22*s22T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 41:  j1112 = C2223 = 0
     * 42:  j1120 = C2321 = 0
     * 43:  j1121 = C2322 = 0
     * 44:  j1122 = C2323 = -1/(2*(alpha**2*deltaT*gammaFTau*s23**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s23*s23T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 45:  j1200 = C2131 = 0
     * 46:  j1201 = C2132 = 0
     * 47:  j1202 = C2133 = 0
     * 48:  j1210 = C2231 = 0
     * 49:  j1211 = C2232 = 0
     * 50:  j1212 = C2233 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s22*s33*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s22*s33T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 51:  j1220 = C2331 = 0
     * 52:  j1221 = C2332 = -1/(2*(alpha**2*deltaT*gammaFTau*s23**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s23*s23T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 53:  j1222 = C2333 = 0
     * 54:  j2000 = C3111 = 0
     * 55:  j2001 = C3112 = 0
     * 56:  j2002 = C3113 = -1/(2*(alpha**2*deltaT*gammaFTau*s13**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s13*s13T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 57:  j2010 = C3211 = 0
     * 58:  j2011 = C3212 = 0
     * 59:  j2012 = C3213 = 0
     * 60:  j2020 = C3311 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s11*s33*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s11T*s33*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 61:  j2021 = C3312 = 0
     * 62:  j2022 = C3313 = 0
     * 63:  j2100 = C3121 = 0
     * 64:  j2101 = C3122 = 0
     * 65:  j2102 = C3123 = 0
     * 66:  j2110 = C3221 = 0
     * 67:  j2111 = C3222 = 0
     * 68:  j2112 = C3223 = -1/(2*(alpha**2*deltaT*gammaFTau*s23**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s23*s23T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 69:  j2120 = C3321 = 0
     * 70:  j2121 = C3322 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s22*s33*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s22T*s33*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 71:  j2122 = C3323 = 0
     * 72:  j2200 = C3131 = -1/(2*(alpha**2*deltaT*gammaFTau*s13**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s13*s13T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 73:  j2201 = C3132 = 0
     * 74:  j2202 = C3133 = 0
     * 75:  j2210 = C3231 = 0
     * 76:  j2211 = C3232 = -1/(2*(alpha**2*deltaT*gammaFTau*s23**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s23*s23T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 77:  j2212 = C3233 = 0
     * 78:  j2220 = C3331 = 0
     * 79:  j2221 = C3332 = 0
     * 80:  j2222 = C3333 = bulkModulus - 2/(3*(alpha**2*deltaT*gammaFTau*s33**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau + alpha*deltaT*gammaFTau*s33*s33T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt) + 1/(2*mu)))
     */

    /* Nonzero Jacobian entries. */
    Jf3[0] -= C1111; /* j0000 */
    Jf3[4] -= C1212; /* j0011 */
    Jf3[8] -= C1313; /* j0022 */
    Jf3[10] -= C1122; /* j0101 */
    Jf3[12] -= C1212; /* j0110 */
    Jf3[20] -= C1133; /* j0202 */
    Jf3[24] -= C1313; /* j0220 */
    Jf3[28] -= C1212; /* j1001 */
    Jf3[30] -= C2211; /* j1010 */
    Jf3[36] -= C1212; /* j1100 */
    Jf3[40] -= C2222; /* j1111 */
    Jf3[44] -= C2323; /* j1122 */
    Jf3[50] -= C2233; /* j1212 */
    Jf3[52] -= C2323; /* j1221 */
    Jf3[56] -= C1313; /* j2002 */
    Jf3[60] -= C3311; /* j2020 */
    Jf3[68] -= C2323; /* j2112 */
    Jf3[70] -= C3322; /* j2121 */
    Jf3[72] -= C1313; /* j2200 */
    Jf3[76] -= C2323; /* j2211 */
    Jf3[80] -= C3333; /* j2222 */

} // Jf3vu


// ---------------------------------------------------------------------------------------------------------------------
/* Jf3_vu entry function for 3-D isotropic power-law viscoelastic material WITH reference stress/strain.
 *
 * stress_ij = C_ijkl strain_kl
 *
 * For reference:
 *
 * Isotropic:
 *  C_ijkl = bulkModulus * delta_ij * delta_kl
 *   + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
 */
void
pylith::fekernels::IsotropicPowerLaw3D::Jf3vu_refstate(const PylithInt dim,
                                                       const PylithInt numS,
                                                       const PylithInt numA,
                                                       const PylithInt sOff[],
                                                       const PylithInt sOff_x[],
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-9;
    const PylithInt i_rstrain = numA-8;
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(constants);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 8; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[8] = {aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate],
                                  aOff[i_powerLawReferenceStress], aOff[i_powerLawExponent], aOff[i_viscousStrain],
                                  aOff[i_stress]};

    // Compute deviatoric stress tensor (9 components).
    PylithScalar devStressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, devStressTensor);
    const PylithScalar devStressTpdt[6] = {devStressTensor[0],
                                           devStressTensor[4],
                                           devStressTensor[8],
                                           devStressTensor[1],
                                           devStressTensor[5],
                                           devStressTensor[2]};

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar powerLawReferenceStrainRate = a[aOff[i_powerLawReferenceStrainRate]];
    const PylithScalar powerLawReferenceStress = a[aOff[i_powerLawReferenceStress]];
    const PylithScalar powerLawExponent = a[aOff[i_powerLawExponent]];
    // ****** I think powerLawAlpha should actually be a kernel constant.
    const PylithScalar powerLawAlpha = 0.5;
    // const PylithScalar powerLawAlpha = a[aOff[i_powerLawAlpha]];
    const PylithScalar dt = constants[0];

    // Compute quantities based on stress at t = T.
    const PylithScalar* stressT = &a[aOff[i_stress]]; // stress_xx, stress_yy, stress_zz, stress_xy, stress_yz,
                                                      // stress_xz at t = T.
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[6] = {stressT[0] - meanStressT,
                                        stressT[1] - meanStressT,
                                        stressT[2] - meanStressT,
                                        stressT[3],
                                        stressT[4],
                                        stressT[5]};
    const PylithScalar devStressProdT = pylith::fekernels::Viscoelasticity::scalarProduct3D(devStressT, devStressT);
    const PylithScalar j2T = sqrt(0.5*devStressProdT);

    // Compute quantities based on stress at t = T + dt.
    const PylithScalar devStressProdTpdt = pylith::fekernels::Viscoelasticity::scalarProduct3D(devStressTpdt, devStressTpdt);
    const PylithScalar j2Tpdt = sqrt(0.5*devStressProdTpdt);

    // Compute quantities at intermediate time tau.
    const PylithScalar j2Tau = powerLawAlpha*j2Tpdt + (1.0 - powerLawAlpha)*j2T;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress),
                                                                  (powerLawExponent - 1.0))/powerLawReferenceStress;

    /* Unique components of Jacobian. */
    const PylithReal C1111 = bulkModulus -
                             2.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[0]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressT[0]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C1122 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[1]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressT[1]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C1133 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[2]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressT[2]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C1212 =
        -1.0/(2.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[3]*devStressTpdt[3]*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) +
                   powerLawAlpha*dt*gammaTau +
                   powerLawAlpha*dt*gammaTau*devStressTpdt[3]*devStressT[3]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C1313 =
        -1.0/(2.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[5]*devStressTpdt[5]*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) +
                   powerLawAlpha*dt*gammaTau +
                   powerLawAlpha*dt*gammaTau*devStressTpdt[5]*devStressT[5]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C2211 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[1]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressT[0]*devStressTpdt[1]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C2222 = bulkModulus -
                             2.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressTpdt[1]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressT[1]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C2233 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressTpdt[2]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressT[2]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C2323 =
        -1.0/(2.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[4]*devStressTpdt[4]*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) +
                   powerLawAlpha*dt*gammaTau +
                   powerLawAlpha*dt*gammaTau*devStressTpdt[4]*devStressT[4]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                   (j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    const PylithReal C3311 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[0]*devStressTpdt[2]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressT[0]*devStressTpdt[2]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C3322 = bulkModulus +
                             1.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[1]*devStressTpdt[2]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau*devStressT[1]*devStressTpdt[2]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt)));
    const PylithReal C3333 = bulkModulus -
                             2.0/(3.0*(powerLawAlpha*powerLawAlpha*dt*gammaTau*devStressTpdt[2]*devStressTpdt[2]*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) +
                                       powerLawAlpha*dt*gammaTau +
                                       powerLawAlpha*dt*gammaTau*devStressTpdt[2]*devStressT[2]*(-powerLawAlpha + 1.0)*(powerLawExponent - 1.0)/
                                       (2.0*j2Tau*j2Tpdt) + 1.0/(2.0*shearModulus)));
    /* j(f,g,df,dg) = C(f,df,g,dg)
     *
     * 0:  j0000 = C1111 = bulkModulus - 2/(3*(alpha**2*deltaT*gammaFTau*s11**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau + alpha*deltaT*gammaFTau*s11*s11T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 1:  j0001 = C1112 = 0
     * 2:  j0002 = C1113 = 0
     * 3:  j0010 = C1211 = 0
     * 4:  j0011 = C1212 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 5:  j0012 = C1213 = 0
     * 6:  j0020 = C1311 = 0
     * 7:  j0021 = C1312 = 0
     * 8:  j0022 = C1313 = -1/(2*(alpha**2*deltaT*gammaFTau*s13**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s13*s13T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 9:  j0100 = C1121 = 0
     * 10:  j0101 = C1122 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s11*s22*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s11*s22T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 11:  j0102 = C1123 = 0
     * 12:  j0110 = C1221 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 13:  j0111 = C1222 = 0
     * 14:  j0112 = C1223 = 0
     * 15:  j0120 = C1321 = 0
     * 16:  j0121 = C1322 = 0
     * 17:  j0122 = C1323 = 0
     * 18:  j0200 = C1131 = 0
     * 19:  j0201 = C1132 = 0
     * 20:  j0202 = C1133 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s11*s33*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s11*s33T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 21:  j0210 = C1231 = 0
     * 22:  j0211 = C1232 = 0
     * 23:  j0212 = C1233 = 0
     * 24:  j0220 = C1331 = -1/(2*(alpha**2*deltaT*gammaFTau*s13**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s13*s13T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 25:  j0221 = C1332 = 0
     * 26:  j0222 = C1333 = 0
     * 27:  j1000 = C2111 = 0
     * 28:  j1001 = C2112 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 29:  j1002 = C2113 = 0
     * 30:  j1010 = C2211 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s11*s22*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s11T*s22*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 31:  j1011 = C2212 = 0
     * 32:  j1012 = C2213 = 0
     * 33:  j1020 = C2311 = 0
     * 34:  j1021 = C2312 = 0
     * 35:  j1022 = C2313 = 0
     * 36:  j1100 = C2121 = -1/(2*(alpha**2*deltaT*gammaFTau*s12**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s12*s12T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 37:  j1101 = C2122 = 0
     * 38:  j1102 = C2123 = 0
     * 39:  j1110 = C2221 = 0
     * 40:  j1111 = C2222 = bulkModulus - 2/(3*(alpha**2*deltaT*gammaFTau*s22**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau + alpha*deltaT*gammaFTau*s22*s22T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 41:  j1112 = C2223 = 0
     * 42:  j1120 = C2321 = 0
     * 43:  j1121 = C2322 = 0
     * 44:  j1122 = C2323 = -1/(2*(alpha**2*deltaT*gammaFTau*s23**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s23*s23T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 45:  j1200 = C2131 = 0
     * 46:  j1201 = C2132 = 0
     * 47:  j1202 = C2133 = 0
     * 48:  j1210 = C2231 = 0
     * 49:  j1211 = C2232 = 0
     * 50:  j1212 = C2233 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s22*s33*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s22*s33T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 51:  j1220 = C2331 = 0
     * 52:  j1221 = C2332 = -1/(2*(alpha**2*deltaT*gammaFTau*s23**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s23*s23T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 53:  j1222 = C2333 = 0
     * 54:  j2000 = C3111 = 0
     * 55:  j2001 = C3112 = 0
     * 56:  j2002 = C3113 = -1/(2*(alpha**2*deltaT*gammaFTau*s13**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s13*s13T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 57:  j2010 = C3211 = 0
     * 58:  j2011 = C3212 = 0
     * 59:  j2012 = C3213 = 0
     * 60:  j2020 = C3311 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s11*s33*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s11T*s33*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 61:  j2021 = C3312 = 0
     * 62:  j2022 = C3313 = 0
     * 63:  j2100 = C3121 = 0
     * 64:  j2101 = C3122 = 0
     * 65:  j2102 = C3123 = 0
     * 66:  j2110 = C3221 = 0
     * 67:  j2111 = C3222 = 0
     * 68:  j2112 = C3223 = -1/(2*(alpha**2*deltaT*gammaFTau*s23**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s23*s23T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 69:  j2120 = C3321 = 0
     * 70:  j2121 = C3322 = bulkModulus + 1/(3*(alpha**2*deltaT*gammaFTau*s22*s33*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau*s22T*s33*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt)))
     * 71:  j2122 = C3323 = 0
     * 72:  j2200 = C3131 = -1/(2*(alpha**2*deltaT*gammaFTau*s13**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s13*s13T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 73:  j2201 = C3132 = 0
     * 74:  j2202 = C3133 = 0
     * 75:  j2210 = C3231 = 0
     * 76:  j2211 = C3232 = -1/(2*(alpha**2*deltaT*gammaFTau*s23**2*(n - 1)/(j2FTau*j2FTplusDt) + alpha*deltaT*gammaFTau
     * + alpha*deltaT*gammaFTau*s23*s23T*(-alpha + 1)*(n - 1)/(j2FTau*j2FTplusDt) + 1/(2*mu)))
     * 77:  j2212 = C3233 = 0
     * 78:  j2220 = C3331 = 0
     * 79:  j2221 = C3332 = 0
     * 80:  j2222 = C3333 = bulkModulus - 2/(3*(alpha**2*deltaT*gammaFTau*s33**2*(n - 1)/(2*j2FTau*j2FTplusDt) +
     * alpha*deltaT*gammaFTau + alpha*deltaT*gammaFTau*s33*s33T*(-alpha + 1)*(n - 1)/(2*j2FTau*j2FTplusDt) + 1/(2*mu)))
     */

    /* Nonzero Jacobian entries. */
    Jf3[0] -= C1111; /* j0000 */
    Jf3[4] -= C1212; /* j0011 */
    Jf3[8] -= C1313; /* j0022 */
    Jf3[10] -= C1122; /* j0101 */
    Jf3[12] -= C1212; /* j0110 */
    Jf3[20] -= C1133; /* j0202 */
    Jf3[24] -= C1313; /* j0220 */
    Jf3[28] -= C1212; /* j1001 */
    Jf3[30] -= C2211; /* j1010 */
    Jf3[36] -= C1212; /* j1100 */
    Jf3[40] -= C2222; /* j1111 */
    Jf3[44] -= C2323; /* j1122 */
    Jf3[50] -= C2233; /* j1212 */
    Jf3[52] -= C2323; /* j1221 */
    Jf3[56] -= C1313; /* j2002 */
    Jf3[60] -= C3311; /* j2020 */
    Jf3[68] -= C2323; /* j2112 */
    Jf3[70] -= C3322; /* j2121 */
    Jf3[72] -= C1313; /* j2200 */
    Jf3[76] -= C2323; /* j2211 */
    Jf3[80] -= C3333; /* j2222 */

} // Jf3vu_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress tensor for 3-D isotropic power-law
 * viscoelastic material WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = 2*shearModulus*visStrain_ij
 */
void
pylith::fekernels::IsotropicPowerLaw3D::deviatoricStress(const PylithInt dim,
                                                         const PylithInt numS,
                                                         const PylithInt numA,
                                                         const PylithInt sOff[],
                                                         const PylithInt sOff_x[],
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
                                                         PylithScalar devStressTpdt[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_powerLawReferenceStrainRate = 1;
    const PylithInt i_powerLawReferenceStress = 2;
    const PylithInt i_powerLawExponent = 3;
    const PylithInt i_viscousStrain = 4;
    const PylithInt i_stress = 5;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 6);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(constants);

    // Constants.
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar powerLawReferenceStrainRate = a[aOff[i_powerLawReferenceStrainRate]];
    const PylithScalar powerLawReferenceStress = a[aOff[i_powerLawReferenceStress]];
    const PylithScalar powerLawExponent = a[aOff[i_powerLawExponent]];
    // ****** I think powerLawAlpha should actually be a kernel constant.
    const PylithScalar powerLawAlpha = 0.5;
    // const PylithScalar powerLawAlpha = a[aOff[i_powerLawAlpha]];
    const PylithScalar dt = constants[0];
    const PylithScalar ae = 1.0/(2.0*shearModulus);
    const PylithScalar timeFac = dt*(1.0 - powerLawAlpha);

    // Compute quantities based on stress at t = T.
    const PylithScalar* stressT = &a[aOff[i_stress]]; // stress_xx, stress_yy, stress_zz, stress_xy, stress_yz,
                                                      // stress_xz at t = T.
    const PylithScalar* visStrainT = &a[aOff[i_viscousStrain]]; // visStrain_xx, visStrain_yy, visStrain_zz,
                                                                // visStrain_xy, visStrain_yz, visStrain_xz at t = T.
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[6] = {stressT[0] - meanStressT,
                                        stressT[1] - meanStressT,
                                        stressT[2] - meanStressT,
                                        stressT[3],
                                        stressT[4],
                                        stressT[5]};
    const PylithScalar devStressProdT = pylith::fekernels::Viscoelasticity::scalarProduct3D(devStressT, devStressT);
    const PylithScalar j2T = sqrt(0.5*devStressProdT);

    // Compute quantities based on strain at t = T + dt.
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar strainTpdt[6] = {disp_x[0*_dim+0],
                                        disp_x[1*_dim+1],
                                        disp_x[2*_dim+2],
                                        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]),
                                        0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]),
                                        0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0])};
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1] + strainTpdt[2])/3.0;
    const PylithScalar devStrainTpdt[6] = {strainTpdt[0] - meanStrainTpdt,
                                           strainTpdt[1] - meanStrainTpdt,
                                           strainTpdt[2] - meanStrainTpdt,
                                           strainTpdt[3],
                                           strainTpdt[4],
                                           strainTpdt[5]};
    const PylithScalar strainPPTpdt[6] = {devStrainTpdt[0] - visStrainT[0],
                                          devStrainTpdt[1] - visStrainT[1],
                                          devStrainTpdt[2] - visStrainT[2],
                                          devStrainTpdt[3] - visStrainT[3],
                                          devStrainTpdt[4] - visStrainT[4],
                                          devStrainTpdt[5] - visStrainT[5]};
    const PylithScalar strainPPInvar2Tpdt = 0.5*pylith::fekernels::Viscoelasticity::scalarProduct3D(strainPPTpdt, strainPPTpdt);
    const PylithScalar strainStressInvar2T = pylith::fekernels::Viscoelasticity::scalarProduct3D(strainPPTpdt, devStressT);

    // Finish defining parameters needed for root-finding algorithm.
    const PylithScalar b = strainPPInvar2Tpdt;
    const PylithScalar c = strainStressInvar2T*timeFac;
    const PylithScalar d = timeFac*j2T;
    PylithScalar j2Tpdt = 0.0;
    if ((b != 0.0) || (c != 0.0) || (d != 0.0)) {
        const PylithScalar j2InitialGuess = j2T;
        const PylithScalar stressScale = shearModulus;
        j2Tpdt = IsotropicPowerLawEffectiveStress::computeEffectiveStress(j2InitialGuess, stressScale, ae, b, c, d, powerLawAlpha,
                                                                          dt, j2T, powerLawExponent, powerLawReferenceStrainRate,
                                                                          powerLawReferenceStress);
    } // if
    // Compute deviatoric stresses from effective stress.
    const PylithScalar j2Tau = (1.0 - powerLawAlpha) * j2T + powerLawAlpha * j2Tpdt;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress),
                                                                  (powerLawExponent - 1.0))/powerLawReferenceStress;
    const PylithScalar factor1 = 1.0/(ae + powerLawAlpha*dt*gammaTau);
    const PylithScalar factor2 = timeFac*gammaTau;
    devStressTpdt[0] += factor1*(strainPPTpdt[0] - factor2*devStressT[0]);
    devStressTpdt[4] += factor1*(strainPPTpdt[1] - factor2*devStressT[1]);
    devStressTpdt[8] += factor1*(strainPPTpdt[2] - factor2*devStressT[2]);
    devStressTpdt[1] += factor1*(strainPPTpdt[3] - factor2*devStressT[3]);
    devStressTpdt[2] += factor1*(strainPPTpdt[5] - factor2*devStressT[5]);
    devStressTpdt[5] += factor1*(strainPPTpdt[4] - factor2*devStressT[4]);
    devStressTpdt[3] += devStressTpdt[1];
    devStressTpdt[6] += devStressTpdt[2];
    devStressTpdt[7] += devStressTpdt[5];

} // deviatoricStress


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate deviatoric stress tensor for 3-D isotropic power-law
 * viscoelastic material WITH reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = 2*shearModulus*visStrain_ij
 */
void
pylith::fekernels::IsotropicPowerLaw3D::deviatoricStress_refstate(const PylithInt dim,
                                                                  const PylithInt numS,
                                                                  const PylithInt numA,
                                                                  const PylithInt sOff[],
                                                                  const PylithInt sOff_x[],
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
                                                                  PylithScalar devStressTpdt[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = 0;
    const PylithInt i_rstrain = 1;
    const PylithInt i_shearModulus = 2;
    const PylithInt i_powerLawReferenceStrainRate = 3;
    const PylithInt i_powerLawReferenceStress = 4;
    const PylithInt i_powerLawExponent = 5;
    const PylithInt i_viscousStrain = 6;
    const PylithInt i_stress = 7;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 8);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(constants);

    // Constants.
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar powerLawReferenceStrainRate = a[aOff[i_powerLawReferenceStrainRate]];
    const PylithScalar powerLawReferenceStress = a[aOff[i_powerLawReferenceStress]];
    const PylithScalar powerLawExponent = a[aOff[i_powerLawExponent]];
    // ****** I think powerLawAlpha should actually be a kernel constant.
    const PylithScalar powerLawAlpha = 0.5;
    // const PylithScalar powerLawAlpha = a[aOff[i_powerLawAlpha]];
    const PylithScalar dt = constants[0];
    const PylithScalar ae = 1.0/(2.0*shearModulus);
    const PylithScalar timeFac = dt*(1.0 - powerLawAlpha);

    // Compute quantities based on reference stress.
    const PylithScalar* stressR = &a[aOff[i_rstress]];
    const PylithScalar meanStressR = (stressR[0] + stressR[1] + stressR[2])/3.0;
    const PylithScalar devStressR[6] = {stressR[0] - meanStressR,
                                        stressR[1] - meanStressR,
                                        stressR[2] - meanStressR,
                                        stressR[3],
                                        stressR[4],
                                        stressR[5]};
    const PylithScalar j2RSquared = 0.5*pylith::fekernels::Viscoelasticity::scalarProduct3D(devStressR, devStressR);

    // Compute quantities based on reference strain.
    const PylithScalar* strainR = &a[aOff[i_rstrain]];
    const PylithScalar meanStrainR = (strainR[0] + strainR[1] + strainR[2])/3.0;

    // Compute quantities based on stress at t = T.
    const PylithScalar* stressT = &a[aOff[i_stress]]; // stress_xx, stress_yy, stress_zz, stress_xy, stress_yz,
                                                      // stress_xz at t = T.
    const PylithScalar* visStrainT = &a[aOff[i_viscousStrain]]; // visStrain_xx, visStrain_yy, visStrain_zz,
                                                                // visStrain_xy at t = T.
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[6] = {stressT[0] - meanStressT,
                                        stressT[1] - meanStressT,
                                        stressT[2] - meanStressT,
                                        stressT[3],
                                        stressT[4],
                                        stressT[5]};
    const PylithScalar devStressProdT = pylith::fekernels::Viscoelasticity::scalarProduct3D(devStressT, devStressT);
    const PylithScalar j2T = sqrt(0.5*devStressProdT);

    // Compute quantities based on strain at t = T + dt.
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar strainTpdt[6] = {disp_x[0*_dim+0],
                                        disp_x[1*_dim+1],
                                        disp_x[2*_dim+2],
                                        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]),
                                        0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]),
                                        0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0])};
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1] + strainTpdt[2])/3.0 - meanStrainR;
    const PylithScalar strainPPTpdt[6] = {strainTpdt[0] - meanStrainTpdt - visStrainT[0] - strainR[0],
                                          strainTpdt[1] - meanStrainTpdt - visStrainT[1] - strainR[1],
                                          strainTpdt[2] - meanStrainTpdt - visStrainT[2] - strainR[2],
                                          strainTpdt[3] - visStrainT[3] - strainR[3],
                                          strainTpdt[4] - visStrainT[4] - strainR[4],
                                          strainTpdt[5] - visStrainT[5] - strainR[5]};
    const PylithScalar strainPPInvar2Tpdt = 0.5*pylith::fekernels::Viscoelasticity::scalarProduct3D(strainPPTpdt, strainPPTpdt);
    const PylithScalar strainStressInvar2T = pylith::fekernels::Viscoelasticity::scalarProduct3D(strainPPTpdt, devStressT);
    const PylithScalar strainStressInvar2R = pylith::fekernels::Viscoelasticity::scalarProduct3D(strainPPTpdt, devStressR);
    const PylithScalar stressInvar2RT = pylith::fekernels::Viscoelasticity::scalarProduct3D(devStressT, devStressR);

    // Finish defining parameters needed for root-finding algorithm.
    const PylithScalar b = strainPPInvar2Tpdt + ae*strainStressInvar2R + ae*ae*j2RSquared;
    const PylithScalar c = strainStressInvar2T*timeFac + ae*stressInvar2RT;
    const PylithScalar d = timeFac*j2T;
    PylithScalar j2Tpdt = 0.0;
    if ((b != 0.0) || (c != 0.0) || (d != 0.0)) {
        const PylithScalar j2InitialGuess = j2T;
        const PylithScalar stressScale = shearModulus;
        j2Tpdt = IsotropicPowerLawEffectiveStress::computeEffectiveStress(j2InitialGuess, stressScale, ae, b, c, d, powerLawAlpha,
                                                                          dt, j2T, powerLawExponent, powerLawReferenceStrainRate,
                                                                          powerLawReferenceStress);
    } // if
    // Compute deviatoric stresses from effective stress.
    const PylithScalar j2Tau = (1.0 - powerLawAlpha) * j2T + powerLawAlpha * j2Tpdt;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress),
                                                                  (powerLawExponent - 1.0))/powerLawReferenceStress;
    const PylithScalar factor1 = 1.0/(ae + powerLawAlpha*dt*gammaTau);
    const PylithScalar factor2 = timeFac*gammaTau;
    devStressTpdt[0] += factor1*(strainPPTpdt[0] - factor2*devStressT[0] + ae*devStressR[0]);
    devStressTpdt[4] += factor1*(strainPPTpdt[1] - factor2*devStressT[1] + ae*devStressR[1]);
    devStressTpdt[8] += factor1*(strainPPTpdt[2] - factor2*devStressT[2] + ae*devStressR[2]);
    devStressTpdt[1] += factor1*(strainPPTpdt[3] - factor2*devStressT[3] + ae*devStressR[3]);
    devStressTpdt[2] += factor1*(strainPPTpdt[5] - factor2*devStressT[5] + ae*devStressR[5]);
    devStressTpdt[5] += factor1*(strainPPTpdt[4] - factor2*devStressT[4] + ae*devStressR[4]);
    devStressTpdt[3] += devStressTpdt[1];
    devStressTpdt[6] += devStressTpdt[2];
    devStressTpdt[7] += devStressTpdt[5];

} // deviatoricStress_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Update stress for a 3D power-law viscoelastic material.
 *
 * IMPORTANT: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 *
 */
void
pylith::fekernels::IsotropicPowerLaw3D::updateStress(const PylithInt dim,
                                                     const PylithInt numS,
                                                     const PylithInt numA,
                                                     const PylithInt sOff[],
                                                     const PylithInt sOff_x[],
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
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

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
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 6; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[6] = {aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate], aOff[i_powerLawReferenceStress],
                                  aOff[i_powerLawExponent], aOff[i_viscousStrain], aOff[i_stress]};

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    IsotropicLinearElasticity3D::meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL, a, a_t,
                                            NULL, t, x, numConstants, constants, stressTensor);

    // Compute deviatoric stress tensor.
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);

    stress[0] += stressTensor[0];
    stress[1] += stressTensor[4];
    stress[2] += stressTensor[8];
    stress[3] += stressTensor[1];
    stress[4] += stressTensor[5];
    stress[5] += stressTensor[2];

} // updateStress


// ---------------------------------------------------------------------------------------------------------------------
/* Update stress for a 3D power-law viscoelastic material WITH reference stress and strain.
 *
 * IMPORTANT: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 *
 */
void
pylith::fekernels::IsotropicPowerLaw3D::updateStress_refstate(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
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
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-9;
    const PylithInt i_rstrain = numA-8;
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 8; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[8] = {aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate],
                                  aOff[i_powerLawReferenceStress], aOff[i_powerLawExponent], aOff[i_viscousStrain],
                                  aOff[i_stress]};

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    IsotropicLinearElasticity3D::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffMean, NULL,
                                                     a, a_t, NULL, t, x, numConstants, constants, stressTensor);

    // Compute deviatoric stress tensor.
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);

    stress[0] += stressTensor[0];
    stress[1] += stressTensor[4];
    stress[2] += stressTensor[8];
    stress[3] += stressTensor[1];
    stress[4] += stressTensor[5];
    stress[5] += stressTensor[2];

} // updateStress_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Update viscous strain for a 3D power-law viscoelastic material.
 *
 * :ATTENTION: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 */
void
pylith::fekernels::IsotropicPowerLaw3D::updateViscousStrain(const PylithInt dim,
                                                            const PylithInt numS,
                                                            const PylithInt numA,
                                                            const PylithInt sOff[],
                                                            const PylithInt sOff_x[],
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
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

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
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(constants);

    // Constants.
    const PylithScalar powerLawReferenceStrainRate = a[aOff[i_powerLawReferenceStrainRate]];
    const PylithScalar powerLawReferenceStress = a[aOff[i_powerLawReferenceStress]];
    const PylithScalar powerLawExponent = a[aOff[i_powerLawExponent]];
    // ****** I think powerLawAlpha should actually be a kernel constant.
    const PylithScalar powerLawAlpha = 0.5;
    // const PylithScalar powerLawAlpha = a[aOff[i_powerLawAlpha]];
    const PylithScalar dt = constants[0];

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 6; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[6] = {aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate], aOff[i_powerLawReferenceStress],
                                  aOff[i_powerLawExponent], aOff[i_viscousStrain], aOff[i_stress]};

    // Compute current deviatoric stress.
    PylithScalar devStressTpdt[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, devStressTpdt);
    const PylithScalar devStressProdTpdt = pylith::fekernels::Viscoelasticity::scalarProduct3D(devStressTpdt, devStressTpdt);
    const PylithScalar j2Tpdt = sqrt(0.5*devStressProdTpdt);

    // Compute stress quantities at time T.
    const PylithScalar* stressT = &a[aOff[i_stress]]; // stress_xx, stress_yy, stress_zz, stress_xy, stress_yz,
                                                      // stress_xz at t = T.
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[6] = {stressT[0] - meanStressT,
                                        stressT[1] - meanStressT,
                                        stressT[2] - meanStressT,
                                        stressT[3],
                                        stressT[4],
                                        stressT[5]};
    const PylithScalar devStressProdT = pylith::fekernels::Viscoelasticity::scalarProduct3D(devStressT, devStressT);
    const PylithScalar j2T = sqrt(0.5*devStressProdT);

    // Compute quantities at intermediate time.
    const PylithScalar devStressTau[6] = {(1.0 - powerLawAlpha)*devStressT[0] + powerLawAlpha*devStressTpdt[0],
                                          (1.0 - powerLawAlpha)*devStressT[1] + powerLawAlpha*devStressTpdt[1],
                                          (1.0 - powerLawAlpha)*devStressT[2] + powerLawAlpha*devStressTpdt[2],
                                          (1.0 - powerLawAlpha)*devStressT[3] + powerLawAlpha*devStressTpdt[3],
                                          (1.0 - powerLawAlpha)*devStressT[4] + powerLawAlpha*devStressTpdt[4],
                                          (1.0 - powerLawAlpha)*devStressT[5] + powerLawAlpha*devStressTpdt[5]};
    const PylithScalar j2Tau = (1.0 - powerLawAlpha)*j2T + powerLawAlpha*j2Tpdt;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress),
                                                                  (powerLawExponent - 1.0))/powerLawReferenceStress;

    // Update viscous strain.
    visStrain[0] += dt*gammaTau*devStressTau[0];
    visStrain[1] += dt*gammaTau*devStressTau[1];
    visStrain[2] += dt*gammaTau*devStressTau[2];
    visStrain[3] += dt*gammaTau*devStressTau[3];
    visStrain[4] += dt*gammaTau*devStressTau[4];
    visStrain[5] += dt*gammaTau*devStressTau[5];

} // updateViscousStrain


// ---------------------------------------------------------------------------------------------------------------------
/* Update viscous strain for a 3D power-law viscoelastic material WITH reference stress and strain.
 *
 * :ATTENTION: The order of the auxiliary field and solution field are reversed compared to the residual and Jacobian
 * kernels.
 */
void
pylith::fekernels::IsotropicPowerLaw3D::updateViscousStrain_refstate(const PylithInt dim,
                                                                     const PylithInt numS,
                                                                     const PylithInt numA,
                                                                     const PylithInt sOff[],
                                                                     const PylithInt sOff_x[],
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
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-9;
    const PylithInt i_rstrain = numA-8;
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(constants);

    // Constants.
    const PylithScalar powerLawReferenceStrainRate = a[aOff[i_powerLawReferenceStrainRate]];
    const PylithScalar powerLawReferenceStress = a[aOff[i_powerLawReferenceStress]];
    const PylithScalar powerLawExponent = a[aOff[i_powerLawExponent]];
    // ****** I think powerLawAlpha should actually be a kernel constant.
    const PylithScalar powerLawAlpha = 0.5;
    // const PylithScalar powerLawAlpha = a[aOff[i_powerLawAlpha]];
    const PylithScalar dt = constants[0];

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 8; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[8] = {aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate],
                                  aOff[i_powerLawReferenceStress], aOff[i_powerLawExponent], aOff[i_viscousStrain],
                                  aOff[i_stress]};

    // Compute current deviatoric stress.
    PylithScalar devStressTpdt[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, devStressTpdt);
    const PylithScalar devStressProdTpdt = pylith::fekernels::Viscoelasticity::scalarProduct3D(devStressTpdt, devStressTpdt);
    const PylithScalar j2Tpdt = sqrt(0.5*devStressProdTpdt);

    // Compute stress quantities at time T.
    const PylithScalar* stressT = &a[aOff[i_stress]]; // stress_xx, stress_yy, stress_zz, stress_xy, stress_yz,
                                                      // stress_xz at t = T.
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[6] = {stressT[0] - meanStressT,
                                        stressT[1] - meanStressT,
                                        stressT[2] - meanStressT,
                                        stressT[3],
                                        stressT[4],
                                        stressT[5]};
    const PylithScalar devStressProdT = pylith::fekernels::Viscoelasticity::scalarProduct3D(devStressT, devStressT);
    const PylithScalar j2T = sqrt(0.5*devStressProdT);

    // Compute quantities at intermediate time.
    const PylithScalar devStressTau[6] = {(1.0 - powerLawAlpha)*devStressT[0] + powerLawAlpha*devStressTpdt[0],
                                          (1.0 - powerLawAlpha)*devStressT[1] + powerLawAlpha*devStressTpdt[1],
                                          (1.0 - powerLawAlpha)*devStressT[2] + powerLawAlpha*devStressTpdt[2],
                                          (1.0 - powerLawAlpha)*devStressT[3] + powerLawAlpha*devStressTpdt[3],
                                          (1.0 - powerLawAlpha)*devStressT[4] + powerLawAlpha*devStressTpdt[4],
                                          (1.0 - powerLawAlpha)*devStressT[5] + powerLawAlpha*devStressTpdt[5]};
    const PylithScalar j2Tau = (1.0 - powerLawAlpha)*j2T + powerLawAlpha*j2Tpdt;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress),
                                                                  (powerLawExponent - 1.0))/powerLawReferenceStress;

    // Update viscous strain.
    visStrain[0] += dt*gammaTau*devStressTau[0];
    visStrain[1] += dt*gammaTau*devStressTau[1];
    visStrain[2] += dt*gammaTau*devStressTau[2];
    visStrain[3] += dt*gammaTau*devStressTau[3];
    visStrain[4] += dt*gammaTau*devStressTau[4];
    visStrain[5] += dt*gammaTau*devStressTau[5];

} // updateViscousStrain_refstate


// ---------------------------------------------------------------------------------------------------------------------
// Calculate stress for 3-D isotropic power-law viscoelastic material WITHOUT a reference stress and strain.
void
pylith::fekernels::IsotropicPowerLaw3D::cauchyStress(const PylithInt dim,
                                                     const PylithInt numS,
                                                     const PylithInt numA,
                                                     const PylithInt sOff[],
                                                     const PylithInt sOff_x[],
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
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

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
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 6; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[6] = {aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate], aOff[i_powerLawReferenceStress],
                                  aOff[i_powerLawExponent], aOff[i_viscousStrain], aOff[i_stress]};

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    IsotropicLinearElasticity3D::meanStress(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                            aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);

    // Compute deviatoric stress (4 components).
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                     t, x, numConstants, constants, stressTensor);

    stressVector[0] = stressTensor[0];
    stressVector[1] = stressTensor[4];
    stressVector[2] = stressTensor[8];
    stressVector[3] = stressTensor[1];
    stressVector[4] = stressTensor[5];
    stressVector[5] = stressTensor[2];

} // stress


// ---------------------------------------------------------------------------------------------------------------------
// Calculate stress for 3-D isotropic power-law viscoelastic material WITH a reference stress/strain.
void
pylith::fekernels::IsotropicPowerLaw3D::cauchyStress_refstate(const PylithInt dim,
                                                              const PylithInt numS,
                                                              const PylithInt numA,
                                                              const PylithInt sOff[],
                                                              const PylithInt sOff_x[],
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
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-9;
    const PylithInt i_rstrain = numA-8;
    const PylithInt i_shearModulus = numA-7;
    const PylithInt i_bulkModulus = numA-6;
    const PylithInt i_powerLawReferenceStrainRate = numA-5;
    const PylithInt i_powerLawReferenceStress = numA-4;
    const PylithInt i_powerLawExponent = numA-3;
    const PylithInt i_viscousStrain = numA-2;
    const PylithInt i_stress = numA-1;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_powerLawReferenceStrainRate] >= 0);
    assert(aOff[i_powerLawReferenceStress] >= 0);
    assert(aOff[i_powerLawExponent] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_stress] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Pass bulk modulus, reference stress, and reference strain.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 8; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[8] = {aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus], aOff[i_powerLawReferenceStrainRate],
                                  aOff[i_powerLawReferenceStress], aOff[i_powerLawExponent], aOff[i_viscousStrain],
                                  aOff[i_stress]};

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    IsotropicLinearElasticity3D::meanStress_refstate(_dim, _numS, numAMean, sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                     aOffMean, NULL, a, a_t, NULL, t, x, numConstants, constants, stressTensor);

    // Compute deviatoric stress tensor.
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                              t, x, numConstants, constants, stressTensor);

    stressVector[0] = stressTensor[0];
    stressVector[1] = stressTensor[4];
    stressVector[2] = stressTensor[8];
    stressVector[3] = stressTensor[1];
    stressVector[4] = stressTensor[5];
    stressVector[5] = stressTensor[2];

} // stress_refstate


// =====================================================================================================================
// Functions for isotropic, power-law viscoelastic material.
// =====================================================================================================================
// ---------------------------------------------------------------------------------------------------------------------
// Get effective stress from initial guess.
PylithScalar
pylith::fekernels::IsotropicPowerLawEffectiveStress::computeEffectiveStress(const PylithScalar j2InitialGuess,
                                                                            const PylithScalar stressScale,
                                                                            const PylithScalar ae,
                                                                            const PylithScalar b,
                                                                            const PylithScalar c,
                                                                            const PylithScalar d,
                                                                            const PylithScalar powerLawAlpha,
                                                                            const PylithScalar dt,
                                                                            const PylithScalar j2T,
                                                                            const PylithScalar powerLawExponent,
                                                                            const PylithScalar powerLawReferenceStrainRate,
                                                                            const PylithScalar powerLawReferenceStress) { //
                                                                                                                          //
                                                                                                                          //
                                                                                                                          //
                                                                                                                          // computeEffectiveStress
    // Check parameters
    assert(j2InitialGuess >= 0.0);
    // If initial guess is too low, use stress scale instead.
    const PylithScalar xMin = 1.0e-10;

    // Bracket the root.
    PylithScalar x1 = 0.0;
    PylithScalar x2 = 0.0;
    if (j2InitialGuess > xMin) {
        x1 = j2InitialGuess - 0.5 * j2InitialGuess;
        x2 = j2InitialGuess + 0.5 * j2InitialGuess;
    } else {
        x1 = stressScale - 0.5 * stressScale;
        x2 = stressScale + 0.5 * stressScale;
    } // else

    _bracket(&x1, &x2, ae, b, c, d, powerLawAlpha, dt, j2T, powerLawExponent, powerLawReferenceStrainRate,
             powerLawReferenceStress);

    // Find effective stress using Newton's method with bisection.
    PylithScalar effStress = _search(x1, x2, ae, b, c, d, powerLawAlpha, dt, j2T, powerLawExponent,
                                     powerLawReferenceStrainRate, powerLawReferenceStress);

    PetscLogFlops(4); // Log flops

    return effStress;
} // computeEffectiveStress


// ----------------------------------------------------------------------
// Calculate effective stress function for power-law material.
PylithScalar
pylith::fekernels::IsotropicPowerLawEffectiveStress::_effStressFunc(const PylithScalar j2Tpdt,
                                                                    const PylithScalar ae,
                                                                    const PylithScalar b,
                                                                    const PylithScalar c,
                                                                    const PylithScalar d,
                                                                    const PylithScalar powerLawAlpha,
                                                                    const PylithScalar dt,
                                                                    const PylithScalar j2T,
                                                                    const PylithScalar powerLawExponent,
                                                                    const PylithScalar powerLawReferenceStrainRate,
                                                                    const PylithScalar powerLawReferenceStress) { // _effStressFunc
    const PylithScalar factor1 = 1.0 - powerLawAlpha;
    const PylithScalar j2Tau = factor1*j2T + powerLawAlpha*j2Tpdt;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress), (powerLawExponent - 1.0))/
                                  powerLawReferenceStress;
    const PylithScalar a = ae + powerLawAlpha*dt*gammaTau;
    const PylithScalar y = a*a*j2Tpdt*j2Tpdt - b + c*gammaTau - d*d*gammaTau*gammaTau;

    return y;
} // _effStressFunc


// ---------------------------------------------------------------------------------------------------------------------
// Calculate effective stress function and its derivative for power-law material.
void
pylith::fekernels::IsotropicPowerLawEffectiveStress::_effStressFuncDerivFunc(PylithScalar* func,
                                                                             PylithScalar* dfunc,
                                                                             const PylithScalar j2Tpdt,
                                                                             const PylithScalar ae,
                                                                             const PylithScalar b,
                                                                             const PylithScalar c,
                                                                             const PylithScalar d,
                                                                             const PylithScalar powerLawAlpha,
                                                                             const PylithScalar dt,
                                                                             const PylithScalar j2T,
                                                                             const PylithScalar powerLawExponent,
                                                                             const PylithScalar powerLawReferenceStrainRate,
                                                                             const PylithScalar powerLawReferenceStress) { //
                                                                                                                           //
                                                                                                                           //
                                                                                                                           //
                                                                                                                           // _effStressFuncDerivFunc
    PylithScalar y = *func;
    PylithScalar dy = *dfunc;

    const PylithScalar factor1 = 1.0 - powerLawAlpha;
    const PylithScalar j2Tau = factor1*j2T + powerLawAlpha*j2Tpdt;
    const PylithScalar gammaTau = powerLawReferenceStrainRate*pow((j2Tau/powerLawReferenceStress), (powerLawExponent - 1.0))/
                                  powerLawReferenceStress;
    const PylithScalar dGammaTau = powerLawReferenceStrainRate*powerLawAlpha*(powerLawExponent - 1.0)*
                                   pow((j2Tau/powerLawReferenceStress), (powerLawExponent - 2.0))/(powerLawReferenceStress*powerLawReferenceStress);
    const PylithScalar a = ae + powerLawAlpha*dt*gammaTau;
    y = a*a*j2Tpdt*j2Tpdt - b + c*gammaTau - d*d*gammaTau*gammaTau;
    dy = 2.0*a*a*j2Tpdt + dGammaTau*(2.0*a*powerLawAlpha*dt*j2Tpdt*j2Tpdt + c - 2.0*d*d*gammaTau);

    *func = y;
    *dfunc = dy;

} // _effStressFuncDerivFunc


// ----------------------------------------------------------------------
// Bracket effective stress.
void
pylith::fekernels::IsotropicPowerLawEffectiveStress::_bracket(PylithScalar* px1,
                                                              PylithScalar* px2,
                                                              const PylithScalar ae,
                                                              const PylithScalar b,
                                                              const PylithScalar c,
                                                              const PylithScalar d,
                                                              const PylithScalar powerLawAlpha,
                                                              const PylithScalar dt,
                                                              const PylithScalar j2T,
                                                              const PylithScalar powerLawExponent,
                                                              const PylithScalar powerLawReferenceStrainRate,
                                                              const PylithScalar powerLawReferenceStress) { // _bracket
    // Arbitrary number of iterations to bracket the root
    const int maxIterations = 50;

    // Arbitrary factor by which to increase the brackets.
    const PylithScalar bracketFactor = 2;
    // Minimum allowed value for effective stress.
    const PylithScalar xMin = 0.0;
    PylithScalar x1 = *px1;
    PylithScalar x2 = *px2;

    PylithScalar funcValue1 = _effStressFunc(x1, ae, b, c, d, powerLawAlpha, dt, j2T, powerLawExponent,
                                             powerLawReferenceStrainRate, powerLawReferenceStress);
    PylithScalar funcValue2 = _effStressFunc(x2, ae, b, c, d, powerLawAlpha, dt, j2T, powerLawExponent,
                                             powerLawReferenceStrainRate, powerLawReferenceStress);

    int iteration = 0;
    bool bracketed = false;
    while (iteration < maxIterations) {
        if ((funcValue1 * funcValue2) < 0.0) {
            bracketed = true;
            break;
        } // if

        if (fabs(funcValue1) < fabs(funcValue2)) {
            x1 += bracketFactor * (x1 - x2);
            x1 = std::max(x1, xMin);
            funcValue1 = _effStressFunc(x1, ae, b, c, d, powerLawAlpha, dt, j2T, powerLawExponent,
                                        powerLawReferenceStrainRate, powerLawReferenceStress);
        } else {
            x2 += bracketFactor * (x1 - x2);
            x2 = std::max(x2, xMin);
            funcValue2 = _effStressFunc(x2, ae, b, c, d, powerLawAlpha, dt, j2T, powerLawExponent,
                                        powerLawReferenceStrainRate, powerLawReferenceStress);
        } // else
        ++iteration;
    } // while

    *px1 = x1;
    *px2 = x2;

    PetscLogFlops(5 * iteration);
    if (!bracketed) {
        throw std::runtime_error("Unable to bracket effective stress.");
    }
} // _bracket


// ----------------------------------------------------------------------
// Find root using Newton's method with bisection.
PylithScalar
pylith::fekernels::IsotropicPowerLawEffectiveStress::_search(const PylithScalar x1,
                                                             const PylithScalar x2,
                                                             const PylithScalar ae,
                                                             const PylithScalar b,
                                                             const PylithScalar c,
                                                             const PylithScalar d,
                                                             const PylithScalar powerLawAlpha,
                                                             const PylithScalar dt,
                                                             const PylithScalar j2T,
                                                             const PylithScalar powerLawExponent,
                                                             const PylithScalar powerLawReferenceStrainRate,
                                                             const PylithScalar powerLawReferenceStress) { // _search
    // Arbitrary number of iterations to find the root
    const int maxIterations = 100;

    // Desired accuracy for root. This is a bit arbitrary for now.
    const PylithScalar accuracy = 1.0e-16;

    // Organize search so that _effStressFunc(xLow) is less than zero.
    PylithScalar funcValueLow = _effStressFunc(x1, ae, b, c, d, powerLawAlpha, dt, j2T, powerLawExponent,
                                               powerLawReferenceStrainRate, powerLawReferenceStress);
    PylithScalar funcValueHigh = _effStressFunc(x2, ae, b, c, d, powerLawAlpha, dt, j2T, powerLawExponent,
                                                powerLawReferenceStrainRate, powerLawReferenceStress);
    assert(funcValueLow * funcValueHigh <= 0.0);

    PylithScalar effStress = 0.0;
    PylithScalar xLow = 0.0;
    PylithScalar xHigh = 0.0;
    bool converged = false;

    if (funcValueLow < 0.0) {
        xLow = x1;
        xHigh = x2;
    } else {
        xLow = x2;
        xHigh = x1;
    } // if/else

    effStress = 0.5 * (x1 + x2);
    PylithScalar dxPrevious = fabs(x2 - x1);
    PylithScalar dx = dxPrevious;
    PylithScalar funcValue = 0.0;
    PylithScalar funcDeriv = 0.0;
    PylithScalar funcXHigh = 0.0;
    PylithScalar funcXLow = 0.0;
    _effStressFuncDerivFunc(&funcValue, &funcDeriv, effStress, ae, b, c, d, powerLawAlpha, dt, j2T, powerLawExponent,
                            powerLawReferenceStrainRate, powerLawReferenceStress);
    int iteration = 0;

    while (iteration < maxIterations) {
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
        _effStressFuncDerivFunc(&funcValue, &funcDeriv, effStress, ae, b, c, d, powerLawAlpha, dt, j2T, powerLawExponent,
                                powerLawReferenceStrainRate, powerLawReferenceStress);
        if (funcValue < 0.0) {
            xLow = effStress;
        } else {
            xHigh = effStress;
        } // else
        ++iteration;
    } // while

    if (converged == false) {
        throw std::runtime_error("Cannot find root of effective stress function.");
    }

    PetscLogFlops(5 + 15 * iteration); // Log flops

    return effStress;
} // _search


// End of file
