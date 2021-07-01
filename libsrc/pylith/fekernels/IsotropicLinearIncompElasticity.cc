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

#include "pylith/fekernels/IsotropicLinearIncompElasticity.hh"

#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity kernels

#include <cassert> // USES assert()

// =====================================================================================================================
// Kernels for isotropic, linear incompressible elasticity (dimension independent).
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
// f0p function for isotropic linear incompressible elasticity plane strain WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticity::f0p(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
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
    // Incoming solution subfields
    const PylithInt i_disp = 0;
    const PylithInt i_pres = 1;

    // Incoming auxiliary subfields
    const PylithInt i_bulkModulus = numA-1;

    assert(numS >= 2);
    assert(sOff);
    assert(sOff[i_pres] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s);
    assert(s_x);
    assert(numA >= 1);
    assert(aOff);
    assert(aOff[i_bulkModulus] >= 0);
    assert(a);
    assert(f0);

    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar pressure = s[sOff[i_pres]];

    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    PylithScalar strainTrace = 0;
    for (PylithInt i = 0; i < dim; ++i) {
        strainTrace += disp_x[i*dim+i];
    } // for
    f0[0] += strainTrace + pressure / bulkModulus;
} // f0p


// ---------------------------------------------------------------------------------------------------------------------
// f1 function for isotropic linear incompressible elasticity plane strain WITH reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticity::f0p_refstate(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
    // Incoming solution subfields
    const PylithInt i_disp = 0;
    const PylithInt i_pres = 1;

    // Incoming auxiliary subfields
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_bulkModulus = numA-1;

    assert(numS >= 2);
    assert(sOff);
    assert(sOff[i_pres] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(s);
    assert(s_x);
    assert(numA >= 4);
    assert(aOff);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(a);
    assert(f0);

    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar pressure = s[sOff[i_pres]];

    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, ...
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, ...

    const PylithScalar meanRefStress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithScalar refStrainTrace = (refstrain[0] + refstrain[1] + refstrain[2]);

    PylithScalar strainTrace = 0;
    for (PylithInt i = 0; i < dim; ++i) {
        strainTrace += disp_x[i*dim+i];
    } // for
    f0[0] += (strainTrace - refStrainTrace) + (pressure + meanRefStress) / bulkModulus;
} // f0p


// ---------------------------------------------------------------------------------------------------------------------
// Jf0_pp function for isotropic linear incompressible elasticity.
void
pylith::fekernels::IsotropicLinearIncompElasticity::Jf0pp(const PylithInt dim,
                                                          const PylithInt numS,
                                                          const PylithInt numA,
                                                          const PylithInt sOff[],
                                                          const PylithInt sOff_x[],
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


// =====================================================================================================================
// Kernels for 2-D plane strain isotropic, linear incompressible elasticity.
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
// f1 function for isotropic linear incompressible elasticity plane strain WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::f1u(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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

    // Incoming solution subfields.
    const PylithInt i_disp = 0;
    const PylithInt i_pres = 1;

    // Incoming auxiliary subfields.
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 2);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff[i_pres] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(sOff_x[i_pres] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(f1);

    const PylithInt numSMean = 2; // Number passed on to mean stress kernel.
    const PylithInt sOffMean[2] = { sOff[i_disp], sOff[i_pres] };
    const PylithInt sOffMean_x[2] = { sOff_x[i_disp], sOff_x[i_pres] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numSDev = 1; // Number passed on to deviatoric stress kernel.
    const PylithInt sOffDev[1] = { sOff[i_disp] };
    const PylithInt sOffDev_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, numSMean, numAMean, sOffMean, sOffMean_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    IsotropicLinearElasticityPlaneStrain::deviatoricStress(_dim, numSDev, numADev, sOffDev, sOffDev_x, s, s_t, s_x,
                                                           aOffDev, NULL, a, a_t, NULL,
                                                           t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1u


// ---------------------------------------------------------------------------------------------------------------------
// f1 function for isotropic linear incompressible elasticity plane strain WITH reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::f1u_refstate(const PylithInt dim,
                                                                            const PylithInt numS,
                                                                            const PylithInt numA,
                                                                            const PylithInt sOff[],
                                                                            const PylithInt sOff_x[],
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

    // Incoming solution subfields.
    const PylithInt i_disp = 0;
    const PylithInt i_pres = 1;

    // Incoming auxiliary subfields.
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 4);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff[i_pres] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(sOff_x[i_pres] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f1);

    const PylithInt numSMean = 2; // Number passed on to mean stress kernel.
    const PylithInt sOffMean[2] = { sOff[i_disp], sOff[i_pres] };
    const PylithInt sOffMean_x[2] = { sOff_x[i_disp], sOff_x[i_pres] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numSDev = 1; // Number passed on to deviatoric stress kernel.
    const PylithInt sOffDev[1] = { sOff[i_disp] };
    const PylithInt sOffDev_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, numSMean, numAMean, sOffMean, sOffMean_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    IsotropicLinearElasticityPlaneStrain::deviatoricStress(_dim, numSDev, numADev, sOffDev, sOffDev_x, s, s_t, s_x,
                                                           aOffDev, NULL, a, a_t, NULL,
                                                           t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1u_refstate


// ---------------------------------------------------------------------------------------------------------------------
// Jf3_uu entry function for 2-D plane strain isotropic linear incompressible elasticity.
void
pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::Jf3uu(const PylithInt dim,
                                                                     const PylithInt numS,
                                                                     const PylithInt numA,
                                                                     const PylithInt sOff[],
                                                                     const PylithInt sOff_x[],
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


// ---------------------------------------------------------------------------------------------------------------------
// Calculate stress for 2-D plane strain isotropic linear incompressible elasticity WITHOUT a reference stress and
void
pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::cauchyStress(const PylithInt dim,
                                                                            const PylithInt numS,
                                                                            const PylithInt numA,
                                                                            const PylithInt sOff[],
                                                                            const PylithInt sOff_x[],
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
    const PylithInt i_pres = 1;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff[i_pres] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(sOff_x[i_pres] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(stressVector);

    const PylithInt numSDev = 1; // Number passed on to deviatoric stress kernel.
    const PylithInt sOffDev[1] = { sOff[i_disp] };
    const PylithInt sOffDev_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
    meanStress(_dim, numS, numAMean, sOff, sOff_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    IsotropicLinearElasticityPlaneStrain::deviatoricStress(_dim, numSDev, numADev, sOffDev, sOffDev_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
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


// ---------------------------------------------------------------------------------------------------------------------
// Calculate stress for 2-D plane strain isotropic linear incompressible elasticity WITH a reference stress/strain.
void
pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::cauchyStress_refstate(const PylithInt dim,
                                                                                     const PylithInt numS,
                                                                                     const PylithInt numA,
                                                                                     const PylithInt sOff[],
                                                                                     const PylithInt sOff_x[],
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
    assert(numS >= 2);
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
    assert(stressVector);

    const PylithInt numSDev = 1; // Number passed on to deviatoric stress kernels.
    const PylithInt sOffDev[1] = { sOff[i_disp] };
    const PylithInt sOffDev_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };
    meanStress_refstate(_dim, numS, numAMean, sOff, sOff_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    IsotropicLinearElasticityPlaneStrain::deviatoricStress_refstate(_dim, numSDev, numADev, sOffDev, sOffDev_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                                                                    t, x, numConstants, constants, stressTensor);

    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar* rstress = &a[aOff[i_rstress]];
    const PylithScalar stress_zz = rstress[2] +
                                   0.5*lambda/(lambda+shearModulus) * (stressTensor[0*_dim+0]-rstress[0] + stressTensor[1*_dim+1]-rstress[1]);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stress_zz;
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
} // cauchyStress_refstate


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate mean stress for 2-D plane strain isotropic, incompressible linear
 * elasticity WITHOUT reference stress and strain.
 *
 * meanStress = bulkModulus * strain_kk
 *
 * stress += meanStress * delta_ij
 */
void
pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::meanStress(const PylithInt dim,
                                                                          const PylithInt numS,
                                                                          const PylithInt numA,
                                                                          const PylithInt sOff[],
                                                                          const PylithInt sOff_x[],
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
                                                                          PylithScalar stressTensor[]) {
    const PylithInt _dim = 2;

    const PylithInt i_pres = 1;

    assert(_dim == dim);
    assert(numS >= 2);
    assert(sOff);
    assert(sOff[i_pres] >= 0);
    assert(s);
    assert(stressTensor);

    const PylithScalar pressure = s[sOff[i_pres]];

    for (PylithInt i = 0; i < _dim; ++i) {
        stressTensor[i*_dim+i] -= pressure;
    } // for
} // meanStress


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate mean stress for 2-D plane strain isotropic, incompressible linear
 * elasticity WITH reference stress and reference strain.
 *
 * We compute the stress relative to a reference stress/strain state.
 *
 * meanStress = meanRefStress + bulkModulus * (strain_kk - refstrain_kk)
 *
 * stress += meanStress * delta_ij
 */
void
pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::meanStress_refstate(const PylithInt dim,
                                                                                   const PylithInt numS,
                                                                                   const PylithInt numA,
                                                                                   const PylithInt sOff[],
                                                                                   const PylithInt sOff_x[],
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
                                                                                   PylithScalar stressTensor[]) {
    const PylithInt _dim = 2;

    // Incoming solution subfields.
    const PylithInt i_pres = 1;

    // Incoming auxiliary subfields.
    const PylithInt i_rstress = numA-4;

    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 4);
    assert(sOff);
    assert(sOff[i_pres] >= 0);
    assert(aOff);
    assert(aOff[i_rstress] >= 0);
    assert(s);
    assert(a);
    assert(stressTensor);

    const PylithScalar pressure = s[sOff[i_pres]];
    const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, ...

    const PylithScalar meanRefStress = (refStress[0] + refStress[1] + refStress[3]) / 3.0;
    const PylithScalar meanStress = meanRefStress - pressure;

    for (PylithInt i = 0; i < _dim; ++i) {
        stressTensor[i*_dim+i] += meanStress;
    } // for
} // meanStress_refstate


// =====================================================================================================================
// Kernels for 3-D isotropic, linear incompressible elasticity.
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
// f1 function for 3-D isotropic linear incompressible elasticity WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticity3D::f1u(const PylithInt dim,
                                                          const PylithInt numS,
                                                          const PylithInt numA,
                                                          const PylithInt sOff[],
                                                          const PylithInt sOff_x[],
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

    // Incoming solution subfields.
    const PylithInt i_disp = 0;
    const PylithInt i_pres = 1;

    // Incoming auxiliary subfields.
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 2);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff[i_pres] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(sOff_x[i_pres] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(f1);

    const PylithInt numSMean = 2; // Number passed on to mean stress kernel.
    const PylithInt sOffMean[2] = { sOff[i_disp], sOff[i_pres] };
    const PylithInt sOffMean_x[2] = { sOff_x[i_disp], sOff_x[i_pres] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numSDev = 1; // Number passed on to deviatoric stress kernel.
    const PylithInt sOffDev[1] = { sOff[i_disp] };
    const PylithInt sOffDev_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, numSMean, numAMean, sOffMean, sOffMean_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    IsotropicLinearElasticity3D::deviatoricStress(_dim, numSDev, numADev, sOffDev, sOffDev_x, s, s_t, s_x,
                                                  aOffDev, NULL, a, a_t, NULL,
                                                  t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1u


// ---------------------------------------------------------------------------------------------------------------------
// f1 function for 3-D isotropic linear incompressible elasticity WITH reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticity3D::f1u_refstate(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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

    // Incoming solution subfields.
    const PylithInt i_disp = 0;
    const PylithInt i_pres = 1;

    // Incoming auxiliary subfields.
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 4);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff[i_pres] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(sOff_x[i_pres] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(f1);

    const PylithInt numSMean = 2; // Number passed on to mean stress kernel.
    const PylithInt sOffMean[2] = { sOff[i_disp], sOff[i_pres] };
    const PylithInt sOffMean_x[2] = { sOff_x[i_disp], sOff_x[i_pres] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_bulkModulus] };

    const PylithInt numSDev = 1; // Number passed on to deviatoric stress kernel.
    const PylithInt sOffDev[1] = { sOff[i_disp] };
    const PylithInt sOffDev_x[1] = { sOff_x[i_disp] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_rstress], aOff[i_rstrain], aOff[i_shearModulus] };

    PylithScalar stressTensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    meanStress(_dim, numSMean, numAMean, sOffMean, sOffMean_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    IsotropicLinearElasticity3D::deviatoricStress(_dim, numSDev, numADev, sOffDev, sOffDev_x, s, s_t, s_x,
                                                  aOffDev, NULL, a, a_t, NULL,
                                                  t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        f1[i] -= stressTensor[i];
    } // for
} // f1u_refstate


// ---------------------------------------------------------------------------------------------------------------------
// Jf3_uu entry function for 3-D isotropic linear incompressible elasticity.
void
pylith::fekernels::IsotropicLinearIncompElasticity3D::Jf3uu(const PylithInt dim,
                                                            const PylithInt numS,
                                                            const PylithInt numA,
                                                            const PylithInt sOff[],
                                                            const PylithInt sOff_x[],
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


// ---------------------------------------------------------------------------------------------------------------------
// Calculate stress for 3-D  isotropic linear incompressible elasticity WITHOUT a reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticity3D::cauchyStress(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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
    const PylithInt i_pres = 1;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;

    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff[i_pres] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(sOff_x[i_pres] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(stressVector);

    const PylithInt numSDev = 1; // Number passed on to deviatoric stress kernel.
    const PylithInt sOffDev[1] = { sOff[i_disp] };
    const PylithInt sOffDev_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    meanStress(_dim, numS, numAMean, sOff, sOff_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
               t, x, numConstants, constants, stressTensor);
    IsotropicLinearElasticity3D::deviatoricStress(_dim, numSDev, numADev, sOffDev, sOffDev_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                                                  t, x, numConstants, constants, stressTensor);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stressTensor[2*_dim+2]; // stress_zz
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
    stressVector[4] = stressTensor[1*_dim+2]; // stress_yz
    stressVector[5] = stressTensor[0*_dim+2]; // stress_xz
} // stress


// ---------------------------------------------------------------------------------------------------------------------
// Calculate stress for 3-D isotropic linear incompressible elasticity WITH a reference stress/strain.
void
pylith::fekernels::IsotropicLinearIncompElasticity3D::cauchyStress_refstate(const PylithInt dim,
                                                                            const PylithInt numS,
                                                                            const PylithInt numA,
                                                                            const PylithInt sOff[],
                                                                            const PylithInt sOff_x[],
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
    const PylithInt i_pres = 1;

    // Incoming auxiliary fields.
    const PylithInt i_rstress = numA-4;
    const PylithInt i_rstrain = numA-3;
    const PylithInt i_shearModulus = numA-2;
    const PylithInt i_bulkModulus = numA-1;

    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 4);
    assert(sOff);
    assert(sOff[i_disp] >= 0);
    assert(sOff[i_pres] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(sOff_x[i_pres] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_bulkModulus] >= 0);
    assert(aOff[i_rstress] >= 0);
    assert(aOff[i_rstrain] >= 0);
    assert(stressVector);

    const PylithInt numSDev = 1; // Number passed on to deviratoric stress kernels.
    const PylithInt sOffDev[1] = { sOff[i_disp] };
    const PylithInt sOffDev_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_bulkModulus], aOff[i_rstress], aOff[i_rstrain] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    meanStress_refstate(_dim, numS, numAMean, sOff, sOff_x, s, s_t, s_x, aOffMean, NULL, a, a_t, NULL,
                        t, x, numConstants, constants, stressTensor);
    IsotropicLinearElasticity3D::deviatoricStress_refstate(_dim, numSDev, numADev, sOffDev, sOffDev_x, s, s_t, s_x, aOffDev, NULL, a, a_t, NULL,
                                                           t, x, numConstants, constants, stressTensor);
    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stressTensor[2*_dim+2]; // stress_zz
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
    stressVector[4] = stressTensor[1*_dim+2]; // stress_yz
    stressVector[5] = stressTensor[0*_dim+2]; // stress_xz
} // stress_refstate


// ---------------------------------------------------------------------------------------------------------------------
// Calculate mean stress for 3-D isotropic, incompressible linear elasticity WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticity3D::meanStress(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
                                                                 PylithScalar stressTensor[]) {
    const PylithInt _dim = 3;

    const PylithInt i_pres = 1;

    assert(_dim == dim);
    assert(numS >= 2);
    assert(sOff);
    assert(sOff[i_pres] >= 0);
    assert(s);
    assert(stressTensor);

    const PylithScalar pressure = s[sOff[i_pres]];

    for (PylithInt i = 0; i < _dim; ++i) {
        stressTensor[i*_dim+i] -= pressure;
    } // for
} // meanStress


// ---------------------------------------------------------------------------------------------------------------------
// Calculate mean stress for 3-D isotropic, incompressible linear elasticity WITH reference stress and reference strain.
void
pylith::fekernels::IsotropicLinearIncompElasticity3D::meanStress_refstate(const PylithInt dim,
                                                                          const PylithInt numS,
                                                                          const PylithInt numA,
                                                                          const PylithInt sOff[],
                                                                          const PylithInt sOff_x[],
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
                                                                          PylithScalar stressTensor[]) {
    const PylithInt _dim = 2;

    // Incoming solution subfields.
    const PylithInt i_pres = 1;

    // Incoming auxiliary subfields.
    const PylithInt i_rstress = numA-4;

    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 4);
    assert(sOff);
    assert(sOff[i_pres] >= 0);
    assert(aOff);
    assert(aOff[i_rstress] >= 0);
    assert(s);
    assert(a);
    assert(stressTensor);

    const PylithScalar pressure = s[sOff[i_pres]];
    const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, ...

    const PylithScalar meanRefStress = (refStress[0] + refStress[1] + refStress[3]) / 3.0;
    const PylithScalar meanStress = meanRefStress - pressure;

    for (PylithInt i = 0; i < _dim; ++i) {
        stressTensor[i*_dim+i] += meanStress;
    } // for
} // meanStress_refstate


// End of file
