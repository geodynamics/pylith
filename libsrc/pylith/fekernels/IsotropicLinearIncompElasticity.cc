/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University of Chicago
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2015 University of California, Davis
 *
 * See COPYING for license information.
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
// g1 function for isotropic linear incompressible elasticity plane strain WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticity::g0p(const PylithInt dim,
                                                        const PylithInt numS,
                                                        const PylithInt numA,
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
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
    assert(g0);

    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar pressure = s[sOff[i_pres]];

    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    PylithScalar strainTrace = 0;
    for (PylithInt i = 0; i < dim; ++i) {
        strainTrace += disp_x[i*dim+i];
    } // for
    g0[0] += strainTrace + pressure / bulkModulus;
} // g0p


// ---------------------------------------------------------------------------------------------------------------------
// g1 function for isotropic linear incompressible elasticity plane strain WITH reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticity::g0p_refstate(const PylithInt dim,
                                                                 const PylithInt numS,
                                                                 const PylithInt numA,
                                                                 const PylithInt sOff[],
                                                                 const PylithInt sOff_x[],
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
    assert(g0);

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
    g0[0] += (strainTrace - refStrainTrace) + (pressure + meanRefStress) / bulkModulus;
} // g0p


// ---------------------------------------------------------------------------------------------------------------------
// Jg0_pp function for isotropic linear incompressible elasticity.
void
pylith::fekernels::IsotropicLinearIncompElasticity::Jg0pp(const PylithInt dim,
                                                          const PylithInt numS,
                                                          const PylithInt numA,
                                                          const PylithInt sOff[],
                                                          const PylithInt sOff_x[],
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
                                                          PylithScalar Jg0[]) {
    // Incoming auxiliary subfields
    const PylithInt i_bulkModulus = numA-1;

    assert(numA >= 1);
    assert(aOff);
    assert(aOff[i_bulkModulus] >= 0);
    assert(a);
    assert(Jg0);

    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    Jg0[0] += 1.0 / bulkModulus;
} // Jg0pp


// =====================================================================================================================
// Kernels for 2-D plane strain isotropic, linear incompressible elasticity.
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
// g1 function for isotropic linear incompressible elasticity plane strain WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::g1u(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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
    assert(g1);

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
        g1[i] -= stressTensor[i];
    } // for
} // g1u


// ---------------------------------------------------------------------------------------------------------------------
// g1 function for isotropic linear incompressible elasticity plane strain WITH reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::g1u_refstate(const PylithInt dim,
                                                                            const PylithInt numS,
                                                                            const PylithInt numA,
                                                                            const PylithInt sOff[],
                                                                            const PylithInt sOff_x[],
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
    assert(g1);

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
        g1[i] -= stressTensor[i];
    } // for
} // g1u_refstate


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
// g1 function for 3-D isotropic linear incompressible elasticity WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticity3D::g1u(const PylithInt dim,
                                                          const PylithInt numS,
                                                          const PylithInt numA,
                                                          const PylithInt sOff[],
                                                          const PylithInt sOff_x[],
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
    assert(g1);

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
    IsotropicLinearElasticityPlaneStrain::deviatoricStress(_dim, numSDev, numADev, sOffDev, sOffDev_x, s, s_t, s_x,
                                                           aOffDev, NULL, a, a_t, NULL,
                                                           t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stressTensor[i];
    } // for
} // g1u


// ---------------------------------------------------------------------------------------------------------------------
// g1 function for 3-D isotropic linear incompressible elasticity WITH reference stress and strain.
void
pylith::fekernels::IsotropicLinearIncompElasticity3D::g1u_refstate(const PylithInt dim,
                                                                   const PylithInt numS,
                                                                   const PylithInt numA,
                                                                   const PylithInt sOff[],
                                                                   const PylithInt sOff_x[],
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
    assert(g1);

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
    IsotropicLinearElasticityPlaneStrain::deviatoricStress(_dim, numSDev, numADev, sOffDev, sOffDev_x, s, s_t, s_x,
                                                           aOffDev, NULL, a, a_t, NULL,
                                                           t, x, numConstants, constants, stressTensor);
    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stressTensor[i];
    } // for
} // g1u_refstate


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
