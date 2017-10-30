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

#include "pylith/fekernels/IsotropicLinearIncompElasticityPlaneStrain.hh"

#include "pylith/fekernels/ElasticityPlaneStrain.hh" // USES ElasticityPlaneStrain kernels

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for isotropic, incompressible plane strain elasticity.
 *
 * Solution fields:  [disp(dim), pressure(1)]
 * Auxiliary fields:
 * - 0: density(1)
 * - 1: shear_modulus(1)
 * - 2: bulk_modulus(1)
 * - 3: gravity_field (2, optional)
 * - 4: body_force(2,optional)
 * - 5: reference_stress(4,optional) (stress_xx, stress_yy, stress_xy, stress_zz)
 * - 6: reference_strain(4,optional) (strain_xx, strain_yy, strain_xy, strain_zz)
 *
 * ======================================================================
 */

// ----------------------------------------------------------------------
/* g1 function for isotropic linear incompressible elasticity plane
 * strain WITHOUT reference stress and strain.
 */
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };
    const PylithInt aOffMean_x[1] = { aOff_x[i_bulkModulus] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };
    const PylithInt aOffDev_x[1] = { aOff_x[i_shearModulus] };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0}; // Full stress tensor
    PylithInt i;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    meanStress(_dim, _numS, numAMean,
               sOffDisp, sOffDisp_x, s, s_t, s_x,
               aOffMean, aOffMean_x, a, a_t, a_x,
               t, x, numConstants, constants, stress);

    pylith::fekernels::ElasticityPlaneStrain::deviatoricStress(_dim, _numS, numADev,
                                                               sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                               aOffDev, aOffDev_x, a, a_t, a_x,
                                                               t, x, numConstants, constants, stress);

    for (i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1u

// ----------------------------------------------------------------------
/** g1 function for isotropic linear incompressible elasticity plane strain
 * with reference stress and strain.
 */
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_pres = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;

    const PylithInt _numS = 2; // Number passed on to stress kernels.
    const PylithInt sOffDisp[2] = { sOff[i_disp], sOff[i_pres] };
    const PylithInt sOffDisp_x[2] = { sOff_x[i_disp], sOff_x[i_pres] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_bulkModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[3] = { aOff_x[i_bulkModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0};
    PylithInt i;

    assert(_dim == dim);
    assert(2 == numS);
    assert(numA >= 5);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    meanStress_refstate(_dim, _numS, numAMean,
                        sOffDisp, sOffDisp_x, s, s_t, s_x,
                        aOffMean, aOffMean_x, a, a_t, a_x,
                        t, x, numConstants, constants, stress);

    pylith::fekernels::ElasticityPlaneStrain::deviatoricStress_refstate(_dim, _numS, numADev,
                                                                        sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                        aOffDev, aOffDev_x, a, a_t, a_x,
                                                                        t, x, numConstants, constants, stress);

    for (i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1u_refstate


// ----------------------------------------------------------------------
/* Jg1_pu entry function for 2-D plane strain incompressible isotropic linear
 * elasticity.
 *
 * Solution fields = [disp(dim), vel(dim), pres]
 * Auxiliary fields = [density, shearModulus, bulkModulus, ...]
 */
void
pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::Jg1pu(const PylithInt dim,
                                                                     const PylithInt numS,
                                                                     const PylithInt numA,
                                                                     const PylithInt sOff[],
                                                                     const PylithInt sOff_x[],
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
                                                                     PylithScalar Jg1[]) {
    const PylithInt _dim = 2;

    assert(_dim == dim);
    assert(2 == numS);
    assert(3 <= numA);

    for (PylithInt i = 0; i < _dim; ++i) {
        Jg1[i*_dim+i] += 1.0;
    } // for
} // Jg1pu


// ----------------------------------------------------------------------
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
                                                                          PylithScalar stress[]) {
    const PylithInt _dim = 2;

    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    const PylithInt i_bulkModulus = 0;
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    assert(_dim == dim);
    assert(1 == numS);
    assert(1 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal meanStress = bulkModulus * strainTrace;

    for (PylithInt i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += meanStress;
    } // for
} // meanStress


// ----------------------------------------------------------------------
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
                                                                                   PylithScalar stress[]) {
    const PylithInt _dim = 2;

    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    const PylithInt i_bulkModulus = 0;
    const PylithInt i_rstress = 1;
    const PylithInt i_rstrain = 2;
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // sigma_11, sigma_22, sigma_33, sigma_12
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // epsilon_11, epsilon_22, sigma_33, epsilon_12

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[3];

    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[3]) / 3.0;
    const PylithReal meanStress = meanrstress + bulkModulus * (strainTrace - refstrainTrace);

    for (PylithInt i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += meanStress;
    } // for
} // meanStress_refstate


// End of file
