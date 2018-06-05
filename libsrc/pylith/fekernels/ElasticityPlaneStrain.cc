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
 * Copyright (c) 2010-2017 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/ElasticityPlaneStrain.hh"
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for elatsicity plane strain.
 * ======================================================================
 */

// ----------------------------------------------------------------------
// f0 function for isotropic linear elasticity plane strain.
void
pylith::fekernels::ElasticityPlaneStrain::f0v(const PylithInt dim,
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
    const PylithInt _dim = 2;

    /* Incoming auxiliary fields. */
    const PylithInt i_density = 0;

    const PylithInt _numS = 2; // Number passed on to f0_inertia.

    const PylithInt _numA = 1; // Number passed on to f0_inertia.
    const PylithInt aOffInertia[1] = { aOff[i_density] };
    const PylithInt aOffInertia_x[1] = { aOff_x[i_density] };

    assert(_dim == dim);
    assert(2 == numS);
    assert(3 >= numA && 5 <= numA);
    assert(aOff);
    assert(aOff_x);

    pylith::fekernels::Elasticity::f0v_inertia(_dim, _numS, _numA,
                                               sOff, sOff_x, s, s_t, s_x,
                                               aOffInertia, aOffInertia_x, a, a_t, a_x,
                                               t, x, numConstants, constants, f0);
} // f0v


// ----------------------------------------------------------------------
// Jf0 function for isotropic linear elasticity plane strain.
void
pylith::fekernels::ElasticityPlaneStrain::Jf0vv(const PylithInt dim,
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
                                                PylithScalar Jf0[]) {
    const PylithInt _dim = 2;

    const PylithInt i_density = 0;
    const PylithScalar density = a[aOff[i_density]];

    PylithInt i;

    assert(_dim == dim);
    assert(2 == numS);
    assert(numA >= 5);
    assert(aOff);
    assert(a);

    for (i = 0; i < _dim; ++i) {
        Jf0[i*_dim+i] += utshift * density;
    } // for
} // Jf0vv


// ----------------------------------------------------------------------
/* Calculate mean stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT reference stress and strain.
 *
 * meanStress = bulkModulus * strain_kk
 *
 * stress += meanStress * delta_ij
 */
void
pylith::fekernels::ElasticityPlaneStrain::meanStress(const PylithInt dim,
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

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    // Incoming auxiliary field.
    const PylithInt i_bulkModulus = 0;
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    PylithInt i;

    assert(_dim == dim);
    assert(1 == numS);
    assert(1 == numA);
    assert(sOff_x);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal meanStress = bulkModulus * strainTrace;

    for (i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += meanStress;
    } // for
} // meanStress


// ----------------------------------------------------------------------
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
pylith::fekernels::ElasticityPlaneStrain::meanStress_refstate(const PylithInt dim,
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

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    // Incoming auxiliary fields.
    const PylithInt i_bulkModulus = 0;
    const PylithInt i_rstress = 1;
    const PylithInt i_rstrain = 2;
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, streass_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    PylithInt i;

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff_x);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[2];

    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithReal meanStress = meanrstress + bulkModulus * (strainTrace - refstrainTrace);

    for (i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += meanStress;
    } // for
} // meanStress_refstate


// ----------------------------------------------------------------------
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
pylith::fekernels::ElasticityPlaneStrain::deviatoricStress(const PylithInt dim,
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

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    // Incoming auxiliary field.
    const PylithInt i_shearModulus = 0;
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    assert(_dim == dim);
    assert(1 == numS);
    assert(1 == numA);
    assert(sOff_x);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal traceTerm = -2.0/3.0*shearModulus * strainTrace;
    const PylithReal twomu = 2.0*shearModulus;

    const PylithScalar stress_xx = twomu*disp_x[0*_dim+0] + traceTerm;
    const PylithScalar stress_yy = twomu*disp_x[1*_dim+1] + traceTerm;
    const PylithScalar stress_xy = shearModulus * (disp_x[0*_dim+1] + disp_x[1*_dim+0]);

    stress[0*_dim+0] += stress_xx;
    stress[1*_dim+1] += stress_yy;
    stress[0*_dim+1] += stress_xy;
    stress[1*_dim+0] += stress_xy;
} // deviatoricStress


// ----------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITH reference stress and reference strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * i==j
 * devStress_ii = refstress_ii - meanRefstress + 2*shearModulus*(strain_ii - refstrain_ii) - 2/3*shearModulus*(strain_kk - refstrain_kk)
 *
 * i!=j
 * devStress_ij = refstress_ij + 2*shearModulus*(strain_ij - refstrain_ij)
 */
void
pylith::fekernels::ElasticityPlaneStrain::deviatoricStress_refstate(const PylithInt dim,
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

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_rstress = 1;
    const PylithInt i_rstrain = 2;
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff_x);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[2];
    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithReal traceTerm = -2.0/3.0*shearModulus * (strainTrace - refstrainTrace);
    const PylithReal twomu = 2.0*shearModulus;

    const PylithScalar stress_xx = refstress[0] - meanrstress + twomu*(disp_x[0*_dim+0]-refstrain[0]) + traceTerm;
    const PylithScalar stress_yy = refstress[1] - meanrstress + twomu*(disp_x[1*_dim+1]-refstrain[1]) + traceTerm;
    const PylithScalar stress_xy = refstress[3] + twomu * (0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]) - refstrain[3]);

    stress[0*_dim+0] += stress_xx;
    stress[1*_dim+1] += stress_yy;
    stress[0*_dim+1] += stress_xy;
    stress[1*_dim+0] += stress_xy;

} // deviatoricStress_refstate


// ----------------------------------------------------------------------
/* Calculate strain for 2-D plane strain isotropic linear
 * elasticity.
 *
 * Used to output the strain field.
 *
 * Solution fields: [disp(dim)]
 */
void
pylith::fekernels::ElasticityPlaneStrain::strain(const PylithInt dim,
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
                                                 PylithScalar strain[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    assert(_dim == dim);
    assert(numS >= 1);
    assert(sOff_x);
    assert(s_x);
    assert(strain);

    const PylithScalar strain_xx = disp_x[0*_dim+0];
    const PylithScalar strain_yy = disp_x[1*_dim+1];
    const PylithScalar strain_zz = 0.0;
    const PylithScalar strain_xy = 0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]);

    strain[0] = strain_xx;
    strain[1] = strain_yy;
    strain[2] = strain_zz;
    strain[3] = strain_xy;
} // strain


// End of file
