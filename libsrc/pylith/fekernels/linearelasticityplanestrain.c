/* -*- C -*-
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

#include "pylith/fekernels/linearelasticityplanestrain.h"

#include "pylith/fekernels/elasticity.h" /* USES Elasticity_f0_inertia, Elasticity_g0_bodyforce */

/* ======================================================================
 * Kernels for isotropic, linear elatsicity plane strain.
 * ======================================================================
 */


/* ---------------------------------------------------------------------- */
/** f0 function for isotropic linear elasticity plane strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_f0v(
						const PylithInt dim,
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
						PylithScalar f0[])
{ /* IsotropicLinearElasticityPlaneStrain_f0v */
    const PylithInt _dim = 2;

    const PylithInt _numS = 2; /* Number passed on to f0_inertia. */

    const PylithInt _numA = 1; /* Number passed on to f0_inertia. */
    const PylithInt i_density = 0;

    assert(_dim == dim);
    assert(2 == numS);
    assert(3 == numA || 4 == numA);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_Elasticity_f0v_inertia(
	_dim, _numS, _numA, sOff, sOff_x, s, s_t, s_x, &aOff[i_density],
	&aOff_x[i_density], a, a_t, a_x, t, x, f0);
} /* IsotropicLinearElasticityPlaneStrain_f0v */


/* ---------------------------------------------------------------------- */
/** g0 function for isotropic linear elasticity plane strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g0v(
						const PylithInt dim,
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
						PylithScalar g0[])
{ /* IsotropicLinearElasticityPlaneStrain_g0v */
    const PylithInt _dim = 2;

    const PylithInt _numS = 2; /* Number passed on to g0_bodyforce. */

    const PylithInt _numA = 1; /* Number passed on to g0_bodyforce. */
    const PylithInt i_bodyforce = 3;

    assert(_dim == dim);
    assert(2 == numS);
    assert(4 == numA);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_Elasticity_g0v_bodyforce(
	_dim, _numS, _numA, sOff, sOff_x, s, s_t, s_x, &aOff[i_bodyforce],
	&aOff_x[i_bodyforce], a, a_t, a_x, t, x, g0);
} /* IsotropicLinearElasticityPlaneStrain_g0v */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear elasticity plane strain WITHOUT reference stress and strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1v(
						const PylithInt dim,
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
						PylithScalar g1[])
{ /* IsotropicLinearElasticityPlaneStrain_g1v */
    const PylithInt _dim = 2;

    const PylithInt _numS = 1; /* Pass displacement field. */
    const PylithInt i_disp = 0;

    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;

    const PylithInt numAMean = 1; /* Pass bulk modulus. */
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };
    const PylithInt aOffMean_x[1] = { aOff_x[i_bulkModulus] };

    const PylithInt numADev = 1; /* Pass shear modulus. */
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };
    const PylithInt aOffDev_x[1] = { aOff_x[i_shearModulus] };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0};
    PylithInt i;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(3 == numA || 4 == numA);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress(
	_dim, _numS, numAMean, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x,
	aOffMean, aOffMean_x, a, a_t, a_x, t, x, stress);
    pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress(
	_dim, _numS, numADev, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x,
	aOffDev, aOffDev_x, a, a_t, a_x, t, x, stress);

    for (i=0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } /* for */
} /* IsotropicLinearElasticityPlaneStrain_g1v */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear elasticity plane strain with reference stress and strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1v_refstate(
						const PylithInt dim,
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
						PylithScalar g1[])
{ /* IsotropicLinearElasticityPlaneStrain_g1v_refstate */
    const PylithInt _dim = 2;

    const PylithInt _numS = 1; /* Pass displacement field. */
    const PylithInt i_disp = 0;

    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;

    const PylithInt numAMean = 3; /* Pass shear modulus, reference stress, and reference strain. */
    const PylithInt aOffMean[3] = { aOff[i_bulkModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[3] = { aOff_x[i_bulkModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 3; /* Pass bulk modulus, reference stress, and reference strain. */
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0};
    PylithInt i;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(5 == numA || 6 == numA);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress_refstate(
	_dim, _numS, numAMean, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x,
	aOffMean, aOffMean_x, a, a_t, a_x, t, x, stress);
    pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress_refstate(
	_dim, _numS, numADev, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x,
	aOffDev, aOffDev_x, a, a_t, a_x, t, x, stress);

    for (i=0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } /* for */
} /* IsotropicLinearElasticityPlaneStrain_g1v_refstate */


/* ---------------------------------------------------------------------- */
/** Jf0 function for isotropic linear elasticity plane strain with implicit time stepping.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [2].
 * @param numA Number of registered subfields in auxiliary field [0].
 * @param sOff Offset of registered subfields in solution field [numS].
 * @param sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param a Auxiliary field with all subfields.
 * @param a_t Time derivative of auxiliary field.
 * @param a_x Gradient of auxiliary field.
 * @param t Time for residual evaluation.
 * @param utshift Coefficient for dF/ds_t term in Jacobian.
 * @param x Coordinates of point evaluation.
 * @param f0 Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0vv_implicit(
						const PylithInt dim,
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
						PylithScalar Jf0[])
{ /* IsotropicLinearElasticityPlaneStrain_Jf0vv_implicit */
    const PylithInt _dim = 2;

    const PylithInt i_density = 0;
    const PylithScalar density = a[aOff[i_density]];
    PylithInt i;

    assert(_dim == dim);
    assert(2 == numS);
    assert(3 <= numA);
    assert(aOff);
    assert(a);

    for (i = 0; i < _dim; ++i) {
        Jf0[i*_dim+i] += utshift * density;
    } /* for */
} /* IsotropicLinearElasticityPlaneStrain_Jf0vv_implicit */


/* ---------------------------------------------------------------------- */
/** Jf0 function for isotropic linear elasticity plane strain with explicit time stepping.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [2].
 * @param numA Number of registered subfields in auxiliary field [0].
 * @param sOff Offset of registered subfields in solution field [numS].
 * @param sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param a Auxiliary field with all subfields.
 * @param a_t Time derivative of auxiliary field.
 * @param a_x Gradient of auxiliary field.
 * @param t Time for residual evaluation.
 * @param utshift Coefficient for dF/ds_t term in Jacobian.
 * @param x Coordinates of point evaluation.
 * @param f0 Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0vv_explicit(
						const PylithInt dim,
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
						PylithScalar Jf0[])
{ /* IsotropicLinearElasticityPlaneStrain_Jf0vv_explicit */
    const PylithInt _dim = 2;

    const PylithInt i_density = 0;
    const PylithScalar density = a[aOff[i_density]];
    PylithInt i;

    assert(_dim == dim);
    assert(2 == numS);
    assert(3 <= numA);
    assert(aOff);
    assert(a);

    for (i = 0; i < _dim; ++i) {
        Jf0[i*_dim+i] += density;
    } /* for */

} /* IsotropicLinearElasticityPlaneStrain_Jf0vv_explicit*/


/* ---------------------------------------------------------------------- */
/* Jg3_vu entry function for 2-D plane strain isotropic linear elasticity.
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
 *  C_ijkl = bulkModulus * delta_ij * delta_kl + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jg3vu(
						const PylithInt dim,
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
						PylithScalar Jg3[])
{ /* IsotropicLinearElasticityPlaneStrain_Jg3vu */
    const PylithInt _dim = 2;

    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar lambda2mu = lambda + 2.0*shearModulus;

    const PylithReal C1111 = lambda2mu;
    const PylithReal C2222 = lambda2mu;
    const PylithReal C1122 = lambda;
    const PylithReal C1212 = shearModulus;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(3 <= numA);
    assert(aOff);
    assert(a);
    assert(Jg3);

    /* j(f,g,df,dg) = C(f,df,g,dg)

       0: j0000 = C1111
       1: j0001 = C1112 = 0
       4: j0100 = C1121, symmetry C1112 = 0
       5: j0101 = C1122

       2: j0010 = C1211 = 0
       3: j0011 = C1212
       6: j0110 = C1221, symmetry C1212
       7: j0111 = C1222 = 0

       8: j1000 = C2111 = 0
       9: j1001 = C2112, symmetry C1212
       12: j1100 = C2121, symmetry C1212
       13: j1101 = C2122, symmetry C1222 = 0

       10: j1010 = C2211, symmetry C1122
       11: j1011 = C2212, symmetry C1222 = 0
       14: j1110 = C2221, symmetry C1222 = 0
       15: j1111 = C2222
     */

    Jg3[ 0] -= C1111; /* j0000 */
    Jg3[ 3] -= C1212; /* j0011 */
    Jg3[ 5] -= C1122; /* j0101 */
    Jg3[ 6] -= C1212; /* j0110, C1221 */
    Jg3[ 9] -= C1212; /* j1001, C2112 */
    Jg3[10] -= C1122; /* j1010, C2211 */
    Jg3[12] -= C1212; /* j1100, C2121 */
    Jg3[15] -= C2222; /* j1111 */

} /* IsotropicLinearElasticityPlaneStrain_Jg3vu */


/* ---------------------------------------------------------------------- */
/* Calculate mean stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT reference stress and strain.
 *
 * meanStress = bulkModulus * strain_kk
 *
 * stress += meanStress * delta_ij
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress(
						const PylithInt dim,
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
						PylithScalar stress[])
{ /* IsotropicLinearElasticityPlaneStrain_meanStress */
    const PylithInt _dim = 2;

    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    const PylithInt i_bulkModulus = 0;
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    PylithInt i;

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

    for (i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += meanStress;
    } /* for */
} /* IsotropicLinearElasticityPlaneStrain_meanStress */


/* ---------------------------------------------------------------------- */
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress_refstate(
						const PylithInt dim,
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
						PylithScalar stress[])
{ /* IsotropicLinearElasticityPlaneStrain_meanStress_refstate */
    const PylithInt _dim = 2;

    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    const PylithInt i_bulkModulus = 0;
    const PylithInt i_rstress = 1;
    const PylithInt i_rstrain = 2;
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; /* sigma_11, sigma_22, sigma_12, sigma_33 */
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; /* epsilon_11, epsilon_22, epsilon_12, epsilon_33 */

    PylithInt i;

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

    for (i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += meanStress;
    } /* for */
} /* IsotropicLinearElasticityPlaneStrain_meanStress_refstate */


/* ---------------------------------------------------------------------- */
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress(
						const PylithInt dim,
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
						PylithScalar stress[])
{ /* IsotropicLinearElasticityPlaneStrain_deviatoricStress */
    const PylithInt _dim = 2;

    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    const PylithInt i_shearModulus = 0;
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    assert(_dim == dim);
    assert(1 == numS);
    assert(1 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal traceTerm = -2.0/3.0*shearModulus * strainTrace;
    const PylithReal twomu = 2.0*shearModulus;

    const PylithScalar sigma_11 = twomu*disp_x[0*_dim+0] + traceTerm;
    const PylithScalar sigma_22 = twomu*disp_x[1*_dim+1] + traceTerm;
    const PylithScalar sigma_12 = shearModulus * (disp_x[0*_dim+1] + disp_x[1*_dim+0]);

    stress[0*_dim+0] += sigma_11;
    stress[1*_dim+1] += sigma_22;
    stress[0*_dim+1] += sigma_12;
    stress[1*_dim+0] += sigma_12;
} /* IsotropicLinearElasticityPlaneStrain_deviatoricStress */


/* ---------------------------------------------------------------------- */
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress_refstate(
	                                   const PylithInt dim,
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
					   PylithScalar stress[])
{ /* deviatoricStress_IsotropicLinearElasticityPlaneStrain_refstate */
    const PylithInt _dim = 2;

    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    const PylithInt i_shearModulus = 0;
    const PylithInt i_rstress = 1;
    const PylithInt i_rstrain = 2;
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; /* sigma_11, sigma_22, sigma_12, sigma_33 */
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; /* epsilon_11, epsilon_22, epsilon_12, epsilon_33 */

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[2];
    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithReal traceTerm = -2.0/3.0*shearModulus * (strainTrace - refstrainTrace);
    const PylithReal twomu = 2.0*shearModulus;

    const PylithScalar sigma_11 = refstress[0] - meanrstress + twomu*(disp_x[0*_dim+0]-refstrain[0]) + traceTerm;
    const PylithScalar sigma_22 = refstress[1] - meanrstress + twomu*(disp_x[1*_dim+1]-refstrain[1]) + traceTerm;
    const PylithScalar sigma_12 = refstress[2] + twomu * (0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]) - refstrain[2]);

    stress[0*_dim+0] += sigma_11;
    stress[1*_dim+1] += sigma_22;
    stress[0*_dim+1] += sigma_12;
    stress[1*_dim+0] += sigma_12;

} /* deviatoricStress_IsotropicLinearElasticityPlaneStrain_refstate */


/* ======================================================================
 * Kernels for incompressible elasticity.
 *
 * Solution fields:  [disp(dim), vel(dim), pres(1)]
 * Auxiliary fields: [density(1), shearModulus(1), bulkModulus(1), bodyforce(dim), initialstrain(4), initialstress(4)]
 *
 * ======================================================================
 */

/* ---------------------------------------------------------------------- */
/** f0 function for isotropic linear incompressible elasticity plane strain.
 * Solution fields:  [disp(dim), vel(dim)]
 * Auxiliary fields: [density(1), ...]
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_f0v(
						const PylithInt dim,
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
						PylithScalar f0[])
{ /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_f0v */
    const PylithInt _dim = 2;

    const PylithInt _numS = 2; /* Passing on displacement and velocity. */
    const PylithInt i_disp = 0;
    const PylithInt i_vel = 1;

    const PylithInt _numA = 1; /* Passing on density. */
    const PylithInt i_density = 0;


    assert(_dim == dim);
    assert(3 == numS);
    assert(1 <= numA);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _sOff[2] = {sOff[i_disp], sOff[i_vel]};
    const PylithInt _sOff_x[2] = {sOff_x[i_disp], sOff_x[i_vel]};

    pylith_fekernels_Elasticity_f0v_inertia(
	_dim, _numS, _numA, _sOff, _sOff_x, s, s_t, s_x, &aOff[i_density],
	&aOff_x[i_density], a, a_t, a_x, t, x, f0);
} /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_f0v */


/* ---------------------------------------------------------------------- */
/** g0 function for isotropic linear incompressible elasticity plane strain.
 *
 * Solution fields:  [disp(dim), vel(dim), pres(1)]
 * Auxiliary fields: [density(1), shearModulus(1), bulkModulus(1), bodyforce(dim), ...]
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g0v(
						const PylithInt dim,
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
						PylithScalar g0[])
{ /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g0 */
    const PylithInt _dim = 2;

    const PylithInt _numS = 2; /* Passing on displacement and velocity. */
    const PylithInt i_disp = 0;
    const PylithInt i_vel = 1;

    const PylithInt _numA = 1; /* Passing on body force. */
    const PylithInt i_bodyforce = 3;

    assert(_dim == dim);
    assert(3 == numS);
    assert(4 == numA || 6 == numA);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _sOff[2] = {sOff[i_disp], sOff[i_vel]};
    const PylithInt _sOff_x[2] = {sOff_x[i_disp], sOff_x[i_vel]};

    pylith_fekernels_Elasticity_g0v_bodyforce(
	_dim, _numS, _numA, _sOff, _sOff_x, s, s_t, s_x, &aOff[i_bodyforce],
	&aOff_x[i_bodyforce], a, a_t, a_x, t, x, g0);
} /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g0 */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear incompressible elasticity plane
 * strain WITHOUT reference stress and strain.
 *
 * Solution fields: [disp(dim), vel(dim), pres(1)]
 * Auxiliary fields: [density(1), shearModulus(1), bulkModulus(1), ...]
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g1v(
						const PylithInt dim,
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
						PylithScalar g1[])
{ /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g1v */
    const PylithInt _dim = 2;

    const PylithInt _numS = 1; /* Passing displacement field. */
    const PylithInt i_disp = 0;
    const PylithInt i_pres = 2;
    const PylithScalar pres = s[sOff[i_pres]];

    const PylithInt i_shearModulus = 1;

    const PylithInt numADev = 1; /* Passing shear modulus */
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };
    const PylithInt aOffDev_x[1] = { aOff_x[i_shearModulus] };

    PylithInt i;
    PylithScalar stress[4] = {0.0,0.0, 0.0, 0.0};

    assert(_dim == dim);
    assert(3 == numS);
    assert(3 == numA || 4 == numA);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress(
	_dim, _numS, numADev, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x,
	aOffDev, aOffDev_x, a, a_t, a_x, t, x, stress);

    for (i=0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } /* for */
    for (i=0; i < _dim; ++i) {
        g1[i*_dim+i] += pres;
    } /* for */
} /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g1v */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear incompressible elasticity plane strain
 * with reference stress and strain.
 *
 * Solution fields:  [disp(dim), vel(dim), pres]
 * Auxiliary fields: [density(1), shearModulus(1), bulkModulus(1), ..., refstress(4), refstrain(4)]
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g1v_refstate(
						const PylithInt dim,
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
						PylithScalar g1[])
{ /* IsotropicLinearIncompElasticityPlaneStrain_g1v_refstate */
    const PylithInt _dim = 2;

    const PylithInt _numS = 1; /* Passing displacement field. */
    const PylithInt i_disp = 0;
    const PylithInt i_pres = 2;

    const PylithInt numADev = 3;
    const PylithInt i_shearModulus = 1;
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;
    const PylithScalar* refstress = &a[aOff[i_rstress]]; /* sigma_11, sigma_22, sigma_12, sigma_33 */

    const PylithScalar pres = s[sOff[i_pres]];

    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[3]) / 3.0;
    const PylithScalar meanstress = meanrstress - pres;

    PylithInt i;
    PylithScalar stress[4] = {0.0,0.0, 0.0, 0.0};

    assert(_dim == dim);
    assert(3 == numS);
    assert(5 == numA || 6 == numA);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress_refstate(
	_dim, _numS, numADev, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x,
	aOffDev, aOffDev_x, a, a_t, a_x, t, x, stress);

    for (i=0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } /* for */
    for (i=0; i < _dim; ++i) {
        g1[i*_dim+i] += pres - meanstress;
    } /* for */
} /* IsotropicLinearIncompElasticityPlaneStrain_g1v_refstate */


/* ---------------------------------------------------------------------- */
/* Jg3_vu entry function for 2-D plane strain incompressible isotropic linear
 * elasticity.
 *
 * Solution fields:  [disp(dim), vel(dim), pres]
 * Auxiliary fields: [density(1), shearModulus(1), bulkModulus[1], ...]
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg3vu(
						const PylithInt dim,
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
						PylithScalar Jg3[])
{ /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg3vu */
    const PylithInt _dim = 2;

    const PylithInt i_shearModulus = 1;
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithReal C1111 = 4.0 / 3.0 * shearModulus;
    const PylithReal C2222 =  C1111;
    const PylithReal C1122 = -2.0 / 3.0 * shearModulus;
    const PylithReal C1212 = shearModulus;

    assert(_dim == dim);
    assert(3 == numS);
    assert(3 <= numA);
    assert(aOff);

    /* g(f,g,df,dg) = C^{dev}(f,df,g,dg)

       0: g0000 = C1111
       1: g0001 = C1112
       4: g0100 = C1121, symmetry C1112
       5: g0101 = C1122

       2: g0010 = C1211
       3: g0011 = C1212
       6: g0110 = C1221, symmetry C1212
       7: g0111 = C1222

       8: g1000 = C2111
       9: g1001 = C2112, symmetry C1212
       12: g1100 = C2121, symmetry C1212
       13: g1101 = C2122, symmetry C1222

       10: g1010 = C2211, symmetry C1122
       11: g1011 = C2212, symmetry C1222
       14: g1110 = C2221, symmetry C1222
       15: g1111 = C2222
     */

    Jg3[ 0] -= C1111; /* g0000 */
    Jg3[ 3] -= C1212; /* g0011 */
    Jg3[ 5] -= C1122; /* g0101 */
    Jg3[ 6] -= C1212; /* g0110, C1221 */
    Jg3[ 9] -= C1212; /* g1001, C2112 */
    Jg3[10] -= C1122; /* g1010, C2211 */
    Jg3[12] -= C1212; /* g1100, C2121 */
    Jg3[15] -= C2222; /* g1111 */
} /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg3vu */


/* ---------------------------------------------------------------------- */
/* Jg1_pu entry function for 2-D plane strain incompressible isotropic linear
 * elasticity.
 *
 * Solution fields = [disp(dim), vel(dim), pres]
 * Auxiliary fields = [density, shearModulus, bulkModulus, ...]
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg1pu(
						const PylithInt dim,
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
						PylithScalar Jg1[])
{ /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg1pu */
    const PylithInt _dim = 2;

    PylithInt i;

    assert(2 == _dim);
    assert(3 == numS);
    assert(3 <= numA);

    for (i=0; i < _dim; ++i) {
        Jg1[i*_dim+i] += 1.0;
    } /* for */
} /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg1pu */


/* ---------------------------------------------------------------------- */
/* Jg2_vp entry function for 2-D plane strain incompressible isotropic linear
 * elasticity.
 *
 * Solution fields = [disp(dim), vel(dim), pres]
 * Auxiliary fields = [density, shearModulus, bulkModulus, ...]
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg2vp(
						const PylithInt dim,
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
						PylithScalar Jg2[])
{ /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg2vp */
    const PylithInt _dim = 2;

    PylithInt i;

    assert(2 == _dim);
    assert(3 == numS);
    assert(3 <= numA);

    for (i=0; i < _dim; ++i) {
        Jg2[i*_dim+i] += 1.0;
    } /* for */
} /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg2vp */


/* End of file */
