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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_f0v(const PylithInt dim,
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

    /* Incoming auxiliary fields. */
    const PylithInt i_density = 0;

    const PylithInt _numS = 2; /* Number passed on to f0_inertia. */

    const PylithInt _numA = 1; /* Number passed on to f0_inertia. */
    const PylithInt aOffInertia[1] = { aOff[i_density] };
    const PylithInt aOffInertia_x[1] = { aOff_x[i_density] };

    assert(_dim == dim);
    assert(2 == numS);
    assert(3 == numA || 4 == numA);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_Elasticity_f0v_inertia(_dim, _numS, _numA,
                                            sOff, sOff_x, s, s_t, s_x,
                                            aOffInertia, aOffInertia_x, a, a_t, a_x,
                                            t, x, f0);
} /* IsotropicLinearElasticityPlaneStrain_f0v */


/* ---------------------------------------------------------------------- */
/** g0 function for isotropic linear elasticity plane strain with both gravity and body forces.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g0v_gravbodyforce(const PylithInt dim,
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
{ /* IsotropicLinearElasticityPlaneStrain_g0v_gravbodyforce */
    const PylithInt _dim = 2;

    /* Incoming auxiliary fields. */
    const PylithInt i_density = 0;
    const PylithInt i_gravityField = 3;
    const PylithInt i_bodyForce = 4;

    const PylithInt _numS = 0; /* Number passed on to g0_bodyforce. */

    const PylithInt numAGrav = 2; /* Number passed on to g0_grav. */
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    const PylithInt numABody = 1; /* Number passed on to g0_bodyforce. */
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 5);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_Elasticity_g0v_grav(_dim, _numS, numAGrav,
                                         NULL, NULL, NULL, NULL, NULL,
                                         aOffGrav, aOffGrav_x, a, a_t, a_x,
                                         t, x, g0);
    pylith_fekernels_Elasticity_g0v_bodyforce(_dim, _numS, numABody,
                                              NULL, NULL, NULL, NULL, NULL,
                                              aOffBody, aOffBody_x, a, a_t, a_x,
                                              t, x, g0);
} /* IsotropicLinearElasticityPlaneStrain_g0v_gravbodyforce */


/* ---------------------------------------------------------------------- */
/** g0 function for isotropic linear elasticity plane strain with both gravity and body forces.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g0v_grav(const PylithInt dim,
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
{ /* IsotropicLinearElasticityPlaneStrain_g0v_grav */
    const PylithInt _dim = 2;

    const PylithInt _numS = 0; /* Number passed on to g0_bodyforce. */

    /* Incoming auxiliary fields. */
    const PylithInt i_density = 0;
    const PylithInt i_gravityField = 3;

    const PylithInt numAGrav = 2; /* Number passed on to g0_grav. */
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 4);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_Elasticity_g0v_grav(_dim, _numS, numAGrav,
                                         NULL, NULL, NULL, NULL, NULL,
                                         aOffGrav, aOffGrav_x, a, a_t, a_x,
                                         t, x, g0);
} /* IsotropicLinearElasticityPlaneStrain_g0v_grav */


/* ---------------------------------------------------------------------- */
/** g0 function for isotropic linear elasticity plane strain with both gravity and body forces.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g0v_bodyforce(const PylithInt dim,
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
{ /* IsotropicLinearElasticityPlaneStrain_g0v_bodyforce */
    const PylithInt _dim = 2;

    const PylithInt _numS = 0; /* Number passed on to g0_bodyforce. */

    /* Incoming auxiliary fields. */
    const PylithInt i_bodyForce = 3;

    const PylithInt numABody = 1; /* Number passed on to g0_bodyforce. */
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 4);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_Elasticity_g0v_bodyforce(_dim, _numS, numABody,
                                              NULL, NULL, NULL, NULL, NULL,
                                              aOffBody, aOffBody_x, a, a_t, a_x,
                                              t, x, g0);
} /* IsotropicLinearElasticityPlaneStrain_g0v_bodyforce */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear elasticity plane strain WITHOUT reference stress and strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1v(const PylithInt dim,
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

    /* Incoming solution fields. */
    const PylithInt i_disp = 0;

    /* Incoming auxiliary fields. */
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;

    const PylithInt _numS = 1; /* Number passed on to stress kernels. */
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; /* Number passed to mean stress kernel. */
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };
    const PylithInt aOffMean_x[1] = { aOff_x[i_bulkModulus] };

    const PylithInt numADev = 1; /* Number passed to deviatoric stress kernel. */
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };
    const PylithInt aOffDev_x[1] = { aOff_x[i_shearModulus] };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0}; /* Full stress tensor */
    PylithInt i;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(3 == numA || 4 == numA);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress(_dim, _numS, numAMean,
                                                                     sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                     aOffMean, aOffMean_x, a, a_t, a_x,
                                                                     t, x, stress);

    pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress(_dim, _numS, numADev,
                                                                           sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                           aOffDev, aOffDev_x, a, a_t, a_x,
                                                                           t, x, stress);

    for (i=0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } /* for */
} /* IsotropicLinearElasticityPlaneStrain_g1v */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear elasticity plane strain with reference stress and strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1v_refstate(const PylithInt dim,
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

    /* Incoming solution fields. */
    const PylithInt i_disp = 0;

    /* Incoming auxiliary fields. */
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;

    const PylithInt _numS = 1; /* Number passed on to stress kernels. */
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; /* Number passed to mean stress kernel. */
    const PylithInt aOffMean[3] = { aOff[i_bulkModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[3] = { aOff_x[i_bulkModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 1; /* Number passed to deviatoric stress kernel. */
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

    pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress_refstate(_dim, _numS, numAMean,
                                                                              sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                              aOffMean, aOffMean_x, a, a_t, a_x,
                                                                              t, x, stress);

    pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress_refstate(_dim, _numS, numADev,
                                                                                    sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                                    aOffDev, aOffDev_x, a, a_t, a_x,
                                                                                    t, x, stress);

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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0vv_implicit(const PylithInt dim,
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0vv_explicit(const PylithInt dim,
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jg3vu(const PylithInt dim,
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

    /* Incoming auxiliary fields. */
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress(const PylithInt dim,
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

    /* Incoming solution field. */
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    /* Incoming auxiliary field. */
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress_refstate(const PylithInt dim,
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

    /* Incoming solution field. */
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    /* Incoming auxiliary fields. */
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
    assert(sOff_x);
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress(const PylithInt dim,
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

    /* Incoming solution field. */
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    /* Incoming auxiliary field. */
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress_refstate(const PylithInt dim,
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

    /* Incoming solution field. */
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    /* Incoming auxiliary fields. */
    const PylithInt i_shearModulus = 0;
    const PylithInt i_rstress = 1;
    const PylithInt i_rstrain = 2;
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; /* sigma_11, sigma_22, sigma_12, sigma_33 */
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; /* epsilon_11, epsilon_22, epsilon_12, epsilon_33 */

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

    const PylithScalar sigma_11 = refstress[0] - meanrstress + twomu*(disp_x[0*_dim+0]-refstrain[0]) + traceTerm;
    const PylithScalar sigma_22 = refstress[1] - meanrstress + twomu*(disp_x[1*_dim+1]-refstrain[1]) + traceTerm;
    const PylithScalar sigma_12 = refstress[2] + twomu * (0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]) - refstrain[2]);

    stress[0*_dim+0] += sigma_11;
    stress[1*_dim+1] += sigma_22;
    stress[0*_dim+1] += sigma_12;
    stress[1*_dim+0] += sigma_12;

} /* deviatoricStress_IsotropicLinearElasticityPlaneStrain_refstate */


/* End of file */
