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

#include "pylith/fekernels/linearmaxwellplanestrain.h"

#include "pylith/fekernels/elasticity.h" /* USES Elasticity_f0_inertia, Elasticity_g0_bodyforce */

/* ======================================================================
 * Kernels for isotropic, linear Maxwell viscoelastic plane strain.
 * ======================================================================
 */


/* ---------------------------------------------------------------------- */
/** f0 function for isotropic linear Maxwell viscoelastic plane strain.
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_f0v(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar f0[])
{ /* IsotropicLinearMaxwellPlaneStrain_f0v */
    const PylithInt _dim = 2;

    /* Incoming auxiliary fields. */
    const PylithInt i_density = 0;

    const PylithInt _numS = 2; /* Number passed on to f0_inertia. */

    const PylithInt _numA = 1; /* Number passed on to f0_inertia. */
    const PylithInt aOffInertia[1] = { aOff[i_density] };
    const PylithInt aOffInertia_x[1] = { aOff_x[i_density] };

    assert(_dim == dim);
    assert(2 == numS);
    assert(6 >= numA && 10 <= numA);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_Elasticity_f0v_inertia(_dim, _numS, _numA,
                                            sOff, sOff_x, s, s_t, s_x,
                                            aOffInertia, aOffInertia_x, a, a_t, a_x,
                                            t, x, numConstants, constants, f0);
} /* IsotropicLinearMaxwellPlaneStrain_f0v */


/* ---------------------------------------------------------------------- */
/** g0 function for isotropic linear Maxwell viscoelastic plane strain with both gravity and body forces.
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_g0v_gravbodyforce(const PylithInt dim,
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
                                                                     PylithScalar g0[])
{ /* IsotropicLinearMaxwellPlaneStrain_g0v_gravbodyforce */
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
    assert(numA >= 6);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_Elasticity_g0v_grav(_dim, _numS, numAGrav,
                                         NULL, NULL, NULL, NULL, NULL,
                                         aOffGrav, aOffGrav_x, a, a_t, a_x,
                                         t, x, numConstants, constants, g0);
    pylith_fekernels_Elasticity_g0v_bodyforce(_dim, _numS, numABody,
                                              NULL, NULL, NULL, NULL, NULL,
                                              aOffBody, aOffBody_x, a, a_t, a_x,
                                              t, x, numConstants, constants, g0);
} /* IsotropicLinearMaxwellPlaneStrain_g0v_gravbodyforce */


/* ---------------------------------------------------------------------- */
/** g0 function for isotropic linear Maxwell viscoelastic plane strain with just gravity.
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_g0v_grav(const PylithInt dim,
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
                                                            PylithScalar g0[])
{ /* IsotropicLinearMaxwellPlaneStrain_g0v_grav */
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
    assert(numA >= 6);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_Elasticity_g0v_grav(_dim, _numS, numAGrav,
                                         NULL, NULL, NULL, NULL, NULL,
                                         aOffGrav, aOffGrav_x, a, a_t, a_x,
                                         t, x, numConstants, constants, g0);
} /* IsotropicLinearMaxwellPlaneStrain_g0v_grav */


/* ---------------------------------------------------------------------- */
/** g0 function for isotropic linear Maxwell viscoelastic plane strain with only body forces.
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_g0v_bodyforce(const PylithInt dim,
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
                                                                 PylithScalar g0[])
{ /* IsotropicLinearMaxwellPlaneStrain_g0v_bodyforce */
    const PylithInt _dim = 2;

    const PylithInt _numS = 0; /* Number passed on to g0_bodyforce. */

    /* Incoming auxiliary fields. */
    const PylithInt i_bodyForce = 3;

    const PylithInt numABody = 1; /* Number passed on to g0_bodyforce. */
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 6);
    assert(aOff);
    assert(aOff_x);

    pylith_fekernels_Elasticity_g0v_bodyforce(_dim, _numS, numABody,
                                              NULL, NULL, NULL, NULL, NULL,
                                              aOffBody, aOffBody_x, a, a_t, a_x,
                                              t, x, numConstants, constants, g0);
} /* IsotropicLinearMaxwellPlaneStrain_g0v_bodyforce */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear Maxwell plane strain WITHOUT reference stress and strain.
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_g1v(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar g1[])
{ /* IsotropicLinearMaxwellPlaneStrain_g1v */
    const PylithInt _dim = 2;

    /* Incoming solution fields. */
    const PylithInt i_disp = 0;

    /* Incoming auxiliary fields. */
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_totalStrain = 4;
    const PylithInt i_viscousStrain = 5;

    const PylithInt _numS = 1; /* Number passed on to stress kernels. */
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; /* Number passed to mean stress kernel. */
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };
    const PylithInt aOffMean_x[1] = { aOff_x[i_bulkModulus] };

    const PylithInt numADev = 4; /* Number passed to deviatoric stress kernel. */
    const PylithInt aOffDev[4] = { aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_totalStrain],
                                   aOff[i_viscousStrain] };
    const PylithInt aOffDev_x[4] = { aOff_x[i_shearModulus], aOff_x[i_maxwellTime],
                                     aOff_x[i_totalStrain], aOff_x[i_viscousStrain] };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0}; /* Full stress tensor */
    PylithInt i;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 6);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);
    assert(numConstants == 1);

    pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_meanStress(_dim, _numS, numAMean,
                                                                  sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                  aOffMean, aOffMean_x, a, a_t, a_x,
                                                                  t, x, numConstants, constants, stress);
    pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_deviatoricStress(_dim, _numS, numADev,
                                                                        sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                        aOffDev, aOffDev_x, a, a_t, a_x,
                                                                        t, x, numConstants, constants, stress);

    for (i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } /* for */
} /* IsotropicLinearMaxwellPlaneStrain_g1v */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear Maxwell viscoelastic plane strain with reference stress and strain.
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_g1v_refstate(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar g1[])
{ /* IsotropicLinearMaxwellPlaneStrain_g1v_refstate */
    const PylithInt _dim = 2;

    /* Incoming solution fields. */
    const PylithInt i_disp = 0;

    /* Incoming auxiliary fields. */
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_totalStrain = 4;
    const PylithInt i_viscousStrain = 5;
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;

    const PylithInt _numS = 1; /* Number passed on to stress kernels. */
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; /* Pass bulk modulus, reference stress, and reference strain. */
    const PylithInt aOffMean[3] = { aOff[i_bulkModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[3] = { aOff_x[i_bulkModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 6; /* Pass shear modulus, Maxwell time, total strain, viscous strain,
                                    reference stress, and reference strain. */
    const PylithInt aOffDev[6] = { aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_totalStrain],
                                   aOff[i_viscousStrain], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[6] = { aOff_x[i_shearModulus], aOff_x[i_maxwellTime], aOff_x[i_totalStrain],
                                     aOff_x[i_viscousStrain], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0};
    PylithInt i;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 6);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);
    assert(numConstants == 1);

    pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_meanStress_refstate(_dim, _numS, numAMean,
                                                                           sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                           aOffMean, aOffMean_x, a, a_t, a_x,
                                                                           t, x, numConstants, constants, stress);
    pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_deviatoricStress_refstate(_dim, _numS, numADev,
                                                                                 sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                                 aOffDev, aOffDev_x, a, a_t, a_x,
                                                                                 t, x, numConstants, constants, stress);

    for (i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } /* for */
} /* IsotropicLinearMaxwellPlaneStrain_g1v_refstate */


/* ---------------------------------------------------------------------- */
/** Jf0 function for isotropic linear Maxwell viscoelastic plane strain with implicit time stepping.
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_Jf0vv_implicit(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar Jf0[])
{ /* IsotropicLinearMaxwellPlaneStrain_Jf0vv_implicit */
    const PylithInt _dim = 2;

    const PylithInt i_density = 0;
    const PylithScalar density = a[aOff[i_density]];

    PylithInt i;

    assert(_dim == dim);
    assert(2 == numS);
    assert(numA >= 6);
    assert(aOff);
    assert(a);

    for (i = 0; i < _dim; ++i) {
        Jf0[i*_dim+i] += utshift * density;
    } /* for */
} /* IsotropicLinearMaxwellPlaneStrain_Jf0vv_implicit */


/* ---------------------------------------------------------------------- */
/** Jf0 function for isotropic linear Maxwell viscoelastic plane strain with explicit time stepping.
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_Jf0vv_explicit(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar Jf0[])
{ /* IsotropicLinearMaxwellPlaneStrain_Jf0vv_explicit */
    const PylithInt _dim = 2;

    const PylithInt i_density = 0;
    const PylithScalar density = a[aOff[i_density]];

    PylithInt i;

    assert(_dim == dim);
    assert(2 == numS);
    assert(numA >= 6);
    assert(aOff);
    assert(a);

    for (i = 0; i < _dim; ++i) {
        Jf0[i*_dim+i] += density;
    } /* for */

} /* IsotropicLinearMaxwellPlaneStrain_Jf0vv_explicit*/


/* ---------------------------------------------------------------------- */
/* Jg3_vu entry function for 2-D plane strain isotropic linear Maxwell viscoelastic.
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
 *  C_ijkl = bulkModulus * delta_ij * delta_kl + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_Jg3vu(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar Jg3[])
{ /* IsotropicLinearMaxwellPlaneStrain_Jg3vu */
    const PylithInt _dim = 2;

    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_maxwellTime = 3;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];
    const PylithScalar dt = constants[0];

    const PylithScalar dq = pylith_fekernels_Elasticity_Maxwell_VisStrain_Coeff(dt, maxwellTime);

    /* Unique components of Jacobian. */
    const PylithReal C1111 = bulkModulus + 4.0 * shearModulus * dq/3.0;
    const PylithReal C1122 = bulkModulus - 2.0 * shearModulus * dq/3.0;
    const PylithReal C1212 = shearModulus * dq;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 6);
    assert(aOff);
    assert(a);
    assert(Jg3);
    assert(numConstants == 1);

    /* j(f,g,df,dg) = C(f,df,g,dg)

       0:  j0000 = C1111 = 1.0*bulkModulus + 1.33333333333333*delHM*shearModulus
       1:  j0001 = C1112 = 0
       2:  j0010 = C1211 = 0
       3:  j0011 = C1212 = 1.0*delHM*shearModulus
       4:  j0100 = C1121 = 0
       5:  j0101 = C1122 = 1.0*bulkModulus - 0.666666666666667*delHM*shearModulus
       6:  j0110 = C1221 = 1.0*delHM*shearModulus
       7:  j0111 = C1222 = 0
       8:  j1000 = C2111 = 0
       9:  j1001 = C2112 = 1.0*delHM*shearModulus
       10:  j1010 = C2211 = 1.0*bulkModulus - 0.666666666666667*delHM*shearModulus
       11:  j1011 = C2212 = 0
       12:  j1100 = C2121 = 1.0*delHM*shearModulus
       13:  j1101 = C2122 = 0
       14:  j1110 = C2221 = 0
       15:  j1111 = C2222 = 1.0*bulkModulus + 1.33333333333333*delHM*shearModulus
     */

    /* Nonzero Jacobian entries. */
    Jg3[0] -= C1111; /* j0000 */
    Jg3[3] -= C1212; /* j0011 */
    Jg3[5] -= C1122; /* j0101 */
    Jg3[6] -= C1212; /* j0110 */
    Jg3[9] -= C1212; /* j1001 */
    Jg3[10] -= C1122; /* j1010 */
    Jg3[12] -= C1212; /* j1100 */
    Jg3[15] -= C1111; /* j1111 */

} /* IsotropicLinearMaxwellPlaneStrain_Jg3vu */


/* ---------------------------------------------------------------------- */
/* Calculate mean stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT reference stress and strain.
 *
 * meanStress = bulkModulus * strain_kk
 *
 * stress += meanStress * delta_ij
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_meanStress(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar stress[])
{ /* IsotropicLinearMaxwellPlaneStrain_meanStress */
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
} /* IsotropicLinearMaxwellPlaneStrain_meanStress */


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
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_meanStress_refstate(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar stress[])
{ /* IsotropicLinearMaxwellPlaneStrain_meanStress_refstate */
    const PylithInt _dim = 2;

    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    const PylithInt i_bulkModulus = 0;
    const PylithInt i_rstress = 1;
    const PylithInt i_rstrain = 2;
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; /* sigma_11, sigma_22, sigma_33, sigma_12 */
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; /* epsilon_11, epsilon_22, epsilon_33, epsilon_12 */

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
    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[2];

    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithReal meanStress = meanrstress + bulkModulus * (strainTrace - refstrainTrace);

    for (i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += meanStress;
    } /* for */
} /* IsotropicLinearMaxwellPlaneStrain_meanStress_refstate */


/* ---------------------------------------------------------------------- */
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = 2*shearModulus*visStrain_ij
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_deviatoricStress(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar stress[])
{ /* IsotropicLinearElasticityPlaneStrain_deviatoricStress */
    const PylithInt _dim = 2;

    const PylithInt i_disp = 0;

    const PylithInt i_shearModulus = 0;
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithInt i_maxwellTime = 1;
    const PylithInt i_totalStrain = 2;
    const PylithInt i_viscousStrain = 3;

    const PylithInt _numS = 1; /* Number passed on to visStrain kernel. */
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; /* Number passed on to visStrain kernel. */
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_totalStrain], aOff[i_viscousStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_totalStrain],
                                     aOff_x[i_viscousStrain] };

    PylithScalar visStrainTpdt[4] = {0.0, 0.0, 0.0, 0.0}; /* Viscous strain tensor. */

    assert(_dim == dim);
    assert(1 == numS);
    assert(4 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);
    assert(numConstants == 1);

    pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_computeVisStrain(_dim, _numS, numAVis,
                                                                        sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                        aOffVis, aOffVis_x,
                                                                        a, a_t, a_x,
                                                                        t, x, numConstants, constants, visStrainTpdt);

    stress[0] += 2.0 * shearModulus * visStrainTpdt[0]; /* sigma_11 */
    stress[1] += 2.0 * shearModulus * visStrainTpdt[3]; /* sigma_12 */
    stress[2] += 2.0 * shearModulus * visStrainTpdt[3]; /* sigma_21 */
    stress[3] += 2.0 * shearModulus * visStrainTpdt[1]; /* sigma_22 */

} /* IsotropicLinearMaxwellPlaneStrain_deviatoricStress */


/* ---------------------------------------------------------------------- */
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITH reference stress and reference strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = devrefstress_ij + 2*shearModulus*(visstrain_ij - devrefstrain_ij)
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_deviatoricStress_refstate(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar stress[])
{ /* deviatoricStress_IsotropicLinearMaxwellPlaneStrain_refstate */
    const PylithInt _dim = 2;

    const PylithInt i_disp = 0;

    const PylithInt i_shearModulus = 0;

    const PylithInt i_maxwellTime = 1;
    const PylithInt i_totalStrain = 2;
    const PylithInt i_viscousStrain = 3;
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; /* sigma_11, sigma_22, sigma_33, sigma_12 */
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; /* epsilon_11, epsilon_22, epsilon_33, epsilon_12 */

    const PylithInt _numS = 1; /* Number passed on to visStrain kernel. */
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; /* Number passed on to visStrain kernel. */
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_totalStrain], aOff[i_viscousStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_totalStrain],
                                     aOff_x[i_viscousStrain] };

    PylithScalar visStrainTpdt[4] = {0.0, 0.0, 0.0, 0.0}; /* Viscous strain tensor. */

    assert(_dim == dim);
    assert(1 == numS);
    assert(numA >= 8);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);
    assert(numConstants == 1);

    // Compute viscous strain for current time step.
    pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_computeVisStrain(_dim, _numS, numAVis,
                                                                        sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                        aOffVis, aOffVis_x,
                                                                        a, a_t, a_x,
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

    // Compute stress components -- note that we are including reference deviatoric stress for now.
    // This may need to be removed after testing.
    const PylithReal twomu = 2.0*shearModulus;

    const PylithScalar sigma_11 = devRefStress[0] + twomu * (visStrainTpdt[0] - devRefStrain[0]);
    const PylithScalar sigma_22 = devRefStress[1] + twomu * (visStrainTpdt[1] - devRefStrain[1]);
    const PylithScalar sigma_12 = devRefStress[3] + twomu * (visStrainTpdt[3] - devRefStrain[3]);

    stress[0*_dim+0] += sigma_11;
    stress[1*_dim+1] += sigma_22;
    stress[0*_dim+1] += sigma_12;
    stress[1*_dim+0] += sigma_12;

} /* deviatoricStress_IsotropicLinearMaxwellPlaneStrain_refstate */


/* ---------------------------------------------------------------------- */
/* Calculate viscous strain for a Maxwell viscoelastic material.
 *
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_computeVisStrain(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar visStrainTpdt[])
{ /* IsotropicLinearMaxwellPlaneStrain_computeVisStrain */
    const PylithInt _dim = 2;

    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    const PylithInt i_maxwellTime = 0;
    const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];

    const PylithScalar dt = constants[0];

    const PylithInt i_totalStrain = 1;
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];

    const PylithInt i_viscousStrain = 2;
    const PylithScalar* viscousStrain = &a[aOff[i_viscousStrain]];

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(visStrainTpdt);
    assert(numConstants == 1);

    const PylithScalar dq = pylith_fekernels_Elasticity_Maxwell_VisStrain_Coeff(dt, maxwellTime);
    const PylithScalar expFac = exp(-dt/maxwellTime);

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
        totalStrain[3] - meanStrainT,
        totalStrain[2]
    };

    for (int iComp = 0; iComp < 4; ++iComp) {
        visStrainTpdt[iComp] = expFac * viscousStrain[iComp] + dq * (devStrainTpdt[iComp] - devStrainT[iComp]);
    } // for

} /* IsotropicLinearMaxwellPlaneStrain_computeVisStrain */


/* ---------------------------------------------------------------------- */
/* Update total strain for a Maxwell viscoelastic material.
 *
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_updateTotalStrain(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar totalStrainTpdt[])
{ /* IsotropicLinearMaxwellPlaneStrain_updateTotalStrain */
    const PylithInt _dim = 2;

    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(6 >= numA && 10 <= numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(totalStrainTpdt);

    totalStrainTpdt[0] = disp_x[0*_dim+0];
    totalStrainTpdt[1] = disp_x[1*_dim+1];
    totalStrainTpdt[2] = 0.0;
    totalStrainTpdt[3] = 0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]);

} /* IsotropicLinearMaxwellPlaneStrain_updateTotalStrain */


/* ---------------------------------------------------------------------- */
/* Update viscous strain for a Maxwell viscoelastic material.
 *
 */
void
pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_updateVisStrain(
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
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar visStrainTpdt[])
{ /* IsotropicLinearElasticityPlaneStrain_updateVisStrain */
    const PylithInt _dim = 2;

    const PylithInt i_disp = 0;

    const PylithInt i_maxwellTime = 3;
    const PylithInt i_totalStrain = 4;
    const PylithInt i_viscousStrain = 5;

    const PylithInt _numS = 1; /* Number passed on to statevars kernel. */
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; /* Number passed on to viscous strain kernel. */
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_totalStrain], aOff[i_viscousStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_totalStrain],
                                     aOff_x[i_viscousStrain] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 6);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(visStrainTpdt);
    assert(numConstants == 1);

    pylith_fekernels_IsotropicLinearMaxwellPlaneStrain_computeVisStrain(_dim, _numS, numAVis,
                                                                        sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                        aOffVis, aOffVis_x,
                                                                        a, a_t, a_x,
                                                                        t, x, numConstants, constants, visStrainTpdt);


} /* IsotropicLinearMaxwellPlaneStrain_updateVisStrain */


/* End of file */
