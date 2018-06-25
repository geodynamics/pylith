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

#include "pylith/fekernels/IsotropicLinearMaxwell3D.hh"
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/Viscoelastic.hh" // USES Viscoelastic kernels
#include "pylith/fekernels/Elasticity3D.hh" // USES Elasticity3D kernels

#include <cassert> // USES assert()
#include <cmath> // USES exp()
#include <iostream> // debugging.

/* ======================================================================
 * Kernels for isotropic, linear Maxwell viscoelastic 3D material.
 * ======================================================================
 */

// ----------------------------------------------------------------------
// g0 function for isotropic linear Maxwell viscoelastic 3D with both gravity
// and body forces.
void
pylith::fekernels::IsotropicLinearMaxwell3D::g0v_gravbodyforce(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_gravityField = 6;
    const PylithInt i_bodyForce = 7;

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 8);
    assert(aOff);
    assert(aOff_x);

    pylith::fekernels::Elasticity::g0v_grav(_dim, _numS, numAGrav,
                                            NULL, NULL, NULL, NULL, NULL,
                                            aOffGrav, aOffGrav_x, a, a_t, a_x,
                                            t, x, numConstants, constants, g0);
    pylith::fekernels::Elasticity::g0v_bodyforce(_dim, _numS, numABody,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 aOffBody, aOffBody_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, g0);
} // g0v_gravbodyforce


// ----------------------------------------------------------------------
// g0 function for isotropic linear Maxwell viscoelastic 3D with just gravity.
void
pylith::fekernels::IsotropicLinearMaxwell3D::g0v_grav(const PylithInt dim,
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
    const PylithInt _dim = 3;

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_gravityField = 6;

    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 7);
    assert(aOff);
    assert(aOff_x);

    pylith::fekernels::Elasticity::g0v_grav(_dim, _numS, numAGrav,
                                            NULL, NULL, NULL, NULL, NULL,
                                            aOffGrav, aOffGrav_x, a, a_t, a_x,
                                            t, x, numConstants, constants, g0);
} // g0v_grav


// ----------------------------------------------------------------------
// g0 function for isotropic linear Maxwell viscoelastic 3D with only body
// forces.
void
pylith::fekernels::IsotropicLinearMaxwell3D::g0v_bodyforce(const PylithInt dim,
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
    const PylithInt _dim = 3;

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    // Incoming auxiliary fields.
    const PylithInt i_bodyForce = 6;

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 7);
    assert(aOff);
    assert(aOff_x);

    pylith::fekernels::Elasticity::g0v_bodyforce(_dim, _numS, numABody,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 aOffBody, aOffBody_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, g0);
} // g0v_bodyforce


// ----------------------------------------------------------------------
// g1 function for isotropic linear Maxwell 3D WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearMaxwell3D::g1v(const PylithInt dim,
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_viscousStrain = 4;
    const PylithInt i_totalStrain = 5;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 6);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };
    const PylithInt aOffMean_x[1] = { aOff_x[i_bulkModulus] };

    const PylithInt numADev = 4; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[4] = {
        aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain]
    };
    const PylithInt aOffDev_x[4] = {
        aOff_x[i_shearModulus], aOff_x[i_maxwellTime], aOff_x[i_viscousStrain], aOff_x[i_totalStrain]
    };

    PylithScalar stress[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Full stress tensor

    pylith::fekernels::Elasticity3D::meanStress(_dim, _numS, numAMean,
												sOffDisp, sOffDisp_x, s, s_t, s_x,
												aOffMean, aOffMean_x, a, a_t, a_x,
												t, x, numConstants, constants, stress);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, aOffDev_x, a, a_t, a_x,
                     t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v


// ----------------------------------------------------------------------
// g1 function for isotropic linear Maxwell viscoelastic 3D with reference
// stress and strain.
void
pylith::fekernels::IsotropicLinearMaxwell3D::g1v_refstate(const PylithInt dim,
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_viscousStrain = 4;
    const PylithInt i_totalStrain = 5;
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 8);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Pass bulk modulus, reference stress, and reference strain.
    const PylithInt aOffMean[3] = { aOff[i_bulkModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[3] = { aOff_x[i_bulkModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 6; // Pass shear modulus, Maxwell time, viscous strain, total strain,
                                 // reference stress, and reference strain.
    const PylithInt aOffDev[6] = { aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_viscousStrain],
                                   aOff[i_totalStrain], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[6] = { aOff_x[i_shearModulus], aOff_x[i_maxwellTime], aOff_x[i_viscousStrain],
                                     aOff_x[i_totalStrain], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stress[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    pylith::fekernels::Elasticity3D::meanStress_refstate(_dim, _numS, numAMean,
														 sOffDisp, sOffDisp_x, s, s_t, s_x,
														 aOffMean, aOffMean_x, a, a_t, a_x,
														 t, x, numConstants, constants, stress);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, aOffDev_x, a, a_t, a_x,
                              t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v_refstate


// ----------------------------------------------------------------------
/* Jg3_vu entry function for 3-D isotropic linear Maxwell viscoelastic.
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
pylith::fekernels::IsotropicLinearMaxwell3D::Jg3vu(const PylithInt dim,
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
												   PylithScalar Jg3[]) {
    const PylithInt _dim = 3;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_maxwellTime = 3;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 6);
    assert(aOff);
    assert(a);
    assert(Jg3);
    assert(numConstants == 1);
    assert(constants);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];
    const PylithScalar dt = constants[0];

    const PylithScalar dq = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime);

	/* Unique components of Jacobian. */
	const PylithReal C1111 = bulkModulus + 4.0*dq*shearModulus/3.0;
	const PylithReal C1122 = bulkModulus - 2.0*dq*shearModulus/3.0;
	const PylithReal C1212 = dq*shearModulus;

    /* j(f,g,df,dg) = C(f,df,g,dg)

	   0:  j0000 = C1111 = bulkModulus + 4*dq*shearModulus/3
	   1:  j0001 = C1112 = 0
	   2:  j0002 = C1113 = 0
	   3:  j0010 = C1211 = 0
	   4:  j0011 = C1212 = dq*shearModulus
	   5:  j0012 = C1213 = 0
	   6:  j0020 = C1311 = 0
	   7:  j0021 = C1312 = 0
	   8:  j0022 = C1313 = dq*shearModulus
	   9:  j0100 = C1121 = 0
	   10:  j0101 = C1122 = bulkModulus - 2*dq*shearModulus/3
	   11:  j0102 = C1123 = 0
	   12:  j0110 = C1221 = dq*shearModulus
	   13:  j0111 = C1222 = 0
	   14:  j0112 = C1223 = 0
	   15:  j0120 = C1321 = 0
	   16:  j0121 = C1322 = 0
	   17:  j0122 = C1323 = 0
	   18:  j0200 = C1131 = 0
	   19:  j0201 = C1132 = 0
	   20:  j0202 = C1133 = bulkModulus - 2*dq*shearModulus/3
	   21:  j0210 = C1231 = 0
	   22:  j0211 = C1232 = 0
	   23:  j0212 = C1233 = 0
	   24:  j0220 = C1331 = dq*shearModulus
	   25:  j0221 = C1332 = 0
	   26:  j0222 = C1333 = 0
	   27:  j1000 = C2111 = 0
	   28:  j1001 = C2112 = dq*shearModulus
	   29:  j1002 = C2113 = 0
	   30:  j1010 = C2211 = bulkModulus - 2*dq*shearModulus/3
	   31:  j1011 = C2212 = 0
	   32:  j1012 = C2213 = 0
	   33:  j1020 = C2311 = 0
	   34:  j1021 = C2312 = 0
	   35:  j1022 = C2313 = 0
	   36:  j1100 = C2121 = dq*shearModulus
	   37:  j1101 = C2122 = 0
	   38:  j1102 = C2123 = 0
	   39:  j1110 = C2221 = 0
	   40:  j1111 = C2222 = bulkModulus + 4*dq*shearModulus/3
	   41:  j1112 = C2223 = 0
	   42:  j1120 = C2321 = 0
	   43:  j1121 = C2322 = 0
	   44:  j1122 = C2323 = dq*shearModulus
	   45:  j1200 = C2131 = 0
	   46:  j1201 = C2132 = 0
	   47:  j1202 = C2133 = 0
	   48:  j1210 = C2231 = 0
	   49:  j1211 = C2232 = 0
	   50:  j1212 = C2233 = bulkModulus - 2*dq*shearModulus/3
	   51:  j1220 = C2331 = 0
	   52:  j1221 = C2332 = dq*shearModulus
	   53:  j1222 = C2333 = 0
	   54:  j2000 = C3111 = 0
	   55:  j2001 = C3112 = 0
	   56:  j2002 = C3113 = dq*shearModulus
	   57:  j2010 = C3211 = 0
	   58:  j2011 = C3212 = 0
	   59:  j2012 = C3213 = 0
	   60:  j2020 = C3311 = bulkModulus - 2*dq*shearModulus/3
	   61:  j2021 = C3312 = 0
	   62:  j2022 = C3313 = 0
	   63:  j2100 = C3121 = 0
	   64:  j2101 = C3122 = 0
	   65:  j2102 = C3123 = 0
	   66:  j2110 = C3221 = 0
	   67:  j2111 = C3222 = 0
	   68:  j2112 = C3223 = dq*shearModulus
	   69:  j2120 = C3321 = 0
	   70:  j2121 = C3322 = bulkModulus - 2*dq*shearModulus/3
	   71:  j2122 = C3323 = 0
	   72:  j2200 = C3131 = dq*shearModulus
	   73:  j2201 = C3132 = 0
	   74:  j2202 = C3133 = 0
	   75:  j2210 = C3231 = 0
	   76:  j2211 = C3232 = dq*shearModulus
	   77:  j2212 = C3233 = 0
	   78:  j2220 = C3331 = 0
	   79:  j2221 = C3332 = 0
	   80:  j2222 = C3333 = bulkModulus + 4*dq*shearModulus/3
	*/

	/* Nonzero Jacobian entries. */
	Jg3[0] -=  C1111; /* j0000 */
	Jg3[4] -=  C1212; /* j0011 */
	Jg3[8] -=  C1212; /* j0022 */
	Jg3[10] -=  C1122; /* j0101 */
	Jg3[12] -=  C1212; /* j0110 */
	Jg3[20] -=  C1122; /* j0202 */
	Jg3[24] -=  C1212; /* j0220 */
	Jg3[28] -=  C1212; /* j1001 */
	Jg3[30] -=  C1122; /* j1010 */
	Jg3[36] -=  C1212; /* j1100 */
	Jg3[40] -=  C1111; /* j1111 */
	Jg3[44] -=  C1212; /* j1122 */
	Jg3[50] -=  C1122; /* j1212 */
	Jg3[52] -=  C1212; /* j1221 */
	Jg3[56] -=  C1212; /* j2002 */
	Jg3[60] -=  C1122; /* j2020 */
	Jg3[68] -=  C1212; /* j2112 */
	Jg3[70] -=  C1122; /* j2121 */
	Jg3[72] -=  C1212; /* j2200 */
	Jg3[76] -=  C1212; /* j2211 */
	Jg3[80] -=  C1111; /* j2222 */

} // Jg3vu


// ----------------------------------------------------------------------
/* Calculate deviatoric stress for 3-D isotropic linear
 * Maxwell viscoelastic WITHOUT reference stress and strain.
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
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithInt _numS = 1; // Number passed on to visStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain],
								   aOff[i_totalStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_viscousStrain],
                                     aOff_x[i_totalStrain] };

    PylithScalar visStrainTpdt[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor (vector).

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x,
						 aOffVis, aOffVis_x, a, a_t, a_x,
                         t, x, numConstants, constants, visStrainTpdt);

    stress[0] += 2.0 * shearModulus * visStrainTpdt[0]; // sigma_11
    stress[1] += 2.0 * shearModulus * visStrainTpdt[3]; // sigma_12
    stress[2] += 2.0 * shearModulus * visStrainTpdt[5]; // sigma_13
    stress[3] += stress[1]; // sigma_21
    stress[4] += 2.0 * shearModulus * visStrainTpdt[1]; // sigma_22
    stress[5] += 2.0 * shearModulus * visStrainTpdt[4]; // sigma_23
    stress[6] += stress[2]; // sigma_31
    stress[7] += stress[5]; // sigma_32
    stress[8] += 2.0 * shearModulus * visStrainTpdt[2]; // sigma_33

} // deviatoricStress


// ----------------------------------------------------------------------
/* Calculate deviatoric stress for 3-D isotropic linear
 * Maxwell viscoelastic WITH reference stress and reference strain.
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
    const PylithInt i_shearModulus = 0;
    const PylithInt i_maxwellTime = 1;
    const PylithInt i_viscousStrain = 2;
    const PylithInt i_totalStrain = 3;
    const PylithInt i_rstress = 4;
    const PylithInt i_rstrain = 5;

    assert(_dim == dim);
    assert(1 == numS);
    assert(6 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // sigma_11, sigma_22, sigma_33, sigma_12, sigma_23, sigma_13
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // epsilon_11, epsilon_22, epsilon_33, epsilon_12, epsilon_23, epsilon_13

    const PylithInt _numS = 1; // Number passed on to visStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_viscousStrain],
                                     aOff_x[i_totalStrain] };

    PylithScalar visStrainTpdt[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor (vector).

    // Compute viscous strain for current time step.
    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x,
						 aOffVis, aOffVis_x, a, a_t, a_x,
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

    const PylithScalar sigma_11 = devRefStress[0] + twomu * (visStrainTpdt[0] - devRefStrain[0]);
    const PylithScalar sigma_22 = devRefStress[1] + twomu * (visStrainTpdt[1] - devRefStrain[1]);
    const PylithScalar sigma_33 = devRefStress[2] + twomu * (visStrainTpdt[2] - devRefStrain[2]);
    const PylithScalar sigma_12 = devRefStress[3] + twomu * (visStrainTpdt[3] - devRefStrain[3]);
    const PylithScalar sigma_23 = devRefStress[4] + twomu * (visStrainTpdt[4] - devRefStrain[4]);
    const PylithScalar sigma_13 = devRefStress[5] + twomu * (visStrainTpdt[5] - devRefStrain[5]);

    stress[0*_dim+0] += sigma_11;
    stress[1*_dim+1] += sigma_22;
    stress[2*_dim+2] += sigma_33;
    stress[0*_dim+1] += sigma_12;
    stress[1*_dim+0] += sigma_12;
    stress[0*_dim+2] += sigma_13;
    stress[2*_dim+0] += sigma_13;
    stress[1*_dim+2] += sigma_23;
    stress[2*_dim+1] += sigma_23;

} // deviatoricStress_refstate


// ----------------------------------------------------------------------
/* Calculate viscous strain for a Maxwell viscoelastic material.
 *
 */
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
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(visStrainTpdt);
    assert(1 == numConstants);
    assert(constants);

    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];
    const PylithScalar* viscousStrain = &a[aOff[i_viscousStrain]];
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];

    const PylithScalar dt = constants[0];

    const PylithScalar dq = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime);
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


// ----------------------------------------------------------------------
/* Update total strain for a Maxwell viscoelastic material.
 * NOTE:  We are assuming right now that solution and auxiliary variables
 * are interchanged for this function.
 *
 */
void
pylith::fekernels::IsotropicLinearMaxwell3D::updateTotalStrain(const PylithInt dim,
															   const PylithInt numA,
															   const PylithInt numS,
															   const PylithInt aOff[],
															   const PylithInt aOff_x[],
															   const PylithScalar a[],
															   const PylithScalar a_t[],
															   const PylithScalar a_x[],
															   const PylithInt sOff[],
															   const PylithInt sOff_x[],
															   const PylithScalar s[],
															   const PylithScalar s_t[],
															   const PylithScalar s_x[],
															   const PylithReal t,
															   const PylithScalar x[],
															   const PylithInt numConstants,
															   const PylithScalar constants[],
															   PylithScalar totalStrainTpdt[]) {
    const PylithInt _dim = 3;

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
    std::cout << "totalStrainTpdt[0]:  " << totalStrainTpdt[0] << std::endl;
#endif

    assert(_dim == dim);
    // assert(1 == numS || 2 == numS);
    assert(6 <= numA && 10 >= numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(totalStrainTpdt);

    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
#if 0 // :DEBUG:
    const PylithInt i_totalStrain = 5;
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];
	std::cout << "fekernels::IsotropicLinearMaxwell3D::updateTotalStrain" << std::endl;
    std::cout << "totalStrain[0]:  " << totalStrain[0] << std::endl;
#endif

    totalStrainTpdt[0] = disp_x[0*_dim+0];
    totalStrainTpdt[1] = disp_x[1*_dim+1];
    totalStrainTpdt[2] = disp_x[2*_dim+2];
    totalStrainTpdt[3] = 0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]);
    totalStrainTpdt[4] = 0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]);
    totalStrainTpdt[5] = 0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0]);
#if 0 // :DEBUG:
    std::cout << "totalStrainTpdt[0]:  " << totalStrainTpdt[0] << std::endl;
#endif

} // updateTotalStrain


// ----------------------------------------------------------------------
/* Update viscous strain for a Maxwell viscoelastic material.
 *
 */
void
pylith::fekernels::IsotropicLinearMaxwell3D::updateViscousStrain(const PylithInt dim,
																 const PylithInt numA,
																 const PylithInt numS,
																 const PylithInt aOff[],
																 const PylithInt aOff_x[],
																 const PylithScalar a[],
																 const PylithScalar a_t[],
																 const PylithScalar a_x[],
																 const PylithInt sOff[],
																 const PylithInt sOff_x[],
																 const PylithScalar s[],
																 const PylithScalar s_t[],
																 const PylithScalar s_x[],
																 const PylithReal t,
																 const PylithScalar x[],
																 const PylithInt numConstants,
																 const PylithScalar constants[],
																 PylithScalar visStrainTpdt[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_viscousStrain = 4;
    const PylithInt i_totalStrain = 5;

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
    std::cout << "visStrainTpdt[0]:  " << visStrainTpdt[0] << std::endl;
#endif

    assert(_dim == dim);
    // assert(1 == numS || 2 == numS);
    assert(numA >= 6);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(visStrainTpdt);
#if 0 // :DEBUG:
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar totalStrainTpdt[6] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        disp_x[2*_dim+2],
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]),
        0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]),
        0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0])
    };
	std::cout << "fekernels::IsotropicLinearMaxwell3D::updateViscousStrain" << std::endl;
    std::cout << "totalStrain[0]:  " << totalStrain[0] << std::endl;
    std::cout << "totalStrainTpdt[0]:  " << totalStrainTpdt[0] << std::endl;
#endif

    const PylithInt _numS = 1; // Number passed on to statevars kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; // Number passed on to viscous strain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain],
								   aOff[i_totalStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_viscousStrain],
                                     aOff_x[i_totalStrain] };

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x,
						 aOffVis, aOffVis_x, a, a_t, a_x,
                         t, x, numConstants, constants, visStrainTpdt);


} // updateViscousStrain


/* End of file */
