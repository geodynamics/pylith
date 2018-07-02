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

#include "pylith/fekernels/IsotropicLinearElasticity3D.hh"
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/Elasticity3D.hh" // USES Elasticity3D kernels

#include <cassert> // USES assert()
#include <iostream> // debugging.

/* ======================================================================
 * Kernels for isotropic, linear 3D elasticity.
 * ======================================================================
 */


// ----------------------------------------------------------------------
// g0 function for isotropic linear elasticity 3D with both gravity and body forces.
void
pylith::fekernels::IsotropicLinearElasticity3D::g0v_gravbodyforce(const PylithInt dim,
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
    const PylithInt i_gravityField = 3;
    const PylithInt i_bodyForce = 4;
#if 0 // :DEBUG:
	std::cout << "IsotropicLinearElasticity3D::g0v_gravbodyforce" << std::endl;
#endif

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 5);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

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
// g0 function for isotropic linear elasticity 3D with just gravity.
void
pylith::fekernels::IsotropicLinearElasticity3D::g0v_grav(const PylithInt dim,
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
    const PylithInt i_gravityField = 3;
#if 0 // :DEBUG:
	std::cout << "IsotropicLinearElasticity3D::g0v_grav" << std::endl;
#endif

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 4);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 0; // Number passed on to g0_grav.

    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    pylith::fekernels::Elasticity::g0v_grav(_dim, _numS, numAGrav,
                                            NULL, NULL, NULL, NULL, NULL,
                                            aOffGrav, aOffGrav_x, a, a_t, a_x,
                                            t, x, numConstants, constants, g0);
} // g0v_grav


// ----------------------------------------------------------------------
// g0 function for isotropic linear elasticity plane strain with just body forces.
void
pylith::fekernels::IsotropicLinearElasticity3D::g0v_bodyforce(const PylithInt dim,
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
    const PylithInt i_bodyForce = 3;
#if 0 // :DEBUG:
	std::cout << "IsotropicLinearElasticity3D::g0v_bodyforce" << std::endl;
#endif

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 4);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    pylith::fekernels::Elasticity::g0v_bodyforce(_dim, _numS, numABody,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 aOffBody, aOffBody_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, g0);
} // 0v_bodyforce


// ----------------------------------------------------------------------
// g1 function for isotropic linear elasticity 3D WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearElasticity3D::g1v(const PylithInt dim,
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
#if 0 // :DEBUG:
	std::cout << "IsotropicLinearElasticity3D::g1v" << std::endl;
#endif

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 3);
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

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };
    const PylithInt aOffDev_x[1] = { aOff_x[i_shearModulus] };

    PylithScalar stress[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Full stress tensor

    pylith::fekernels::Elasticity3D::meanStress(_dim, _numS, numAMean,
												sOffDisp, sOffDisp_x, s, s_t, s_x,
												aOffMean, aOffMean_x, a, a_t, a_x,
												t, x, numConstants, constants, stress);

    pylith::fekernels::Elasticity3D::deviatoricStress(_dim, _numS, numADev,
													  sOffDisp, sOffDisp_x, s, s_t, s_x,
													  aOffDev, aOffDev_x, a, a_t, a_x,
													  t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v

// ----------------------------------------------------------------------
// g1 function for isotropic linear elasticity 3D with reference stress and strain.
void
pylith::fekernels::IsotropicLinearElasticity3D::g1v_refstate(const PylithInt dim,
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
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;
#if 0 // :DEBUG:
	std::cout << "IsotropicLinearElasticity3D::g1v_refstate" << std::endl;
#endif

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 5);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_bulkModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[3] = { aOff_x[i_bulkModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stress[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    pylith::fekernels::Elasticity3D::meanStress_refstate(_dim, _numS, numAMean,
														 sOffDisp, sOffDisp_x, s, s_t, s_x,
														 aOffMean, aOffMean_x, a, a_t, a_x,
														 t, x, numConstants, constants, stress);

    pylith::fekernels::Elasticity3D::deviatoricStress_refstate(_dim, _numS, numADev,
															   sOffDisp, sOffDisp_x, s, s_t, s_x,
															   aOffDev, aOffDev_x, a, a_t, a_x,
															   t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v_refstate


// ----------------------------------------------------------------------
/* Jg3_vu entry function for 3-D isotropic linear elasticity.
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
pylith::fekernels::IsotropicLinearElasticity3D::Jg3vu(const PylithInt dim,
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
#if 0 // :DEBUG:
	std::cout << "IsotropicLinearElasticity3D::Jg3vu" << std::endl;
#endif

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(3 <= numA);
    assert(aOff);
    assert(a);
    assert(Jg3);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar lambda2mu = lambda + 2.0*shearModulus;

    // All other values are either zero or equal to one of these.
    const PylithReal C1111 = lambda2mu;
    const PylithReal C1122 = lambda;
    const PylithReal C1212 = shearModulus;
    /* j(f,g,df,dg) = C(f,df,g,dg)

	   0:  j0000 = C1111 = 1.0*lambda + 2.0*shearModulus
	   1:  j0001 = C1112 = 0
	   2:  j0002 = C1113 = 0
	   3:  j0010 = C1211 = 0
	   4:  j0011 = C1212 = 1.0*shearModulus
	   5:  j0012 = C1213 = 0
	   6:  j0020 = C1311 = 0
	   7:  j0021 = C1312 = 0
	   8:  j0022 = C1313 = 1.0*shearModulus
	   9:  j0100 = C1121 = 0
	   10:  j0101 = C1122 = 1.0*lambda
	   11:  j0102 = C1123 = 0
	   12:  j0110 = C1221 = 1.0*shearModulus
	   13:  j0111 = C1222 = 0
	   14:  j0112 = C1223 = 0
	   15:  j0120 = C1321 = 0
	   16:  j0121 = C1322 = 0
	   17:  j0122 = C1323 = 0
	   18:  j0200 = C1131 = 0
	   19:  j0201 = C1132 = 0
	   20:  j0202 = C1133 = 1.0*lambda
	   21:  j0210 = C1231 = 0
	   22:  j0211 = C1232 = 0
	   23:  j0212 = C1233 = 0
	   24:  j0220 = C1331 = 1.0*shearModulus
	   25:  j0221 = C1332 = 0
	   26:  j0222 = C1333 = 0
	   27:  j1000 = C2111 = 0
	   28:  j1001 = C2112 = 1.0*shearModulus
	   29:  j1002 = C2113 = 0
	   30:  j1010 = C2211 = 1.0*lambda
	   31:  j1011 = C2212 = 0
	   32:  j1012 = C2213 = 0
	   33:  j1020 = C2311 = 0
	   34:  j1021 = C2312 = 0
	   35:  j1022 = C2313 = 0
	   36:  j1100 = C2121 = 1.0*shearModulus
	   37:  j1101 = C2122 = 0
	   38:  j1102 = C2123 = 0
	   39:  j1110 = C2221 = 0
	   40:  j1111 = C2222 = 1.0*lambda + 2.0*shearModulus
	   41:  j1112 = C2223 = 0
	   42:  j1120 = C2321 = 0
	   43:  j1121 = C2322 = 0
	   44:  j1122 = C2323 = 1.0*shearModulus
	   45:  j1200 = C2131 = 0
	   46:  j1201 = C2132 = 0
	   47:  j1202 = C2133 = 0
	   48:  j1210 = C2231 = 0
	   49:  j1211 = C2232 = 0
	   50:  j1212 = C2233 = 1.0*lambda
	   51:  j1220 = C2331 = 0
	   52:  j1221 = C2332 = 1.0*shearModulus
	   53:  j1222 = C2333 = 0
	   54:  j2000 = C3111 = 0
	   55:  j2001 = C3112 = 0
	   56:  j2002 = C3113 = 1.0*shearModulus
	   57:  j2010 = C3211 = 0
	   58:  j2011 = C3212 = 0
	   59:  j2012 = C3213 = 0
	   60:  j2020 = C3311 = 1.0*lambda
	   61:  j2021 = C3312 = 0
	   62:  j2022 = C3313 = 0
	   63:  j2100 = C3121 = 0
	   64:  j2101 = C3122 = 0
	   65:  j2102 = C3123 = 0
	   66:  j2110 = C3221 = 0
	   67:  j2111 = C3222 = 0
	   68:  j2112 = C3223 = 1.0*shearModulus
	   69:  j2120 = C3321 = 0
	   70:  j2121 = C3322 = 1.0*lambda
	   71:  j2122 = C3323 = 0
	   72:  j2200 = C3131 = 1.0*shearModulus
	   73:  j2201 = C3132 = 0
	   74:  j2202 = C3133 = 0
	   75:  j2210 = C3231 = 0
	   76:  j2211 = C3232 = 1.0*shearModulus
	   77:  j2212 = C3233 = 0
	   78:  j2220 = C3331 = 0
	   79:  j2221 = C3332 = 0
	   80:  j2222 = C3333 = 1.0*lambda + 2.0*shearModulus
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
/* Calculate stress for 3-D isotropic linear
 * elasticity WITHOUT a reference stress and strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
 */
void
pylith::fekernels::IsotropicLinearElasticity3D::stress(const PylithInt dim,
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
{
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
#if 0 // :DEBUG:
	std::cout << "IsotropicLinearElasticity3D::stress" << std::endl;
#endif

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 3);
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

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };
    const PylithInt aOffDev_x[1] = { aOff_x[i_shearModulus] };

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    pylith::fekernels::Elasticity3D::meanStress(_dim, _numS, numAMean,
												sOffDisp, sOffDisp_x, s, s_t, s_x,
												aOffMean, aOffMean_x, a, a_t, a_x,
												t, x, numConstants, constants, stressTensor);

    pylith::fekernels::Elasticity3D::deviatoricStress(_dim, _numS, numADev,
													  sOffDisp, sOffDisp_x, s, s_t, s_x,
													  aOffDev, aOffDev_x, a, a_t, a_x,
													  t, x, numConstants, constants, stressTensor);

    stress[0] = stressTensor[0*_dim+0]; // stress_xx
    stress[1] = stressTensor[1*_dim+1]; // stress_yy
    stress[2] = stressTensor[2*_dim+2]; // stress_zz
    stress[3] = stressTensor[0*_dim+1]; // stress_xy
    stress[4] = stressTensor[1*_dim+2]; // stress_yz
    stress[5] = stressTensor[0*_dim+2]; // stress_xz
} // stress


// ----------------------------------------------------------------------
/* Calculate stress for 3-D isotropic linear
 * elasticity WITH a reference stress/strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
 */
void
pylith::fekernels::IsotropicLinearElasticity3D::stress_refstate(const PylithInt dim,
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
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;
#if 0 // :DEBUG:
	std::cout << "IsotropicLinearElasticity3D::stress_refstate" << std::endl;
#endif

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 5);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_bulkModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[3] = { aOff_x[i_bulkModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    pylith::fekernels::Elasticity3D::meanStress_refstate(_dim, _numS, numAMean,
														 sOffDisp, sOffDisp_x, s, s_t, s_x,
														 aOffMean, aOffMean_x, a, a_t, a_x,
														 t, x, numConstants, constants, stressTensor);

    pylith::fekernels::Elasticity3D::deviatoricStress_refstate(_dim, _numS, numADev,
															   sOffDisp, sOffDisp_x, s, s_t, s_x,
															   aOffDev, aOffDev_x, a, a_t, a_x,
															   t, x, numConstants, constants, stressTensor);

    stress[0] = stressTensor[0*_dim+0]; // stress_xx
    stress[1] = stressTensor[1*_dim+1]; // stress_yy
    stress[2] = stressTensor[2*_dim+2]; // stress_zz
    stress[3] = stressTensor[0*_dim+1]; // stress_xy
    stress[4] = stressTensor[1*_dim+2]; // stress_yz
    stress[5] = stressTensor[0*_dim+2]; // stress_xz

} // stress_refstate


// End of file
