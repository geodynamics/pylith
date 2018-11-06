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

#include "pylith/fekernels/IsotropicLinearGenMaxwell3D.hh"
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/Viscoelastic.hh" // USES Viscoelastic kernels
#include "pylith/fekernels/Elasticity3D.hh" // USES Elasticity3D kernels

#include <cassert> // USES assert()
#include <cmath> // USES exp()
#include <iostream> // debugging.

/* ======================================================================
 * Kernels for isotropic, linear Generalized Maxwell viscoelastic 3D.
 * ======================================================================
 */

// ----------------------------------------------------------------------
// g0 function for isotropic linear Generalized Maxwell viscoelastic 3D with
// both gravity and body forces.
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::g0v_gravbodyforce(const PylithInt dim,
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
    const PylithInt i_gravityField = 7;
    const PylithInt i_bodyForce = 8;

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 9);
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
// g0 function for isotropic linear generalized Maxwell viscoelastic 3D with
// just gravity.
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::g0v_grav(const PylithInt dim,
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
    const PylithInt i_gravityField = 7;

    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 8);
    assert(aOff);
    assert(aOff_x);

    pylith::fekernels::Elasticity::g0v_grav(_dim, _numS, numAGrav,
                                            NULL, NULL, NULL, NULL, NULL,
                                            aOffGrav, aOffGrav_x, a, a_t, a_x,
                                            t, x, numConstants, constants, g0);
} // g0v_grav


// ----------------------------------------------------------------------
// g0 function for isotropic linear generalized Maxwell viscoelastic 3D with
// only body forces.
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::g0v_bodyforce(const PylithInt dim,
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
    const PylithInt i_bodyForce = 7;

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 8);
    assert(aOff);
    assert(aOff_x);

    pylith::fekernels::Elasticity::g0v_bodyforce(_dim, _numS, numABody,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 aOffBody, aOffBody_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, g0);
} // g0v_bodyforce


// ----------------------------------------------------------------------
// g1 function for isotropic linear generalized Maxwell 3D WITHOUT reference
// stress and strain.
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::g1v(const PylithInt dim,
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
    const PylithInt i_shearModulusRatio = 4;
    const PylithInt i_viscousStrain = 5;
    const PylithInt i_totalStrain = 6;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 7);
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

    const PylithInt numADev = 5; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[5] = {
		aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_shearModulusRatio],
		aOff[i_viscousStrain], aOff[i_totalStrain]
    };
    const PylithInt aOffDev_x[5] = {
      aOff_x[i_shearModulus], aOff_x[i_maxwellTime], aOff_x[i_shearModulusRatio],
      aOff_x[i_viscousStrain], aOff_x[i_totalStrain]
    };

    PylithScalar stress[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Full stress tensor

    pylith::fekernels::Elasticity3D::meanStress(_dim, _numS, numAMean,
												sOffDisp, sOffDisp_x, s, s_t, s_x,
												aOffMean, aOffMean_x, a, a_t, a_x,
												t, x, numConstants, constants, stress);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, aOffDev_x,
					 a, a_t, a_x, t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v


// ----------------------------------------------------------------------
// g1 function for isotropic linear generalized Maxwell viscoelastic 3D with
// reference stress and strain.
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::g1v_refstate(const PylithInt dim,
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
    const PylithInt i_shearModulusRatio = 4;
    const PylithInt i_viscousStrain = 5;
    const PylithInt i_totalStrain = 6;
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 14);
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

    const PylithInt numADev = 7; // Pass shear modulus, Maxwell times, shear modulus ratios,
                                  // total strain, viscous strains,
                                  // reference stress, and reference strain.
    const PylithInt aOffDev[7] = {
		aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_shearModulusRatio],
		aOff[i_viscousStrain], aOff[i_totalStrain], aOff[i_rstress], aOff[i_rstrain]
    };
    const PylithInt aOffDev_x[7] = {
		aOff_x[i_shearModulus], aOff_x[i_maxwellTime], aOff_x[i_shearModulusRatio],
		aOff_x[i_viscousStrain], aOff_x[i_totalStrain], aOff_x[i_rstress],
		aOff_x[i_rstrain]
    };

    PylithScalar stress[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    pylith::fekernels::Elasticity3D::meanStress_refstate(_dim, _numS, numAMean,
														 sOffDisp, sOffDisp_x, s, s_t, s_x,
														 aOffMean, aOffMean_x, a, a_t, a_x,
														 t, x, numConstants, constants,
														 stress);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev,
							  aOffDev_x, a, a_t, a_x, t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v_refstate


// ----------------------------------------------------------------------
/* Jg3_vu entry function for 3-D isotropic linear generalized Maxwell
 * viscoelastic.
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
pylith::fekernels::IsotropicLinearGenMaxwell3D::Jg3vu(const PylithInt dim,
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
    const PylithInt i_shearModulusRatio = 4;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 7);
    assert(aOff);
    assert(a);
    assert(Jg3);
    assert(numConstants == 1);
    assert(constants);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar maxwellTime_1 = a[aOff[i_maxwellTime]];
    const PylithScalar maxwellTime_2 = a[aOff[i_maxwellTime] + 1];
    const PylithScalar maxwellTime_3 = a[aOff[i_maxwellTime] + 2];
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio] + 1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio] + 2];
    const PylithScalar dt = constants[0];

    const PylithScalar dq_1 = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime_1);
    const PylithScalar dq_2 = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime_2);
    const PylithScalar dq_3 = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime_3);

    // Unique components of Jacobian.
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2 - shearModulusRatio_3;
    const PylithScalar shearFactor = shearModulus *
      (dq_1 * shearModulusRatio_1 + dq_2 * shearModulusRatio_2 + dq_3 * shearModulusRatio_3 +
       shearModulusRatio_0);

    const PylithReal C1111 = bulkModulus + 4.0 * shearFactor/3.0;
    const PylithReal C1122 = bulkModulus - 2.0 * shearFactor/3.0;
    const PylithReal C1212 = shearFactor;

	/* j(f,g,df,dg) = C(f,df,g,dg)

	   0:  j0000 = C1111 = bulkModulus + 2*shearModulus*(2*dq_1*shearModulusRatio_1/3 + 2*dq_2*shearModulusRatio_2/3 + 2*dq_3*shearModulusRatio_3/3 + 2*shearModulusRatio_0/3)
	   1:  j0001 = C1112 = 0
	   2:  j0002 = C1113 = 0
	   3:  j0010 = C1211 = 0
	   4:  j0011 = C1212 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 + dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
	   5:  j0012 = C1213 = 0
	   6:  j0020 = C1311 = 0
	   7:  j0021 = C1312 = 0
	   8:  j0022 = C1313 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 + dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
	   9:  j0100 = C1121 = 0
	   10:  j0101 = C1122 = bulkModulus + 2*shearModulus*(-dq_1*shearModulusRatio_1/3 - dq_2*shearModulusRatio_2/3 - dq_3*shearModulusRatio_3/3 - shearModulusRatio_0/3)
	   11:  j0102 = C1123 = 0
	   12:  j0110 = C1221 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 + dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
	   13:  j0111 = C1222 = 0
	   14:  j0112 = C1223 = 0
	   15:  j0120 = C1321 = 0
	   16:  j0121 = C1322 = 0
	   17:  j0122 = C1323 = 0
	   18:  j0200 = C1131 = 0
	   19:  j0201 = C1132 = 0
	   20:  j0202 = C1133 = bulkModulus + 2*shearModulus*(-dq_1*shearModulusRatio_1/3 - dq_2*shearModulusRatio_2/3 - dq_3*shearModulusRatio_3/3 - shearModulusRatio_0/3)
	   21:  j0210 = C1231 = 0
	   22:  j0211 = C1232 = 0
	   23:  j0212 = C1233 = 0
	   24:  j0220 = C1331 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 + dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
	   25:  j0221 = C1332 = 0
	   26:  j0222 = C1333 = 0
	   27:  j1000 = C2111 = 0
	   28:  j1001 = C2112 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 + dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
	   29:  j1002 = C2113 = 0
	   30:  j1010 = C2211 = bulkModulus + 2*shearModulus*(-dq_1*shearModulusRatio_1/3 - dq_2*shearModulusRatio_2/3 - dq_3*shearModulusRatio_3/3 - shearModulusRatio_0/3)
	   31:  j1011 = C2212 = 0
	   32:  j1012 = C2213 = 0
	   33:  j1020 = C2311 = 0
	   34:  j1021 = C2312 = 0
	   35:  j1022 = C2313 = 0
	   36:  j1100 = C2121 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 + dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
	   37:  j1101 = C2122 = 0
	   38:  j1102 = C2123 = 0
	   39:  j1110 = C2221 = 0
	   40:  j1111 = C2222 = bulkModulus + 2*shearModulus*(2*dq_1*shearModulusRatio_1/3 + 2*dq_2*shearModulusRatio_2/3 + 2*dq_3*shearModulusRatio_3/3 + 2*shearModulusRatio_0/3)
	   41:  j1112 = C2223 = 0
	   42:  j1120 = C2321 = 0
	   43:  j1121 = C2322 = 0
	   44:  j1122 = C2323 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 + dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
	   45:  j1200 = C2131 = 0
	   46:  j1201 = C2132 = 0
	   47:  j1202 = C2133 = 0
	   48:  j1210 = C2231 = 0
	   49:  j1211 = C2232 = 0
	   50:  j1212 = C2233 = bulkModulus + 2*shearModulus*(-dq_1*shearModulusRatio_1/3 - dq_2*shearModulusRatio_2/3 - dq_3*shearModulusRatio_3/3 - shearModulusRatio_0/3)
	   51:  j1220 = C2331 = 0
	   52:  j1221 = C2332 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 + dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
	   53:  j1222 = C2333 = 0
	   54:  j2000 = C3111 = 0
	   55:  j2001 = C3112 = 0
	   56:  j2002 = C3113 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 + dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
	   57:  j2010 = C3211 = 0
	   58:  j2011 = C3212 = 0
	   59:  j2012 = C3213 = 0
	   60:  j2020 = C3311 = bulkModulus + 2*shearModulus*(-dq_1*shearModulusRatio_1/3 - dq_2*shearModulusRatio_2/3 - dq_3*shearModulusRatio_3/3 - shearModulusRatio_0/3)
	   61:  j2021 = C3312 = 0
	   62:  j2022 = C3313 = 0
	   63:  j2100 = C3121 = 0
	   64:  j2101 = C3122 = 0
	   65:  j2102 = C3123 = 0
	   66:  j2110 = C3221 = 0
	   67:  j2111 = C3222 = 0
	   68:  j2112 = C3223 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 + dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
	   69:  j2120 = C3321 = 0
	   70:  j2121 = C3322 = bulkModulus + 2*shearModulus*(-dq_1*shearModulusRatio_1/3 - dq_2*shearModulusRatio_2/3 - dq_3*shearModulusRatio_3/3 - shearModulusRatio_0/3)
	   71:  j2122 = C3323 = 0
	   72:  j2200 = C3131 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 + dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
	   73:  j2201 = C3132 = 0
	   74:  j2202 = C3133 = 0
	   75:  j2210 = C3231 = 0
	   76:  j2211 = C3232 = 2*shearModulus*(dq_1*shearModulusRatio_1/2 + dq_2*shearModulusRatio_2/2 + dq_3*shearModulusRatio_3/2 + shearModulusRatio_0/2)
	   77:  j2212 = C3233 = 0
	   78:  j2220 = C3331 = 0
	   79:  j2221 = C3332 = 0
	   80:  j2222 = C3333 = bulkModulus + 2*shearModulus*(2*dq_1*shearModulusRatio_1/3 + 2*dq_2*shearModulusRatio_2/3 + 2*dq_3*shearModulusRatio_3/3 + 2*shearModulusRatio_0/3)
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
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * generalized Maxwell viscoelastic WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = 2*shearModulus*visStrain_ij
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::deviatoricStress(const PylithInt dim,
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
	const PylithInt _strainSize = 6;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_maxwellTime = 1;
    const PylithInt i_shearModulusRatio = 2;
    const PylithInt i_viscousStrain = 3;
    const PylithInt i_totalStrain = 4;

    assert(_dim == dim);
    assert(1 == numS);
    assert(5 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithInt _numS = 1; // Number passed on to visStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    // Viscous strain.
    const PylithInt numAVis = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain],
								   aOff[i_totalStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_viscousStrain],
									 aOff_x[i_totalStrain] };

    PylithScalar visStrainTpdt[18] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
									  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
									  0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x,
						 aOffVis, aOffVis_x, a, a_t, a_x,
						 t, x, numConstants, constants, visStrainTpdt);

    // Shear modulus ratio factors.
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio] + 1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio] + 2];
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 -
		shearModulusRatio_2 - shearModulusRatio_3;

    // Deviatoric strains.
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar strainTpdt[6] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        disp_x[2*_dim+2],
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]),
        0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]),
        0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0])
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1] + strainTpdt[2])/3.0;

    const PylithScalar devStrainTpdt[6] = {
        strainTpdt[0] - meanStrainTpdt,
        strainTpdt[1] - meanStrainTpdt,
        strainTpdt[2] - meanStrainTpdt,
        strainTpdt[3],
        strainTpdt[4],
        strainTpdt[5]
    };

    // Stresses.
    stress[0] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[0] +
				       shearModulusRatio_1 * visStrainTpdt[0] + 
				       shearModulusRatio_2 * visStrainTpdt[_strainSize] + 
				       shearModulusRatio_3 * visStrainTpdt[2 * _strainSize]); // sigma_11
    stress[1] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[3] +
				       shearModulusRatio_1 * visStrainTpdt[3] + 
				       shearModulusRatio_2 * visStrainTpdt[3 + _strainSize] + 
				       shearModulusRatio_3 * visStrainTpdt[3 + 2 * _strainSize]); // sigma_12
    stress[2] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[5] +
				       shearModulusRatio_1 * visStrainTpdt[5] + 
				       shearModulusRatio_2 * visStrainTpdt[5 + _strainSize] + 
				       shearModulusRatio_3 * visStrainTpdt[5 + 2 * _strainSize]); // sigma_13
    stress[3] += stress[1];                                       // sigma_21
    stress[4] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[1] +
				       shearModulusRatio_1 * visStrainTpdt[1] + 
				       shearModulusRatio_2 * visStrainTpdt[1 + _strainSize] + 
				       shearModulusRatio_3 * visStrainTpdt[1 + 2 * _strainSize]); // sigma_22
    stress[5] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[4] +
				       shearModulusRatio_1 * visStrainTpdt[4] + 
				       shearModulusRatio_2 * visStrainTpdt[4 + _strainSize] + 
				       shearModulusRatio_3 * visStrainTpdt[4 + 2 * _strainSize]); // sigma_23
    stress[6] += stress[2];                                       // sigma_31
    stress[7] += stress[5];                                       // sigma_32
    stress[8] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[2] +
				       shearModulusRatio_1 * visStrainTpdt[2] + 
				       shearModulusRatio_2 * visStrainTpdt[2 + _strainSize] + 
				       shearModulusRatio_3 * visStrainTpdt[2 + 2 * _strainSize]); // sigma_33
	
} // deviatoricStress


// ----------------------------------------------------------------------
/* Calculate deviatoric stress for 3-D isotropic linear
 * generalized Maxwell viscoelastic WITH reference stress and reference strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = devrefstress_ij + 2*shearModulus*(visstrain_ij - devrefstrain_ij)
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::deviatoricStress_refstate(const PylithInt dim,
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
	const PylithInt _strainSize = 6;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_maxwellTime = 1;
    const PylithInt i_shearModulusRatio = 2;
    const PylithInt i_viscousStrain = 3;
    const PylithInt i_totalStrain = 4;
    const PylithInt i_rstress = 5;
    const PylithInt i_rstrain = 6;

    assert(_dim == dim);
    assert(1 == numS);
    assert(7 == numA);
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

    // Viscous strain.
    const PylithInt numAVis = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain],
								   aOff[i_totalStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_viscousStrain],
									 aOff_x[i_totalStrain] };

    PylithScalar visStrainTpdt[18] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
									  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
									  0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

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

    // Shear modulus ratio factors.
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio] + 1];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio] + 2];
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 -
		shearModulusRatio_2 - shearModulusRatio_3;

    // Deviatoric strains.
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar strainTpdt[6] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        disp_x[2*_dim+2],
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]),
        0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]),
        0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0])
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1] + strainTpdt[2])/3.0;

    const PylithScalar devStrainTpdt[6] = {
        strainTpdt[0] - meanStrainTpdt,
        strainTpdt[1] - meanStrainTpdt,
        strainTpdt[2] - meanStrainTpdt,
        strainTpdt[3],
        strainTpdt[4],
        strainTpdt[5]
    };

    // Compute stress components -- note that we are including reference deviatoric
	// stress for now. This may need to be removed after testing.
    const PylithReal twomu = 2.0*shearModulus;
    stress[0] += devRefStress[0] + twomu * (shearModulusRatio_0 * devStrainTpdt[0] +
					    shearModulusRatio_1 * visStrainTpdt[0] + 
					    shearModulusRatio_2 * visStrainTpdt[_strainSize] + 
					    shearModulusRatio_3 * visStrainTpdt[2 * _strainSize] -
					    devRefStrain[0]); // sigma_11
    stress[1] += devRefStress[3] + twomu * (shearModulusRatio_0 * devStrainTpdt[3] +
					    shearModulusRatio_1 * visStrainTpdt[3] + 
					    shearModulusRatio_2 * visStrainTpdt[3 + _strainSize] + 
					    shearModulusRatio_3 * visStrainTpdt[3 + 2 * _strainSize] -
					    devRefStrain[3]); // sigma_12
    stress[2] += devRefStress[5] + twomu * (shearModulusRatio_0 * devStrainTpdt[5] +
					    shearModulusRatio_1 * visStrainTpdt[5] + 
					    shearModulusRatio_2 * visStrainTpdt[5 + _strainSize] + 
					    shearModulusRatio_3 * visStrainTpdt[5 + 2 * _strainSize] -
					    devRefStrain[5]); // sigma_13
    stress[3] += stress[1];               // sigma_21
    stress[4] += devRefStress[1] + twomu * (shearModulusRatio_0 * devStrainTpdt[1] +
					    shearModulusRatio_1 * visStrainTpdt[1] + 
					    shearModulusRatio_2 * visStrainTpdt[1 + _strainSize] + 
					    shearModulusRatio_3 * visStrainTpdt[1 + 2 * _strainSize] -
					    devRefStrain[1]); // sigma_22
    stress[5] += devRefStress[4] + twomu * (shearModulusRatio_0 * devStrainTpdt[4] +
					    shearModulusRatio_1 * visStrainTpdt[4] + 
					    shearModulusRatio_2 * visStrainTpdt[4 + _strainSize] + 
					    shearModulusRatio_3 * visStrainTpdt[4 + 2 * _strainSize] -
					    devRefStrain[4]); // sigma_23
    stress[6] += stress[2];               // sigma_31
    stress[7] += stress[5];               // sigma_32
    stress[8] += devRefStress[2] + twomu * (shearModulusRatio_0 * devStrainTpdt[2] +
					    shearModulusRatio_1 * visStrainTpdt[2] + 
					    shearModulusRatio_2 * visStrainTpdt[2 + _strainSize] + 
					    shearModulusRatio_3 * visStrainTpdt[2 + 2 * _strainSize] -
					    devRefStrain[2]); // sigma_33

} // deviatoricStress_refstate


// ----------------------------------------------------------------------
/* Calculate viscous strain for a generalized Maxwell viscoelastic material.
 *
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::computeViscousStrain(const PylithInt dim,
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
	const PylithInt _strainSize = 6;

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
    assert(numConstants == 1);
    assert(constants);

    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    const PylithScalar maxwellTime_1 = a[aOff[i_maxwellTime]];
    const PylithScalar maxwellTime_2 = a[aOff[i_maxwellTime] + 1];
    const PylithScalar maxwellTime_3 = a[aOff[i_maxwellTime] + 2];
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];
    const PylithScalar* viscousStrain_1 = &a[aOff[i_viscousStrain]];
    const PylithScalar* viscousStrain_2 = &a[aOff[i_viscousStrain] + _strainSize];
    const PylithScalar* viscousStrain_3 = &a[aOff[i_viscousStrain] + 2 * _strainSize];

    const PylithScalar dt = constants[0];

    const PylithScalar dq_1 = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime_1);
    const PylithScalar expFac_1 = exp(-dt/maxwellTime_1);

    const PylithScalar dq_2 = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime_2);
    const PylithScalar expFac_2 = exp(-dt/maxwellTime_2);

    const PylithScalar dq_3 = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime_3);
    const PylithScalar expFac_3 = exp(-dt/maxwellTime_3);

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

	PylithScalar strainDiff = 0.0;
    for (int iComp = 0; iComp < _strainSize; ++iComp) {
		strainDiff = devStrainTpdt[iComp] - devStrainT[iComp];
        visStrainTpdt[iComp] = expFac_1 * viscousStrain_1[iComp] + dq_1 * strainDiff;
		visStrainTpdt[iComp + _strainSize] = expFac_2 * viscousStrain_2[iComp] + dq_2 * strainDiff;
		visStrainTpdt[iComp + 2 * _strainSize] = expFac_3 * viscousStrain_3[iComp] + dq_3 * strainDiff;
    } // for

} // computeViscousStrain


// ----------------------------------------------------------------------
/* Update total strain for a Maxwell viscoelastic material.
 * NOTE:  We are assuming right now that solution and auxiliary variables
 * are interchanged for this function.
 *
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::updateTotalStrain(const PylithInt dim,
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
																  PylithScalar totalStrain[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_disp = 2;

	// Assertions.
    assert(_dim == dim);
    assert(3 <= numS);
    assert(7 <= numA && 11 >= numA);
    assert(sOff_x);
    assert(sOff_x[i_disp] >= 0);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(totalStrain);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    totalStrain[0] = disp_x[0*_dim+0];
    totalStrain[1] = disp_x[1*_dim+1];
    totalStrain[2] = disp_x[2*_dim+2];
    totalStrain[3] = 0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]);
    totalStrain[4] = 0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]);
    totalStrain[5] = 0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0]);

#if 0 // :DEBUG:
	std::cout << "fekernels::IsotropicLinearGenMaxwell3D::updateTotalStrain" << std::endl;
	std::cout << "dim:  " << dim << std::endl;
	std::cout << "numS:  " << numS << std::endl;
	std::cout << "numA:  " << numA << std::endl;
	std::cout << "t:  " << t << std::endl;
	std::cout << "x[0]:  " << x[0] << std::endl;
	std::cout << "x[1]:  " << x[1] << std::endl;
	std::cout << "x[2]:  " << x[2] << std::endl;
	const PylithScalar* disp = &s[sOff[i_disp]];
	const PylithInt i_totalStrainPrevious = 1;
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_shearModulusRatio = 4;
	const PylithScalar maxwellTime_1 = a[aOff[i_maxwellTime]];
	const PylithScalar maxwellTime_2 = a[aOff[i_maxwellTime] + 1];
	const PylithScalar maxwellTime_3 = a[aOff[i_maxwellTime] + 2];
	const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]];
	const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio] + 1];
	const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio] + 2];
	const PylithScalar* totalStrainPrevious = &s[sOff[i_totalStrainPrevious]];
	const PylithScalar aa = 1.0e-4;
	const PylithScalar b = 2.5e-4;
	const PylithScalar c = 3.0e-4;
	const PylithScalar d = 3.5e-4;
	const PylithScalar e = 4.0e-4;
	const PylithScalar f = 4.5e-4;
	const PylithScalar g = 9.0e-8;
	const PylithScalar dt = constants[0];
	const PylithScalar dispxPredPrevious =
		(aa*x[0]*x[0] + 2.0*b*x[0]*x[1] + c*x[1]*x[1] +
		 2.0*d*x[0]*x[2] + 2.0*e*x[1]*x[2] + f*x[2]*x[2])*
		(shearModulusRatio_1*exp(t*(1.0/maxwellTime_3 + 1.0/maxwellTime_2)) +
		 shearModulusRatio_2*exp(t*(1.0/maxwellTime_3 + 1.0/maxwellTime_1)) +
		 shearModulusRatio_3*exp(t*(1.0/maxwellTime_2 + 1.0/maxwellTime_1)))*
		exp(-t*(1.0/maxwellTime_3 + 1.0/maxwellTime_2 + 1.0/maxwellTime_1));
	const PylithScalar dispxPred = dispxPredPrevious + g*x[0];
	const PylithScalar totalStrainxxPredPrevious = (2.0*aa*x[0] + 2.0*b*x[1] + 2.0*d*x[2]) *
		(shearModulusRatio_1*exp(-t/maxwellTime_1) +
		 shearModulusRatio_2*exp(-t/maxwellTime_2) +
		 shearModulusRatio_3*exp(-t/maxwellTime_3));
	const PylithScalar totalStrainxxPred = totalStrainxxPredPrevious + g;
	std::cout << "dispx:  " << disp[0] << std::endl;
	std::cout << "dispxPred:  " << dispxPred << std::endl;
	std::cout << "totalStrainxxPrevious:  " << totalStrainPrevious[0] << std::endl;
	std::cout << "totalStrainxxPredPrevious:  " << totalStrainxxPredPrevious << std::endl;
	std::cout << "totalStrainxx:  " << totalStrain[0] << std::endl;
	std::cout << "totalStrainxxPred:  " << totalStrainxxPred << std::endl;
#endif

} // updateTotalStrain


// ----------------------------------------------------------------------
/* Update viscous strain for generalized Maxwell.
 *
 */
void
pylith::fekernels::IsotropicLinearGenMaxwell3D::updateViscousStrain(const PylithInt dim,
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
																	PylithScalar visStrain[]) {
    const PylithInt _dim = 3;
    const PylithInt _strainSize = 6;

    // Incoming solution fields.
    const PylithInt i_viscousStrainPrevious = 0;
    const PylithInt i_totalStrainPrevious = 1;
    const PylithInt i_disp = 2;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_viscousStrain = 5;
    const PylithInt i_totalStrain = 6;

	// Assertions.
    assert(_dim == dim);
    assert(3 == numS);
    assert(numA >= 7);
    assert(sOff);
	assert(sOff[i_viscousStrainPrevious] >= 0);
	assert(sOff[i_totalStrainPrevious] >= 0);
	assert(sOff[i_disp] >= 0);
    assert(aOff);
	assert(aOff[i_maxwellTime] >= 0);
    assert(aOff[i_viscousStrain] >= 0);
    assert(aOff[i_totalStrain] >= 0);
    assert(s_x);
    assert(a);
    assert(visStrain);

	// Compute strain, deviatoric strain, etc.
	const PylithScalar maxwellTime_1 = a[aOff[i_maxwellTime]];
	const PylithScalar maxwellTime_2 = a[aOff[i_maxwellTime] + 1];
	const PylithScalar maxwellTime_3 = a[aOff[i_maxwellTime] + 2];
    const PylithScalar* viscousStrainPrevious_1 = &s[sOff[i_viscousStrainPrevious]];
    const PylithScalar* viscousStrainPrevious_2 = &s[sOff[i_viscousStrainPrevious] + _strainSize];
    const PylithScalar* viscousStrainPrevious_3 = &s[sOff[i_viscousStrainPrevious] + 2*_strainSize];
    const PylithScalar* totalStrainPrevious = &s[sOff[i_totalStrainPrevious]];
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

	const PylithScalar dt = constants[0];
	
	const PylithScalar dq_1 = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime_1);
    const PylithScalar expFac_1 = exp(-dt/maxwellTime_1);
	
	const PylithScalar dq_2 = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime_2);
    const PylithScalar expFac_2 = exp(-dt/maxwellTime_2);
	
	const PylithScalar dq_3 = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime_3);
    const PylithScalar expFac_3 = exp(-dt/maxwellTime_3);

	const PylithScalar strain[6] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        disp_x[2*_dim+2],
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]),
		0.5 * (disp_x[1*_dim+2] + disp_x[2*_dim+1]),
        0.5 * (disp_x[0*_dim+2] + disp_x[2*_dim+0])
    };
    const PylithReal meanStrain = (strain[0] + strain[1] + strain[2])/3.0;

	const PylithScalar devStrain[6] = {
        strain[0] - meanStrain,
        strain[1] - meanStrain,
        strain[2] - meanStrain,
        strain[3],
        strain[4],
        strain[5]
    };

    const PylithReal meanStrainPrevious = (totalStrainPrevious[0] + totalStrainPrevious[1] + totalStrainPrevious[2])/3.0;

	const PylithScalar devStrainPrevious[6] = {
        totalStrainPrevious[0] - meanStrainPrevious,
        totalStrainPrevious[1] - meanStrainPrevious,
        totalStrainPrevious[2] - meanStrainPrevious,
        totalStrainPrevious[3],
        totalStrainPrevious[4],
        totalStrainPrevious[5]
    };

	PylithScalar strainDiff = 0.0;
    for (int iComp = 0; iComp < _strainSize; ++iComp) {
		strainDiff = devStrain[iComp] - devStrainPrevious[iComp];
        visStrain[iComp] = expFac_1*viscousStrainPrevious_1[iComp] + dq_1*strainDiff;
		visStrain[iComp + _strainSize] = expFac_2*viscousStrainPrevious_2[iComp] + dq_2*strainDiff;
		visStrain[iComp + 2*_strainSize] = expFac_3*viscousStrainPrevious_3[iComp] + dq_3*strainDiff;
    } // for
	 
#if 0 // :DEBUG:
	std::cout << "fekernels::IsotropicLinearGenMaxwell3D::updateViscousStrain" << std::endl;
	std::cout << "dim:  " << dim << std::endl;
	std::cout << "numS:  " << numS << std::endl;
	std::cout << "numA:  " << numA << std::endl;
	std::cout << "t:  " << t << std::endl;
	std::cout << "x[0]:  " << x[0] << std::endl;
	std::cout << "x[1]:  " << x[1] << std::endl;
	std::cout << "x[2]:  " << x[2] << std::endl;
	const PylithScalar* disp = &s[sOff[i_disp]];
	const PylithInt i_shearModulusRatio = 4;
	const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio]];
	const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio] + 1];
	const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio] + 2];
	const PylithScalar aa = 1.0e-4;
	const PylithScalar b = 2.5e-4;
	const PylithScalar c = 3.0e-4;
	const PylithScalar d = 3.5e-4;
	const PylithScalar e = 4.0e-4;
	const PylithScalar f = 4.5e-4;
	const PylithScalar g = 9.0e-8;
	const PylithScalar eps = 1.0e-10;
	const PylithScalar viscousStrain_1_xxPredPrevious = 2.0*maxwellTime_1*(exp(t/maxwellTime_1) - 1.0)*
		(shearModulusRatio_1*exp(t*(maxwellTime_2 + maxwellTime_3)/(maxwellTime_2*maxwellTime_3)) +
		 shearModulusRatio_2*exp(t*(maxwellTime_1 + maxwellTime_3)/(maxwellTime_1*maxwellTime_3)) +
		 shearModulusRatio_3*exp(t*(maxwellTime_1 + maxwellTime_2)/(maxwellTime_1*maxwellTime_2)))*
		(aa*(2.0*x[0] - x[1] - x[2]) + b*(2.0*x[1] - x[0]) + d*(2.0*x[2] - x[0]) - e*(x[1] + x[2]))*
		exp(-t*(maxwellTime_1*maxwellTime_2 + maxwellTime_1*maxwellTime_3 + 2.0*maxwellTime_2*maxwellTime_3)/
			(maxwellTime_1*maxwellTime_2*maxwellTime_3))/(3.0*t);
	const PylithScalar viscousStrain_1_xyPredPrevious = maxwellTime_1*(exp(t/maxwellTime_1) - 1.0)*
		(shearModulusRatio_1*exp(t*(maxwellTime_2 + maxwellTime_3)/(maxwellTime_2*maxwellTime_3)) +
		 shearModulusRatio_2*exp(t*(maxwellTime_1 + maxwellTime_3)/(maxwellTime_1*maxwellTime_3)) +
		 shearModulusRatio_3*exp(t*(maxwellTime_1 + maxwellTime_2)/(maxwellTime_1*maxwellTime_2)))*
		(b*x[0] + b*x[1] + c*x[0] + c*x[1] + d*x[2] + e*x[2])*
		exp(-t*(maxwellTime_1*maxwellTime_2 + maxwellTime_1*maxwellTime_3 + 2.0*maxwellTime_2*maxwellTime_3)/
			(maxwellTime_1*maxwellTime_2*maxwellTime_3))/t;
	const PylithScalar viscousStrain_1_xzPredPrevious = maxwellTime_1*(exp(t/maxwellTime_1) - 1.0)*
		(shearModulusRatio_1*exp(t*(maxwellTime_2 + maxwellTime_3)/(maxwellTime_2*maxwellTime_3)) +
		 shearModulusRatio_2*exp(t*(maxwellTime_1 + maxwellTime_3)/(maxwellTime_1*maxwellTime_3)) +
		 shearModulusRatio_3*exp(t*(maxwellTime_1 + maxwellTime_2)/(maxwellTime_1*maxwellTime_2)))*
		(b*x[1] + d*x[0] + d*x[2] + e*x[1] + f*x[0] + f*x[2])*
		exp(-t*(maxwellTime_1*maxwellTime_2 + maxwellTime_1*maxwellTime_3 + 2.0*maxwellTime_2*maxwellTime_3)/
			(maxwellTime_1*maxwellTime_2*maxwellTime_3))/t;
	const PylithScalar viscousStrain_2_xxPredPrevious = 2.0*maxwellTime_2*(exp(t/maxwellTime_2) - 1.0)*
		(shearModulusRatio_1*exp(t*(maxwellTime_2 + maxwellTime_3)/(maxwellTime_2*maxwellTime_3)) +
		 shearModulusRatio_2*exp(t*(maxwellTime_1 + maxwellTime_3)/(maxwellTime_1*maxwellTime_3)) +
		 shearModulusRatio_3*exp(t*(maxwellTime_1 + maxwellTime_2)/(maxwellTime_1*maxwellTime_2)))*
		(aa*(2.0*x[0] - x[1] - x[2]) + b*(2.0*x[1] - x[0]) + d*(2.0*x[2] - x[0]) - e*(x[1] + x[2]))*
		exp(-t*(maxwellTime_1*maxwellTime_2 + 2.0*maxwellTime_1*maxwellTime_3 + maxwellTime_2*maxwellTime_3)/
			(maxwellTime_1*maxwellTime_2*maxwellTime_3))/(3.0*t);
	const PylithScalar viscousStrain_2_xyPredPrevious = maxwellTime_2*(exp(t/maxwellTime_2) - 1.0)*
		(shearModulusRatio_1*exp(t*(maxwellTime_2 + maxwellTime_3)/(maxwellTime_2*maxwellTime_3)) +
		 shearModulusRatio_2*exp(t*(maxwellTime_1 + maxwellTime_3)/(maxwellTime_1*maxwellTime_3)) +
		 shearModulusRatio_3*exp(t*(maxwellTime_1 + maxwellTime_2)/(maxwellTime_1*maxwellTime_2)))*
		(b*x[0] + b*x[1] + c*x[0] + c*x[1] + d*x[2] + e*x[2])*
		exp(-t*(maxwellTime_1*maxwellTime_2 + 2.0*maxwellTime_1*maxwellTime_3 + maxwellTime_2*maxwellTime_3)/
			(maxwellTime_1*maxwellTime_2*maxwellTime_3))/t;
	const PylithScalar viscousStrain_2_xzPredPrevious = maxwellTime_2*(exp(t/maxwellTime_2) - 1.0)*
		(shearModulusRatio_1*exp(t*(maxwellTime_2 + maxwellTime_3)/(maxwellTime_2*maxwellTime_3)) +
		 shearModulusRatio_2*exp(t*(maxwellTime_1 + maxwellTime_3)/(maxwellTime_1*maxwellTime_3)) +
		 shearModulusRatio_3*exp(t*(maxwellTime_1 + maxwellTime_2)/(maxwellTime_1*maxwellTime_2)))*
		(b*x[1] + d*x[0] + d*x[2] + e*x[1] + f*x[0] + f*x[2])*
		exp(-t*(maxwellTime_1*maxwellTime_2 + 2.0*maxwellTime_1*maxwellTime_3 + maxwellTime_2*maxwellTime_3)/
			(maxwellTime_1*maxwellTime_2*maxwellTime_3))/t;
	const PylithScalar totalStrainxxPredPrevious = (2.0*aa*x[0] + 2.0*b*x[1] + 2.0*d*x[2]) *
		(shearModulusRatio_1*exp(-t/maxwellTime_1) +
		 shearModulusRatio_2*exp(-t/maxwellTime_2) +
		 shearModulusRatio_3*exp(-t/maxwellTime_3));
	const PylithScalar totalStrainxyPredPrevious =
		(shearModulusRatio_1*exp(t*(maxwellTime_2 + maxwellTime_3)/(maxwellTime_2*maxwellTime_3)) +
		 shearModulusRatio_2*exp(t*(maxwellTime_1 + maxwellTime_3)/(maxwellTime_1*maxwellTime_3)) +
		 shearModulusRatio_3*exp(t*(maxwellTime_1 + maxwellTime_2)/(maxwellTime_1*maxwellTime_2)))*
		(b*x[0] + b*x[1] + c*x[0] + c*x[1] + d*x[2] + e*x[2])*
		exp(-t*(maxwellTime_1*maxwellTime_2 + maxwellTime_1*maxwellTime_3 + maxwellTime_2*maxwellTime_3)/
			(maxwellTime_1*maxwellTime_2*maxwellTime_3));
	const PylithScalar totalStrainxzPredPrevious =
		(shearModulusRatio_1*exp(t*(maxwellTime_2 + maxwellTime_3)/(maxwellTime_2*maxwellTime_3)) +
		 shearModulusRatio_2*exp(t*(maxwellTime_1 + maxwellTime_3)/(maxwellTime_1*maxwellTime_3)) +
		 shearModulusRatio_3*exp(t*(maxwellTime_1 + maxwellTime_2)/(maxwellTime_1*maxwellTime_2)))*
		(b*x[1] + d*x[0] + d*x[2] + e*x[1] + f*x[0] + f*x[2])*
		exp(-t*(maxwellTime_1*maxwellTime_2 + maxwellTime_1*maxwellTime_3 + maxwellTime_2*maxwellTime_3)/
			(maxwellTime_1*maxwellTime_2*maxwellTime_3));
	const PylithScalar totalStrainxxPred = totalStrainxxPredPrevious + g;
	const PylithScalar totalStrainxyPred =
		(g*exp(t*(1.0/maxwellTime_3 + 1.0/maxwellTime_2 + 1.0/maxwellTime_1))/2.0 +
		 (b*x[0] + c*x[1] + e*x[2])*(shearModulusRatio_1*exp(t*(1.0/maxwellTime_3 + 1.0/maxwellTime_2)) +
									 shearModulusRatio_2*exp(t*(1.0/maxwellTime_3 + 1.0/maxwellTime_1)) +
									 shearModulusRatio_3*exp(t*(1.0/maxwellTime_2 + 1.0/maxwellTime_1))) +
		 (b*x[1] + c*x[0] + d*x[2])*(shearModulusRatio_1*exp(t*(1.0/maxwellTime_3 + 1.0/maxwellTime_2)) +
									 shearModulusRatio_2*exp(t*(1.0/maxwellTime_3 + 1.0/maxwellTime_1)) +
									 shearModulusRatio_3*exp(t*(1.0/maxwellTime_2 + 1.0/maxwellTime_1))))*
		exp(-t*(1.0/maxwellTime_3 + 1.0/maxwellTime_2 + 1.0/maxwellTime_1));
	const PylithScalar totalStrainxzPred =
		(g*exp(t*(1.0/maxwellTime_3 + 1.0/maxwellTime_2 + 1.0/maxwellTime_1))/2.0 +
		 (b*x[1] + d*x[2] + f*x[0])*(shearModulusRatio_1*exp(t*(1.0/maxwellTime_3 + 1.0/maxwellTime_2)) +
									 shearModulusRatio_2*exp(t*(1.0/maxwellTime_3 + 1.0/maxwellTime_1)) +
									 shearModulusRatio_3*exp(t*(1.0/maxwellTime_2 + 1.0/maxwellTime_1))) +
		 (d*x[0] + e*x[1] + f*x[2])*(shearModulusRatio_1*exp(t*(1.0/maxwellTime_3 + 1.0/maxwellTime_2)) +
									 shearModulusRatio_2*exp(t*(1.0/maxwellTime_3 + 1.0/maxwellTime_1)) +
									 shearModulusRatio_3*exp(t*(1.0/maxwellTime_2 + 1.0/maxwellTime_1))))*
		exp(-t*(1.0/maxwellTime_3 + 1.0/maxwellTime_2 + 1.0/maxwellTime_1));
	const PylithScalar viscousStrainxxPrevious_1_diff = viscousStrainPrevious_1[0] - viscousStrain_1_xxPredPrevious;
	const PylithScalar viscousStrainxxPrevious_2_diff = viscousStrainPrevious_2[0] - viscousStrain_2_xxPredPrevious;
	const PylithScalar viscousStrainxyPrevious_1_diff = viscousStrainPrevious_1[3] - viscousStrain_1_xyPredPrevious;
	const PylithScalar viscousStrainxyPrevious_2_diff = viscousStrainPrevious_2[3] - viscousStrain_2_xyPredPrevious;
	const PylithScalar viscousStrainxzPrevious_1_diff = viscousStrainPrevious_1[5] - viscousStrain_1_xzPredPrevious;
	const PylithScalar viscousStrainxzPrevious_2_diff = viscousStrainPrevious_2[5] - viscousStrain_2_xzPredPrevious;
	const PylithScalar totalStrainxxPrevious_diff = totalStrainPrevious[0] - totalStrainxxPredPrevious;
	const PylithScalar totalStrainxyPrevious_diff = totalStrainPrevious[3] - totalStrainxyPredPrevious;
	const PylithScalar totalStrainxzPrevious_diff = totalStrainPrevious[5] - totalStrainxzPredPrevious;
	const PylithScalar totalStrainxx_diff = strain[0] - totalStrainxxPred;
	const PylithScalar totalStrainxy_diff = strain[3] - totalStrainxyPred;
	const PylithScalar totalStrainxz_diff = strain[5] - totalStrainxzPred;
	std::cout << "viscousStrainxxPrevious_1_diff:  " << viscousStrainxxPrevious_1_diff << std::endl;
	std::cout << "viscousStrainxxPrevious_2_diff:  " << viscousStrainxxPrevious_2_diff << std::endl;
	std::cout << "viscousStrainxyPrevious_1_diff:  " << viscousStrainxyPrevious_1_diff << std::endl;
	std::cout << "viscousStrainxyPrevious_2_diff:  " << viscousStrainxyPrevious_2_diff << std::endl;
	std::cout << "viscousStrainxzPrevious_1_diff:  " << viscousStrainxzPrevious_1_diff << std::endl;
	std::cout << "viscousStrainxzPrevious_2_diff:  " << viscousStrainxzPrevious_2_diff << std::endl;
	std::cout << "totalStrainxxPrevious_diff:  " << totalStrainxxPrevious_diff << std::endl;
	std::cout << "totalStrainxyPrevious_diff:  " << totalStrainxyPrevious_diff << std::endl;
	std::cout << "totalStrainxzPrevious_diff:  " << totalStrainxzPrevious_diff << std::endl;
	std::cout << "totalStrainxx_diff:  " << totalStrainxx_diff << std::endl;
	std::cout << "totalStrainxy_diff:  " << totalStrainxy_diff << std::endl;
	std::cout << "totalStrainxz_diff:  " << totalStrainxz_diff << std::endl;
	assert(abs(viscousStrainxxPrevious_1_diff) < eps);
	assert(abs(viscousStrainxxPrevious_2_diff) < eps);
	assert(abs(viscousStrainxyPrevious_1_diff) < eps);
	assert(abs(viscousStrainxyPrevious_2_diff) < eps);
	assert(abs(viscousStrainxzPrevious_1_diff) < eps);
	assert(abs(viscousStrainxzPrevious_2_diff) < eps);
	assert(abs(totalStrainxxPrevious_diff) < eps);
	assert(abs(totalStrainxyPrevious_diff) < eps);
	assert(abs(totalStrainxzPrevious_diff) < eps);
	assert(abs(totalStrainxx_diff) < eps);
	assert(abs(totalStrainxy_diff) < eps);
	assert(abs(totalStrainxz_diff) < eps);
	// std::cout << "viscousStrainxxPrevious_1:  " << viscousStrainPrevious_1[0] << std::endl;
	// std::cout << "viscousStrain_1_xxPredPrevious:  " << viscousStrain_1_xxPredPrevious << std::endl;
	// std::cout << "viscousStrainxxPrevious_2:  " << viscousStrainPrevious_2[0] << std::endl;
	// std::cout << "viscousStrain_2_xxPredPrevious:  " << viscousStrain_2_xxPredPrevious << std::endl;
	// std::cout << "viscousStrainxyPrevious_1:  " << viscousStrainPrevious_1[3] << std::endl;
	// std::cout << "viscousStrain_1_xyPredPrevious:  " << viscousStrain_1_xyPredPrevious << std::endl;
	// std::cout << "viscousStrainxyPrevious_2:  " << viscousStrainPrevious_2[3] << std::endl;
	// std::cout << "viscousStrain_2_xyPredPrevious:  " << viscousStrain_2_xyPredPrevious << std::endl;
	std::cout << "visStrainxx_1:  " << visStrain[0] << std::endl;
	std::cout << "visStrainxx_2:  " << visStrain[0 + _strainSize] << std::endl;
	std::cout << "visStrainxy_1:  " << visStrain[3] << std::endl;
	std::cout << "visStrainxy_2:  " << visStrain[3 + _strainSize] << std::endl;
	//std::cout << "totalStrainxx:  " << disp_x[0] << std::endl;
	//std::cout << "totalStrainxxPred  " << totalStrainxxPred << std::endl;
#endif
	
} // updateViscousStrain


/* End of file */
