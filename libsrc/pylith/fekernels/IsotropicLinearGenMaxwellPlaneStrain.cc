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

#include "pylith/fekernels/IsotropicLinearGenMaxwellPlaneStrain.hh"
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/Viscoelastic.hh" // USES Viscoelastic kernels
#include "pylith/fekernels/ElasticityPlaneStrain.hh" // USES ElasticityPlaneStrain kernels

#include <cassert> // USES assert()
#include <cmath> // USES exp()

/* ======================================================================
 * Kernels for isotropic, linear Generalized Maxwell viscoelastic plane strain.
 * ======================================================================
 */

// ----------------------------------------------------------------------
// g0 function for isotropic linear Generalized Maxwell viscoelastic plane strain with both
// gravity and body forces.
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::g0v_gravbodyforce(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_gravityField = 13;
    const PylithInt i_bodyForce = 14;

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 15);
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
// g0 function for isotropic linear Maxwell viscoelastic plane strain with just gravity.
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::g0v_grav(const PylithInt dim,
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
    const PylithInt _dim = 2;

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_gravityField = 13;

    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 14);
    assert(aOff);
    assert(aOff_x);

    pylith::fekernels::Elasticity::g0v_grav(_dim, _numS, numAGrav,
                                            NULL, NULL, NULL, NULL, NULL,
                                            aOffGrav, aOffGrav_x, a, a_t, a_x,
                                            t, x, numConstants, constants, g0);
} // g0v_grav


// ----------------------------------------------------------------------
// g0 function for isotropic linear Maxwell viscoelastic plane strain with only body forces.
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::g0v_bodyforce(const PylithInt dim,
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
    const PylithInt _dim = 2;

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    // Incoming auxiliary fields.
    const PylithInt i_bodyForce = 13;

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 13);
    assert(aOff);
    assert(aOff_x);

    pylith::fekernels::Elasticity::g0v_bodyforce(_dim, _numS, numABody,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 aOffBody, aOffBody_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, g0);
} // g0v_bodyforce


// ----------------------------------------------------------------------
// g1 function for isotropic linear Maxwell plane strain WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::g1v(const PylithInt dim,
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
    const PylithInt i_maxwellTime_1 = 3;
    const PylithInt i_maxwellTime_2 = 4;
    const PylithInt i_maxwellTime_3 = 5;
    const PylithInt i_shearModulusRatio_1 = 6;
    const PylithInt i_shearModulusRatio_2 = 7;
    const PylithInt i_shearModulusRatio_3 = 8;
    const PylithInt i_totalStrain = 9;
    const PylithInt i_viscousStrain_1 = 10;
    const PylithInt i_viscousStrain_2 = 11;
    const PylithInt i_viscousStrain_3 = 12;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 13);
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

    const PylithInt numADev = 11; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[11] = {
      aOff[i_shearModulus], aOff[i_maxwellTime_1], aOff[i_maxwellTime_2], aOff[i_maxwellTime_3],
      aOff[i_shearModulusRatio_1], aOff[i_shearModulusRatio_2], aOff[i_shearModulusRatio_3],
      aOff[i_totalStrain], aOff[i_viscousStrain_1], aOff[i_viscousStrain_2], aOff[i_viscousStrain_3]
    };
    const PylithInt aOffDev_x[11] = {
      aOff_x[i_shearModulus], aOff_x[i_maxwellTime_1], aOff_x[i_maxwellTime_2],
      aOff_x[i_maxwellTime_3], aOff_x[i_shearModulusRatio_1], aOff_x[i_shearModulusRatio_2],
      aOff_x[i_shearModulusRatio_3], aOff_x[i_totalStrain], aOff_x[i_viscousStrain_1],
      aOff_x[i_viscousStrain_2], aOff_x[i_viscousStrain_3]
    };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0}; // Full stress tensor

    pylith::fekernels::ElasticityPlaneStrain::meanStress(_dim, _numS, numAMean,
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
// g1 function for isotropic linear Maxwell viscoelastic plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::g1v_refstate(const PylithInt dim,
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
    const PylithInt i_maxwellTime_1 = 3;
    const PylithInt i_maxwellTime_2 = 4;
    const PylithInt i_maxwellTime_3 = 5;
    const PylithInt i_shearModulusRatio_1 = 6;
    const PylithInt i_shearModulusRatio_2 = 7;
    const PylithInt i_shearModulusRatio_3 = 8;
    const PylithInt i_totalStrain = 9;
    const PylithInt i_viscousStrain_1 = 10;
    const PylithInt i_viscousStrain_2 = 11;
    const PylithInt i_viscousStrain_3 = 12;
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

    const PylithInt numADev = 13; // Pass shear modulus, Maxwell times, shear modulus ratios,
                                  // total strain, viscous strains,
                                  // reference stress, and reference strain.
    const PylithInt aOffDev[13] = {
      aOff[i_shearModulus], aOff[i_maxwellTime_1], aOff[i_maxwellTime_2], aOff[i_maxwellTime_3],
      aOff[i_shearModulusRatio_1], aOff[i_shearModulusRatio_2], aOff[i_shearModulusRatio_3],
      aOff[i_totalStrain], aOff[i_viscousStrain_1], aOff[i_viscousStrain_2],
      aOff[i_viscousStrain_3], aOff[i_rstress], aOff[i_rstrain]
    };
    const PylithInt aOffDev_x[13] = {
      aOff_x[i_shearModulus], aOff_x[i_maxwellTime_1], aOff_x[i_maxwellTime_2],
      aOff_x[i_maxwellTime_3], aOff_x[i_shearModulusRatio_1], aOff_x[i_shearModulusRatio_2],
      aOff_x[i_shearModulusRatio_3], aOff_x[i_totalStrain], aOff_x[i_viscousStrain_1],
      aOff_x[i_viscousStrain_2], aOff_x[i_viscousStrain_3], aOff_x[i_rstress], aOff_x[i_rstrain]
    };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0};

    pylith::fekernels::ElasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean,
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
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::Jg3vu(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_maxwellTime_1 = 3;
    const PylithInt i_maxwellTime_2 = 4;
    const PylithInt i_maxwellTime_3 = 5;
    const PylithInt i_shearModulusRatio_1 = 6;
    const PylithInt i_shearModulusRatio_2 = 7;
    const PylithInt i_shearModulusRatio_3 = 8;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 13);
    assert(aOff);
    assert(a);
    assert(Jg3);
    assert(numConstants == 1);
    assert(constants);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar maxwellTime_1 = a[aOff[i_maxwellTime_1]];
    const PylithScalar maxwellTime_2 = a[aOff[i_maxwellTime_2]];
    const PylithScalar maxwellTime_3 = a[aOff[i_maxwellTime_3]];
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio_1]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio_2]];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio_3]];
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

       0:  j0000 = C1111 = 1.0*bulkModulus + 2.0*shearModulus*(0.666666666666667*dq_1*shearModulusRatio_1 + 0.666666666666667*dq_2*shearModulusRatio_2 + 0.666666666666667*dq_3*shearModulusRatio_3 + 0.666666666666667*shearModulusRatio_0)
       1:  j0001 = C1112 = 0
       2:  j0010 = C1211 = 0
       3:  j0011 = C1212 = 2.0*shearModulus*(0.5*dq_1*shearModulusRatio_1 + 0.5*dq_2*shearModulusRatio_2 + 0.5*dq_3*shearModulusRatio_3 + 0.5*shearModulusRatio_0)
       4:  j0100 = C1121 = 0
       5:  j0101 = C1122 = 1.0*bulkModulus + 2.0*shearModulus*(-0.333333333333333*dq_1*shearModulusRatio_1 - 0.333333333333333*dq_2*shearModulusRatio_2 - 0.333333333333333*dq_3*shearModulusRatio_3 - 0.333333333333333*shearModulusRatio_0)
       6:  j0110 = C1221 = 2.0*shearModulus*(0.5*dq_1*shearModulusRatio_1 + 0.5*dq_2*shearModulusRatio_2 + 0.5*dq_3*shearModulusRatio_3 + 0.5*shearModulusRatio_0)
       7:  j0111 = C1222 = 0
       8:  j1000 = C2111 = 0
       9:  j1001 = C2112 = 2.0*shearModulus*(0.5*dq_1*shearModulusRatio_1 + 0.5*dq_2*shearModulusRatio_2 + 0.5*dq_3*shearModulusRatio_3 + 0.5*shearModulusRatio_0)
       10:  j1010 = C2211 = 1.0*bulkModulus + 2.0*shearModulus*(-0.333333333333333*dq_1*shearModulusRatio_1 - 0.333333333333333*dq_2*shearModulusRatio_2 - 0.333333333333333*dq_3*shearModulusRatio_3 - 0.333333333333333*shearModulusRatio_0)
       11:  j1011 = C2212 = 0
       12:  j1100 = C2121 = 2.0*shearModulus*(0.5*dq_1*shearModulusRatio_1 + 0.5*dq_2*shearModulusRatio_2 + 0.5*dq_3*shearModulusRatio_3 + 0.5*shearModulusRatio_0)
       13:  j1101 = C2122 = 0
       14:  j1110 = C2221 = 0
       15:  j1111 = C2222 = 1.0*bulkModulus + 2.0*shearModulus*(0.666666666666667*dq_1*shearModulusRatio_1 + 0.666666666666667*dq_2*shearModulusRatio_2 + 0.666666666666667*dq_3*shearModulusRatio_3 + 0.666666666666667*shearModulusRatio_0)
    */

    /* Nonzero Jacobian entries. */
    Jg3[0] -=  C1111; /* j0000 */
    Jg3[3] -=  C1212; /* j0011 */
    Jg3[5] -=  C1122; /* j0101 */
    Jg3[6] -=  C1212; /* j0110 */
    Jg3[9] -=  C1212; /* j1001 */
    Jg3[10] -=  C1122; /* j1010 */
    Jg3[12] -=  C1212; /* j1100 */
    Jg3[15] -=  C1111; /* j1111 */

} // Jg3vu


// ----------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = 2*shearModulus*visStrain_ij
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::deviatoricStress(const PylithInt dim,
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_maxwellTime_1 = 1;
    const PylithInt i_maxwellTime_2 = 2;
    const PylithInt i_maxwellTime_3 = 3;
    const PylithInt i_shearModulusRatio_1 = 4;
    const PylithInt i_shearModulusRatio_2 = 5;
    const PylithInt i_shearModulusRatio_3 = 6;
    const PylithInt i_totalStrain = 7;
    const PylithInt i_viscousStrain_1 = 8;
    const PylithInt i_viscousStrain_2 = 9;
    const PylithInt i_viscousStrain_3 = 10;

    assert(_dim == dim);
    assert(1 == numS);
    assert(11 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithInt _numS = 1; // Number passed on to visStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    // Viscous strain for first Maxwell element.
    const PylithInt numAVis_1 = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis_1[3] = { aOff[i_maxwellTime_1], aOff[i_totalStrain],
				     aOff[i_viscousStrain_1] };
    const PylithInt aOffVis_1_x[3] = { aOff_x[i_maxwellTime_1], aOff_x[i_totalStrain],
				       aOff_x[i_viscousStrain_1] };

    PylithScalar visStrainTpdt_1[4] = {0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    computeViscousStrain(_dim, _numS, numAVis_1, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis_1,
			 aOffVis_1_x, a, a_t, a_x, t, x, numConstants, constants, visStrainTpdt_1);

    // Viscous strain for second Maxwell element.
    const PylithInt numAVis_2 = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis_2[3] = { aOff[i_maxwellTime_2], aOff[i_totalStrain],
				     aOff[i_viscousStrain_2] };
    const PylithInt aOffVis_2_x[3] = { aOff_x[i_maxwellTime_2], aOff_x[i_totalStrain],
				       aOff_x[i_viscousStrain_2] };

    PylithScalar visStrainTpdt_2[4] = {0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    computeViscousStrain(_dim, _numS, numAVis_2, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis_2,
			 aOffVis_2_x, a, a_t, a_x, t, x, numConstants, constants, visStrainTpdt_2);

    // Viscous strain for third Maxwell element.
    const PylithInt numAVis_3 = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis_3[3] = { aOff[i_maxwellTime_3], aOff[i_totalStrain],
				     aOff[i_viscousStrain_3] };
    const PylithInt aOffVis_3_x[3] = { aOff_x[i_maxwellTime_3], aOff_x[i_totalStrain],
				       aOff_x[i_viscousStrain_3] };

    PylithScalar visStrainTpdt_3[4] = {0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    computeViscousStrain(_dim, _numS, numAVis_3, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis_3,
			 aOffVis_3_x, a, a_t, a_x, t, x, numConstants, constants, visStrainTpdt_3);

    // Shear modulus ratio factors.
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio_1]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio_2]];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio_3]];
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2
      - shearModulusRatio_3;

    // Deviatoric strains.
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar strainTpdt[4] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        0.0,
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0])
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0;

    const PylithScalar devStrainTpdt[4] = {
        strainTpdt[0] - meanStrainTpdt,
        strainTpdt[1] - meanStrainTpdt,
        strainTpdt[2] - meanStrainTpdt,
        strainTpdt[3]
    };

    // Stresses.
    stress[0] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[0] +
				       shearModulusRatio_1 * visStrainTpdt_1[0] + 
				       shearModulusRatio_2 * visStrainTpdt_2[0] + 
				       shearModulusRatio_3 * visStrainTpdt_3[0]); // sigma_11
    stress[1] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[3] +
				       shearModulusRatio_1 * visStrainTpdt_1[3] + 
				       shearModulusRatio_2 * visStrainTpdt_2[3] + 
				       shearModulusRatio_3 * visStrainTpdt_3[3]); // sigma_12
    stress[2] += stress[1];                                                       //sigma_21
    stress[3] += 2.0 * shearModulus * (shearModulusRatio_0 * devStrainTpdt[1] +
				       shearModulusRatio_1 * visStrainTpdt_1[1] + 
				       shearModulusRatio_2 * visStrainTpdt_2[1] + 
				       shearModulusRatio_3 * visStrainTpdt_3[1]); // sigma_22
	
} // deviatoricStress


// ----------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITH reference stress and reference strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = devrefstress_ij + 2*shearModulus*(visstrain_ij - devrefstrain_ij)
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::deviatoricStress_refstate(const PylithInt dim,
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_maxwellTime_1 = 1;
    const PylithInt i_maxwellTime_2 = 2;
    const PylithInt i_maxwellTime_3 = 3;
    const PylithInt i_shearModulusRatio_1 = 4;
    const PylithInt i_shearModulusRatio_2 = 5;
    const PylithInt i_shearModulusRatio_3 = 6;
    const PylithInt i_totalStrain = 7;
    const PylithInt i_viscousStrain_1 = 8;
    const PylithInt i_viscousStrain_2 = 9;
    const PylithInt i_viscousStrain_3 = 10;
    const PylithInt i_rstress = 11;
    const PylithInt i_rstrain = 12;

    assert(_dim == dim);
    assert(1 == numS);
    assert(13 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // sigma_11, sigma_22, sigma_33, sigma_12
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // epsilon_11, epsilon_22, epsilon_33, epsilon_12

    const PylithInt _numS = 1; // Number passed on to visStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    // Viscous strain for first Maxwell element.
    const PylithInt numAVis_1 = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis_1[3] = { aOff[i_maxwellTime_1], aOff[i_totalStrain],
				     aOff[i_viscousStrain_1] };
    const PylithInt aOffVis_1_x[3] = { aOff_x[i_maxwellTime_1], aOff_x[i_totalStrain],
				       aOff_x[i_viscousStrain_1] };

    PylithScalar visStrainTpdt_1[4] = {0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    computeViscousStrain(_dim, _numS, numAVis_1, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis_1,
			 aOffVis_1_x, a, a_t, a_x, t, x, numConstants, constants, visStrainTpdt_1);

    // Viscous strain for second Maxwell element.
    const PylithInt numAVis_2 = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis_2[3] = { aOff[i_maxwellTime_2], aOff[i_totalStrain],
				     aOff[i_viscousStrain_2] };
    const PylithInt aOffVis_2_x[3] = { aOff_x[i_maxwellTime_2], aOff_x[i_totalStrain],
				       aOff_x[i_viscousStrain_2] };

    PylithScalar visStrainTpdt_2[4] = {0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    computeViscousStrain(_dim, _numS, numAVis_2, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis_2,
			 aOffVis_2_x, a, a_t, a_x, t, x, numConstants, constants, visStrainTpdt_2);

    // Viscous strain for third Maxwell element.
    const PylithInt numAVis_3 = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis_3[3] = { aOff[i_maxwellTime_3], aOff[i_totalStrain],
				     aOff[i_viscousStrain_3] };
    const PylithInt aOffVis_3_x[3] = { aOff_x[i_maxwellTime_3], aOff_x[i_totalStrain],
				       aOff_x[i_viscousStrain_3] };

    PylithScalar visStrainTpdt_3[4] = {0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    computeViscousStrain(_dim, _numS, numAVis_3, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis_3,
			 aOffVis_3_x, a, a_t, a_x, t, x, numConstants, constants, visStrainTpdt_3);

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

    // Shear modulus ratio factors.
    const PylithScalar shearModulusRatio_1 = a[aOff[i_shearModulusRatio_1]];
    const PylithScalar shearModulusRatio_2 = a[aOff[i_shearModulusRatio_2]];
    const PylithScalar shearModulusRatio_3 = a[aOff[i_shearModulusRatio_3]];
    const PylithScalar shearModulusRatio_0 = 1.0 - shearModulusRatio_1 - shearModulusRatio_2
      - shearModulusRatio_3;

    // Deviatoric strains.
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar strainTpdt[4] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        0.0,
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0])
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0;

    const PylithScalar devStrainTpdt[4] = {
        strainTpdt[0] - meanStrainTpdt,
        strainTpdt[1] - meanStrainTpdt,
        strainTpdt[2] - meanStrainTpdt,
        strainTpdt[3]
    };

    // Compute stress components -- note that we are including reference deviatoric stress for now.
    // This may need to be removed after testing.
    const PylithReal twomu = 2.0*shearModulus;
    stress[0] += devRefStress[0] + twomu * (shearModulusRatio_0 * devStrainTpdt[0] +
					    shearModulusRatio_1 * visStrainTpdt_1[0] + 
					    shearModulusRatio_2 * visStrainTpdt_2[0] + 
					    shearModulusRatio_3 * visStrainTpdt_3[0] -
					    devRefStrain[0]); // sigma_11
    stress[1] += devRefStress[3] + twomu * (shearModulusRatio_0 * devStrainTpdt[3] +
					    shearModulusRatio_1 * visStrainTpdt_1[3] + 
					    shearModulusRatio_2 * visStrainTpdt_2[3] + 
					    shearModulusRatio_3 * visStrainTpdt_3[3] -
					    devRefStrain[3]); // sigma_12
    stress[2] += stress[1];                                   // sigma_21
    stress[3] += devRefStress[1] + twomu * (shearModulusRatio_0 * devStrainTpdt[1] +
					    shearModulusRatio_1 * visStrainTpdt_1[1] + 
					    shearModulusRatio_2 * visStrainTpdt_2[1] + 
					    shearModulusRatio_3 * visStrainTpdt_3[1] -
					    devRefStrain[3]); // sigma_22

} // deviatoricStress_refstate


// ----------------------------------------------------------------------
/* Calculate viscous strain for a Maxwell viscoelastic material.
 *
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::computeViscousStrain(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = 0;
    const PylithInt i_totalStrain = 1;
    const PylithInt i_viscousStrain = 2;

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

    const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];
    const PylithScalar* viscousStrain = &a[aOff[i_viscousStrain]];

    const PylithScalar dt = constants[0];

    const PylithScalar dq = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime);
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

} // computeViscousStrain


// ----------------------------------------------------------------------
/* Update total strain for a Maxwell viscoelastic material.
 *
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::updateTotalStrain(const PylithInt dim,
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
									   PylithScalar totalStrainTpdt[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(13 >= numA && 17 <= numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(totalStrainTpdt);

    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    totalStrainTpdt[0] = disp_x[0*_dim+0];
    totalStrainTpdt[1] = disp_x[1*_dim+1];
    totalStrainTpdt[2] = 0.0;
    totalStrainTpdt[3] = 0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]);

} // updateTotalStrain


// ----------------------------------------------------------------------
/* Update viscous strain for first Maxwell element.
 *
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::updateViscousStrain_1(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_totalStrain = 9;
    const PylithInt i_viscousStrain = 10;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 13);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(visStrainTpdt);

    const PylithInt _numS = 1; // Number passed on to statevars kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; // Number passed on to viscous strain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_totalStrain], aOff[i_viscousStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_totalStrain],
                                     aOff_x[i_viscousStrain] };

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis,
			 aOffVis_x, a, a_t, a_x, t, x, numConstants, constants, visStrainTpdt);


} // updateViscousStrain_1


// ----------------------------------------------------------------------
/* Update viscous strain for a Maxwell viscoelastic material.
 *
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::updateViscousStrain_2(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = 4;
    const PylithInt i_totalStrain = 9;
    const PylithInt i_viscousStrain = 11;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 13);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(visStrainTpdt);

    const PylithInt _numS = 1; // Number passed on to statevars kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; // Number passed on to viscous strain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_totalStrain], aOff[i_viscousStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_totalStrain],
                                     aOff_x[i_viscousStrain] };

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis,
			 aOffVis_x, a, a_t, a_x, t, x, numConstants, constants, visStrainTpdt);


} // updateViscousStrain_2


// ----------------------------------------------------------------------
/* Update viscous strain for a Maxwell viscoelastic material.
 *
 */
void
pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::updateViscousStrain_3(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = 5;
    const PylithInt i_totalStrain = 9;
    const PylithInt i_viscousStrain = 12;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 13);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(visStrainTpdt);

    const PylithInt _numS = 1; // Number passed on to statevars kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; // Number passed on to viscous strain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_totalStrain], aOff[i_viscousStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_totalStrain],
                                     aOff_x[i_viscousStrain] };

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffVis,
			 aOffVis_x, a, a_t, a_x, t, x, numConstants, constants, visStrainTpdt);


} // updateViscousStrain_3


/* End of file */
