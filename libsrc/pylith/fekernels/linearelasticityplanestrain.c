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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_f0(const PylithInt dim,
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
{ /* IsotropicLinearElasticityPlaneStrain_f0 */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2;

  const PylithInt _numA = 1;
  const PylithInt i_density = 0;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA <= numA);
  assert(aOff);
  assert(aOff_x);

  pylith_fekernels_Elasticity_f0_inertia(_dim, _numS, _numA, sOff, sOff_x, s, s_t, s_x, &aOff[i_density], &aOff_x[i_density], a, a_t, a_x, t, x, f0);
} /* IsotropicLinearElasticityPlaneStrain_f0 */


/* ---------------------------------------------------------------------- */
/** g0 function for isotropic linear elasticity plane strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g0(const PylithInt dim,
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
{ /* IsotropicLinearElasticityPlaneStrain_g0 */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2;

  const PylithInt _numA = 1;
  const PylithInt i_bodyforce = 3;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA <= numA);
  assert(aOff);
  assert(aOff_x);

  pylith_fekernels_Elasticity_g0_bodyforce(_dim, _numS, _numA, sOff, sOff_x, s, s_t, s_x, &aOff[i_bodyforce], &aOff_x[i_bodyforce], a, a_t, a_x, t, x, g0);
} /* IsotropicLinearElasticityPlaneStrain_g0 */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear elasticity plane strain WITHOUT initial stress and strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1(const PylithInt dim,
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
{ /* IsotropicLinearElasticityPlaneStrain_g1 */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2;
  const PylithInt i_disp = 0;

  const PylithInt _numA = 3;
  const PylithInt i_mu = 1;
  const PylithInt i_lambda = 2;

  const PylithInt numAVol = 1;
  const PylithInt aOffVol[1] = { aOff[i_lambda] };
  const PylithInt aOffVol_x[1] = { aOff_x[i_lambda] };

  const PylithInt numADev = 1;
  const PylithInt aOffDev[1] = { aOff[i_mu] };
  const PylithInt aOffDev_x[1] = { aOff_x[i_mu] };

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA <= numA);
  assert(sOff);
  assert(sOff_x);
  assert(aOff);
  assert(aOff_x);

  pylith_fekernels_IsotropicLinearElasticityPlaneStrain_volumetricStress(_dim, 1, numAVol, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x, aOffVol, aOffVol_x, a, a_t, a_x, t, x, g1);
  pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress(_dim, 1, numADev, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x, aOffDev, aOffDev_x, a, a_t, a_x, t, x, g1);
} /* IsotropicLinearElasticityPlaneStrain_g1 */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear elasticity plane strain with initial stress and strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1_initstate(const PylithInt dim,
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
{ /* IsotropicLinearElasticityPlaneStrain_g1_initstate */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2;
  const PylithInt i_disp = 0;

  const PylithInt _numA = 5;
  const PylithInt i_mu = 1;
  const PylithInt i_lambda = 2;
  const PylithInt i_istress = numA-2;
  const PylithInt i_istrain = numA-1;

  const PylithInt numAVol = 3;
  const PylithInt aOffVol[3] = { aOff[i_lambda], aOff[i_istress], aOff[i_istrain] };
  const PylithInt aOffVol_x[3] = { aOff_x[i_lambda], aOff_x[i_istress], aOff_x[i_istrain] };

  const PylithInt numADev = 3;
  const PylithInt aOffDev[3] = { aOff[i_mu], aOff[i_istress], aOff[i_istrain] };
  const PylithInt aOffDev_x[3] = { aOff_x[i_mu], aOff_x[i_istress], aOff_x[i_istrain] };

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA <= numA);
  assert(sOff);
  assert(sOff_x);
  assert(aOff);
  assert(aOff_x);

  pylith_fekernels_IsotropicLinearElasticityPlaneStrain_volumetricStress_initstate(_dim, 1, numAVol, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x, aOffVol, aOffVol_x, a, a_t, a_x, t, x, g1);
  pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress_initstate(_dim, 1, numADev, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x, aOffDev, aOffDev_x, a, a_t, a_x, t, x, g1);
} /* IsotropicLinearElasticityPlaneStrain_g1_initstate */


/* ---------------------------------------------------------------------- */
/** Jf0 function for isotropoc linear elasticity plane strain with implicit time stepping.
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0_uv_implicit(const PylithInt dim,
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
{ /* IsotropicLinearElasticityPlaneStrain_Jf0_uv_implicit */
  assert(0);
} /* IsotropicLinearElasticityPlaneStrain_Jf0_uv_implicit */


/** Jf0 function for isotropoc linear elasticity plane strain with explicit time stepping.
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0_uv_explicit(const PylithInt dim,
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
{ /* IsotropicLinearElasticityPlaneStrain_Jf0_uv_explicit */
  assert(0);
} /* IsotropicLinearElasticityPlaneStrain_Jf0_uv_explicit*/


/* ---------------------------------------------------------------------- */
/* g3_uu entry function for 2-D plane strain isotropic linear elasticity.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jg3_uu(const PylithInt dim,
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
{ /* IsotropicLinearElasticityPlaneStrain_g3_uu */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2;

  const PylithInt _numA = 4;
  const PylithInt i_lambda = 1;
  const PylithInt i_mu = 0;

  const PylithScalar lambda = a[aOff[i_lambda]];
  const PylithScalar mu = a[aOff[i_mu]];

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar lambda2mu = lambda + mu2;
   
  const PylithReal C1111 = lambda2mu;
  const PylithReal C2222 = lambda2mu;
  const PylithReal C1122 = lambda;
  const PylithReal C1212 = mu2;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);

  /* j(f,g,df,dg) = C(f,df,g,dg)

     0: j0000 = C1111
     1: j0001 = C1112
     4: j0100 = C1121, symmetry C1112
     5: j0101 = C1122

     2: j0010 = C1211
     3: j0011 = C1212
     6: j0110 = C1221, symmetry C1212
     7: j0111 = C1222
  
     8: j1000 = C2111
     9: j1001 = C2112, symmetry C1212
    12: j1100 = C2121, symmetry C1212
    13: j1101 = C2122, symmetry C1222

    10: j1010 = C2211, symmetry C1122
    11: j1011 = C2212, symmetry C1222
    14: j1110 = C2221, symmetry C1222
    15: j1111 = C2222
  */

  Jg3[ 0] += C1111; /* j0000 */
  Jg3[ 3] += C1212; /* j0011 */
  Jg3[ 5] += C1122; /* j0101 */
  Jg3[ 6] += C1212; /* j0110, C1221 */
  Jg3[ 9] += C1212; /* j1001, C2112 */
  Jg3[10] += C1122; /* j1010, C2211 */
  Jg3[12] += C1212; /* j1100, C2121 */
  Jg3[15] += C2222; /* j1111 */
} /* IsotropicLinearElasticityPlaneStrain_Jg3_uu */


/* ---------------------------------------------------------------------- */
/* Calculate volumetric stress for 2-D plane strain isotropic linear elasticity WITHOUT initial stress and initial strain.
 *
 * -volStress_ij = -lambda * strain_kk * delta_ij
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_volumetricStress(const PylithInt dim,
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
{ /* IsotropicLinearElasticityPlaneStrain_volumetricStress */
  const PylithInt _dim = 2;

  const PylithInt _numS = 1;
  const PylithInt i_disp = 0;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];

  const PylithInt _numA = 1;
  const PylithInt i_lambda = 0;
  const PylithScalar lambda = a[aOff[i_lambda]];

  PylithInt i;
  PylithScalar strainTrace = 0;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(stress);

  for (i = 0; i < _dim; ++i) {
    strainTrace += disp_x[i*_dim+i];
  } /* for */
  for (i = 0; i < _dim; ++i) {
    stress[i*_dim+i] -= lambda * strainTrace;
  } /* for */
} /* volumetricStress_IsotropicLinearElasticityPlaneStrain */


/* ---------------------------------------------------------------------- */
/* Calculate volumetric stress for 2-D plane strain isotropic linear elasticity WITH initial stress and strain.
 *
 * We compute the stress relative to a reference stress/strain state.
 *
 * -volStress_ij + meanInitialStress_ij = -lambda * (strain_kk + initialstrain_kk) * delta_ij
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_volumetricStress_initstate(const PylithInt dim,
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
{ /* IsotropicLinearElasticityPlaneStrain_initState_volumetricStress */
  const PylithInt _dim = 2;

  const PylithInt _numS = 1;
  const PylithInt i_disp = 0;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];

  const PylithInt _numA = 3;
  const PylithInt i_lambda = 0;
  const PylithInt i_istress = 1;
  const PylithInt i_istrain = 2;
  const PylithScalar lambda = a[aOff[i_lambda]];
  const PylithScalar* initialstress = &a[aOff[i_istress]];
  const PylithScalar* initialstrain = &a[aOff[i_istrain]];

  PylithInt i;
  PylithScalar strainTrace = 0;
  PylithScalar meanistress = 0;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(stress);

  for (i=0; i < _dim; ++i) {
    strainTrace += disp_x[i*_dim+i] - initialstrain[i*_dim+i];
    meanistress += initialstress[i*_dim+i];
  } /* for */
  meanistress /= (PylithScalar)_dim;
  for (i = 0; i < _dim; ++i) {
    stress[i*_dim+i] -= lambda * strainTrace + meanistress;
  } /* for */
} /* IsotropicLinearElasticityPlaneStrain_initState_volumetricStress */


/* ---------------------------------------------------------------------- */
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT initial stress and strain.
 *
 * -devStress_ij = -2.0*mu * strain_ij
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

  const PylithInt _numS = 1;
  const PylithInt i_disp = 0;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];

  const PylithInt _numA = 1;
  const PylithInt i_mu = 0;
  const PylithScalar mu = a[aOff[i_mu]];

  PylithInt i, j;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(stress);

  for (i=0; i < _dim; ++i) {
    for (j=0; j < _dim; ++j) {
      stress[i*_dim+j] -= mu * (disp_x[i*_dim+j] + disp_x[j*_dim+i]);
    } /* for */
  } /* for */
} /* IsotropicLinearElasticityPlaneStrain_deviatoricStress */


/* ---------------------------------------------------------------------- */
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITH initial stress and initial strain.
 *
 * -devStress_ij + (initialStress_ij + meanInitialStress_ij * delta_ij) = -2.0*mu * (strain_ij + initialstrain_ij)
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress_initstate(const PylithInt dim,
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
{ /* deviatoricStress_IsotropicLinearElasticityPlaneStrain_initstate */
  const PylithInt _dim = 2;

  const PylithInt _numS = 1;
  const PylithInt i_disp = 0;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];

  const PylithInt _numA = 3;
  const PylithInt i_mu = 0;
  const PylithInt i_istress = 1;
  const PylithInt i_istrain = 2;
  const PylithScalar mu = a[aOff[i_mu]];
  const PylithScalar* initialstress = &a[aOff[i_istress]];
  const PylithScalar* initialstrain = &a[aOff[i_istrain]];

  PylithInt i, j;
  PylithScalar meanistress = 0;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(stress);

  for (i=0; i < _dim; ++i) {
    meanistress += initialstress[i*_dim+i];
  } /* for */
  meanistress /= (PylithScalar)_dim;
  for (i=0; i < _dim; ++i) {
    for (j=0; j < _dim; ++j) {
      stress[i*_dim+j] -= mu * (disp_x[i*_dim+j] + disp_x[j*_dim+i] - initialstrain[i*_dim+j]) + initialstress[i*_dim+j];
    } /* for */
    stress[i*_dim+i] += meanistress;
  } /* for */
} /* deviatoricStress_IsotropicLinearElasticityPlaneStrain_initstate */


/* ====================================================================== 
 * Kernels for incompressible elasticity volume integral.
 *
 * \int_V \tensor{S}:\nabla \vec{\phi}_u - p \tensor{I}:\vec{\nabla} 
 * \vec{\phi}_u - \vec{\phi}_u \cdot \vec{f} \, dV
 * ====================================================================== 
 */

/* ---------------------------------------------------------------------- */
/* f1 entry function for 2-D plane strain incompressible isotropic linear
 * elasticity.
 *
 * Solution fields = [disp(dim), pres]
 * Auxiliary fields = [lambda(1), mu(1), initialstress(dim*dim),
 *                     initialstrain(dim*dim)]
 */
void
pylith_fekernels_f1_IncomprUIntegralPlaneStrain(const PylithInt dim,
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
						PylithScalar f1[])
{ /* f1_IncomprUIntegralPlaneStrain */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2;
  const PylithInt i_disp = 0;
  const PylithInt i_pres = 1;
  const PylithScalar pres = s[sOff[i_pres]];

  const PylithInt _numA = 4;

  const PylithInt i_mu = 1;
  const PylithInt i_istress = 2;
  const PylithInt i_istrain = 3;

  const PylithInt numADev = 3;
  const PylithInt aOffDev[3] = { aOff[i_mu], aOff[i_istress], aOff[i_istrain] };
  const PylithInt aOffDev_x[3] = { aOff_x[i_mu], aOff_x[i_istress],
				   aOff_x[i_istrain] };

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(sOff_x);
  assert(aOff);
  assert(aOff_x);

  // NOTE:  At present, only the deviatoric initial strains and stresses are
  // being used. Not sure at present how to incorporate the volumetric parts.
  pylith_fekernels_deviatoricStress_IsotropicLinearElasticityIncomprPlaneStrain(
					dim, 1, numADev, &sOff[i_disp],
					&sOff_x[i_disp], s, s_t, s_x, aOffDev,
					aOffDev_x, a, a_t, a_x, t, x, f1);
  PylithInt i;

  for (i=0; i < _dim; ++i) {
    f1[i*_dim+i] -= pres;
  } /* for */
} /* f1_IncomprUIntegralPlaneStrain */


/* ---------------------------------------------------------------------- */
/* Calculate deviatoric stress for 2-D plane strain incompressible isotropic
 * linear elasticity.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_deviatoricStress_IsotropicLinearElasticityIncomprPlaneStrain(
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
{ /* deviatoricStress_IsotropicLinearElasticityIncomprPlaneStrain */
  const PylithInt _dim = 2;

  const PylithInt _numS = 1;
  const PylithInt i_disp = 0;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];

  const PylithInt _numA = 3;
  const PylithInt i_mu = 0;
  const PylithInt i_istress = 1;
  const PylithInt i_istrain = 2;
  const PylithScalar mu = a[aOff[i_mu]];
  const PylithScalar* initialstress = &a[aOff[i_istress]];
  const PylithScalar* initialstrain = &a[aOff[i_istrain]];

  PylithInt i, j;
  PylithScalar meanistress = 0;
  PylithScalar meanistrain = 0;
  PylithScalar meanstrain = 0;
  PylithScalar devistrain[_dim*_dim];
  PylithScalar devstrain[_dim*_dim];

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(stress);

  // Need to make sure whether all the decomposition into deviatoric parts is
  // necessary.
  for (i=0; i < _dim; ++i) {
    meanistress += initialstress[i*_dim+i];
    meanistrain += initialstrain[i*_dim+i];
    meanstrain += disp_x[i];
  } /* for */
  meanistress /= (PylithScalar)_dim;
  meanistrain /= (PylithScalar)_dim;
  meanstrain /= (PylithScalar)_dim;
  for (i=0; i < _dim; ++i) {
    for (j=0; j < _dim; ++j) {
      devistrain[i*_dim+j] = initialstrain[i*_dim+j];
      devstrain[i*_dim+j] = 0.5 * (disp_x[i*_dim+j] + disp_x[j*dim+i]);
    } /* for */
    devistrain[i*_dim+i] -= meanistrain;
    devstrain[i*_dim+i] -= meanstrain;
  } /* for */
  
  for (i=0; i < _dim; ++i) {
    for (j=0; j < _dim; ++j) {
      stress[i*_dim+j] += mu * (devstrain[i*_dim+j] - devistrain[i*_dim+j]) +
	initialstress[i*_dim+j];
    } /* for */
    stress[i*_dim+i] -= meanistress;
  } /* for */
} /* deviatoricStress_IsotropicLinearElasticityIncomprPlaneStrain */

/* ---------------------------------------------------------------------- */
/* g3_uu entry function for 2-D plane strain incompressible isotropic linear
 * elasticity.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [lambda(1), mu(1), initialstress(dim*dim),
 *                     initialstrain(dim*dim)]
 */
void
pylith_fekernels_g3_uu_IncomprIsotropicLinearElasticityPlaneStrain(
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
					PylithScalar g3[])
{ /* g3_uu_IncomprIsotropicLinearElasticityPlaneStrain */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2;

  const PylithInt _numA = 4;
  const PylithInt i_mu = 1;

  const PylithScalar mu = a[aOff[i_mu]];

  const PylithReal C1111 = 5.0 * mu/3.0;
  const PylithReal C2222 = C1111;
  const PylithReal C1122 = -mu/3.0;
  const PylithReal C1212 = mu;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);

  /* g(f,g,df,dg) = C^{\prime}(f,df,g,dg) - \frac{1}{6}C^{\prime}(f,df,h,h) 
                                          \delta(g,dg)

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

  g3[ 0] += C1111; /* g0000 */
  g3[ 3] += C1212; /* g0011 */
  g3[ 5] += C1122; /* g0101 */
  g3[ 6] += C1212; /* g0110, C1221 */
  g3[ 9] += C1212; /* g1001, C2112 */
  g3[10] += C1122; /* g1010, C2211 */
  g3[12] += C1212; /* g1100, C2121 */
  g3[15] += C2222; /* g1111 */
} /* g3_uu_IncomprIsotropicLinearElasticityPlaneStrain */
					      

/* ---------------------------------------------------------------------- */
/* g2_uv entry function for 2-D plane strain incompressible isotropic linear
 * elasticity.
 *
 * Solution fields = [disp(dim), pres]
 * Auxiliary fields = None
 */
void
pylith_fekernels_g2_uv_IncomprIsotropicLinearElasticityPlaneStrain(
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
					PylithScalar g2[])
{ /* g2_uv_IncomprIsotropicLinearElasticityPlaneStrain */
  const PylithInt _numS = 2;

  const PylithInt _numA = 0;

  PylithInt i;

  assert(_numS == numS);
  assert(_numA == numA);

  for (i=0; i < dim; ++i) {
    g2[i*dim+i] += +1.0;
  } /* for */
} /* g2_uv_IncomprIsotropicLinearElasticityPlaneStrain */

