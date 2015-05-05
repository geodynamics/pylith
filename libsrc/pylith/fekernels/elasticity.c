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

#include "pylith/fekernels/elasticity.h"

/* ====================================================================== 
 * Kernels for displacement/velocity.
 *
 * \int_V \vec{\phi}_v \cdot \left(  \vec{v} - \frac{\partial \vec{u}}{\partial t} \right) \, dV
 * ====================================================================== 
 */

/* ---------------------------------------------------------------------- */
/* f0 entry function for time evolution of elasticity.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = None
 */
PetscErrorCode
pylith_fekernels_f0_DispVel(const PylithInt dim,
			    const PylithInt numS,
			    const PylithInt numA,
			    const PylithInt sOff[],
			    const PylithInt aOff[],
			    const PylithScalar s[],
			    const PylithScalar s_t[],
			    const PylithScalar s_x[],
			    const PylithScalar a[],
			    const PylithScalar a_t[],
			    const PylithScalar a_x[],
			    const PylithReal t,
			    const PylithScalar x[],
			    PylithScalar f0[])
{ /* f0_DispVel */
  const PylithInt _numS = 2;
  const PylithInt i_disp = 0;
  const PylithInt i_vel = 1;
  const PylithScalar* disp_t = &s_t[sOff[i_disp]];
  const PylithScalar* vel = &s[sOff[i_vel]];

  const PylithInt _numA = 0;

  PylithInt i;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(s);
  assert(s_t);

  for (i=0; i < dim; ++i) {
    f0[i] += vel[i] - disp_t[i];
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* f0_DispVel */
					      

/* ---------------------------------------------------------------------- */
/* g0_vv entry function for time evolution of elasticity.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = None
 */
PetscErrorCode
pylith_fekernels_g0_vv_DispVel(const PylithInt dim,
			       const PylithInt numS,
			       const PylithInt numA,
			       const PylithInt sOff[],
			       const PylithInt aOff[],
			       const PylithScalar s[],
			       const PylithScalar s_t[],
			       const PylithScalar s_x[],
			       const PylithScalar a[],
			       const PylithScalar a_t[],
			       const PylithScalar a_x[],
			       const PylithReal t,
			       const PylithReal utshift,
			       const PylithScalar x[],
			       PylithScalar g0[])
{ /* g0_vv_DispVel */
  const PylithInt _numS = 2;

  const PylithInt _numA = 0;

  PylithInt i;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);

  for (i=0; i < dim; ++i) {
    g0[i*dim+i] += +1.0;
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* g0_vv_DispVel */
					      

/* ---------------------------------------------------------------------- */
/* g0_vu entry function for time evolution of elasticity.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = None
 */
PetscErrorCode
pylith_fekernels_g0_vu_DispVel(const PylithInt dim,
			       const PylithInt numS,
			       const PylithInt numA,
			       const PylithInt sOff[],
			       const PylithInt aOff[],
			       const PylithScalar s[],
			       const PylithScalar s_t[],
			       const PylithScalar s_x[],
			       const PylithScalar a[],
			       const PylithScalar a_t[],
			       const PylithScalar a_x[],
			       const PylithReal t,
			       const PylithReal utshift,
			       const PylithScalar x[],
			       PylithScalar g0[])
{ /* g0_uv_DispVel */
  const PylithInt _numS = 2;

  const PylithInt _numA = 0;

  PylithInt i;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);

  for (i=0; i < dim; ++i) {
    g0[i*dim+i] += -utshift;
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* g0_vu_DispVel */
					      

/* ====================================================================== 
 * Kernels for inertia and body foces.
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}}{\partial t} - \vec{f} \right) \, dV
 * ====================================================================== 
 */

/* ---------------------------------------------------------------------- */
/* f0 entry function for inertia and body forces.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [density(1), body force(dim)]
 */
PetscErrorCode
pylith_fekernels_f0_ElasticityInertiaBodyForce(const PylithInt dim,
					       const PylithInt numS,
					       const PylithInt numA,
					       const PylithInt sOff[],
					       const PylithInt aOff[],
					       const PylithScalar s[],
					       const PylithScalar s_t[],
					       const PylithScalar s_x[],
					       const PylithScalar a[],
					       const PylithScalar a_t[],
					       const PylithScalar a_x[],
					       const PylithReal t,
					       const PylithScalar x[],
					       PylithScalar f0[])
{ /* f0_ElasticityInertia */
  const PylithInt _numS = 2;
  const PylithInt i_vel = 1;

  const PylithInt _numA = 2;
  const PylithInt i_density = 0;
  const PylithInt i_bodyforce = 1;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);

  pylith_fekernels_Inertia(dim, 1, 1, &sOff[i_vel], &aOff[i_density], s, s_t, s_x, a, a_t, a_x, t, x, f0);
  pylith_fekernels_BodyForce(dim, 0, 1, NULL, &aOff[i_bodyforce], s, s_t, s_x, a, a_t, a_x, t, x, f0);
  
  PYLITH_METHOD_RETURN(0);
} /* f0_ElasticityInertia */
					      

/* ---------------------------------------------------------------------- */
/* f0 entry function for inertia (no body force).
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [density(1)]
 */
PetscErrorCode
pylith_fekernels_f0_ElasticityInertia(const PylithInt dim,
				      const PylithInt numS,
				      const PylithInt numA,
				      const PylithInt sOff[],
				      const PylithInt aOff[],
				      const PylithScalar s[],
				      const PylithScalar s_t[],
				      const PylithScalar s_x[],
				      const PylithScalar a[],
				      const PylithScalar a_t[],
				      const PylithScalar a_x[],
				      const PylithReal t,
				      const PylithScalar x[],
				      PylithScalar f0[])
{ /* f0_ElasticityInertia */
  const PylithInt _numS = 2;
  const PylithInt i_vel = 1;

  const PylithInt _numA = 1;
  const PylithInt i_density = 0;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);

  pylith_fekernels_Inertia(dim, 1, 1, &sOff[i_vel], &aOff[i_density], s, s_t, s_x, a, a_t, a_x, t, x, f0);
  
  PYLITH_METHOD_RETURN(0);
} /* f0_ElasticityInertia */
					      

/* ---------------------------------------------------------------------- */
/* f0 entry function for body forces (no inertia).
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [body force(dim)]
 */
PetscErrorCode
pylith_fekernels_f0_ElasticityBodyForce(const PylithInt dim,
					const PylithInt numS,
					const PylithInt numA,
					const PylithInt sOff[],
					const PylithInt aOff[],
					const PylithScalar s[],
					const PylithScalar s_t[],
					const PylithScalar s_x[],
					const PylithScalar a[],
					const PylithScalar a_t[],
					const PylithScalar a_x[],
					const PylithReal t,
					const PylithScalar x[],
					PylithScalar f0[])
{ /* f0_ElasticityBodyForce */
  const PylithInt _numS = 2;

  const PylithInt _numA = 1;
  const PylithInt i_bodyforce = 0;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);

  pylith_fekernels_BodyForce(dim, 0, 1, NULL, &aOff[i_bodyforce], s, s_t, s_x, a, a_t, a_x, t, x, f0);
  
  PYLITH_METHOD_RETURN(0);
} /* f0_ElasticityBodyForce */
					      

/* ---------------------------------------------------------------------- */
/* g0_uv entry function for inertia.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [density(1)]
 */
PetscErrorCode
pylith_fekernels_g0_uv_ElasticityInertia(const PylithInt dim,
					 const PylithInt numS,
					 const PylithInt numA,
					 const PylithInt sOff[],
					 const PylithInt aOff[],
					 const PylithScalar s[],
					 const PylithScalar s_t[],
					 const PylithScalar s_x[],
					 const PylithScalar a[],
					 const PylithScalar a_t[],
					 const PylithScalar a_x[],
					 const PylithReal t,
					 const PylithReal utshift,
					 const PylithScalar x[],
					 PylithScalar g0[])
{ /* g0_uv_ElasticityInertia */
  const PylithInt _numS = 2;

  const PylithInt _numA = 0;
  const PylithInt i_density = 0;
  const PylithScalar density = a[aOff[i_density]];

  PylithInt i;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);

  for (i=0; i < dim; ++i) {
    g0[i*dim+i] += density*utshift;
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* g0_uv_ElasticityInertia */
					      

/* ---------------------------------------------------------------------- */
/* f0 function for inertia.
 *
 * Solution fields = [vel(dim)]
 * Auxiliary fields = [density]
 */
PetscErrorCode
pylith_fekernels_Inertia(const PylithInt dim,
			 const PylithInt numS,
			 const PylithInt numA,
			 const PylithInt sOff[],
			 const PylithInt aOff[],
			 const PylithScalar s[],
			 const PylithScalar s_t[],
			 const PylithScalar s_tt[],
			 const PylithScalar s_x[],
			 const PylithScalar a[],
			 const PylithScalar a_x[],
			 const PylithReal t,
			 const PylithScalar x[],
			 PylithScalar f0[])
{ /* Inertia */
  const PylithInt _numS = 1;
  const PylithInt i_vel = 0;
  const PylithScalar* acc = &s_t[sOff[i_vel]];

  const PylithInt _numA = 1;
  const PylithInt i_density = 0;
  const PylithScalar density = a[aOff[i_density]];

  PylithInt i;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);
  assert(s_t);
  assert(a);
  assert(f0);

  for (i=0; i < dim; ++i) {
    f0[i] += acc[i]*density;
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* Inertia */
					      

/* ---------------------------------------------------------------------- */
/* f0 function for body force.
 *
 * Solution fields = NONE
 * Auxiliary fields = [body force(dim)]
 */
PetscErrorCode
pylith_fekernels_BodyForce(const PylithInt dim,
			   const PylithInt numS,
			   const PylithInt numA,
			   const PylithInt sOff[],
			   const PylithInt aOff[],
			   const PylithScalar s[],
			   const PylithScalar s_t[],
			   const PylithScalar s_x[],
			   const PylithScalar a[],
			   const PylithScalar a_t[],
			   const PylithScalar a_x[],
			   const PylithReal t,
			   const PylithScalar x[],
			   PylithScalar f0[])
{ /* BodyForce */
  const PylithInt _numS = 0;

  const PylithInt _numA = 1;
  const PylithInt i_bodyforce = 0;
  const PylithScalar* bodyforce = &a[aOff[i_bodyforce]];

  PylithInt i;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);
  assert(a);
  assert(f0);

  for (i=0; i < dim; ++i) {
    f0[i] -= bodyforce[i];
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* BodyForce */


/* ====================================================================== 
 * Kernels for stress.
 *
 * \int_V \nabla \vec{\phi}_u : \tensor{\sigma} \, dV
 * ====================================================================== 
 */

/* ---------------------------------------------------------------------- */
/* f1 entry function for isotropic linear elasticity in 3-D.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [lambda(1), mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
PetscErrorCode
pylith_fekernels_f1_IsotropicLinearElasticity3D(const PylithInt dim,
						const PylithInt numS,
						const PylithInt numA,
						const PylithInt sOff[],
						const PylithInt aOff[],
						const PylithScalar s[],
						const PylithScalar s_t[],
						const PylithScalar s_x[],
						const PylithScalar a[],
						const PylithScalar a_t[],
						const PylithScalar a_x[],
						const PylithReal t,
						const PylithScalar x[],
						PylithScalar f1[])
{ /* f1_IsotropicLinearElasticity3D */
  const PylithInt _dim = 3;

  const PylithInt _numS = 2;
  const PylithInt i_disp = 0;

  const PylithInt _numA = 4;
  const PylithInt i_lambda = 0;
  const PylithInt i_mu = 1;
  const PylithInt i_istress = 2;
  const PylithInt i_istrain = 3;

  const PylithInt numAVol = 3;
  const PylithInt aOffVol[3] = { aOff[i_lambda], aOff[i_istress], aOff[i_istrain] };

  const PylithInt numADev = 3;
  const PylithInt aOffDev[3] = { aOff[i_mu], aOff[i_istress], aOff[i_istrain] };

  PYLITH_METHOD_BEGIN;
  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);

  pylith_fekernels_volumetricStress_IsotropicLinearElasticity3D(dim, 1, numAVol, &sOff[i_disp], aOffVol, s, s_t, s_x, a, a_t, a_x, t, x, f1);
  pylith_fekernels_deviatoricStress_IsotropicLinearElasticity3D(dim, 1, numADev, &sOff[i_disp], aOffDev, s, s_t, s_x, a, a_t, a_x, t, x, f1);
  
  PYLITH_METHOD_RETURN(0);
} /* f1_IsotropicLinearElasticity3D */


/* ---------------------------------------------------------------------- */
/* g3_uu entry function for isotropic linear elasticity in 3-D.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [lambda(1), mu(1)]
 */
PetscErrorCode
pylith_fekernels_g3_uu_IsotropicLinearElasticity3D(const PylithInt dim,
						   const PylithInt numS,
						   const PylithInt numA,
						   const PylithInt sOff[],
						   const PylithInt aOff[],
						   const PylithScalar s[],
						   const PylithScalar s_t[],
						   const PylithScalar s_x[],
						   const PylithScalar a[],
						   const PylithScalar a_t[],
						   const PylithScalar a_x[],
						   const PylithReal t,
						   const PylithReal utshift,
						   const PylithScalar x[],
						   PylithScalar g3[])
{ /* g3_uu_IsotropicLinearElasticity3D */
  const PylithInt _dim = 3;

  const PylithInt _numS = 2;

  const PylithInt _numA = 2;
  const PylithInt i_lambda = 0;
  const PylithInt i_mu = 1;

  const PylithScalar lambda = a[aOff[i_lambda]];
  const PylithScalar mu = a[aOff[i_mu]];
  const PylithScalar mu2 = 2.0*mu;
  const PylithScalar lambda2mu = lambda + mu2;

  PYLITH_METHOD_BEGIN;
  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);

  g3[ 0] = lambda2mu; // C1111
  g3[ 1] = lambda; // C1122
  g3[ 2] = lambda; // C1133
  g3[ 3] = 0; // C1112
  g3[ 4] = 0; // C1123
  g3[ 5] = 0; // C1113

  g3[ 6] = lambda; // C2211
  g3[ 7] = lambda2mu; // C2222
  g3[ 8] = lambda; // C2233
  g3[ 9] = 0; // C2212
  g3[10] = 0; // C2223
  g3[11] = 0; // C2213

  g3[12] = lambda; // C3311
  g3[13] = lambda; // C3322
  g3[14] = lambda2mu; // C3333
  g3[15] = 0; // C3312
  g3[16] = 0; // C3323
  g3[17] = 0; // C3313

  g3[18] = 0; // C1211
  g3[19] = 0; // C1222
  g3[20] = 0; // C1233
  g3[21] = mu2; // C1212
  g3[22] = 0; // C1223
  g3[23] = 0; // C1213

  g3[24] = 0; // C2311
  g3[25] = 0; // C2322
  g3[26] = 0; // C2333
  g3[27] = 0; // C2312
  g3[28] = mu2; // C2323
  g3[29] = 0; // C2313

  g3[30] = 0; // C1311
  g3[31] = 0; // C1322
  g3[32] = 0; // C1333
  g3[33] = 0; // C1312
  g3[34] = 0; // C1323
  g3[35] = mu2; // C1313
  
  PYLITH_METHOD_RETURN(0);
} /* g3_uu_IsotropicLinearElasticity3D */


/* ---------------------------------------------------------------------- */
/* Calculate volumetic stress for isotropic linear elasticity.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [lambda(1), initialstress(dim*dim), initialstrain(dim*dim)]
 *
 * We compute the stress relative to a reference stress/strain state.
 *
 * stress_ij - initialstress_ij = lambda * (strain_kk - initialstrain_kk) * delta_ij + 2*mu * (strain_ij - initialstrain_ij)
 */
PetscErrorCode
pylith_fekernels_volumetricStress_IsotropicLinearElasticity3D(const PylithInt dim,
							      const PylithInt numS,
							      const PylithInt numA,
							      const PylithInt sOff[],
							      const PylithInt aOff[],
							      const PylithScalar s[],
							      const PylithScalar s_t[],
							      const PylithScalar s_x[],
							      const PylithScalar a[],
							      const PylithScalar a_t[],
							      const PylithScalar a_x[],
							      const PylithReal t,
							      const PylithScalar x[],
							      PylithScalar stress[])
{ /* volumetricStress_IsotropicLinearElasticity3D */
  const PylithInt _dim = 3;

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
  PylithScalar trace = 0;
  PylithScalar meanistress = 0;

  PYLITH_METHOD_BEGIN;
  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(stress);

  for (i=0; i < _dim; ++i) {
    trace += disp_x[i] - initialstrain[i*_dim+i];
    meanistress += initialstress[i*_dim+i];
  } /* for */
  meanistress /= (PylithScalar)_dim;
  for (i = 0; i < _dim; ++i) {
    stress[i*_dim+i] += lambda * trace + meanistress;
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* volumetricStress_IsotropicLinearElasticity */


/* ---------------------------------------------------------------------- */
/* Calculate deviatoric stress for isotropic linear elasticity.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
PetscErrorCode
pylith_fekernels_deviatoricStress_IsotropicLinearElasticity3D(const PylithInt dim,
							      const PylithInt numS,
							      const PylithInt numA,
							      const PylithInt sOff[],
							      const PylithInt aOff[],
							      const PylithScalar s[],
							      const PylithScalar s_t[],
							      const PylithScalar s_x[],
							      const PylithScalar a[],
							      const PylithScalar a_t[],
							      const PylithScalar a_x[],
							      const PylithReal t,
							      const PylithScalar x[],
							      PylithScalar stress[])
{ /* deviatoricStress_IsotropicLinearElasticity3D */
  const PylithInt _dim = 3;

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

  PYLITH_METHOD_BEGIN;
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
      stress[i*_dim+j] += mu * (disp_x[i*_dim+j] + disp_x[j*_dim+i] - initialstrain[i*_dim+j]) + initialstress[i*_dim+j];
    } /* for */
    stress[i*_dim+i] -= meanistress;
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* deviatoricStress_IsotropicLinearElasticity */


/* ---------------------------------------------------------------------- */
/* f1 entry function for 2-D plane strain isotropic linear elasticity.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [lambda(1), mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
PetscErrorCode
pylith_fekernels_f1_IsotropicLinearElasticityPlaneStrain(const PylithInt dim,
							 const PylithInt numS,
							 const PylithInt numA,
							 const PylithInt sOff[],
							 const PylithInt aOff[],
							 const PylithScalar s[],
							 const PylithScalar s_t[],
							 const PylithScalar s_x[],
							 const PylithScalar a[],
							 const PylithScalar a_t[],
							 const PylithScalar a_x[],
							 const PylithReal t,
							 const PylithScalar x[],
							 PylithScalar f1[])
{ /* f1_IsotropicLinearElasticityPlaneStrain */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2;
  const PylithInt i_disp = 0;

  const PylithInt _numA = 4;
  const PylithInt i_lambda = 0;
  const PylithInt i_mu = 1;
  const PylithInt i_istress = 2;
  const PylithInt i_istrain = 3;

  const PylithInt numAVol = 3;
  const PylithInt aOffVol[3] = { aOff[i_lambda], aOff[i_istress], aOff[i_istrain] };

  const PylithInt numADev = 3;
  const PylithInt aOffDev[3] = { aOff[i_mu], aOff[i_istress], aOff[i_istrain] };

  PYLITH_METHOD_BEGIN;
  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);

  pylith_fekernels_volumetricStress_IsotropicLinearElasticityPlaneStrain(dim, 1, numAVol, &sOff[i_disp], aOffVol, s, s_t, s_x, a, a_t, a_x, t, x, f1);
  pylith_fekernels_deviatoricStress_IsotropicLinearElasticityPlaneStrain(dim, 1, numADev, &sOff[i_disp], aOffDev, s, s_t, s_x, a, a_t, a_x, t, x, f1);
  
  PYLITH_METHOD_RETURN(0);
} /* f1_IsotropicLinearElasticityPlaneStrain */


/* ---------------------------------------------------------------------- */
/* g3_uu entry function for 2-D plane strain isotropic linear elasticity.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [lambda(1), mu(1)]
 */
PetscErrorCode
pylith_fekernels_g3_uu_IsotropicLinearElasticityPlaneStrain(const PylithInt dim,
							    const PylithInt numS,
							    const PylithInt numA,
							    const PylithInt sOff[],
							    const PylithInt aOff[],
							    const PylithScalar s[],
							    const PylithScalar s_t[],
							    const PylithScalar s_x[],
							    const PylithScalar a[],
							    const PylithScalar a_t[],
							    const PylithScalar a_x[],
							    const PylithReal t,
							    const PylithReal utshift,
							    const PylithScalar x[],
							    PylithScalar g3[])
{ /* g3_uu_IsotropicLinearElasticityPlaneStrain */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2;

  const PylithInt _numA = 4;
  const PylithInt i_lambda = 0;
  const PylithInt i_mu = 1;

  const PylithScalar lambda = a[aOff[i_lambda]];
  const PylithScalar mu = a[aOff[i_mu]];

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar lambda2mu = lambda + mu2;
   
  PYLITH_METHOD_BEGIN;
  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);

  g3[0] = lambda2mu; // C1111
  g3[1] = lambda; // C1122
  g3[2] = 0; // C1112

  g3[3] = lambda; // C2211
  g3[4] = lambda2mu; // C2222
  g3[5] = 0; // C2212

  g3[6] = 0; // C1211
  g3[7] = 0; // C1222
  g3[8] = mu2; // C1212
  
  PYLITH_METHOD_RETURN(0);
} /* g3_uu_IsotropicLinearElasticityPlaneStrain */


/* ---------------------------------------------------------------------- */
/* Calculate volumetic stress for 2-D plane strain isotropic linear elasticity.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [lambda(1), initialstress(dim*dim), initialstrain(dim*dim)]
 *
 * We compute the stress relative to a reference stress/strain state.
 *
 * stress_ij - initialstress_ij = lambda * (strain_kk - initialstrain_kk) * delta_ij + 2*mu * (strain_ij - initialstrain_ij)
 */
PetscErrorCode
pylith_fekernels_volumetricStress_IsotropicLinearElasticityPlaneStrain(const PylithInt dim,
								       const PylithInt numS,
								       const PylithInt numA,
								       const PylithInt sOff[],
								       const PylithInt aOff[],
								       const PylithScalar s[],
								       const PylithScalar s_t[],
								       const PylithScalar s_x[],
								       const PylithScalar a[],
								       const PylithScalar a_t[],
								       const PylithScalar a_x[],
								       const PylithReal t,
								       const PylithScalar x[],
								       PylithScalar stress[])
{ /* volumetricStress_IsotropicLinearElasticityPlaneStrain */
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
  PylithScalar trace = 0;
  PylithScalar meanistress = 0;

  PYLITH_METHOD_BEGIN;
  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(stress);

  for (i=0; i < _dim; ++i) {
    trace += disp_x[i] - initialstrain[i*_dim+i];
    meanistress += initialstress[i*_dim+i];
  } /* for */
  meanistress /= (PylithScalar)_dim;
  for (i = 0; i < _dim; ++i) {
    stress[i*_dim+i] += lambda * trace + meanistress;
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* volumetricStress_IsotropicLinearElasticityPlaneStrain */


/* ---------------------------------------------------------------------- */
/* Calculate deviatoric stress for 2-D plane strain isotropic linear elasticity.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
PetscErrorCode
pylith_fekernels_deviatoricStress_IsotropicLinearElasticityPlaneStrain(const PylithInt dim,
								       const PylithInt numS,
								       const PylithInt numA,
								       const PylithInt sOff[],
								       const PylithInt aOff[],
								       const PylithScalar s[],
								       const PylithScalar s_t[],
								       const PylithScalar s_x[],
								       const PylithScalar a[],
								       const PylithScalar a_t[],
								       const PylithScalar a_x[],
								       const PylithReal t,
								       const PylithScalar x[],
								       PylithScalar stress[])
{ /* deviatoricStress_IsotropicLinearElasticityPlaneStrain */
  const PylithInt _dim = 3;

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

  PYLITH_METHOD_BEGIN;
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
      stress[i*_dim+j] += mu * (disp_x[i*_dim+j] + disp_x[j*_dim+i] - initialstrain[i*_dim+j]) + initialstress[i*_dim+j];
    } /* for */
    stress[i*_dim+i] -= meanistress;
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* deviatoricStress_IsotropicLinearElasticityPlaneStrain */


/* End of file */
