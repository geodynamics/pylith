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

/* ---------------------------------------------------------------------- */
/* f0 entry function for time evolution of elasticity.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = None
 */
PetscErrorCode
pylith_fekernels_f0_EvolutionDispVel(const PylithInt dim,
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
				     const PylithScalar x[],
				     PylithScalar f0[])
{ /* f0_EvolutionDispVel */
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
} /* f0_EvolutionDispVel */
					      

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

  pylith_fekernels_Inertia(dim, 1, 1, &sOff[i_vel], &aOff[i_density], s, s_t, s_x, a, a_t, a_x, x, f0);
  pylith_fekernels_BodyForce(dim, 0, 1, NULL, &aOff[i_bodyforce], s, s_t, s_x, a, a_t, a_x, x, f0);
  
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

  pylith_fekernels_Inertia(dim, 1, 1, &sOff[i_vel], &aOff[i_density], s, s_t, s_x, a, a_t, a_x, x, f0);
  
  PYLITH_METHOD_RETURN(0);
} /* f0_ElasticityInertia */
					      

/* ---------------------------------------------------------------------- */
/* f0 entry function for body forces (no inertia).
 *
 * QUESTION: Make a separate f0 entry function if we don't have body
 * forces? Most dynamic simulations don't use body forces.
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

  pylith_fekernels_BodyForce(dim, 0, 1, NULL, &aOff[i_bodyforce], s, s_t, s_x, a, a_t, a_x, x, f0);
  
  PYLITH_METHOD_RETURN(0);
} /* f0_ElasticityBodyForce */
					      

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

  pylith_fekernels_volumetricStress_IsotropicLinearElasticity3D(dim, 1, numAVol, &sOff[i_disp], aOffVol, s, s_t, s_x, a, a_t, a_x, x, f1);
  pylith_fekernels_deviatoricStress_IsotropicLinearElasticity3D(dim, 1, numADev, &sOff[i_disp], aOffDev, s, s_t, s_x, a, a_t, a_x, x, f1);
  
  PYLITH_METHOD_RETURN(0);
} /* f1_IsotropicLinearElasticity3D */


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
/* f1 entry function for isotropic linear elasticity in 3-D.
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

  pylith_fekernels_volumetricStress_IsotropicLinearElasticityPlaneStrain(dim, 1, numAVol, &sOff[i_disp], aOffVol, s, s_t, s_x, a, a_t, a_x, x, f1);
  pylith_fekernels_deviatoricStress_IsotropicLinearElasticityPlaneStrain(dim, 1, numADev, &sOff[i_disp], aOffDev, s, s_t, s_x, a, a_t, a_x, x, f1);
  
  PYLITH_METHOD_RETURN(0);
} /* f1_IsotropicLinearElasticityPlaneStrain */


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


// End of file 
