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
				     const PylithInt indicesS[],
				     const PylithInt numA,
				     const PylithInt indicesA[],
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
  const PylithInt i_disp = 1;
  const PylithInt i_vel = 1;
  const PylithInt f_disp = indicesS[i_disp];
  const PylithInt f_vel = indicesS[i_vel];
  const PylithScalar* disp_t = &s_t[sOff[f_disp]];
  const PylithScalar* vel = &s[sOff[f_vel]];

  const PylithInt _numA = 0;

  PylithInt i;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);
  assert(indicesS);
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
 * QUESTION: Make a separate f0 entry function if we don't have body
 * forces? Most dynamic simulations don't use body forces.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [density(1), body force(dim)]
 */
PetscErrorCode
pylith_fekernels_f0_ElasticityInertia(const PylithInt dim,
				      const PylithInt numS,
				      const PylithInt indicesS[],
				      const PylithInt numA,
				      const PylithInt indicesA[],
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
  assert(indicesS);
  assert(indicesA);

  pylith_fekernels_f0_Inertia(dim, 1, &indicesS[i_vel], 1, &indicesA[i_density], sOff, aOff, s, s_t, s_x, a, a_t, a_x, x, f0);
  pylith_fekernels_f0_BodyForce(dim, 0, NULL, 1, &indicesA[i_bodyforce], sOff, aOff, s, s_t, s_x, a, a_t, a_x, x, f0);
  
  PYLITH_METHOD_RETURN(0);
} /* f0_ElasticityInertia */
					      

/* ---------------------------------------------------------------------- */
/* f0 function for inertia.
 *
 * Solution fields = [vel(dim)]
 * Auxiliary fields = [density]
 */
PetscErrorCode
pylith_fekernels_f0_Inertia(const PylithInt dim,
			    const PylithInt numS,
			    const PylithInt indicesS[],
			    const PylithInt numA,
			    const PylithInt indicesA[],
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
{ /* f0_Inertia */
  const PylithInt _numS = 1;
  const PylithInt i_vel = 0;
  const PylithInt f_vel = indicesS[i_vel];
  const PylithScalar* acc = &s_t[sOff[f_vel]];

  const PylithInt _numA = 1;
  const PylithInt i_density = 0;
  const PylithInt f_density = indicesA[i_density];
  const PylithScalar density = a[aOff[f_density]];

  PylithInt i;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);
  assert(indicesS);
  assert(indicesA);
  assert(sOff);
  assert(aOff);
  assert(s_t);
  assert(a);
  assert(f0);

  for (i=0; i < dim; ++i) {
    f0[i] += acc[i]*density;
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* f0_Inertia */
					      

/* ---------------------------------------------------------------------- */
/* f0 function for body force.
 *
 * Solution fields = NONE
 * Auxiliary fields = [body force(dim)]
 */
PetscErrorCode
pylith_fekernels_f0_BodyForce(const PylithInt dim,
			      const PylithInt numS,
			      const PylithInt indicesS[],
			      const PylithInt numA,
			      const PylithInt indicesA[],
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
{ /* f0_BodyForce */
  const PylithInt _numS = 0;

  const PylithInt _numA = 1;
  const PylithInt i_bodyforce = 0;
  const PylithInt f_bodyforce = indicesA[i_bodyforce];
  const PylithScalar* bodyforce = &a[aOff[f_bodyforce]];

  PylithInt i;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);
  assert(indicesA);
  assert(aOff);
  assert(a);
  assert(f0);

  for (i=0; i < dim; ++i) {
    f0[i] -= bodyforce[i];
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* f0_BodyForce */


/* ---------------------------------------------------------------------- */
/* f1 entry function for isotropic linear elasticity in 3-D.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = [lambda(1), mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
PetscErrorCode
pylith_fekernels_f1_IsotropicLinearElasticity3D(const PylithInt dim,
						const PylithInt numS,
						const PylithInt indicesS[],
						const PylithInt numA,
						const PylithInt indicesA[],
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
  const PylithInt indicesAVol[3] = { i_lambda, i_istress, i_istrain };

  const PylithInt numADev = 3;
  const PylithInt indicesADev[3] = { i_mu, i_istress, i_istrain };

  PYLITH_METHOD_BEGIN;
  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(indicesS);
  assert(indicesA);

  pylith_fekernels_volumetricStress_IsotropicLinearElasticity3D(dim, 1, &indicesS[i_disp], numAVol, indicesAVol, sOff, aOff, s, s_t, s_x, a, a_t, a_x, x, f1);
  pylith_fekernels_deviatoricStress_IsotropicLinearElasticity3D(dim, 1, &indicesS[i_disp], numADev, indicesADev, sOff, aOff, s, s_t, s_x, a, a_t, a_x, x, f1);
  
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
							      const PylithInt indicesS[],
							      const PylithInt numA,
							      const PylithInt indicesA[],
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
  const PylithInt f_disp = indicesS[i_disp];
  const PylithScalar* disp_x = &s_x[sOff[f_disp]];

  const PylithInt _numA = 3;
  const PylithInt i_lambda = 0;
  const PylithInt i_istress = 1;
  const PylithInt i_istrain = 2;
  const PylithInt f_lambda = indicesA[i_lambda];
  const PylithInt f_istress = indicesA[i_istress];
  const PylithInt f_istrain = indicesA[i_istrain];
  const PylithScalar lambda = a[aOff[f_lambda]];
  const PylithScalar* initialstress = &a[aOff[f_istress]];
  const PylithScalar* initialstrain = &a[aOff[f_istrain]];

  PylithInt i;
  PylithScalar trace = 0;
  PylithScalar meanistress = 0;

  PYLITH_METHOD_BEGIN;
  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(indicesS);
  assert(indicesA);
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
							      const PylithInt indicesS[],
							      const PylithInt numA,
							      const PylithInt indicesA[],
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
  const PylithInt f_disp = indicesS[i_disp];
  const PylithScalar* disp_x = &s_x[sOff[f_disp]];

  const PylithInt _numA = 3;
  const PylithInt i_mu = 0;
  const PylithInt i_istress = 1;
  const PylithInt i_istrain = 2;
  const PylithInt f_mu = indicesA[i_mu];
  const PylithInt f_istress = indicesA[i_istress];
  const PylithInt f_istrain = indicesA[i_istrain];
  const PylithScalar mu = a[aOff[f_mu]];
  const PylithScalar* initialstress = &a[aOff[f_istress]];
  const PylithScalar* initialstrain = &a[aOff[f_istrain]];

  PylithInt i, j;
  PylithScalar meanistress = 0;

  PYLITH_METHOD_BEGIN;
  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(indicesS);
  assert(indicesA);
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


// End of file 
