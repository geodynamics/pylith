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
/* Main f0 function for with inertia and body forces in 3-D.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [density(1), body force(dim)]
 */
PetscErrorCode
pylith_f0_ElasticityInertia(const PetscInt dim,
			    const PetscInt numS,
			    const PetscInt indicesS[],
			    const PetscInt numA,
			    const PetscInt indicesA[],
			    const PetscInt uOff[],
			    const PetscInt aOff[],
			    const PetscScalar s[],
			    const PetscScalar s_t[],
			    const PetscScalar s_tt[],
			    const PetscScalar s_x[],
			    const PetscScalar a[],
			    const PetscScalar a_x[],
			    const PetscScalar x[],
			    PetscScalar f0[])
{ /* f0_ElasticityInertia */
  const PetscInt _numS = 1;
  const PetscInt _numA = 2;

  const PetscInt i_disp = 0;
  const PetscInt i_density = 0;
  const PetscInt i_bodyforce = 1;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);

  pylith_f0_Inertia(dim, 1, &indicesS[i_disp], 1, &indicesA[i_density], uOff, aOff, s, s_t, s_tt, s_x, a, a_x, x, f0);
  pylith_f0_BodyForce(dim, 1, &indicesS[i_disp], 1, &indicesA[i_bodyforce], uOff, aOff, s, s_t, s_tt, s_x, a, a_x, x, f0);
  
  PYLITH_METHOD_RETURN(0);
} /* f0_ElasticityInertia */
					      

/* ---------------------------------------------------------------------- */
/* f0 function for inertia.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [density]
 */
PetscErrorCode
pylith_f0_Inertia(const PetscInt dim,
		  const PetscInt numS,
		  const PetscInt indicesS[],
		  const PetscInt numA,
		  const PetscInt indicesA[],
		  const PetscInt uOff[],
		  const PetscInt aOff[],
		  const PetscScalar s[],
		  const PetscScalar s_t[],
		  const PetscScalar s_tt[],
		  const PetscScalar s_x[],
		  const PetscScalar a[],
		  const PetscScalar a_x[],
		  const PetscScalar x[],
		  PetscScalar f0[])
{ /* f0_Inertia */
  const PetscInt _numS = 1;
  const PetscInt _numA = 1;
  const PetscInt i_disp = 0;
  const PetscInt i_density = 0;
  const PetscInt f_disp = indicesS[i_disp];
  const PetscInt f_density = indicesA[i_density];

  const PetscScalar acc[] = &u_tt[uOff[f_disp]];
  const PetscScalar density = a[aOff[f_density]];

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);

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
pylith_f0_BodyForce(const PetscInt dim,
		    const PetscInt numS,
		    const PetscInt indicesS[],
		    const PetscInt numA,
		    const PetscInt indicesA[],
		    const PetscInt uOff[],
		    const PetscInt aOff[],
		    const PetscScalar s[],
		    const PetscScalar s_t[],
		    const PetscScalar s_tt[],
		    const PetscScalar s_x[],
		    const PetscScalar a[],
		    const PetscScalar a_x[],
		    const PetscScalar x[],
		    PetscScalar f0[])
{ /* f0_BodyForce */
  const PetscInt _numS = 0;
  const PetscInt _numA = 1;
  const PetscInt i_bodyforce = 0;
  const PetscInt f_bodyforce = indicesA[i_bodyforce];

  const PetscScalar bodyforce[] = a[aOff[f_bodyforce]];

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);

  for (i=0; i < dim; ++i) {
    f0[i] -= bodyforce[i];
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* f0_BodyForce */


/* ---------------------------------------------------------------------- */
/* f1 function for isotropic linear elasticity in 3-D.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [lambda(1), mu(1)]
 */
PetscErrorCode
pylith_f1_IsotropicLinearElasticity3D(const PetscInt dim,
				      const PetscInt numS,
				      const PetscInt indicesS[],
				      const PetscInt numA,
				      const PetscInt indicesA[],
				      const PetscInt uOff[],
				      const PetscInt aOff[],
				      const PetscScalar s[],
				      const PetscScalar s_t[],
				      const PetscScalar s_tt[],
				      const PetscScalar s_x[],
				      const PetscScalar a[],
				      const PetscScalar a_x[],
				      const PetscScalar x[],
				      PetscScalar f1[])
{ /* f1_IsotropicLinearElasticity3D */
  const PetscInt _dim = 3;
  const PetscInt _numS = 1;
  const PetscInt _numA = 2;

  const PetscInt i_disp = 0;
  const PetscInt i_lambda = 0;
  const PetscInt i_mu = 1;

  PYLITH_METHOD_BEGIN;
  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);

  pylith_f1_IsotropicLinearElasticityVolumetricStress(dim, 1, &indicesS[i_disp], 1, &indicesA[i_lambda], uOff, aOff, s, s_t, s_tt, s_x, a, a_x, x, f1);
  pylith_f1_IsotropicLinearElasticityDeviatoricStress(dim, 1, &indicesS[i_disp], 1, &indicesA[i_mu], uOff, aOff, s, s_t, s_tt, s_x, a, a_x, x, f1);
  
  PYLITH_METHOD_RETURN(0);
} /* f1_IsotropicLinearElasticity3D */


/* ---------------------------------------------------------------------- */
/* Calculate volumetic stress for isotropic linear elasticity.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [lambda(1)]
 */
PetscErrorCode
pylith_volumetricStress_IsotropicLinearElasticity(const PetscInt dim,
						  const PetscInt numS,
						  const PetscInt indicesS[],
						  const PetscInt numA,
						  const PetscInt indicesA[],
						  const PetscInt uOff[],
						  const PetscInt aOff[],
						  const PetscScalar s[],
						  const PetscScalar s_t[],
						  const PetscScalar s_tt[],
						  const PetscScalar s_x[],
						  const PetscScalar a[],
						  const PetscScalar a_x[],
						  const PetscScalar x[],
						  PetscScalar stress[])
{ /* volumetricStress_IsotropicLinearElasticity */
  const PetscInt _numS = 1;
  const PetscInt _numA = 1;
  const PetscInt i_disp = 0;
  const PetscInt i_lambda = 0;
  const PetscInt f_disp = indicesS[i_disp];
  const PetscInt f_lambda = indicesA[i_lambda];

  const PetscScalar disp_x[] = &u_x[uOff[f_disp]];
  const PetscScalar lambda = a[aOff[f_lambda]];

  PetscScalar trace = 0;

  PYLITH_METHOD_BEGIN;
  assert(_numS == numS);
  assert(_numA == numA);

  for (i=0; i < dim; ++i) {
    trace += disp_x[i];
  } /* for */
  for (i = 0; i < dim; ++i) {
    stress[i*dim+i] += lambda * trace;
  } /* for */

  PYLITH_METHOD_RETURN(0);
} /* volumetricStress_IsotropicLinearElasticity */


/* ---------------------------------------------------------------------- */
/* Calculate deviatoric stress for isotropic linear elasticity.
 *
 * Solution fields = [disp(dim)]
 * Auxiliary fields = [mu(1)]
 */
PetscErrorCode
pylith_deviatoricStress_IsotropicLinearElasticity(const PetscInt dim,
						  const PetscInt numS,
						  const PetscInt indicesS[],
						  const PetscInt numA,
						  const PetscInt indicesA[],
						  const PetscInt uOff[],
						  const PetscInt aOff[],
						  const PetscScalar s[],
						  const PetscScalar s_t[],
						  const PetscScalar s_tt[],
						  const PetscScalar s_x[],
						  const PetscScalar a[],
						  const PetscScalar a_x[],
						  const PetscScalar x[],
						  PetscScalar stress[])
{ /* deviatoricStress_IsotropicLinearElasticity */
} /* deviatoricStress_IsotropicLinearElasticity */

#endif // pylith_fekernels_elasticity_h


// End of file 
