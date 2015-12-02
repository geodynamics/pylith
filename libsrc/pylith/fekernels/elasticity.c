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
 * Generic elasticity kernels for inertia and body forces.
 *
 * Solution fields: [disp(dim), vel(dim)]
 *
 * ====================================================================== 
 */

/* ---------------------------------------------------------------------- */
/* f0 function for generic elasticity terms (inertia).
 */
void
pylith_fekernels_Elasticity_f0_inertia(const PylithInt dim,
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
{ /* Elasticity_f0_inertia */
  const PylithInt _numS = 2;
  const PylithInt i_vel = 1;
  const PylithScalar* vel_t = &s_t[sOff[i_vel]]; /* acceleration */

  const PylithInt _numA = 1;
  const PylithInt i_density = 0;
  const PylithScalar density = a[aOff[i_density]];

  PylithInt i;

  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(s_t);
  assert(aOff);
  assert(a);
  
  for (i=0; i < dim; ++i) {
    f0[i] += vel_t[i] * density;
  } /* for */
} /* Elasticity_f0_inertia */
					      

/* ---------------------------------------------------------------------- */
/* g0 function for generic elasticity terms (body forces).
 */
void
pylith_fekernels_Elasticity_g0_bodyforce(const PylithInt dim,
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
{ /* Elasticity_g0_bodyforce */
  const PylithInt _numS = 2;

  const PylithInt _numA = 1;
  const PylithInt i_bodyforce = 0;
  const PylithScalar* bodyforce = &a[aOff[i_bodyforce]];

  PylithInt i;

  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);
  assert(aOff[i_bodyforce] >= 0);
  assert(a);

  for (i=0; i < dim; ++i) {
    g0[i] += bodyforce[i];
  } /* for */
} /* Elasticity_g0_bodyforce */
					      

/** Jf0 function for generic elasticity terms (inertia) with implicit time stepping.
 */
void
pylith_fekernels_Elasticity_Jf0_uv_inertiaimplicit(const PylithInt dim,
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
{ /* Elasticity_Jf0_uv_inertiaimplicit */
  const PylithInt _numS = 2;

  const PylithInt _numA = 1;
  const PylithInt i_density = 0;
  const PylithScalar* density = &a[aOff[i_density]];

  PylithInt i, j;

  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);
  assert(aOff[i_density] >= 0);
  assert(a);

  for (i=0; i < dim; ++i) {
    for (j=0; j < dim; ++j) {
      Jf0[i*dim+j] += utshift * density[i];
    } /* for */
  } /* for */
} /* Elasticity_Jf0_uv_inertiaimplicit */




/** Jf0 function for generic elasticity terms (inertia) with explicit time stepping.
 */
void
pylith_fekernels_Elasticity_Jf0_uv_inertiaexplicit(const PylithInt dim,
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
{ /* Elasticity_Jf0_uv_inertiaexplicit */
  const PylithInt _numS = 2;
  const PylithInt i_disp = 0;
  const PylithInt i_vel = 1;

  const PylithInt _numA = 1;
  const PylithInt i_density = 0;
  const PylithScalar density = a[aOff[i_density]];


  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);
  assert(a);

  Jf0[i_disp*_numS+i_vel] += density;
} /* Elasticity_Jf0_uv_inertiaexplicit */




/* ====================================================================== 
 * Kernels for incompressibility volume integral.
 *
 * \int_V \phi_p \left( \vec{\nabla} \cdot \vec{u} + \frac{p}{K} \right) \, dV
 * ====================================================================== 
 */

/* ---------------------------------------------------------------------- */
/* f0 entry function for incompressibility volume integral.
 *
 * Solution fields = [disp(dim), pres]
 * Auxiliary fields = [lambda, mu]
 */
void
pylith_fekernels_IncompElasticity_f0(const PylithInt dim,
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
{ /* IncompElasticity_f0 */
  const PylithInt _numS = 2;
  const PylithInt i_disp = 0;
  const PylithInt i_pres = 1;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];
  const PylithScalar pres = s[sOff[i_pres]];

  const PylithInt _numA = 2;
  const PylithInt i_lambda = 0;
  const PylithInt i_mu = 1;
  const PylithScalar lambda = a[aOff[i_lambda]];
  const PylithScalar mu = a[aOff[i_mu]];

  const PylithScalar bulkModulus = lambda + 2.0 * mu/3.0;

  PylithInt i;
  PylithScalar volStrain = 0;

  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(s);
  assert(s_x);

  for (i=0; i < dim; ++i) {
    volStrain += disp_x[i];
  } /* for */

  f0[0] += volStrain + pres/bulkModulus;
} /* f0_IncomprPIntegral */
					      

/* ---------------------------------------------------------------------- */
/* g0_vv entry function for incompressibility volume integral.
 *
 * Solution fields = [disp(dim), pres]
 * Auxiliary fields = [lambda, mu]
 */
void
pylith_fekernels_g0_vv_IncomprPIntegral(const PylithInt dim,
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
					PylithScalar g0[])
{ /* g0_vv_IncomprPIntegral */
  const PylithInt _numS = 2;

  const PylithInt _numA = 2;
  const PylithInt i_lambda = 0;
  const PylithInt i_mu = 1;
  const PylithScalar lambda = a[aOff[i_lambda]];
  const PylithScalar mu = a[aOff[i_mu]];

  const PylithScalar bulkModulus = lambda + 2.0 * mu/3.0;

  assert(_numS == numS);
  assert(_numA == numA);

  g0[0] += 1.0/bulkModulus;
} /* g0_vv_IncomprPIntegral */
					      

/* ---------------------------------------------------------------------- */
/* g2_vu entry function for incompressibility volume integral.
 *
 * Solution fields = [disp(dim), pres]
 * Auxiliary fields = None
 */
void
pylith_fekernels_g2_vu_IncomprPIntegral(const PylithInt dim,
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
{ /* g2_vu_IncomprPIntegral */
  const PylithInt _numS = 2;

  const PylithInt _numA = 0;

  PylithInt i;

  assert(_numS == numS);
  assert(_numA == numA);

  for (i=0; i < dim; ++i) {
    g2[i] += 1.0;
  } /* for */
} /* g2_vu_IncomprPIntegral */
					      
					      


/* End of file */
