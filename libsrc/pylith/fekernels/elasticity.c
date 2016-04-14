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
pylith_fekernels_Elasticity_f0v_inertia(const PylithInt dim,
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
{ /* Elasticity_f0v_inertia */
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
} /* Elasticity_f0v_inertia */
					      

/* ---------------------------------------------------------------------- */
/* g0 function for generic elasticity terms (body forces).
 */
void
pylith_fekernels_Elasticity_g0v_bodyforce(const PylithInt dim,
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
{ /* Elasticity_g0v_bodyforce */
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
} /* Elasticity_g0v_bodyforce */
					      

/** Jf0 function for generic elasticity terms (inertia) with implicit time stepping.
 */
void
pylith_fekernels_Elasticity_Jf0vv_inertiaimplicit(const PylithInt dim,
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
{ /* Elasticity_Jf0vv_inertiaimplicit */
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
} /* Elasticity_Jf0vv_inertiaimplicit */




/** Jf0 function for generic elasticity terms (inertia) with explicit time stepping.
 */
void
pylith_fekernels_Elasticity_Jf0vv_inertiaexplicit(const PylithInt dim,
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
{ /* Elasticity_Jf0vv_inertiaexplicit */
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
} /* Elasticity_Jf0vv_inertiaexplicit */


/* End of file */
