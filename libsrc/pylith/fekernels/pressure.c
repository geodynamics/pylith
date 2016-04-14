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

#include "pylith/fekernels/pressure.h"

/* ====================================================================== 
 * Kernels for pressure volume integral.
 *
 * Solution fields = [disp(dim), pres]
 * Auxiliary fields = bulkModulus
 *
 * 0 = \int_V \phi_p \cdot
 *  \left( \vec {\nabla} \cdot \vec{u} + \frac{p}{\kappa} \right) \, dV.
 *
 * ====================================================================== 
 */

/* ---------------------------------------------------------------------- */
/* g0 function for pressure equation.
 */
void
pylith_fekernels_Pressure_g0p(const PylithInt dim,
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
{ /* g0p */
  const PylithInt _numS = 2;
  const PylithInt i_disp = 0;
  const PylithInt i_pres = 2;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];
  const PylithScalar pres = s[sOff[i_pres]];

  const PylithInt i_bulkModulus = 2;

  const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

  PylithInt i;

  assert(_numS <= numS);
  assert(sOff);
  assert(s);
  assert(g0);

  PylithScalar strainTrace = 0;
  
  for (i=0; i < dim; ++i) {
    strainTrace += disp_x[i];
  } /* for */
  g0[0] += strainTrace + pres/bulkModulus;
} /* g0p */
					      

/* ---------------------------------------------------------------------- */
/* Jg0 function for pressure equation.
 */
void
pylith_fekernels_Pressure_Jg0pp(const PylithInt dim,
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
				PylithScalar Jg0[])
{ /* Jg0pp */
  const PylithInt i_bulkModulus = 1;
  const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
  
  Jg0[0] += 1.0/bulkModulus;
} /* Jg0pp_implicit */
					      
/* End of file */

