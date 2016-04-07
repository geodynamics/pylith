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

/** @file libsrc/fekernels/pressure.h
 *
 * Kernels for pressure volume integral for incompressible elasticity.
 *
 * Solution fields: [disp(dim), pres]
 *
 * Auxiliary fields: [bulkModulus]
 *
 * 0 = \int_V \phi_p \cdot
 *  \left( \vec {\nabla} \cdot \vec{u} + \frac{p}{\kappa} \right) \, dV.
 *
 * RHS Residual
 *
 * g0_Pressure: g0 = \phi_p \cdot
 *        \left( \vec {\nabla} \cdot \vec{u} + \frac{p}{\kappa} \right)
 *
 * RHS Jacobian
 *
 * Jg0_pp_Pressure: +1.0/bulkModulus
 *
 * ====================================================================== 
 */

#if !defined(pylith_fekernels_pressure_h)
#define pylith_fekernels_pressure_h

/* Include directives ---------------------------------------------------
 */
#include <portinfo>

#include "pylith/utils/types.hh" 
#include "pylith/utils/error.h" 

/** g0 function for pressure equation.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [3].
 * @param numA Number of registered subfields in auxiliary field [2].
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
 * @param x Coordinates of point evaluation.
 * @param g0 Result [dim].
 */
void
pylith_fekernels_Pressure_g0(const PylithInt dim,
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
			     PylithScalar g0[]);


/** Jg0 function for pressure equation.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [3].
 * @param numA Number of registered subfields in auxiliary field [2].
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
 * @param g0 Result [dim*dim].
 */
void
pylith_fekernels_Pressure_Jg0_pp(const PylithInt dim,
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
				 PylithScalar Jg0[]);


#endif /* pylith_fekernels_pressure_h */


/* End of file */
