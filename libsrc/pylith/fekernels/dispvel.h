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

/** @file libsrc/fekernels/dispvel.h
 *
 * Kernels for time evolution equation with displacement and velocity
 * solution fields.
 *
 * Solution fields: [disp(dim), vel(dim)]
 *
 * Auxiliary fields: None
 *
 * \int_V \vec{\phi}_v \cdot \left( \frac{\partial \vec{u}(t)}{\partial t} \right) \, dV = 
 *   \int_V \vec{\phi}_v \cdot \vec{v}(t) \, dV.
 *
 * LHS Residual
 * 
 * f0_DispVel: \vec{f0} = \frac{\partial \vec{u}(t)}{\partial t}
 *
 * RHS Residual
 *
 * g0_DispVel: \vec{g0} = \vec{v}(t)
 *
 * LHS Jacobian
 *
 * Jf0_veldisp_DispVelImplicit: tshift
 *
 * Jf0_veldisp_DispVelExplicit: 0
 *
 * RHS Jacobian
 *
 * Jg0_velvel_DispVel: +1.0
 *
 * ====================================================================== 
 */

#if !defined(pylith_fekernels_dispvel_h)
#define pylith_fekernels_dispvel_h

/* Include directives ---------------------------------------------------
 */
#include <portinfo>

#include "pylith/utils/types.hh" 
#include "pylith/utils/error.h" 

/** f0 function for disp/vel equation.
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
 * @param x Coordinates of point evaluation.
 * @param f0 Result [dim].
 */
void
pylith_fekernels_DispVel_f0(const PylithInt dim,
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
			    PylithScalar f0[]);


/** g0 function for disp/vel equation.
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
 * @param x Coordinates of point evaluation.
 * @param f0 Result [dim].
 */
void
pylith_fekernels_DispVel_g0(const PylithInt dim,
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


/** Jf0 function for disp/velocity equation with implicit time-stepping.
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
 * @param g0 Result [dim*dim].
 */
void
pylith_fekernels_DispVel_Jf0_vu_implicit(const PylithInt dim,
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
					 PylithScalar Jf0[]);

/** Jf0 function for disp/velocity equation with explicit time-stepping.
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
 * @param g0 Result [dim*dim].
 */
void
pylith_fekernels_DispVel_Jf0_vu_explicit(const PylithInt dim,
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
					 PylithScalar Jf0[]);


/** Jg0 function for disp/velocity equation.
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
 * @param g0 Result [dim*dim].
 */
void
pylith_fekernels_DispVel_Jg0_vv(const PylithInt dim,
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


#endif /* pylith_fekernels_dispvel_h */


/* End of file */
