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

/** @file libsrc/fekernels/elasticity.h
 *
 * We cast the governing equation in the form:
 *
 * F(t,s,\dot{s}) = G(t,s).
 *
 * Using the finite-element method, we manipulate the weak form into
 * integrals over the domain with the form,
 *
 * \int_\Omega \vec{\phi} \cdot \vec{f}_0(t,s,\dot{s}) + 
 *   \nabla \vec{\phi} : \tensor{f}_1(t,s,\dot{s}) \, d\Omega, =
 *   \int_\Omega \vec{\phi} \cdot \vec{g}_0(t,s) + 
 *   \nabla \vec{\phi} : \tensor{g}_1(t,s) \, d\Omega,
 */

#if !defined(pylith_fekernels_elasticity_h)
#define pylith_fekernels_elasticity_h

/* Include directives ---------------------------------------------------
 */
#include <portinfo>

#include "pylith/utils/types.hh" 
#include "pylith/utils/error.h" 

/* ====================================================================== 
 * Kernels for displacement/velocity equation.
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
pylith_fekernels_f0_DispVel(const PylithInt dim,
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
pylith_fekernels_g0_DispVel(const PylithInt dim,
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
pylith_fekernels_Jf0_veldisp_DispVelImplicit(const PylithInt dim,
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
pylith_fekernels_Jf0_veldisp_DispVelExplicit(const PylithInt dim,
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
pylith_fekernels_Jg0_velvel_DispVel(const PylithInt dim,
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


/* ====================================================================== 
 * Generic elasticity kernels for inertia and body forces.
 *
 * Solution fields: [disp(dim), vel(dim)]
 *
 * Auxiliary fields: [density(1), ..., body_force(dim,optional)]
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV + 
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * Note: Ignore stress term because it depends on the constitutive model.
 *
 * LHS Residual
 * 
 * f0_Elasticity: \vec{f0} = \rho \frac{\partial \vec{v}(t)}{\partial t}
 *
 * RHS Residual
 *
 * g0_Elasticity: \vec{g0} = \vec{f}(t)
 *

Jf0_dispdisp_ElasticityExplicit
Jf0_dispdisp_ElasticityImplicit
 * ====================================================================== 
 */

/** f0 function for generic elasticity terms (inertia and body force).
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [2].
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
 * @param f0 Result [dim].
 */
void
pylith_fekernels_f0_Elasticity(const PylithInt dim,
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


/** g0 function for generic elasticity terms (inertia and body force).
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [2].
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
 * @param f0 Result [dim].
 *
 * Solution fields: [disp(dim), vel(dim)]
 *
 * Auxiliary fields: [density(1), mu(1), lambda(1), body force(dim)]
 */
void
pylith_fekernels_g0_Elasticity(const PylithInt dim,
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


/* ====================================================================== 
 * Kernels for stress.
 *
 * \int_V \nabla \vec{\phi}_u : \tensor{\sigma} \, dV

g1_IsotropicLinearElasticityPlaneStrain

Jg3_dispdisp_IsotropicLinearElasticityPlaneStrain

 * ====================================================================== 
 */

/** f1 entry function for isotropic linear elasticity in 3-D.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
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
 * @param f1 Result [dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), mu(1)]
 */
void
pylith_fekernels_f1_IsotropicLinearElasticity3D(const PylithInt dim,
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
						PylithScalar f1[]);
					      

/** g3_uu entry function for isotropic linear elasticity in 3-D.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
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
 * @param g3 Result [dim*dim*dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_g3_uu_IsotropicLinearElasticity3D(const PylithInt dim,
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
						   PylithScalar g3[]);
					      

/** Calculate volumetric stress for 3-D isotropic linear elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [1].
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
 * @param f1 Result [dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_volumetricStress_IsotropicLinearElasticity3D(const PylithInt dim,
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
							      PylithScalar stress[]);


/** Calculate deviatoric stress for 3-D isotropic linear elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [1].
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
 * @param f1 Result [dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_deviatoricStress_IsotropicLinearElasticity3D(const PylithInt dim,
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
							      PylithScalar stress[]);

/** f1 entry function for 2-D plane strain isotropic linear elasticity.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
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
 * @param f1 Result [dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_f1_IsotropicLinearElasticityPlaneStrain(const PylithInt dim,
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
							 PylithScalar f1[]);
					      

/** g3_uu entry function for 2-D plane strain isotropic linear elasticity.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
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
 * @param g3 Result [dim*dim*dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), mu(1)]
 */
void
pylith_fekernels_g3_uu_IsotropicLinearElasticityPlaneStrain(const PylithInt dim,
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
							    PylithScalar g3[]);
					      

/** Calculate volumetric stress for 2-D plane strain isotropic linear elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [1].
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
 * @param f1 Result [dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_volumetricStress_IsotropicLinearElasticityPlaneStrain(const PylithInt dim,
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
								       PylithScalar stress[]);


/** Calculate deviatoric stress for 2-D plane strain isotropic linear elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [1].
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
 * @param f1 Result [dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_deviatoricStress_IsotropicLinearElasticityPlaneStrain(const PylithInt dim,
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
								       PylithScalar stress[]);

/* ====================================================================== 
 * Kernels for incompressibility volume integral.
 *
 * \int_V \phi_p \left( \vec{\nabla} \cdot \vec{u} + \frac{p}{K} \right) \, dV
 * ====================================================================== 
 */

/** f0 entry function for incompressibility volume integral.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [2].
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
 * @param f0 Result [dim].
 *
 * Solution fields: [disp(dim), pres]
 *
 * Auxiliary fields: [lambda, mu]
 */
void
pylith_fekernels_f0_IncomprPIntegral(const PylithInt dim,
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

/** g0_vv entry function for incompressibility volume integral.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [2].
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
 *
 * Solution fields: [disp(dim), pres]
 *
 * Auxiliary fields: [lambda, mu]
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
					PylithScalar g0[]);


/** g2_vu entry function for incompressibility volume integral.
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
 *
 * Solution fields: [disp(dim), pres]
 *
 * Auxiliary fields: None
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
					PylithScalar g0[]);


/* ====================================================================== 
 * Kernels for incompressible elasticity volume integral.
 *
 * \int_V \tensor{S}:\nabla \vec{\phi}_u -
 * p \tensor{I}:\vec{\nabla} \vec{\phi}_u - \vec{\phi}_u \cdot \vec{f} \, dV
 * ====================================================================== 
 */

/** f0 entry function for body forces.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [2].
 * @param numA Number of registered subfields in auxiliary field [1].
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
 *
 * Solution fields: [disp(dim), pres]
 *
 * Auxiliary fields: [bodyforce(dim)]
 */
void
pylith_fekernels_f0_IncomprUIntegral(const PylithInt dim,
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

/** f1 entry function for 2-D plane strain incompressible isotropic linear
 *  elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [2].
 * @param numA Number of registered subfields in auxiliary field [4].
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
 *
 * Solution fields: [disp(dim), pres]
 *
 * Auxiliary fields: [lambda(1), mu(1),
 *                    initialstress(dim*dim), initialstrain(dim*dim)]
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
						PylithScalar f1[]);

/** f1 entry function for 3-D incompressible isotropic linear elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [2].
 * @param numA Number of registered subfields in auxiliary field [4].
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
 *
 * Solution fields: [disp(dim), pres]
 *
 * Auxiliary fields: [lambda(1), mu(1),
 *                    initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_f1_IncomprUIntegral3D(const PylithInt dim,
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
				       PylithScalar f1[]);

/** Calculate deviatoric stress for 2-D plane strain incompressible isotropic
 * linear elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [3].
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
 * @param stress Result [dim*dim].
 *
 * Solution fields: [disp(dim), pres]
 *
 * Auxiliary fields: [mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
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
						PylithScalar stress[]);

/** Calculate deviatoric stress for 3-D incompressible isotropic linear
 *  elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [1].
 * @param numA Number of registered subfields in auxiliary field [3].
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
 * @param stress Result [dim*dim].
 *
 * Solution fields: [disp(dim), pres]
 *
 * Auxiliary fields: [mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_deviatoricStress_IsotropicLinearElasticityIncompr3D(
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
						PylithScalar stress[]);


/** g3_uu entry function for 2-D plane strain incompressible isotropic linear
 * elasticity.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
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
 * @param g3 Result [dim*dim*dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), mu(1),
 *                    initialstress(dim*dim), initialstrain(dim*dim)]
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
						PylithScalar g3[]);

/** g2_uv entry function for 2-D plane strain incompressible isotropic linear
 * elasticity.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
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
 * @param g3 Result [dim*dim*dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), mu(1),
 *                    initialstress(dim*dim), initialstrain(dim*dim)]
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
						PylithScalar g2[]);

/** g3_uu entry function for 3-D incompressible isotropic linear elasticity.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
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
 * @param g3 Result [dim*dim*dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), mu(1),
 *                    initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_g3_uu_IncomprIsotropicLinearElasticity3D(
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
						PylithScalar g3[]);

/** g2_uv entry function for 3-D incompressible isotropic linear elasticity.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
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
 * @param g3 Result [dim*dim*dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [lambda(1), mu(1),
 *                    initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_g2_uv_IncomprIsotropicLinearElasticity3D(
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
						PylithScalar g2[]);

#endif /* pylith_fekernels_elasticity_h */


/* End of file */
