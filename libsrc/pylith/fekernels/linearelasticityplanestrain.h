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

/** @file libsrc/fekernels/linearelasticityplanestrain.h
 *
 * Kernels for linear elasticity plane strain.
 *
 * Solution fields: [disp(dim), vel(dim)]
 *
 * Isotropic, linear elasticity plane strain without initial stress/strain.
 * Auxiliary fields:
 * - 0: density(1)
 * - 1: mu(1)
 * - 2: lambda(1)
 * - 3: bodyforce(2,optional)
 * - 4: initialstress(4,optional)
 * - 5: initialstrain(4,optional)
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV + 
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * ====================================================================== 
 */

#if !defined(pylith_fekernels_linearelasticityplanestrain_h)
#define pylith_fekernels_linearelasticityplanestrain_h

/* Include directives ---------------------------------------------------
 */
#include <portinfo>

#include "pylith/utils/types.hh" 
#include "pylith/utils/error.h" 

/** f0 function for isotropic linear elasticity plane strain.
 *
 * Solution fields: [disp(dim), vel(dim)]
 * Auxiliary fields: [density(1), ...]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [2].
 * @param[in] numA Number of registered subfields in auxiliary field [2].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[out] f0 [dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_f0(const PylithInt dim,
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


/** g0 function for isotropic linear elasticity plane strain.
 *
 * Solution fields: [disp(dim), vel(dim)]
 * Auxiliary fields: [density(1), mu(1), lambda(1), bodyforce(dim), ...]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [2].
 * @param[in] numA Number of registered subfields in auxiliary field [2].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[out] g0 [dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g0(const PylithInt dim,
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


/** g1 function for isotropic linear elasticity plane strain WITHOUT initial stress and initial strain.
 *
 * Solution fields: [disp(dim), vel(dim)]
 * Auxiliary fields: [density(1), mu(1), lambda(1), ...]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [2].
 * @param[in] numA Number of registered subfields in auxiliary field [2].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[out] g1 [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1(const PylithInt dim,
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
							 PylithScalar g1[]);


/** g1 function for isotropic linear elasticity plane strain WITH initial stress and initial strain.
 *
 * Solution fields: [disp(dim), vel(dim)]
 * Auxiliary fields: [density(1), mu(1), lambda(1), ..., initialstress(dim*dim), initialstrain(dim*dim)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [2].
 * @param[in] numA Number of registered subfields in auxiliary field [2].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[out] g1 [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1_initstate(const PylithInt dim,
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
								   PylithScalar g1[]);


/** Jf0 function for isotropoc linear elasticity plane strain with implicit time stepping.
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [2].
 * @param[in] numA Number of registered subfields in auxiliary field [0].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] utshift Coefficient for dF/ds_t term in Jacobian.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] f0 Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0_uv_implicit(const PylithInt dim,
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


/** Jf0 function for isotropoc linear elasticity plane strain with explicit time stepping.
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [2].
 * @param[in] numA Number of registered subfields in auxiliary field [0].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] utshift Coefficient for dF/ds_t term in Jacobian.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] f0 Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0_uv_explicit(const PylithInt dim,
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

/** g3_uu entry function for 2-D plane strain isotropic linear elasticity.
 *
 * @param[in] dim Spatial dimension [3].
 * @param[in] numS Number of registered subfields in solution field [1].
 * @param[in] numA Number of registered subfields in auxiliary field [2].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] utshift Coefficient for dF/ds_t term in Jacobian.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] g3 Result [dim*dim*dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jg3_uu(const PylithInt dim,
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
							     PylithScalar Jg3[]);
					      

/** Calculate volumetric stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT initial stress and strain.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [lambda(1)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [1].
 * @param[in] numA Number of registered subfields in auxiliary field [1].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] f1 Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_volumetricStress(const PylithInt dim,
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


/** Calculate volumetric stress for 2-D plane strain isotropic linear
 * elasticity WITH initial stress and initial strain.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [lambda(1), initialstress(dim*dim), initialstrain(dim*dim)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [1].
 * @param[in] numA Number of registered subfields in auxiliary field [1].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] f1 Result [dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_volumetricStress_initstate(const PylithInt dim,
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
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [mu(1)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [1].
 * @param[in] numA Number of registered subfields in auxiliary field [1].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] f1 Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress(const PylithInt dim,
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
 * Solution fields: disp(dim)
 * Auxiliary fields: mu(1)
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [1].
 * @param[in] numA Number of registered subfields in auxiliary field [1].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] f1 Result [dim*dim].
 *
 * Solution fields: [disp(dim)]
 *
 * Auxiliary fields: [mu(1), initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress_initstate(const PylithInt dim,
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


/** f0 function for isotropic linear incompressible elasticity plane strain.
 *
 * Solution fields: [disp(dim), vel(dim)]
 * Auxiliary fields: [density(1), ...]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [3].
 * @param[in] numA Number of registered subfields in auxiliary field [5].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[out] f0 [dim].
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_f0(
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
						PylithScalar f0[]);


/** g0 function for isotropic linear incompressible elasticity plane strain.
 *
 * Solution fields: [disp(dim), vel(dim)]
 * Auxiliary fields: [bodyforce(dim), ...]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [3].
 * @param[in] numA Number of registered subfields in auxiliary field [5].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[out] f0 [dim].
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g0(
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
						PylithScalar g0[]);


/** g1 function for isotropic linear incompressible elasticity plane strain.
 *
 * Solution fields: [disp(dim), vel(dim), pres]
 * Auxiliary fields: [mu(1), ...]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [3].
 * @param[in] numA Number of registered subfields in auxiliary field [5].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[out] g1 [dim].
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g1(
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
						PylithScalar g1[]);


/** g1 function for isotropic linear incompressible elasticity plane strain.
/** with initial stress and strain.
 *
 * Solution fields: [disp(dim), vel(dim)]
 * Auxiliary fields: [mu(1), initstress(dim*dim), initstrain(dim*dim)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [3].
 * @param[in] numA Number of registered subfields in auxiliary field [5].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[out] f0 [dim].
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g1_initstate(
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
						PylithScalar g1[]);


/** g3_uu entry function for isotropic linear incompressible elasticity plane
/** strain.
 *
 * Solution fields: [disp(dim), vel(dim), pres]
 * Auxiliary fields: [mu(1)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [3].
 * @param[in] numA Number of registered subfields in auxiliary field [5].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[out] f0 [dim].
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg3_uu(
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
						PylithScalar Jg3[]);


/** g2_up entry function for isotropic linear incompressible elasticity plane
/** strain.
 *
 * Solution fields: [disp(dim), vel(dim), pres]
 * Auxiliary fields: None
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [3].
 * @param[in] numA Number of registered subfields in auxiliary field [5].
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[out] Jg2 [dim].
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_Jg2_up(
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
						PylithScalar Jg2[]);


#endif /* pylith_fekernels_linearelasticityplanestrain_h */


/* End of file */
