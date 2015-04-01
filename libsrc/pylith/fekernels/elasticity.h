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
 */

#if !defined(pylith_fekernels_elasticity_h)
#define pylith_fekernels_elasticity_h

// Include directives ---------------------------------------------------
#include <portinfo>

/** Main f0 function for elasticity with inertia and body forces in 3-D.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
 * @param indicesS Indices of subfields in solution field [numS].
 * @param numA Number of registered subfields in auxiliary field [2].
 * @param indicesA Indices of subfields in auxiliary field [numA].
 * @param uOff Offset of registered subfields in solution field [numS].
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_tt Second time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param a Auxiliary field with all subfields.
 * @param a_x Gradient of auxiliary field.
 * @param a_t Time derivative of auxiliary field. DO WE NEED THIS?
 * @param x Position for point evaluation.
 * @param f0 Result [dim].
 *
 * @returns 0 if no errors.
 * 
 * Solution fields = [disp(dim)]
 *
 * Auxiliary fields = [density(1), body force(dim)]
 */
PetscErrorCode
pylith_f0_ElasticityInertia3D(const PetscInt dim,
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
			      PetscScalar f0[]);


/** Main f1 function for isotropic linear elasticity in 3-D.
 *
 * @param dim Spatial dimension [3].
 * @param numS Number of registered subfields in solution field [1].
 * @param indicesS Indices of subfields in solution field [numS].
 * @param numA Number of registered subfields in auxiliary field [2].
 * @param indicesA Indices of subfields in auxiliary field [numA].
 * @param uOff Offset of registered subfields in solution field [numS].
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_tt Second time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param a Auxiliary field with all subfields.
 * @param a_x Gradient of auxiliary field.
 * @param a_t Time derivative of auxiliary field. DO WE NEED THIS?
 * @param x Position for point evaluation.
 * @param f1 Result [dim].
 *
 * @returns 0 if no errors.
 * 
 * Solution fields = [disp(dim)]
 *
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
				      PetscScalar f1[]);
					      

/** f0 function for inertia.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [1].
 * @param indicesS Indices of subfields in solution field [numS].
 * @param numA Number of registered subfields in auxiliary field [1].
 * @param indicesA Indices of subfields in auxiliary field [numA].
 * @param uOff Offset of registered subfields in solution field [numS].
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_tt Second time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param a Auxiliary field with all subfields.
 * @param a_x Gradient of auxiliary field.
 * @param a_t Time derivative of auxiliary field. DO WE NEED THIS?
 * @param x Position for point evaluation.
 * @param f0 Result [dim].
 *
 * @returns 0 if no errors.
 * 
 * Solution fields = [disp(dim)]
 *
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
		  PetscScalar f0[]);
					      

/** f0 function for body force.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [0].
 * @param indicesS Indices of subfields in solution field [numS].
 * @param numA Number of registered subfields in auxiliary field [1].
 * @param indicesA Indices of subfields in auxiliary field [numA].
 * @param uOff Offset of registered subfields in solution field [numS].
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_tt Second time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param a Auxiliary field with all subfields.
 * @param a_x Gradient of auxiliary field.
 * @param a_t Time derivative of auxiliary field. DO WE NEED THIS?
 * @param x Position for point evaluation.
 * @param f0 Result [dim].
 *
 * @returns 0 if no errors.
 * 
 * Solution fields = NONE
 *
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
		    PetscScalar f0[]);


/** f1 function for volumetic stress for isotropic linear elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [1].
 * @param indicesS Indices of subfields in solution field [numS].
 * @param numA Number of registered subfields in auxiliary field [1].
 * @param indicesA Indices of subfields in auxiliary field [numA].
 * @param uOff Offset of registered subfields in solution field [numS].
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_tt Second time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param a Auxiliary field with all subfields.
 * @param a_x Gradient of auxiliary field.
 * @param a_t Time derivative of auxiliary field. DO WE NEED THIS?
 * @param x Position for point evaluation.
 * @param f1 Result [dim*dim].
 *
 * @returns 0 if no errors.
 * 
 * Solution fields = [disp(dim)]
 *
 * Auxiliary fields = [lambda(1)]
 */
PetscErrorCode
pylith_f1_IsotropicLinearElasticityVolumetricStress(const PetscInt dim,
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
						    PetscScalar f0[]);


/** f1 function for deviatoric stress for isotropic linear elasticity.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [1].
 * @param indicesS Indices of subfields in solution field [numS].
 * @param numA Number of registered subfields in auxiliary field [1].
 * @param indicesA Indices of subfields in auxiliary field [numA].
 * @param uOff Offset of registered subfields in solution field [numS].
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_tt Second time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param a Auxiliary field with all subfields.
 * @param a_x Gradient of auxiliary field.
 * @param a_t Time derivative of auxiliary field. DO WE NEED THIS?
 * @param x Position for point evaluation.
 * @param f1 Result [dim*dim].
 *
 * @returns 0 if no errors.
 * 
 * Solution fields = [disp(dim)]
 *
 * Auxiliary fields = [mu(1)]
 */
PetscErrorCode
pylith_f1_IsotropicLinearElasticityDeviatoricStress(const PetscInt dim,
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
						    PetscScalar f0[]);

#endif // pylith_fekernels_elasticity_h


// End of file 
