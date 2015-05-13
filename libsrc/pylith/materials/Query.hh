// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/** @file libsrc/materials/Quuery.hh
 *
 */

#if !defined(pylith_materials_query_hh)
#define pylith_materials_query_hh

// Include directives ---------------------------------------------------
#include "pylith/materials/materialsfwd.hh" // forward declarations

#include "pylith/utils/types.hh" 

// Query ----------------------------------------------------------------
/** @brief Query functions associated with material properties.
 */
class pylith::materials::Query
{ // Query
  friend class TestQuery; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Query for density in 2-D.
   *
   * @param x Coordinates (nondimensioned) of point location for query.
   * @param values Array of values to be returned.
   * @param context Query context.
   */
  static
  void dbQueryDensity2D(const PetscReal x[],
			PetscScalar* values,
			void* context);

  /** Query for shear modulus, $\mu$, in 2-D.
   *
   * @param x Coordinates (nondimensioned) of point location for query.
   * @param values Array of values to be returned.
   * @param context Query context.
   */
  static
  void dbQueryMu2D(const PetscReal x[],
		   PetscScalar* values,
		   void* context);

  /** Query for Lame's constant, $\lambda, in 2-D.
   *
   * @param x Coordinates (nondimensioned) of point location for query.
   * @param values Array of values to be returned.
   * @param context Query context.
   */
  static
  void dbQueryLambda2D(const PetscReal x[],
		       PetscScalar* values,
		       void* context);

  /** Query for body force vector in 2-D.
   *
   * @param x Coordinates (nondimensioned) of point location for query.
   * @param values Array of values to be returned.
   * @param context Query context.
   */
  static
  void dbQueryBodyForce2D(const PetscReal x[],
			  PetscScalar* values,
			  void* context);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Query(void); ///< Not implemented
  ~Query(void); ///< Not implemented
  Query(const Query&); ///< Not implemented
  const Query& operator=(const Query&); ///< Not implemented

}; // Query

#endif // pylith_materials_query_hh */


// End of file
