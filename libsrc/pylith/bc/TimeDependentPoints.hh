// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/TimeDependentPoints.hh
 *
 * @brief C++ implementation of time dependent boundary conditions
 * applied to a set of vertices.
 */

#if !defined(pylith_bc_timedependentpoints_hh)
#define pylith_bc_timedependentpoints_hh

// Include directives ---------------------------------------------------
#include "BoundaryConditionPoints.hh" // ISA BoundaryConditionPoints
#include "TimeDependent.hh" // ISA TimeDependent

#include "pylith/utils/array.hh" // HASA int_array

// TimeDependentPoints ------------------------------------------------------
/// Time dependent boundary conditions applied to a set of vertices.
class pylith::bc::TimeDependentPoints : public BoundaryConditionPoints, 
					public TimeDependent
{ // class TimeDependentPoints
  friend class TestTimeDependentPoints; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  TimeDependentPoints(void);

  /// Destructor.
  ~TimeDependentPoints(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set indices of degrees of freedom associated with BC.
   *
   * Note: Forces at all points are applied to the same degrees of freedom.
   *
   * Example: [0, 1] to apply forces to x and y degrees of freedom in
   * Cartesian system.
   *
   * @param flags Array of indices for degrees of freedom for forces.
   * @param size Size of array
   */
  void bcDOF(const int* flags,
	     const int size);  

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Get label of boundary condition surface.
   *
   * @returns Label of surface (from mesh generator).
   */
  const char* _getLabel(void) const;

  /** Get manager of scales used to nondimensionalize problem.
   *
   * @returns Nondimensionalizer.
   */
  virtual
  const spatialdata::units::Nondimensional& _getNormalizer(void) const = 0;

  /** Query databases for time dependent parameters.
   *
   * @param mesh Finite-element mesh.
   * @param valueScale Dimension scale for value.
   * @param fieldName Name of field associated with value.
   */
  void _queryDatabases(const topology::Mesh& mesh,
		       const PylithScalar valueScale,
		       const char* fieldName);

  /** Query database for values.
   *
   * @param name Name of field in which to store values.
   * @param db Spatial database with values.
   * @param querySize Number of values at each location.
   * @param scale Dimension scale associated with values.
   */
  void _queryDB(const char* name,
		spatialdata::spatialdb::SpatialDB* const db,
		const int querySize,
		const PylithScalar scale);

  /** Calculate spatial and temporal variation of value over the list
   *  of points.
   *
   * @param t Current time.
   */
  void _calculateValue(const PylithScalar t);

  /** Calculate increment in spatial and temporal variation of value
   *  over the list of points.
   *
   * @param t0 Time when increment begins.
   * @param t1 Time when increment ends.
   */
  void _calculateValueIncr(const PylithScalar t0,
			   const PylithScalar t1);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  int_array _bcDOF; ///< Degrees of freedom associated with BC.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  TimeDependentPoints(const TimeDependentPoints&); ///< Not implemented.
  const TimeDependentPoints& operator=(const TimeDependentPoints&); ///< Not implemented.

}; // class TimeDependentPoints

#include "TimeDependentPoints.icc"

#endif // pylith_bc_timedependentpoints_hh


// End of file 
