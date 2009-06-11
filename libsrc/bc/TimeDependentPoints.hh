// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
  
  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Get label of boundary condition surface.
   *
   * @returns Label of surface (from mesh generator).
   */
  const char* _getLabel(void) const;

  /** Query databases for time dependent parameters.
   *
   * @param mesh Finite-element mesh.
   * @param valueScale Dimension scale for value.
   * @param fieldName Name of field associated with value.
   */
  void _queryDatabases(const topology::Mesh& mesh,
		       const double valueScale,
		       const char* fieldName);

  /** Query database for values.
   *
   * @param field Field in which to store values.
   * @param db Spatial database with values.
   * @param querySize Number of values at each location.
   * @param scale Dimension scale associated with values.
   */
  void _queryDB(topology::Field<topology::Mesh>* field,
		spatialdata::spatialdb::SpatialDB* const db,
		const int querySize,
		const double scale);

  /** Calculate spatial and temporal variation of value over the list
   *  of points.
   *
   * @param t Current time.
   */
  void _calculateValue(const double t);

  /** Calculate increment in spatial and temporal variation of value
   *  over the list of points.
   *
   * @param t0 Time when increment begins.
   * @param t1 Time when increment ends.
   */
  void _calculateValueIncr(const double t0,
			   const double t1);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  TimeDependentPoints(const TimeDependentPoints&); ///< Not implemented.
  const TimeDependentPoints& operator=(const TimeDependentPoints&); ///< Not implemented.

}; // class TimeDependentPoints

#include "TimeDependentPoints.icc"

#endif // pylith_bc_timedependentpoints_hh


// End of file 
