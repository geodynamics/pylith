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

/** @file libsrc/bc/PointForce.hh
 *
 * @brief C++ implementation of point force on vertices.
 */

#if !defined(pylith_bc_pointforce_hh)
#define pylith_bc_pointforce_hh

// Include directives ---------------------------------------------------
#include "BoundaryCondition.hh" // ISA BoundaryCondition
#include "pylith/topology/Mesh.hh" // ISA Integrator<Quadrature<Mesh> >
#include "pylith/feassemble/Quadrature.hh" // ISA Integrator<Quadrature<Mesh >
#include "pylith/feassemble/Integrator.hh" // ISA Integrator

#include "pylith/topology/topologyfwd.hh" // USES Fields<Mesh>
#include "pylith/utils/array.hh" // HASA int_array

// PointForce ------------------------------------------------------
class pylith::bc::PointForce : public BoundaryCondition, 
			       public feassemble::Integrator<feassemble::Quadrature<topology::Mesh> >
{ // class PointForce
  friend class TestPointForce; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  PointForce(void);

  /// Destructor.
  ~PointForce(void);

  /** Set indices of vertices with point forces.
   *
   * Note: forces at all points are applied to the same degrees of freedom.
   *
   * Example: [0, 1] to apply forces to x and y degrees of freedom in
   * Cartesian system.
   *
   * @param flags Array of indices for degrees of freedom for forces.
   * @param size Size of array
   */
  void forceDOF(const int* flags,
		const int size);

  /** Set database for initial values.
   *
   * @param db Spatial database
   */
  void dbInitial(spatialdata::spatialdb::SpatialDB* const db);

  /** Set database for rate of change of values.
   *
   * @param db Spatial database
   */
  void dbRate(spatialdata::spatialdb::SpatialDB* const db);

  /** Set database for start time of rate change.
   *
   * @param db Spatial database
   */
  void dbRateTime(spatialdata::spatialdb::SpatialDB* const db);

  /** Set database for change in values.
   *
   * @param db Spatial database
   */
  void dbChange(spatialdata::spatialdb::SpatialDB* const db);

  /** Set database for start time of change in values.
   *
   * @param db Spatial database
   */
  void dbChangeTime(spatialdata::spatialdb::SpatialDB* const db);

  /** Set database for temporal evolution of change in value.
   *
   * @param db Time history database.
   */
  void dbTimeHistory(spatialdata::spatialdb::TimeHistory* const db);


  /** Initialize boundary condition.
   *
   * @param mesh PETSc mesh
   * @param upDir Vertical direction (somtimes used in 3-D problems).
   */
  void initialize(const topology::Mesh& mesh,
		  const double upDir[3]);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidualAssembled(topology::Field<topology::Mesh>* residual,
				  const double t,
				  topology::SolutionFields* const fields);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Get mesh labels for points associated with applied forces.
   *
   * @param mesh Finite-element mesh.
   */
  void _getPoints(const topology::Mesh& mesh);

  /** Query databases for time dependent parameters.
   *
   * @param mesh Finite-element mesh.
   * @param valueScale Dimension scale for value.
   * @param fieldName Name of field associated with value.
   */
  void _queryDatabases(const topology::Mesh& mesh,
		       const double valueScale,
		       const char* fieldName);

  /** Wuery database for values.
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

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Parameters for point forces.
  topology::Fields<topology::Field<topology::Mesh> >* _parameters;

  /// Spatial database for initial values.
  spatialdata::spatialdb::SpatialDB* _dbInitial;

  /// Spatial database for rate of change of values.
  spatialdata::spatialdb::SpatialDB* _dbRate;

  /// Spatial database for start time of rate change.
  spatialdata::spatialdb::SpatialDB* _dbRateTime;

  /// Spatial database for change in value.
  spatialdata::spatialdb::SpatialDB* _dbChange;

  /// Spatial database for start time of change in value.
  spatialdata::spatialdb::SpatialDB* _dbChangeTime;

  /// Temporal evolution of amplitude for change in value;
  spatialdata::spatialdb::TimeHistory* _dbTimeHistory;
  
  int_array _points; ///< Points for forces.
  int_array _bcDOF; ///< Indices of degrees of freedom with boundary condition.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  PointForce(const PointForce&); ///< Not implemented.
  const PointForce& operator=(const PointForce&); ///< Not implemented.

}; // class PointForce

#include "PointForce.icc" // inline methods

#endif // pylith_bc_pointforce_hh


// End of file 
