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
#include "pylith/feassemble/Integrator.hh" // ISA Integrator

#include "pylith/topology/topologyfwd.hh" // USES Fields<Mesh>
#include "pylith/utils/array.hh" // HASA int_array

// PointForce ------------------------------------------------------
class pylith::bc::PointForce : public BoundaryCondition, 
			       public feassemble::Integrator
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
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
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

  /// Setup databases for querying for point forces.
  void _setupQueryDatabases(void);

  /** Query databases for point forces.
   *
   * @param mesh Finite-element mesh.
   */
  void _queryDatabases(const topology::Mesh& mesh);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  PointForce(const PointForce& m);

  /// Not implemented
  const PointForce& operator=(const PointForce& m);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Parameters for point forces.
  topology::Fields<topology::Field<topology::Mesh> >* _parameters;

  /// Start time for point forces.
  spatialdata::spatialdb::SpatialDB* _dbStartTime;

  /// Temporal evolution for amplitude of point forces.
  spatialdata::spatialdb::SpatialDB* _dbTimeAmp;
  

  int_array _points; ///< Points for forces.
  int_array _forceDOF; ///< Indices of degrees of freedom with forces.

}; // class PointForce

#include "PointForce.icc" // inline methods

#endif // pylith_bc_pointforce_hh


// End of file 
