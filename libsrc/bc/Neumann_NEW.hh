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

/** @file libsrc/bc/Neumann_NEW.hh
 *
 * @brief C++ implementation of time dependent Neumann_NEW (traction)
 * boundary conditions applied to a simply-connected boundary.
 */

#if !defined(pylith_bc_neumann_hh)
#define pylith_bc_neumann_hh

// Include directives ---------------------------------------------------
#include "BCIntegratorSubMesh.hh" // ISA BCIntegratorSubMesh
#include "TimeDependent.hh" // ISA TimeDependent

// Neumann_NEW ------------------------------------------------------
class pylith::bc::Neumann_NEW : public BCIntegratorSubMesh, 
				public TimeDependent
{ // class Neumann_NEW
  friend class TestNeumann_NEW; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  Neumann_NEW(void);

  /// Destructor.
  ~Neumann_NEW(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Initialize boundary condition.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Direction perpendicular to horizontal surface tangent 
   *   direction that is not collinear with surface normal.
   */
  void initialize(const topology::Mesh& mesh,
		  const double upDir[3]);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual.
   * @param t Current time.
   * @param fields Solution fields.
   */
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobian(topology::Jacobian* jacobian,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  /** Get cell field with BC information.
   *
   * @param fieldType Type of field.
   * @param name Name of field.
   * @param mesh Finite-element mesh.
   * @param fields Solution fields.
   *
   * @returns Traction vector at integration points.
   */
  const topology::Field<topology::SubMesh>&
  cellField(const char* name,
	    topology::SolutionFields* const fields);

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
  const spatialdata::units::Nondimensional& _getNormalizer(void) const;

  /// Query databases for time dependent parameters.
  void _queryDatabases(void);

  /** Query database for values.
   *
   * @param field Field in which to store values.
   * @param db Spatial database with values.
   * @param querySize Number of values at each location.
   * @param scale Dimension scale associated with values.
   */
  void _queryDB(topology::Field<topology::SubMesh>* field,
		spatialdata::spatialdb::SpatialDB* const db,
		const int querySize,
		const double scale);

  /** Convert parameters in local coordinates to global coordinates.
   *
   * @param upDir Direction perpendicular to horizontal surface tangent 
   *   direction that is not collinear with surface normal.
   */
  void _paramsLocalToGlobal(const double upDir[3]);

  /** Calculate spatial and temporal variation of value over the list
   *  of submesh.
   *
   * @param t Current time.
   */
  void _calculateValue(const double t);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  Neumann_NEW(const Neumann_NEW&); ///< Not implemented.
  const Neumann_NEW& operator=(const Neumann_NEW&); ///< Not implemented.

}; // class Neumann_NEW

#include "Neumann_NEW.icc"

#endif // pylith_bc_neumann_hh


// End of file 
