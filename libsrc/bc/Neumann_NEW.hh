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

/** @file libsrc/bc/Neumann.hh
 *
 * @brief C++ implementation of time dependent Neumann (traction)
 * boundary conditions applied to a simply-connected boundary.
 */

#if !defined(pylith_bc_neumann_hh)
#define pylith_bc_neumann_hh

// Include directives ---------------------------------------------------
#include "BCIntegratorSubmesh.hh" // ISA BCIntegratorSubmesh
#include "TimeDependent.hh" // ISA TimeDependent

// Neumann ------------------------------------------------------
class pylith::bc::Neumann : public BCIntegratorSubmesh, 
			    public TimeDependent
{ // class Neumann
  friend class TestNeumann; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  Neumann(void);

  /// Destructor.
  ~Neumann(void);

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

  /** Query databases for time dependent parameters.
   *
   * @param valueScale Dimension scale for value.
   * @param fieldName Name of field associated with value.
   */
  void _queryDatabases(const double valueScale,
		       const char* fieldName);

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

  /// Convert parameters in local coordinates to global coordinates.
  void _paramsLocalToGlobal(void);

  /** Calculate spatial and temporal variation of value over the list
   *  of submesh.
   *
   * @param t Current time.
   */
  void _calculateValue(const double t);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  Neumann(const Neumann&); ///< Not implemented.
  const Neumann& operator=(const Neumann&); ///< Not implemented.

}; // class Neumann

#include "Neumann.icc"

#endif // pylith_bc_neumann_hh


// End of file 
