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
 * @brief C++ implementation of Neumann (prescribed tractions
 * on a surface) boundary conditions.
 */

#if !defined(pylith_bc_neumann_hh)
#define pylith_bc_neumann_hh

// Include directives ---------------------------------------------------
#include "BCIntegratorSubMesh.hh" // ISA BCIntegratorSubMesh

// Neumann --------------------------------------------------------------
class pylith::bc::Neumann : public BCIntegratorSubMesh
{ // class Neumann
  friend class TestNeumann; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  Neumann(void);

  /// Destructor.
  ~Neumann(void);

  /** Set database for boundary condition parameters.
   *
   * @param db Spatial database
   */
  void db(spatialdata::spatialdb::SpatialDB* const db);

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
	    topology::SolutionFields* const fields =0);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  spatialdata::spatialdb::SpatialDB* _db; ///< Spatial database w/parameters

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  Neumann(const Neumann&); ///< Not implemented
  const Neumann& operator=(const Neumann&); ///< Not implemented

}; // class Neumann

#include "Neumann.icc" // inline methods

#endif // pylith_bc_neumann_hh


// End of file 
