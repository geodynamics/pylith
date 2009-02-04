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
#include "BoundaryCondition.hh" // ISA BoundaryCondition
#include "pylith/feassemble/Integrator.hh" // ISA Integrator

#include "pylith/utils/array.hh" // USES std::vector, double_array, int_array

// Neumann --------------------------------------------------------------
class pylith::bc::Neumann : public BoundaryCondition, 
			    public feassemble::Integrator
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
  void integrateResidual(const topology::Field& residual,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobian(PetscMat* jacobian,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  /** Get boundary mesh.
   *
   * @returns Boundary mesh.
   */
  const topology::SubMesh& boundaryMesh(void) const;

  /** Get cell field with BC information.
   *
   * @param fieldType Type of field.
   * @param name Name of field.
   * @param mesh Finite-element mesh.
   * @param fields Solution fields.
   *
   * @returns Traction vector at integration points.
   */
  const topology::Field&
  cellField(const char* name,
	    topology::SolutionFields* const fields);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  Neumann(const Neumann&);

  /// Not implemented
  const Neumann& operator=(const Neumann&);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Mesh over which tractions are applied
  topology::SubMesh* _boundaryMesh;

  /// Traction vector in global coordinates at integration points.
  topology::FieldSubMesh* _tractions;

}; // class Neumann

#endif // pylith_bc_neumann_hh


// End of file 
