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

#include "BoundaryCondition.hh" // ISA BoundaryCondition
#include "pylith/feassemble/Integrator.hh" // ISA Integrator

#include "pylith/utils/array.hh" // USES std::vector, double_array, int_array
#include "pylith/utils/sievetypes.hh" // USES real_section_type

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class Neumann;
    class TestNeumann; // unit testing
  } // bc
} // pylith


/// C++ implementation of Neumann boundary conditions.
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
   * @param mesh PETSc mesh
   * @param cs Coordinate system for mesh
   * @param upDir Direction perpendicular to horizontal surface tangent 
   *   direction that is not collinear with surface normal.
   */
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const double_array& upDir);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param mesh Finite-element mesh
   */
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 const ALE::Obj<Mesh>& mesh);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param mat Sparse matrix
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void integrateJacobian(PetscMat* mat,
			 const double t,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const ALE::Obj<Mesh>& mesh);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  Neumann(const Neumann&);

  /// Not implemented
  const Neumann& operator=(const Neumann&);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Traction vector in global coordinates at integration points.
  ALE::Obj<real_section_type> _tractionGlobal;

}; // class Neumann

// #include "Neumann.icc" // inline methods

#endif // pylith_bc_neumann_hh


// End of file 
