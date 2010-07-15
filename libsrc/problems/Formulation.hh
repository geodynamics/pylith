// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/**
 * @file libsrc/problems/Formulation.hh
 *
 * @brief C++ Object that manages reforming the Jacobian and residual for
 * the problem.
 */

#if !defined(pylith_problems_formulation_hh)
#define pylith_problems_formulation_hh

// Include directives ---------------------------------------------------
#include "problemsfwd.hh" // forward declarations

#include "pylith/feassemble/feassemblefwd.hh" // USES Integrator
#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field, SolutionFields

#include "pylith/utils/petscfwd.h" // USES PetscVec, PetscMat

#include "pylith/utils/array.hh" // HASA std::vector


// Formulation ----------------------------------------------------------
/// Reform the Jacobian and residual for the problem.
class pylith::problems::Formulation
{ // Formulation
  friend class TestFormulation; // unit testing

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private :

  typedef feassemble::Integrator<feassemble::Quadrature<topology::Mesh> > IntegratorMesh;
  typedef feassemble::Integrator<feassemble::Quadrature<topology::SubMesh> > IntegratorSubMesh;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Formulation(void);

  /// Destructor
  ~Formulation(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Set flag for splitting fields.
   *
   * @param flag True if splitting fields, false otherwise.
   */
  void splitFields(const bool flag);

  /** Get flag for splitting fields.
   *
   * @returns flag True if splitting fields, false otherwise.
   */
  bool splitFields(void) const;

  /** Set flag for using custom preconditioner for Lagrange constraints.
   *
   * @param flag True if using custom fault preconditioner, false otherwise.
   */
  void useCustomConstraintPC(const bool flag);

  /** Get flag indicating use of custom conditioner for Lagrange constraints.
   *
   * @returns True if using custom fault preconditioner, false otherwise.
   */
  bool useCustomConstraintPC(void) const;

  /** Get solution fields.
   *
   * @returns solution fields.
   */
  const topology::SolutionFields& fields(void) const;

  /** Get flag indicating whether Jacobian is symmetric.
   *
   * @returns True if Jacobian is symmetric, otherwise false.
   */
  bool isJacobianSymmetric(void) const;
  
  /** Set handles to integrators over the mesh.
   *
   * @param integrators Integrators over the mesh.
   * @param numIntegrators Number of integrators.
   */
  void meshIntegrators(IntegratorMesh** integrators,
		       const int numIntegrators);
  
  /** Set handles to integrators over lower-dimension meshes.
   *
   * @param integrators Integrators over lower-dimension meshes.
   * @param numIntegrators Number of integrators.
   */
  void submeshIntegrators(IntegratorSubMesh** integrators,
			  const int numIntegrators);

  /** Set handle to preconditioner.
   *
   * @param pc PETSc preconditioner.
   */
  void customPCMatrix(PetscMat& mat);

  /** Update handles and parameters for reforming the Jacobian and
   *  residual.
   *
   * @param jacobian Handle to sparse matrix for Jacobian of system.
   * @param fields Handle to solution fields.
   * @param t Current time (nondimensional).
   * @param dt Time step (nondimension).
   */
  void updateSettings(topology::Jacobian* jacobian,
		      topology::SolutionFields* fields,
		      const double t,
		      const double dt);

  /** Update handles and parameters for reforming the Jacobian and
   *  residual.
   *
   * @param jacobian Handle to diagonal matrix (as Field) for system Jacobian.
   * @param fields Handle to solution fields.
   * @param t Current time (nondimensional).
   * @param dt Time step (nondimension).
   */
  void updateSettings(topology::Field<topology::Mesh>* jacobian,
		      topology::SolutionFields* fields,
		      const double t,
		      const double dt);

  /** Reform system residual.
   *
   * @param tmpResidualVec Temporary PETSc vector for residual.
   * @param tmpSolutionVec Temporary PETSc vector for solution.
   */
  void reformResidual(const PetscVec* tmpResidualVec =0,
		      const PetscVec* tmpSolutionVec =0);
  
  /** Reform system residual.
   *
   * @param tmpResidualVec Temporary PETSc vector for residual.
   * @param tmpSolutionVec Temporary PETSc vector for solution.
   */
  void reformResidualLumped(const PetscVec* tmpResidualVec =0,
          const PetscVec* tmpSolutionVec =0);

  /* Reform system Jacobian.
   *
   * @param tmpSolveSolnVec Temporary PETSc vector for solution.
   */
  void reformJacobian(const PetscVec* tmpSolveSolnVec =0);

  /* Reform system Jacobian.
   */
  void reformJacobianLumped(void);

  /** Constrain solution space.
   *
   * @param tmpSolutionVec Temporary PETSc vector for solution.
   */
  void constrainSolnSpace(const PetscVec* tmpSolutionVec);

  /** Adjust solution from solver with lumped Jacobian to match Lagrange
   *  multiplier constraints.
   */
  void adjustSolnLumped(void);

  /// Compute rate fields (velocity and/or acceleration) at time t.
  virtual
  void calcRateFields(void) = 0;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  double _t; ///< Current time (nondimensional).
  double _dt; ///< Current time step (nondimensional).
  topology::Jacobian* _jacobian; ///< Handle to Jacobian of system.
  PetscMat _precondMatrix; ///< Custom PETSc preconditioning matrix.
  topology::Field<topology::Mesh>* _jacobianLumped; ///< Handle to lumped Jacobian of system.
  topology::SolutionFields* _fields; ///< Handle to solution fields for system.

  /// Integrators over subdomains of the mesh.
  std::vector<IntegratorMesh*> _meshIntegrators;

  ///< Integrators over lower-dimensional subdomains of the mesh.
  std::vector<IntegratorSubMesh*> _submeshIntegrators;

  bool _isJacobianSymmetric; ///< Is system Jacobian symmetric?
  bool _splitFields; ///< True if splitting fields.

  /// True if using custom preconditioner for Lagrange constraints.
  bool _useCustomConstraintPC;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Formulation(const Formulation&); ///< Not implemented
  const Formulation& operator=(const Formulation&); ///< Not implemented

}; // Formulation

#endif // pylith_problems_formulation_hh


// End of file 
