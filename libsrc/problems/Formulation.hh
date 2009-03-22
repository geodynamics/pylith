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
 * @file pylith/problems/Formulation.hh
 *
 * @brief Object for managing forming Jacobian and residual for the problem.
 */

#if !defined(pylith_problems_formulation_hh)
#define pylith_problems_formulation_hh

// Include directives ---------------------------------------------------
#include "problemsfwd.hh" // forward declarations

#include "pylith/feassemble/feassemblefwd.hh" // USES Integrator
#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field, SolutionFields

#include "pylith/utils/petscfwd.h" // USES PetscVec, PetscMat, PetscSNES

#include "pylith/utils/array.hh" // HASA std::vector


// Formulation ----------------------------------------------------------
class pylith::problems::Formulation
{ // Integrator
  friend class TestFormulation; // unit testing

// PUBLIC STRUCTS ///////////////////////////////////////////////////////
public :

  struct ArgsResidual {
    Formulation* object;
    topology::Field<topology::Mesh>* const residual;
    topology::SolutionFields* const fields;
    double t;
    double dt;
  }; // ArgsResidual
  
  struct ArgsJacobian {
    Formulation* object;
    topology::Jacobian* jacobian;
    topology::SolutionFields* const fields;
    double t;
    double dt;
  }; // ArgsJacobian

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

  /** Generic C interface for reformResidual for integration with
   * PETSc SNES solvers.
   *
   * @param snes PETSc scalable nonlinear equation solver.
   * @param solutionVec PETSc vector for solution.
   * @param residualVec PETSc vector for residual.
   * @param context ArgsResidual structure with arguments.
   */
  static
  void reformResidual(PetscSNES snes,
		      PetscVec solutionVec,
		      PetscVec residualVec,
		      void* context);

  /** Generic C interface for reformJacobian for integration with
   * PETSc SNES solvers.
   *
   * @param snes PETSc scalable nonlinear equation solver.
   * @param solutionVec PETSc vector for solution.
   * @param jacobianMat PETSc sparse matrix for system Jacobian.
   * @param preconditionerMat PETSc sparse matrix for preconditioner.
   * @param Flag indicating layout of preconditioner matrix.
   * @param context ArgsJacobian structure with arguments.
   */
  static
  void reformJacobian(PetscSNES snes,
		      PetscVec solutionVec,
		      PetscMat jacobianMat,
		      PetscMat preconditionerMat,
		      int* preconditionerLayout,
		      void* context);

  /** Set integrators over the mesh.
   *
   * @param integrators Integrators over the mesh.
   * @param numIntegrators Number of integrators.
   */
  void meshIntegrators(IntegratorMesh** integrators,
		       const int numIntegrators);
  
  /** Set integrators over lower-dimension meshes.
   *
   * @param integrators Integrators over lower-dimension meshes.
   * @param numIntegrators Number of integrators.
   */
  void submeshIntegrators(IntegratorSubMesh** integrators,
			  const int numIntegrators);

  /** Initialize solver for formulation.
   *
   * @param 

  /** Reform system residual.
   *
   * @param residual Field containing values for residual
   * @param fields Solution fields
   * @param t Current time.
   * @param dt Current time step (t -> t+dt).
   */
  void reformResidual(topology::Field<topology::Mesh>* const residual,
		      topology::SolutionFields* const fields,
		      const double t,
		      const double dt);
  
  /** Reform system Jacobian.
   *
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param fields Solution fields
   * @param t Current time.
   * @param dt Current time step (t -> t+dt).
   */
  void reformJacobian(topology::Jacobian* jacobian,
		      topology::SolutionFields* const fields,
		      const double t,
		      const double dt);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::vector<IntegratorMesh*> _meshIntegrators;
  std::vector<IntegratorSubMesh*> _submeshIntegrators;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Formulation(const Formulation&); ///< Not implemented
  const Formulation& operator=(const Formulation&); ///< Not implemented

}; // Formulation

#endif // pylith_problems_formulation_hh


// End of file 
