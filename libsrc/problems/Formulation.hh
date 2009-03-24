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

  /** Set handle to solution fields.
   *
   * @param fields Solution fields.
   */
  void solutionFields(SolutionFields* fields);

  /** Update current time and time step for advancing from t to t+dt.
   *
   * @param t Current time (nondimensional).
   * @param dt Time step (nondimension).
   */
  void updateTime(const double t,
		  const double dt);

  /** Initialize solver for formulation.
   *
   * @param solver Solver for system.
   */
  void initializeSolver(Solver* solver);

  /** Reform system residual.
   *
   * @param residual Field containing values for residual
   * @param fields Solution fields
   * @param t Current time.
   * @param dt Current time step (t -> t+dt).
   */
  void reformResidual(void);
  
  /** Reform system Jacobian.
   *
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param fields Solution fields
   * @param t Current time.
   * @param dt Current time step (t -> t+dt).
   */
  void reformJacobian(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  double t; ///< Current time (nondimensional).
  double dt; ///< Current time step (nondimensional).

  topology::Jacobian* _jacobian; ///< Jacobian of system.
  topology::SolutionFields* _fields; ///< Solution fields for system.

  /// Integrators over subdomains of the mesh.
  std::vector<IntegratorMesh*> _meshIntegrators;

  ///< Integrators over lower-dimensional subdomains of the mesh.
  std::vector<IntegratorSubMesh*> _submeshIntegrators;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Formulation(const Formulation&); ///< Not implemented
  const Formulation& operator=(const Formulation&); ///< Not implemented

}; // Formulation

#endif // pylith_problems_formulation_hh


// End of file 
