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

#include "pylith/utils/petscfwd.h" // USES PetscVec, PetscMat, PetscSNES

#include "pylith/utils/array.hh" // HASA std::vector


// Formulation ----------------------------------------------------------
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

  /** Reform system residual.
   *
   * @param tmpResidualVec Temporary PETSc vector for residual.
   * @param tmpSolveSolnVec Temporary PETSc vector for solution.
   */
  void reformResidual(const PetscVec* tmpResidualVec =0,
		      const PetscVec* tmpSolveSolnVec =0);
  
  /* Reform system Jacobian.
   *
   * @param tmpSolveSolnVec Temporary PETSc vector for solution.
   */
  void reformJacobian(const PetscVec* tmpSolveSolnVec =0);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  double _t; ///< Current time (nondimensional).
  double _dt; ///< Current time step (nondimensional).
  topology::Jacobian* _jacobian; ///< Handle to Jacobian of system.
  topology::SolutionFields* _fields; ///< Handle to solution fields for system.

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
