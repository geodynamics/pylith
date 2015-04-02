// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/BCIntegratorSubMesh.hh
 *
 * @brief C++ abstract base class for BoundaryCondition object with
 * boundary condition applied at a simply-connected boundary (submesh).
 *
 * Interface definition for boundary conditions applied to a
 * simply-connected boundary (submesh).
 */

#if !defined(pylith_bc_bcintegratorsubmesh_hh)
#define pylith_bc_bcintegratorsubmesh_hh

// Include directives ---------------------------------------------------
#include "BoundaryCondition.hh" // ISA BoundaryCondition

#include "pylith/topology/topologyfwd.hh" // HOLDSA Fields, Mesh

#include "pylith/feassemble/Integrator.hh" // ISA Integrator

// BCIntegratorSubMesh ----------------------------------------------
/// @brief Abstract base classs for BoundaryCondition object with
/// boundary condition applied at a simple-connect surface (submesh).
class pylith::bc::BCIntegratorSubMesh : public BoundaryCondition,
					public feassemble::Integrator
{ // class BCIntegratorSubMesh
  friend class TestBCIntegratorSubMesh; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  BCIntegratorSubMesh(void);

  /// Destructor.
  virtual
  ~BCIntegratorSubMesh(void);
  
  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Get parameter fields.
   *
   * @returns Parameter fields.
   */
  const topology::Fields* parameterFields(void) const;

  /** Get boundary mesh.
   *
   * @return Boundary mesh.
   */
  const topology::Mesh& boundaryMesh(void) const;

  /** Get mesh labels for submesh associated with applied forces.
   *
   * @param mesh Finite-element mesh.
   */
  void createSubMesh(const topology::Mesh& mesh);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  topology::Mesh* _boundaryMesh; ///< Boundary mesh.
  topology::SubMeshIS* _submeshIS; ///< Cache index set for submesh.
  topology::VecVisitorSubMesh* _residualVisitor; ///< Cache residual field visitor.
  topology::MatVisitorSubMesh* _jacobianMatVisitor; ///< Cache jacobian  matrix visitor.
  topology::VecVisitorSubMesh* _jacobianVecVisitor; ///< Cache jacobian field visitor.

  /// Parameters for boundary condition.
  topology::Fields* _parameters;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  BCIntegratorSubMesh(const BCIntegratorSubMesh&);

  /// Not implemented
  const BCIntegratorSubMesh& operator=(const BCIntegratorSubMesh&);

}; // class BCIntegratorSubMesh

#include "BCIntegratorSubMesh.icc" // inline methods

#endif // pylith_bc_bcintegratorsubmesh_hh


// End of file 
