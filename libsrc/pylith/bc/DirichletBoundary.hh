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

/** @file libsrc/bc/DirichletBoundary.hh
 *
 * @brief C++ implementation of Dirichlet (prescribed values at
 * degrees of freedom) boundary conditions with points on a boundary.
 */

#if !defined(pylith_bc_dirichletboundary_hh)
#define pylith_bc_dirichletboundary_hh

// Include directives ---------------------------------------------------
#include "DirichletBC.hh" // ISA DirichletBC

// DirichletBoundary ----------------------------------------------------
/// @brief Dirichlet (prescribed values at degrees of freedom) boundary
/// conditions with points on a boundary.
class pylith::bc::DirichletBoundary : public DirichletBC
{ // class DirichletBoundary
  friend class TestDirichletBoundary; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  DirichletBoundary(void);

  /// Destructor.
  ~DirichletBoundary(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Initialize boundary condition.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Vertical direction (somtimes used in 3-D problems).
   */
  void initialize(const topology::Mesh& mesh,
		  const PylithScalar upDir[3]);

  /** Get boundary mesh.
   *
   * @return Boundary mesh.
   */
  const topology::Mesh& boundaryMesh(void) const;

  /** Get vertex field with BC information.
   *
   * @param name Name of field.
   * @param fields Solution fields.
   *
   * @returns Field over vertices.
   */
  const topology::Field& vertexField(const char* name,
				     const topology::SolutionFields& fields);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Get vertex vector field with BC information.
   *
   * @param name Name of parameter field.
   * @param label Label for buffer field.
   * @param scale Scale used to dimensionalize field.
   *
   * @returns Field over vertices.
   */
  const topology::Field& _bufferVector(const char* name,
				       const char* label,
				       const PylithScalar scale);

  /** Get vertex scalar field with BC information.
   *
   * @param name Name of parameter field.
   * @param label Label for buffer field.
   * @param scale Scale used to dimensionalize field.
   *
   * @returns Field over vertices.
   */
  const topology::Field& _bufferScalar(const char* name,
				       const char* label,
				       const PylithScalar scale);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  topology::Mesh* _boundaryMesh; ///< Boundary mesh.

  /// Fields manager (holds temporary field for output).
  topology::Fields* _outputFields;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  DirichletBoundary(const DirichletBoundary&);

  /// Not implemented
  const DirichletBoundary& operator=(const DirichletBoundary&);

}; // class DirichletBoundary

#include "DirichletBoundary.icc" // inline methods

#endif // pylith_bc_dirichletboundary_hh


// End of file 
