// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file libsrc/topology/Jacobian.hh
 *
 * @brief Jacobian of the system as a PETSc sparse matrix.
 */

#if !defined(pylith_topology_jacobian_hh)
#define pylith_topology_jacobian_hh


// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // HOLDSA PetscMat

#include <string> // USES std::string

// Jacobian -------------------------------------------------------------
/// Jacobian of the system as a PETSc sparse matrix.
class pylith::topology::Jacobian
{ // Jacobian
  friend class TestJacobian; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param field Field associated with mesh and solution of the problem.
   * @param matrixType Type of PETSc sparse matrix.
   * @param blockOkay True if okay to use block size equal to fiberDim
   * (all or none of the DOF at each point are constrained).
   */
  Jacobian(const Field& field,
           const char* matrixType ="aij",
           const bool blockOkay =false);

  /// Destructor.
  ~Jacobian(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Get PETSc matrix.
   *
   * @returns PETSc sparse matrix.
   */
  const PetscMat matrix(void) const;

  /** Get PETSc matrix.
   *
   * @returns PETSc sparse matrix.
   */
  PetscMat matrix(void);

  /** Get matrix type.
   *
   * @returns Matrix type.
   */
  const char* matrixType(void) const;

  /** Assemble matrix.
   *
   * @param mode Assembly mode.
   */
  void assemble(const char* mode);

  /// Set entries in matrix to zero (retain structure).
  void zero(void);

  /// View matrix to stdout.
  void view(void) const;

  /** Write matrix to binary file.
   *
   * @param filename Name of file.
   * @param comm MPI communicator.
   */
  void write(const char* filename,
             const MPI_Comm comm);

  /// Verify symmetry of matrix. For debugger purposes only.
  void verifySymmetry(void) const;

  /** Get flag indicating if sparse matrix values have been
   * updated. This pertains to the values being changed, NOT the
   * pattern of nonzero entries.
   *
   * @returns True if values have been updated/altered.
   */
  bool valuesChanged(void) const;

  /// Reset flag indicating if sparse matrix values have been updated.
  void resetValuesChanged(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  PetscMat _matrix; ///< Sparse matrix for Jacobian of problem.

  bool _valuesChanged; ///< Sparse matrix values have been updated.

  std::string _type; ///< String associated with matrix type.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Jacobian(const Jacobian&); ///< Not implemented
  const Jacobian& operator=(const Jacobian&); ///< Not implemented

}; // Jacobian

#endif // pylith_topology_jacobian_hh


// End of file 
