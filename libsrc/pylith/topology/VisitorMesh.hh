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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/topology/VisitorMesh.hh
 *
 * @brief C++ helper class for accessing field and matrix values at
 * points in a finite-element mesh.
 */

#if !defined(pylith_topology_visitormesh_hh)
#define pylith_topology_visitormesh_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // HASA PetscVec, PetscSection
#include "pylith/utils/arrayfwd.hh" // USES scalar_array

// VecVisitorMesh ----------------------------------------------------------
/** @brief Helper class for accessing field values at points in a
 *  finite-element mesh.
 */
class pylith::topology::VecVisitorMesh
{ // VecVisitorMesh
  friend class TestVecVisitorMesh; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /** Constructor with field over a mesh.
   *
   * @param field Field over a mesh.
   */
  VecVisitorMesh(const Field& field);

  /// Default destructor
  ~VecVisitorMesh(void);

  /** Initialize using field over a mesh or submesh.
   *
   * @param field Field over a mesh/submesh.
   */
  void initialize(const Field& field);

  /// Clear cached data.
  void clear(void);
  
  /** Get the array of values associated with the local PETSc Vec.
   * 
   * @returns Array of values.
   */
  PetscScalar* localArray(void) const;

  /** Get the PETSc section.
   * 
   * @returns PETSc section.
   */
  PetscSection petscSection(void) const;

  /** Get the local PETSc Vec.
   * 
   * @returns PETSc Vec.
   */
  PetscVec localVec(void) const;

  /** Get fiber dimension for values at point.
   *
   * @param point Point in mesh.
   * @returns Fiber dimension.
   */
  PetscInt sectionDof(const PetscInt point) const;

  /** Get fiber dimension for constraints at point.
   *
   * @param point Point in mesh.
   * @returns Fiber dimension.
   */
  PetscInt sectionConstraintDof(const PetscInt point) const;

  /** Get offset into values array for point.
   *
   * @param point Point in mesh.
   * @returns Offset.
   */
  PetscInt sectionOffset(const PetscInt point) const;

  /** Get array of values associated with closure.
   *
   * @pre Must be followed by call to restoreClosure().
   *
   * @param valuesCell Array of values for cell.
   * @param valuesSize Size of values array.
   * @param cell Finite-element cell.
   */
  void getClosure(PetscScalar** valuesCell,
		  PetscInt* valuesSize,
		  const PetscInt cell) const;

  /** Get array of values associated with closure.
   *
   * @param values Array of values for cell.
   * @param cell Finite-element cell.
   */
  void getClosure(scalar_array* values,
		  const PetscInt cell) const;

  /** Restore array of values associated with closure.
   *
   * @pre Must be preceded by call to getClosure().
   *
   * @param valuesCell Array of values for cell.
   * @param valuesSize Size of values array.
   * @param cell Finite-element cell.
   */
  void restoreClosure(PetscScalar** valuesCell,
		      PetscInt* valuesSize,
		      const PetscInt cell) const;

  /** Set values associated with closure.
   *
   * @param valuesCell Array of values for cell.
   * @param valuesSize Size of values array.
   * @param cell Finite-element cell.
   * @param mode Mode for inserting values.
   */
  void setClosure(const PetscScalar* valuesCell,
		  const PetscInt valuesSize,
		  const PetscInt cell,
		  const InsertMode mode) const;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  PetscDM _dm; ///< Cached PETSc dm for mesh.
  PetscVec _localVec; ///< Cached local PETSc Vec.
  PetscSection _section; ///< Cached PETSc section.
  PetscScalar* _localArray; ///< Cached local array

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  VecVisitorMesh(const VecVisitorMesh&); ///< Not implemented
  const VecVisitorMesh& operator=(const VecVisitorMesh&); ///< Not implemented

}; // VecVisitorMesh


// MatVisitorMesh ----------------------------------------------------------
/** @brief Helper class for accessing field values at points in a
 *  finite-element mesh.
 */
class pylith::topology::MatVisitorMesh
{ // MatVisitorMesh
  friend class TestMatVisitorMesh; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param mat PETSc matrix.
   * @param field Field associated with matrix layout.
   */
  MatVisitorMesh(const PetscMat mat,
		 const Field& field);

  /// Default destructor
  ~MatVisitorMesh(void);

  // Initialize.
  void initialize(void);

  /// Clear cached data.
  void clear(void);
  
  /** Set values associated with closure.
   *
   * @param valuesCell Array of values for cell.
   * @param valuesSize Size of values array.
   * @param cell Finite-element cell.
   * @param mode Mode for inserting values.
   */
  void setClosure(const PetscScalar* valuesCell,
		  const PetscInt valuesSize,
		  const PetscInt cell,
		  const InsertMode mode) const;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  const PetscMat _mat; ///< Cached PETSc matrix.
  PetscDM _dm; ///< Cached PETSc dm for mesh.
  PetscSection _section; ///< Cached PETSc section.
  PetscSection _globalSection; ///< Cached PETSc global section.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  MatVisitorMesh(const MatVisitorMesh&); ///< Not implemented
  const MatVisitorMesh& operator=(const MatVisitorMesh&); ///< Not implemented

}; // MatVisitorMesh

#include "VisitorMesh.icc"

#endif // pylith_topology_visitormesh_hh


// End of file
