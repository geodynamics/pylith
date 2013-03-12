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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/topology/MeshVisitor.hh
 *
 * @brief C++ helper class for accessing field values at points in a
 * finite-element mesh.
 */

#if !defined(pylith_topology_meshvisitor_hh)
#define pylith_topology_meshvisitor_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // HASA PetscVec, PetscSection

#include <petscdmmesh.hh>

// MeshVisitor ----------------------------------------------------------
/** @brief Helper class for accessing field values at points in a
 *  finite-element mesh.
 */
template<typename field_type>
class pylith::topology::MeshVisitor
{ // MeshVisitor
  friend class TestMeshVisitor; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  MeshVisitor(const field_type& field);

  /// Default destructor
  ~MeshVisitor(void);

  /// Initialize cached data.
  void initialize(void);

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

  /** Get fiber dimension of coordinates for point.
   *
   * @param point Point in mesh.
   * @returns Fiber dimension.
   */
  PetscInt sectionDof(const PetscInt point) const;

  /** Get offset into coordinates array for point.
   *
   * @param point Point in mesh.
   * @returns Offset.
   */
  PetscInt sectionOffset(const PetscInt point) const;

  /** Get array of values associated with closure.
   *
   * @pre Must be followed by call to getClosure().
   *
   * @param valuesCell Array of values for cell.
   * @param valuesSize Size of values array.
   * @param cell Finite-element cell.
   */
  void getClosure(PetscScalar** valuesCell,
		  PetscInt* valuesSize,
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

  const field_type& _field;

  PetscDM _dm; ///< Cached PETSc dm for mesh.
  PetscVec _localVec; ///< Cached local PETSc Vec.
  PetscSection _section; ///< Cached PETSc section.
  PetscScalar* _localArray; ///< Cached local array

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  MeshVisitor(const MeshVisitor&); ///< Not implemented
  const MeshVisitor& operator=(const MeshVisitor&); ///< Not implemented

}; // MeshVisitor

#include "MeshVisitor.icc"

#endif // pylith_topology_meshvisitor_hh


// End of file
