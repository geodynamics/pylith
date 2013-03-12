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
 * @file libsrc/topology/SubMeshVisitor.hh
 *
 * @brief C++ helper class for accessing field values at points in a
 * submesh within a finite-element mesh.
 *
 * This visitor is used to access values associated with a submesh
 * when the field is defined over the entire mesh. This is why a
 * submesh and index set are passed to the constructor.
 *
 * Use the PointVisitorMesh object when the field and mesh/submesh are
 * defined over the same set of points (i.e., field over a submesh or field
 * of a mesh).
 */

#if !defined(pylith_topology_submeshvisitor_hh)
#define pylith_topology_submeshvisitor_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // HASA PetscVec, PetscSection

#include <petscdmmesh.hh>

// SubMeshVisitor -------------------------------------------------------
/** @brief Helper class for accessing field values at points in a
 *  finite-element mesh.
 */
class pylith::topology::SubMeshVisitor
{ // SubMeshVisitor
  friend class TestSubMeshVisitor; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param field Field associated with visitor.
   * @param submeshIS Submesh index set associated with visitor.
   */
  SubMeshVisitor(const Field<Mesh>& field,
		 const SubMeshIS& submeshIS);

  /// Default destructor
  ~SubMeshVisitor(void);

  /* Initialize cached data.
   *
   * @param submeshIS Submesh index set associated with visitor.
   */
  void initialize(const SubMeshIS& submeshIS);

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
   * @param valuesCell Array of values for cell.
   * @param valuesSize Size of values array.
   * @param cell Finite-element cell.
   */
  void getClosure(PetscScalar** valuesCell,
		  PetscInt* valuesSize,
		  const PetscInt cell) const;

  /** Restore array of values associated with closure.
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

  const Field<Mesh>& _field;

  PetscDM _dm; ///< Cached PETSc dm for mesh.
  PetscVec _localVec; ///< Cached local PETSc Vec.
  PetscSection _section; ///< Cached PETSc subsection.
  PetscScalar* _localArray; ///< Cached local array

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  SubMeshVisitor(const SubMeshVisitor&); ///< Not implemented
  const SubMeshVisitor& operator=(const SubMeshVisitor&); ///< Not implemented

}; // SubMeshVisitor


// SubMeshIS ------------------------------------------------------------
/// Index set associated with submesh.
class pylith::topology::SubMeshIS
{ // SubMeshIS
  friend class TestSubMeshIS; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   * 
   * @param submesh Submesh associated with index set.
   */
  SubMeshIS(const SubMesh& submesh);

  /// Default destructor.
  ~SubMeshIS(void);

  /// Deallocate.
  void deallocate(void);

  /** Get the submesh.
   *
   * @returns Submesh.
   */
  const SubMesh& submesh(void) const;

  /** Get PETSc index set.
   *
   * @returns PETSc index set.
   */
  PetscIS indexSet(void) const;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  const SubMesh& _submesh;
  PetscIS _indexSet; ///< PETSc index set.

}; // SubMeshIS


#include "SubMeshVisitor.icc"

#endif // pylith_topology_submeshvisitor_hh


// End of file
