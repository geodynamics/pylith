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
 * @file libsrc/topology/MeshCoords.hh
 *
 * @brief C++ helper class for accessing coordinates in a finite-element mesh.
 */

#if !defined(pylith_topology_coordinates_hh)
#define pylith_topology_coordinates_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // HASA PetscVec, PetscSection

#include <petscdmmesh.hh>

// Coordinates ----------------------------------------------------------
/** @brief Helper class for accessing coordinates in a finite-element mesh.
 */
template<typename mesh_type>
class pylith::topology::Coordinates
{ // MeshCoords
  friend class TestCoordinates; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  Coordinates(const mesh_type& mesh);

  /// Default destructor
  ~Coordinates(void);

  /** Get the local coordinates array associated with the local PETSc Vec.
   *
   * Must call restoreLocalArray() afterwards.
   * 
   * @returns Local array.
   */
  PetscScalar* getLocalArray(void) const;

  /** Restore local coordinates array associated with the local PETSc Vec.
   *
   * @preq Must be preceded by call to getLocalArray().
   *
   * @param a Local array.
   */
  void restoreLocalArray(PetscScalar** a) const;

  /** Get fiber dimension of coordinates for point.
   *
   * @preq Must call cacheSection().
   *
   * @param point Point in mesh.
   * @returns Fiber dimension.
   */
  PetscInt sectionDof(const PetscInt point) const;

  /** Get offset into coordinates array for point.
   *
   * @preq Must call cacheSection().
   *
   * @param point Point in mesh.
   * @returns Offset.
   */
  PetscInt sectionOffset(const PetscInt point) const;

  /** Get coordinates array associated with closure.
   *
   * @param coordsCell Array of coordinates for cell.
   * @param coordsSize Size of coordinates array.
   * @param cell Finite-element cell.
   */
  void getClosure(PetscScalar** coordsCell,
		  PetscInt* coordsSize,
		  const PetscInt cell) const;

  /** Restore coordinates array associated with closure.
   *
   * @param coordsCell Array of coordinates for cell.
   * @param coordsSize Size of coordinates array.
   * @param cell Finite-element cell.
   */
  void restoreClosure(PetscScalar** coordsCell,
		      PetscInt* coordsSize,
		      const PetscInt cell) const;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  const mesh_type& _mesh;

  PetscDM _dm; ///< Cached PETSc dm for mesh.
  PetscVec _localVec; ///< Cached local PETSc Vec.
  PetscSection _section; ///< Cached PETSc section.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Coordinates(const Coordinates&); ///< Not implemented
  const Coordinates& operator=(const Coordinates&); ///< Not implemented

}; // Coordinates

#include "Coordinates.icc"

#endif // pylith_topology_coordinates_hh


// End of file
