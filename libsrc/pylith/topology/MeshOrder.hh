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
 * @file libsrc/topology/MeshOrder.hh
 *
 * @brief Object for managing order of mesh entities.
 *
 * Entities are stored in contiguous ranges.
 */

#if !defined(pylith_topology_meshorder_hh)
#define pylith_topology_meshorder_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include <petscdmmesh.hh> // HASA ALE::IMesh

// MeshOrder ------------------------------------------------------------
/// Object for managing order of mesh entities.
class ALE::MeshOrder
{ // MeshOrder
  typedef int point_type;
  typedef ALE::IMesh<PetscInt,PetscScalar> mesh_type;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  MeshOrder(void);

  /// Destructor
  ~MeshOrder(void);

  /** Determine order from pre-existing mesh.
   *
   * @param mesh Finite-element mesh.
   */
  void initialize(const ALE::Obj<mesh_type>& mesh);

  /** Set range for normal cells.
   *
   * @param min Minimum cell label.
   * @param max Maximum cell label.
   */
  void cellsNormal(const point_type min, const point_type max); 

  /** Get range for normal cells.
   *
   * @returns Interval with range of normal cells.
   */
  const ALE::Interval<point_type>& cellsNormal(void) const;

  /** Set range for normal vertices.
   *
   * @param min Minimum vertex label.
   * @param max Maximum vertex label.
   */
  void verticesNormal(const point_type min, const point_type max); 

  /** Get range for normal vertices.
   *
   * @returns Interval with range of normal vertices.
   */
  const ALE::Interval<point_type>& verticesNormal(void) const;

  /** Set range for censored cells.
   *
   * @param min Minimum cell label.
   * @param max Maximum cell label.
   */
  void cellsCensored(const point_type min, const point_type max); 

  /** Get range for censored cells.
   *
   * @returns Interval with range of censored cells.
   */
  const ALE::Interval<point_type>& cellsCensored(void) const;

  /** Set range for censored vertices.
   *
   * @param min Minimum vertex label.
   * @param max Maximum vertex label.
   */
  void verticesCensored(const point_type min, const point_type max); 

  /** Get range for censored vertices.
   *
   * @returns Interval with range of censored vertices.
   */
  const ALE::Interval<point_type>& verticesCensored(void) const;


// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  ALE::Interval<point_type> _cellsNormal;
  ALE::Interval<point_type> _verticesNormal;
  ALE::Interval<point_type> _verticesCensored;
  ALE::Interval<point_type> _cellsCensored;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  MeshOrder(const MeshOrder&); ///< Not implemented
  const MeshOrder& operator=(const MeshOrder&); ///< Not implemented

}; // MeshOrder

#endif // pylith_topology_meshorder_hh

 
// End of file 
