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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/topology/MeshRefiner.hh
 *
 * @brief Object for tri3 refinement of cells.
 */

#if !defined(pylith_topology_meshrefiner_hh)
#define pylith_topology_meshrefiner_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

// RefineTri3 --------------------------------------------------------
/// Object for tri3 refinement of cells.
class ALE::MeshRefiner
{ // MeshRefiner
  typedef IMesh<> mesh_type;
  typedef mesh_type::point_type point_type;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor
   *
   * @param mesh Finite-element mesh.
   */
  MeshRefiner(void);

  /// Destructor
  ~MeshRefiner(void);

  /** Refine mesh.
   *
   * @param newMesh New mesh.
   * @param mesh Current mesh.
   * @param refiner Cell refiner.
   */
  void refine(const Obj<mesh_type>& newMesh, 
	      const Obj<mesh_type>& mesh, 
	      CellRefinerTri3& refiner);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Refine mesh.
   *
   * @param newMesh New mesh.
   * @param mesh Current mesh.
   * @param refiner Cell refiner.
   */
  void _refine(const Obj<mesh_type>& newMesh, 
	       const Obj<mesh_type>& mesh, 
	       CellRefinerTri3& refiner);
  
  /** Refine mesh with a censored depth.
   *
   * @param newMesh New mesh.
   * @param mesh Current mesh.
   * @param refiner Cell refiner.
   */
  void _refineCensored(const Obj<mesh_type>& newMesh, 
		       const Obj<mesh_type>& mesh, 
		       CellRefinerTri3& refiner);

  /** Stratify mesh.
   *
   * @param mesh Mesh to stratify.
   */
  void _stratify(const Obj<mesh_type>& mesh);

  /** Calculate new overlap.
   *
   * @param newMesh New (refined) mesh.
   * @param mesh Current (unrefined) mesh with overlap.
   */
  void _calcNewOverlap(const Obj<mesh_type>& newMesh,
		       const Obj<mesh_type>& mesh);
  
  /** Create integer sections in new mesh.
   *
   * :WARNING: Only implemented for integer sections containing vertices.
   *
   * @param newMesh New (refined) mesh.
   * @param mesh Current (unrefined) mesh with integer sections.
   */
  void _createIntSections(const Obj<mesh_type>& newMesh,
			  const Obj<mesh_type>& mesh,
			  CellRefinerTri3& refiner);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  MeshOrder* _orderOldMesh; ///< Order of entities in old mesh.
  MeshOrder* _orderNewMesh; ///< Order of entities in new mesh.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  MeshRefiner(const MeshRefiner&); ///< Not implemented
  const MeshRefiner& operator=(const MeshRefiner&); ///< Not implemented

}; // MeshRefiner

#endif // pylith_topology_meshrefiner_hh

 
// End of file 
