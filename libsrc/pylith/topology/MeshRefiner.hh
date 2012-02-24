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
 * @file libsrc/topology/MeshRefiner.hh
 *
 * @brief Object for refinement of cells.
 */

#if !defined(pylith_topology_meshrefiner_hh)
#define pylith_topology_meshrefiner_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/utils/sievetypes.hh" // USE Obj

// RefineTri3 --------------------------------------------------------
/// Object for refinement of cells.
template<typename cellrefiner_type>
class ALE::MeshRefiner
{ // MeshRefiner
  typedef IMesh<PetscInt,PetscScalar> mesh_type;
  typedef mesh_type::point_type point_type;
  typedef pylith::SieveFlexMesh SieveFlexMesh;

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
	      cellrefiner_type& refiner);

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
	       cellrefiner_type& refiner);
  
  /** Refine mesh with a censored depth.
   *
   * @param newMesh New mesh.
   * @param mesh Current mesh.
   * @param refiner Cell refiner.
   */
  void _refineCensored(const Obj<mesh_type>& newMesh, 
		       const Obj<mesh_type>& mesh, 
		       cellrefiner_type& refiner);

  /** Stratify mesh.
   *
   * @param mesh Mesh to stratify.
   */
  void _stratify(const Obj<mesh_type>& mesh);

  /** Create integer sections in new mesh.
   *
   * :WARNING: Only implemented for integer sections containing vertices.
   *
   * @param newMesh New (refined) mesh.
   * @param mesh Current (unrefined) mesh with integer sections.
   * @param refiner Cell refiner.
   */
  void _createIntSections(const Obj<mesh_type>& newMesh,
			  const Obj<mesh_type>& mesh,
			  cellrefiner_type& refiner);

  /** Create labels in new mesh.
   *
   * :WARNING: Only implemented for integer sections containing vertices.
   *
   * @param newMesh New (refined) mesh.
   * @param mesh Current (unrefined) mesh with integer sections.
   * @param refiner Cell refiner.
   */
  void _createLabels(const Obj<mesh_type>& newMesh,
		     const Obj<mesh_type>& mesh,
		     cellrefiner_type& refiner);

  /** Calculate new overlap.
   *
   * @param newMesh New (refined) mesh.
   * @param mesh Current (unrefined) mesh with overlap.
   * @param refiner Cell refiner.
   */
  void _calcNewOverlap(const Obj<mesh_type>& newMesh,
		       const Obj<mesh_type>& mesh,
		       cellrefiner_type& refiner);
  
// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  MeshOrder* _orderOldMesh; ///< Order of entities in old mesh.
  MeshOrder* _orderNewMesh; ///< Order of entities in new mesh.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  MeshRefiner(const MeshRefiner&); ///< Not implemented
  const MeshRefiner& operator=(const MeshRefiner&); ///< Not implemented

}; // MeshRefiner

#include "MeshRefiner.cc"

#endif // pylith_topology_meshrefiner_hh

 
// End of file 
