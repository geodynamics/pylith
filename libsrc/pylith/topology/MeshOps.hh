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
 * @file libsrc/topology/MeshOps.hh
 *
 * @brief Simple operations on a Mesh object.
 */

#if !defined(pylith_topology_meshops_hh)
#define pylith_topology_meshops_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "spatialdata/units/unitsfwd.hh" // forward declarations

// MeshOps --------------------------------------------------------------
/// Simple operations on a Mesh object.
class pylith::topology::MeshOps
{ // MeshOps
  friend class TestMeshOps; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Create DMPlex mesh.
   *
   * @param mesh Finite-element mesh.
   * @param dim Dimension associated with mesh cells.
   * @param comm MPI communicator for mesh.
   * @param label Label for mesh.
   */
  static
  void createDMMesh(Mesh* const mesh,
		    const int dim =3,
		    const MPI_Comm& comm =PETSC_COMM_WORLD,
		    const char* label ="domain");

  /** Nondimensionalize the finite-element mesh.
   *
   * @param mesh Finite-element mesh.
   * @param normalizer Nondimensionalizer.
   */
  static
  void nondimensionalize(Mesh* const mesh,
			 const spatialdata::units::Nondimensional& normalizer);

  /** Check topology of mesh.
   *
   * @param mesh Finite-element mesh.
   */
  static
  void checkTopology(const Mesh& mesh);

  /** Check to make sure material id of every cell matches the id of
   *  one of the materials.
   *
   * @param mesh Finite-element mesh.
   * @param materialIds Array of ids for all materials and cohesive
   * cell interfaces.
   * @param numMaterials Size of array.
   */
  static
  void checkMaterialIds(const Mesh& mesh,
			int* const materialIds,
			const int numMaterials);

  /** Get number of cells associated with material.
   *
   * @param mesh Finite-element mesh.
   * @param materialId Id of material.
   * @returns Number of cells.
   */
  static
  int numMaterialCells(const Mesh& mesh,
		       int materialId);
  

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  MeshOps(void); ///< Not Implemented
  MeshOps(const MeshOps&); ///< Not implemented
  const MeshOps& operator=(const MeshOps&); ///< Not implemented


}; // MeshOps

#endif // pylith_topology_meshops_hh


// End of file 
