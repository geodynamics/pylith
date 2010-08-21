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
 * @file modulesrc/topology/MeshOps.i
 *
 * @brief Python interface to C++ MeshOps.
 */

%apply(int* IN_ARRAY1, int DIM1) {
  (int* const materialIds, const int numMaterials)
  };
%inline %{
  /** Check to make sure material id of every cell matches the id of
   *  one of the materials.
   *
   * @param mesh PETSc mesh.
   * @param materialIds Array of ids for all materials and cohesive
   * cell interfaces.
   * @param numMaterials Size of array.
   */
  void
  MeshOps_checkMaterialIds(const pylith::topology::Mesh& mesh,
			   int* const materialIds,
			   const int numMaterials) {
    pylith::topology::MeshOps::checkMaterialIds(mesh, 
						materialIds, numMaterials);
  } // checkMaterialIds
%}
%clear(int* const materialIds, const int numMaterials);

%inline %{
  int
  MeshOps_numMaterialCells(const pylith::topology::Mesh& mesh, int materialId) {
    pylith::topology::MeshOps::numMaterialCells(mesh,materialId);
  } // numMaterialCells
%}

// End of file 
