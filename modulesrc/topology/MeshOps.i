// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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


// End of file 
