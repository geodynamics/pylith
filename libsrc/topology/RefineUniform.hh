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
 * @file pylith/topology/RefineUniform.hh
 *
 * @brief Object for managing uniform global mesh refinement.
 */

#if !defined(pylith_topology_refineuniform_hh)
#define pylith_topology_refineuniform_hh

// Include directives ---------------------------------------------------
#include "MeshRefiner.hh" // ISA MeshRefiner

// RefineUniform --------------------------------------------------------
class pylith::topology::RefineUniform : public MeshRefiner
{ // RefineUniform
  friend class TestRefineUniform; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  RefineUniform(void);

  /// Destructor
  ~RefineUniform(void);

  /** Refine mesh.
   *
   * @param newMesh Refined mesh (result).
   * @param mesh Mesh to refine.
   * @param levels Number of levels to refine.
   * @param fields Solution fields.
   */
  void refine(Mesh* const newMesh,
	      const Mesh& mesh,
	      const int levels =1,
	      const SolutionFields* fields =0);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /** Refine tet4 mesh.
   *
   * @param newMesh Refined mesh (result).
   * @param mesh Mesh to refine.
   * @param levels Number of levels to refine.
   */
  void _refineTet4(Mesh* const newMesh,
		   const Mesh& mesh,
		   const int levels =1);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  RefineUniform(const RefineUniform&); ///< Not implemented
  const RefineUniform& operator=(const RefineUniform&); ///< Not implemented

}; // RefineUniform

#endif // pylith_topology_refineuniform_hh


// End of file 
