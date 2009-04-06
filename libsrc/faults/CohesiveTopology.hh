// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/CohesiveTopology.hh
 *
 * @brief C++ object to manage creation of cohesive cells.
 */

#if !defined(pylith_faults_cohesivetopology_hh)
#define pylith_faults_cohesivetopology_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/topology/Mesh.hh" // USES Mesh::IntSection
#include "pylith/utils/sievetypes.hh" // USE ALE::Obj

// CohesiveTopology -----------------------------------------------------
class pylith::faults::CohesiveTopology
{ // class CohesiveTopology

private :
  typedef pylith::topology::Mesh::SieveMesh::point_type point_type;

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Create the fault mesh.
   *
   * @param faultMesh Finite-element mesh of fault (output).
   * @param faultBoundary Finite-element mesh of fault boundary (output).
   * @param mesh Finite-element mesh of domain.
   * @param faultVertices Vertices assocated with faces of cells defining 
   *   fault surface
   */
  static
  void createFault(topology::SubMesh* faultMesh,
		   ALE::Obj<ALE::Mesh>& faultBoundary,
		   const topology::Mesh& mesh,
		   const ALE::Obj<topology::Mesh::IntSection>& groupField,
		   const bool flipFault =false);

  /** Create cohesive cells.
   *
   * @param fault Finite-element mesh of fault (output)
   * @param mesh Finite-element mesh
   * @param materialId Material id for cohesive elements.
   * @param constraintCell True if creating cells constrained with 
   *   Lagrange multipliers that require extra vertices, false otherwise
   */
  static
  void create(topology::Mesh* mesh,
	      const topology::SubMesh& faultMesh,
              const ALE::Obj<ALE::Mesh>& faultBoundary,
              const ALE::Obj<topology::Mesh::IntSection>& groupField,
              const int materialId,
              const bool constraintCell =false);

  /** Create (distributed) fault mesh from cohesive cells.
   *
   * @param faultMesh Finite-element mesh of fault (output).
   * @param cohesiveToFault Mapping of cohesive cell to fault mesh
   *   cell (output).
   * @param mesh Finite-element mesh.
   * @param materialId Material id for cohesive elements.
   * @param constraintCell True if creating cells constrained with 
   *   Lagrange multipliers that require extra vertices, false otherwise.
   */
  static
  void createFaultParallel(topology::SubMesh* faultMesh,
			   std::map<point_type, point_type>* cohesiveToFault,
			   const topology::Mesh& mesh,
			   const int materialId,
			   const bool constraintCell =false);

}; // class CohesiveTopology

#endif // pylith_faults_cohesivetopology_hh


// End of file 
