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

// CohesiveTopology -----------------------------------------------------
class pylith::faults::CohesiveTopology
{ // class CohesiveTopology
public :
  typedef std::set<SieveMesh::point_type> PointSet;
  typedef std::vector<sieve_type::point_type> PointArray;
  typedef std::pair<sieve_type::point_type, int> oPoint_type;
  typedef std::vector<oPoint_type>  oPointArray;

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :
  /** Create the fault mesh.
   *
   * @param fault Finite-element mesh of fault (output)
   * @param mesh Finite-element mesh
   * @param faultVertices Vertices assocated with faces of cells defining 
   *   fault surface
   */
  static
  void createFault(topology::SubMesh* ifault,
                   ALE::Obj<ALE::Mesh>& faultBd,
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
  void create(topology::SubMesh* ifault,
              const ALE::Obj<ALE::Mesh>& faultBd,
              const topology::Mesh& mesh,
              const ALE::Obj<topology::Mesh::IntSection>& groupField,
              const int materialId,
              const bool constraintCell =false);

  /** Create (distributed) fault mesh from cohesive cells.
   *
   * @param fault Finite-element mesh of fault (output).
   * @param cohesiveToFault Mapping of cohesive cell to fault mesh cell.
   * @param mesh Finite-element mesh.
   * @param materialId Material id for cohesive elements.
   * @param constraintCell True if creating cells constrained with 
   *   Lagrange multipliers that require extra vertices, false otherwise.
   */
  static
  void createParallel(topology::SubMesh* ifault,
		      std::map<Mesh::point_type, Mesh::point_type>* cohesiveToFault,
		      const topology::Mesh& mesh,
		      const int materialId,
		      const bool constraintCell =false);

}; // class CohesiveTopology

#endif // pylith_faults_cohesivetopology_hh


// End of file 
