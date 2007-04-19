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

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class CohesiveTopology;
  } // faults
} // pylith

/// C++ object to manage creation of cohesive cells.
class pylith::faults::CohesiveTopology
{ // class Fault

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Create cohesive cells.
   *
   * @param mesh Finite-element mesh
   * @param faultVertices Vertices assocated with faces of cells defining 
   *   fault surface
   */
  static
  void create(const ALE::Obj<Mesh>& mesh,
              const ALE::Obj<Mesh::int_section_type>& groupField);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Get number of vertices on face.
   *
   * @param cell Finite-element cell
   * @param mesh Finite-element mesh
   *
   * @returns Number of vertices on cell face
   */
  static
  unsigned int _numFaceVertices(const Mesh::point_type& cell,
		       const ALE::Obj<Mesh>& mesh);

}; // class CohesiveTopology

#endif // pylith_faults_cohesivetopology_hh


// End of file 
