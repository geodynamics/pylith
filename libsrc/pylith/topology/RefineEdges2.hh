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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/topology/RefineEdges2.hh
 *
 * @brief Object for refinement of cells via refinement of edges
 * comprised of two vertices.
 */

#if !defined(pylith_topology_refineedges2_hh)
#define pylith_topology_refineedges2_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include <list> // USES std::pair

// RefineEdges2 ------------------------------------------------------
/// Object for tri3 refinement of cells.
class ALE::RefineEdges2
{ // RefineEdges2
protected:

  typedef ALE::IMesh<PetscInt,PetscScalar> mesh_type;
  typedef mesh_type::point_type point_type;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor
   *
   * @param mesh Finite-element mesh.
   */
  RefineEdges2(const mesh_type& mesh);

  /// Destructor
  ~RefineEdges2(void);

  /** Set coordinates of new vertices.
   *
   * @param newCoordsSection Coordinates of vertices in new mesh.
   * @param oldCoordsSection Coordinates of vertices in original mesh.
   */
  void setCoordsNewVertices(const ALE::Obj<mesh_type::real_section_type>& newCoordsSection,
			    const ALE::Obj<mesh_type::real_section_type>& oldCoordsSection);

  /** Add space for new vertices in group.
   *
   * @param newGroup Group in refine mesh.
   * @param oldGroup Group in original mesh.
   */
  void groupAddNewVertices(const ALE::Obj<mesh_type::int_section_type>& newGroup,
			   const ALE::Obj<mesh_type::int_section_type>& oldGroup);

  /** Set new vertices in group.
   *
   * @param newGroup Group in refine mesh.
   * @param oldGroup Group in original mesh.
   */
  void groupSetNewVertices(const ALE::Obj<mesh_type::int_section_type>& newGroup,
			   const ALE::Obj<mesh_type::int_section_type>& oldGroup);

  /** Add new vertices to label.
   *
   * @param newMesh Mesh with refined cells.
   * @param oldMesh Original mesh.
   * @param labelName Name of label.
   */
  void labelAddNewVertices(const ALE::Obj<mesh_type>& newMesh,
			   const ALE::Obj<mesh_type>& oldMesh,
			   const char* labelName);

  /** Calculate new overlap.
   *
   * @param newMesh New (refined) mesh.
   * @param orderNewMesh Order in new mesh.
   * @param oldMesh Current (unrefined) mesh with overlap.
   * @param orderOldMesh Order in old mesh.
   */
  void overlapAddNewVertices(const Obj<mesh_type>& newMesh,
			     const MeshOrder& orderNewMesh,
			     const Obj<mesh_type>& oldMesh,
			     const MeshOrder& orderOldMesh);
  
// PROTECTED TYPEDEFS ///////////////////////////////////////////////////
protected :

  template<typename Point>
  class Edge : public std::pair<Point, Point> {
  public:
    Edge() : std::pair<Point, Point>() {};
    Edge(const Point l) : std::pair<Point, Point>(l, l) {};
    Edge(const Point l, const Point r) : std::pair<Point, Point>(l, r) {};
    ~Edge() {};
    friend std::ostream& operator<<(std::ostream& stream, const Edge& edge) {
      stream << "(" << edge.first << ", " << edge.second << ")";
      return stream;
    };
  };

  typedef Edge<point_type> EdgeType;
  typedef std::map<EdgeType, point_type> edge_map_type;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  const mesh_type& _mesh;
  edge_map_type _edgeToVertex;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  RefineEdges2(const RefineEdges2&); ///< Not implemented
  const RefineEdges2& operator=(const RefineEdges2&); ///< Not implemented

}; // RefineEdges2

#endif // pylith_topology_refineedges2_hh

 
// End of file 
