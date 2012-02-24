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
 * @file libsrc/topology/RefineFace4Edges2.hh
 *
 * @brief Object for refinement of cells via refinement of edges
 * comprised of two vertices.
 */

#if !defined(pylith_topology_refineface4edges2_hh)
#define pylith_topology_refineface4edges2_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include <list> // USES std::pair

// RefineFace4Edges2 ------------------------------------------------------
/// Object for tri3 refinement of cells.
class ALE::RefineFace4Edges2
{ // RefineFace4Edges2
protected:

  typedef IMesh<PetscInt,PetscScalar> mesh_type;
  typedef mesh_type::point_type point_type;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor
   *
   * @param mesh Finite-element mesh.
   */
  RefineFace4Edges2(const mesh_type& mesh);

  /// Destructor
  ~RefineFace4Edges2(void);

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
    Edge(void) : std::pair<Point, Point>() {};
    Edge(const Point l) : std::pair<Point, Point>(l, l) {};
    Edge(const Point l, const Point r) : std::pair<Point, Point>(l, r) {};
    ~Edge(void) {};
    friend std::ostream& operator<<(std::ostream& stream, const Edge& edge) {
      stream << "(" << edge.first << ", " << edge.second << ")";
      return stream;
    };
  };
  typedef Edge<point_type> EdgeType;
  typedef std::map<EdgeType, point_type> edge_map_type;

  template<typename Point>
  class Face {
  public:
    Face(void) {};
    Face(const Point p) {
      points[0] = p;
      points[1] = p;
      points[2] = p;
      points[3] = p;
    };
    Face(const Point p0,	 
	 const Point p1,
	 const Point p2,
	 const Point p3) {
      points[0] = p0;
      points[1] = p1;
      points[2] = p2;
      points[3] = p3;      
    };
    ~Face(void) {};
    friend bool operator==(const Face<Point>& a, 
			   const Face<Point>& b) {
      const bool result = 
	a.points[0] == b.points[0] &&
	a.points[1] == b.points[1] &&
	a.points[2] == b.points[2] &&
	a.points[3] == b.points[3];
      return result;
    };
    friend bool operator<(const Face<Point>& a, 
			  const Face<Point>& b) {
      if (a.points[0] < b.points[0]) {
	return true;
      } else if (a.points[0] == b.points[0]) {
	if (a.points[1] < b.points[1]) {
	  return true;
	} else if (a.points[1] == b.points[1]) {
	  if (a.points[2] < b.points[2]) {
	    return true;
	  } else if (a.points[2] == b.points[2]) {
	    if (a.points[3] < b.points[3]) {
	      return true;
	    } // if
	  } // if/else
	} // if/else
      } // if/else
    
      return false;
    };
    friend std::ostream& operator<<(std::ostream& stream, 
				    const Face<Point>& face) {
      stream << "(" << face.points[0]
	     << ", " << face.points[1]
	     << ", " << face.points[2]
	     << ", " << face.points[3]
	     << ")";
      return stream;
    };
  public:
    int points[4];
  };
  template<typename Point>
  class FaceCmp { // for compatibility with gcc 3.4.6
  public :
    bool operator() (const Face<Point>& a, 
		     const Face<Point>& b) const {
      if (a.points[0] < b.points[0]) {
	return true;
      } else if (a.points[0] == b.points[0]) {
	if (a.points[1] < b.points[1]) {
	  return true;
	} else if (a.points[1] == b.points[1]) {
	  if (a.points[2] < b.points[2]) {
	    return true;
	  } else if (a.points[2] == b.points[2]) {
	    if (a.points[3] < b.points[3]) {
	      return true;
	    } // if
	  } // if/else
	} // if/else
      } // if/else
      
      return false;
    };
  };
  typedef Face<point_type> FaceType;
  typedef std::map<FaceType, point_type, FaceCmp<point_type> > face_map_type;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  const mesh_type& _mesh;
  edge_map_type _edgeToVertex;
  face_map_type _faceToVertex;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  RefineFace4Edges2(const RefineFace4Edges2&); ///< Not implemented
  const RefineFace4Edges2& operator=(const RefineFace4Edges2&); ///< Not implemented

}; // RefineFace4Edges2

#endif // pylith_topology_refineface4edges2_hh

 
// End of file 
