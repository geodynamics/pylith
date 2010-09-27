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
 * @file libsrc/topology/CellRefinerTet4.hh
 *
 * @brief Object for tet4 refinement of cells.
 */

#if !defined(pylith_topology_cellrefinertet4_hh)
#define pylith_topology_cellrefinertet4_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include <list> // USES std::pair

// CellRefinerTet4 ------------------------------------------------------
/// Object for tet4 refinement of cells.
class ALE::CellRefinerTet4
{ // CellRefinerTet4
  typedef IMesh<> mesh_type;
  typedef mesh_type::point_type point_type;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor
   *
   * @param mesh Finite-element mesh.
   */
  CellRefinerTet4(const mesh_type& mesh);

  /// Destructor
  ~CellRefinerTet4(void);

  /** Get number of refined cells for each original cell.
   *
   * @param cell Original cell.
   *
   * @returns Number of refined cells.
   */
  int numNewCells(const point_type cell);

  /** Split cell into smaller cells of same type.
   *
   * @param cell Original cell.
   * @param cone Vertices in cell (original mesh).
   * @param coneSize Number of vertices in cell.
   * @param curNewVertex Value for next new vertex.
   */
  void splitCell(const point_type cell,
		 const point_type cone[],
		 const int coneSize,
		 point_type* curNewVertex);

  /** Get refined cells.
   *
   * @param cells Vertices in refined cells (refined mesh).
   * @param numCells Number of refined cells.
   * @param cone Vertices in cell (original mesh).
   * @param coneSize Number of vertices in cell.
   * @param orderOldMesh Order in old mesh.
   * @param orderNewMesh Order in new mesh.
   */
  void getNewCells(const point_type** cells,
		   int* numCells,
		   const point_type cell,
		   const point_type cone[],
		   const int coneSize,
		   const MeshOrder& orderOldMesh,
		   const MeshOrder& orderNewMesh);

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

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private :

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

// PRIVATE ENUMS ////////////////////////////////////////////////////////
private :

  enum CellEnum { 
    TETRAHEDRON, // Normal tetrahedral cell
    TRIANGLE_COHESIVE_LAGRANGE, // Cohesive cell with Lagrange multiplier vertices
  };

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Get cell type.
   *
   * @param cell Cell in original mesh.
   * @returns Cell type.
   */
  CellEnum _cellType(const point_type cell);
  
  /** Get edges of triangular cell.
   *
   * @param edges Edges of cell.
   * @param numEdges Number of edges.
   * @param cone Vertices in cell (original mesh).
   * @param coneSize Number of vertices in cell.
   */
  void _edges_TETRAHEDRON(const EdgeType** edges,
		       int* numEdges,
		       const point_type cone[],
		       const int coneSize);
  
  /** Get edges of line cohesive cell with Lagrange multipler vertices.
   *
   * @param edges Edges of cell.
   * @param numEdges Number of edges.
   * @param cone Vertices in cell (original mesh).
   * @param coneSize Number of vertices in cell.
   */
  void _edges_TRIANGLE_COHESIVE_LAGRANGE(const EdgeType** edges,
					 int* numEdges,
					 const point_type cone[],
					 const int coneSize);
  
  /** Get new cells from refinement of a triangular cell.
   *
   * @param cells Vertices in refined cells (refined mesh).
   * @param numCells Number of refined cells.
   * @param cone Vertices in cell (original mesh).
   * @param coneSize Number of vertices in cell.
   * @param coneVertexOffset Offset for cone vertices.
   */
  void _newCells_TETRAHEDRON(const point_type** cells,
			  int *numCells,
			  const point_type cone[],
			  const int coneSize,
			  const int coneVertexOffset);
  
  /** Get new cells from refinement of a line cohseive cell with
   * Lagrange multiplier vertices.
   *
   * @param cells Vertices in refined cells (refined mesh).
   * @param numCells Number of refined cells.
   * @param cone Vertices in cell (original mesh).
   * @param coneSize Number of vertices in cell.
   * @param coneVertexOffsetNormal Offset for normal cone vertices.
   * @param coneVertexOffset Offset for censored cone vertices.
   */
  void _newCells_TRIANGLE_COHESIVE_LAGRANGE(const point_type** cells,
					    int *numCells,
					    const point_type cone[],
					    const int coneSize,
					    const int coneVertexOffsetNormal,
					    const int coneVertexOffsetCensored);
  
// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  const mesh_type& _mesh;
  edge_map_type _edgeToVertex;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  CellRefinerTet4(const CellRefinerTet4&); ///< Not implemented
  const CellRefinerTet4& operator=(const CellRefinerTet4&); ///< Not implemented

}; // CellRefinerTet4

#endif // pylith_topology_cellrefinertet4_hh

 
// End of file 
