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
 * @file libsrc/topology/CellRefinerQuad4.hh
 *
 * @brief Object for quad4 refinement of cells.
 */

#if !defined(pylith_topology_cellrefinerquad4_hh)
#define pylith_topology_cellrefinerquad4_hh

// Include directives ---------------------------------------------------
#include "RefineFace4Edges2.hh" // ISA RefineFace4Edges2

// CellRefinerQuad4 ------------------------------------------------------
/// Object for quad4 refinement of cells.
class ALE::CellRefinerQuad4 : public RefineFace4Edges2
{ // CellRefinerQuad4
// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor
   *
   * @param mesh Finite-element mesh.
   */
  CellRefinerQuad4(const mesh_type& mesh);

  /// Destructor
  ~CellRefinerQuad4(void);

  /** Get number of refined cells for each original cell.
   *
   * @param cell Original cell.
   *
   * @returns Number of refined cells.
   */
  int numNewCells(const point_type cell);

  /** Split cell into smaller cells of same type. Do not create
   * censored vertices on censored cells.
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

  /** Split cell into smaller cells of same type. Create only censored
   * vertices on censored cells.
   *
   * @param cell Original cell.
   * @param cone Vertices in cell (original mesh).
   * @param coneSize Number of vertices in cell.
   * @param curNewVertex Value for next new vertex.
   */
  void splitCellUncensored(const point_type cell,
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

// PRIVATE ENUMS ////////////////////////////////////////////////////////
private :

  enum CellEnum { 
    QUADRILATERAL, // Normal quadrilateral cell
    LINE_COHESIVE_LAGRANGE, // Cohesive cell with Lagrange multiplier vertices
  };

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Get cell type.
   *
   * @param cell Cell in original mesh.
   * @returns Cell type.
   */
  CellEnum _cellType(const point_type cell);
  
  /** Get edges of quadrilateral cell.
   *
   * @param edges Edges of cell.
   * @param numEdges Number of edges.
   * @param cone Vertices in cell (original mesh).
   * @param coneSize Number of vertices in cell.
   */
  void _edges_QUADRILATERAL(const EdgeType** edges,
			    int* numEdges,
			    const point_type cone[],
			    const int coneSize);
  
  /** Get edges of line cohesive cell with Lagrange multipler vertices.
   *
   * @param edges Edges of cell.
   * @param numEdges Number of edges.
   * @param cone Vertices in cell (original mesh).
   * @param coneSize Number of vertices in cell.
   * @param uncensored True if including edges with censored vertices.
   */
  void _edges_LINE_COHESIVE_LAGRANGE(const EdgeType** edges,
				     int* numEdges,
				     const point_type cone[],
				     const int coneSize,
				     const bool uncensored =false);
  
  /** Get faces of quadrilateral cell.
   *
   * @param faces Faces of cell.
   * @param numFaces Number of faces.
   * @param cone Vertices in cell (original mesh).
   * @param coneSize Number of vertices in cell.
   */
  void _faces_QUADRILATERAL(const FaceType** faces,
			    int* numFaces,
			    const point_type cone[],
			    const int coneSize);
  
  /** Get new cells from refinement of a quadrilateral cell.
   *
   * @param cells Vertices in refined cells (refined mesh).
   * @param numCells Number of refined cells.
   * @param cone Vertices in cell (original mesh).
   * @param coneSize Number of vertices in cell.
   * @param coneVertexOffset Offset for cone vertices.
   */
  void _newCells_QUADRILATERAL(const point_type** cells,
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
  void _newCells_LINE_COHESIVE_LAGRANGE(const point_type** cells,
					int *numCells,
					const point_type cone[],
					const int coneSize,
					const int coneVertexOffsetNormal,
					const int coneVertexOffsetCensored);
  
// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  static const int _edgesSize;
  static const int _facesSize;
  static const int _cellsSize;
  EdgeType _edges[4]; ///< Buffer for edges
  FaceType _faces[1]; ///< Buffer for faces
  point_type _cells[16]; ///< Buffer for cells
  

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  CellRefinerQuad4(const CellRefinerQuad4&); ///< Not implemented
  const CellRefinerQuad4& operator=(const CellRefinerQuad4&); ///< Not implemented

}; // CellRefinerQuad4

#endif // pylith_topology_cellrefinerquad4_hh

 
// End of file 
