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
 * @file libsrc/feassemble/GeometryPoint1D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 0-D
 * point cell.
 */

#if !defined(pylith_feassemble_geometrypoint1d_hh)
#define pylith_feassemble_geometrypoint1d_hh

// Include directives ---------------------------------------------------
#include "CellGeometry.hh" // ISA CellGeometry

// GeometryPoint1D ------------------------------------------------------
/// Cell geometry calculations for 0-D point cell.
class pylith::feassemble::GeometryPoint1D : public CellGeometry
{ // GeometryPoint1D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryPoint1D(void);

  /// Default destructor.
  ~GeometryPoint1D(void);

  /** Create a copy of geometry.
   *
   * @returns Copy of geometry.
   */
  CellGeometry* clone(void) const;

  /** Get cell geometry for lower dimension cell.
   *
   * @returns Pointer to cell geometry object corresponding to next
   * lower dimension, NULL if there is no lower dimension object.
   */
  CellGeometry* geometryLowerDim(void) const;

  /** Transform coordinates in reference cell to global coordinates.
   *
   * @param ptsGlobal Array of points in global coordinate system.
   * @param ptsRef Array of points in reference cell.
   * @param vertices Array of cell vertices in global coordinates.
   * @param dim Dimension of global coordinate system.
   * @param npts Number of points to transform.
   */
  void ptsRefToGlobal(PylithScalar* ptsGlobal,
		      const PylithScalar* ptsRef,
		      const PylithScalar* vertices,
		      const int dim,
		      const int npts =1) const;

  /** Compute Jacobian at location in cell.
   *
   * @param jacobian Jacobian at location.
   * @param det Determinant of Jacobian at location.
   * @param vertices Coordinates of vertices of cell.
   * @param numVertices Number of vertices in cell.
   * @param spaceDim Spatial dimension of coordinates.
   * @param location Location in reference cell at which to compute Jacobian.
   * @param cellDim Dimension of reference cell.
   */
  void jacobian(scalar_array* jacobian,
		PylithScalar* det,
		const PylithScalar* vertices,
		const int numVertices,
		const int spaceDim,
		const PylithScalar* location,
		const int cellDim) const;

  /** Compute Jacobian at location in cell.
   *
   * @param jacobian Jacobian at location.
   * @param det Determinant of Jacobian at location.
   * @param vertices Coordinates of vertices of cell.
   * @param ptsRef Points in reference cell at which to compute Jacobian.
   * @param dim Dimension of coordinate system.
   * @param npts Number of points to transform.
   */
  void jacobian(PylithScalar* jacobian,
		PylithScalar* det,
		const PylithScalar* vertices,
		const PylithScalar* ptsRef,
		const int dim,
		const int npts =1) const;

  /** Compute minimum width across cell.
   *
   * @param coordinatesCell Coordinates of vertices in cell.
   * @param numVertices Number of vertices in cell.
   * @param spaceDim Coordinate dimension.
   *
   * @returns Minimum width across cell.
   */
  PylithScalar minCellWidth(const PylithScalar* coordinatesCell,
			    const int numVertices,
			    const int spaceDim) const;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  GeometryPoint1D(const GeometryPoint1D&); ///< Not implemented
  const GeometryPoint1D& operator=(const GeometryPoint1D&); ///< Not implemented

}; // GeometryPoint1D

#endif // pylith_feassemble_geometrypoint1d_hh


// End of file
