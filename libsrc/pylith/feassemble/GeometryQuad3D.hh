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
 * @file libsrc/feassemble/GeometryQuad3D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 2-D
 * quadrilateral cell in 3-D.
 */

#if !defined(pylith_feassemble_geometryquad3d_hh)
#define pylith_feassemble_geometryquad3d_hh

// Include directives ---------------------------------------------------
#include "CellGeometry.hh" // ISA CellGeometry

// GeometryQuad3D -------------------------------------------------------
/** @brief C++ implementation of cell geometry calculations for 2-D
 * quadrilateral cell in 3-D.
 *
 * Reference cell:
@verbatim
3 -- 2
|    |
|    |
0 -- 1

Vertex   x     y
   0   -1.0  -1.0
   1   +1.0  -1.0
   2   +1.0  +1.0
   3   -1.0  +1.0
@endverbatim
 */
class pylith::feassemble::GeometryQuad3D : public CellGeometry
{ // GeometryQuad3D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryQuad3D(void);

  /// Default destructor.
  ~GeometryQuad3D(void);

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
  void ptsRefToGlobal(double* ptsGlobal,
		      const double* ptsRef,
		      const double* vertices,
		      const int dim,
		      const int npts =1) const;

  /** Compute Jacobian at location in cell.
   *
   * @param jacobian Jacobian at location.
   * @param det Determinant of Jacobian at location.
   * @param vertices Coordinates of vertices of cell.
   * @param location Location in reference cell at which to compute Jacobian.
   */
  void jacobian(double_array* jacobian,
		double* det,
		const double_array& vertices,
		const double_array& location) const;

  /** Compute Jacobian at location in cell.
   *
   * @param jacobian Jacobian at location.
   * @param det Determinant of Jacobian at location.
   * @param vertices Coordinates of vertices of cell.
   * @param ptsRef Points in reference cell at which to compute Jacobian.
   * @param dim Dimension of coordinate system.
   * @param npts Number of points to transform.
   */
  void jacobian(double* jacobian,
		double* det,
		const double* vertices,
		const double* ptsRef,
		const int dim,
		const int npts =1) const;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  GeometryQuad3D(const GeometryQuad3D&); ///< Not implemented
  const GeometryQuad3D& operator=(const GeometryQuad3D&); ///< Not implemented

}; // GeometryQuad3D

#endif // pylith_feassemble_geometryquad3d_hh


// End of file
