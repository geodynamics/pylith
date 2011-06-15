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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/GeometryPoint2D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 0-D
 * point cell in 2-D space.
 */

#if !defined(pylith_feassemble_geometrypoint2d_hh)
#define pylith_feassemble_geometrypoint2d_hh

// Include directives ---------------------------------------------------
#include "CellGeometry.hh" // ISA CellGeometry

// GeometryPoint2D ------------------------------------------------------
/// Cell geometry calculations for 0-D line cell in 2-D space.
class pylith::feassemble::GeometryPoint2D : public CellGeometry
{ // GeometryPoint2D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryPoint2D(void);

  /// Default destructor.
  ~GeometryPoint2D(void);

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

  GeometryPoint2D(const GeometryPoint2D&); ///< Not implemented
  const GeometryPoint2D& operator=(const GeometryPoint2D&); ///< Not implemented

}; // GeometryPoint2D

#endif // pylith_feassemble_geometrypoint2d_hh


// End of file
