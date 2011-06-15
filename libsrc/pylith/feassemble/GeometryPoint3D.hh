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
 * @file libsrc/feassemble/GeometryPoint3D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 0-D
 * point cell in 3-D space.
 */

#if !defined(pylith_feassemble_geometrypoint3d_hh)
#define pylith_feassemble_geometrypoint3d_hh

// Include directives ---------------------------------------------------
#include "CellGeometry.hh" // ISA CellGeometry

// GeometryPoint3D ------------------------------------------------------
/// Cell geometry calculations for 0-D line cell in 3-D space.
class pylith::feassemble::GeometryPoint3D : public CellGeometry
{ // GeometryPoint3D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryPoint3D(void);

  /// Default destructor.
  ~GeometryPoint3D(void);

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

  GeometryPoint3D(const GeometryPoint3D&); ///< Not implemented
  const GeometryPoint3D& operator=(const GeometryPoint3D&); ///< Not implemented

}; // GeometryPoint3D

#endif // pylith_feassemble_geometrypoint3d_hh


// End of file
