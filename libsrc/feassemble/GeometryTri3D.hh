// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/**
 * @file pylith/feassemble/GeometryTri3D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 2-D
 * triangular cell.
 */

#if !defined(pylith_feassemble_geometrytri3d_hh)
#define pylith_feassemble_geometrytri3d_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace feassemble {
    class GeometryTri3D;
  } // feassemble
} // pylith

class pylith::feassemble::GeometryTri3D : public CellGeometry
{ // GeometryTri3D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryTri3D(void);

  /// Default destructor.
  ~GeometryTri3D(void);

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

// PROTECTED ////////////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param g Geometry to copy.
   */
  GeometryTri3D(const GeometryTri3D& g);

}; // GeometryTri3D

#endif // pylith_feassemble_geometrytri3d_hh


// End of file
