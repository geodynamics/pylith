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
 * @file pylith/feassemble/GeometryQuad3D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 2-D
 * quadrilateral cell in 3-D.
 *
 * Reference cell:
 *
 * 3 -- 2
 * |    |
 * |    |
 * 0 -- 1
 *
 * Vertex   x     y
 *    0   -1.0  -1.0
 *    1   +1.0  -1.0
 *    2   +1.0  +1.0
 *    3   -1.0  +1.0
 */

#if !defined(pylith_feassemble_geometryquad3d_hh)
#define pylith_feassemble_geometryquad3d_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace feassemble {
    class GeometryQuad3D;
  } // feassemble
} // pylith

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
  GeometryQuad3D(const GeometryQuad3D& g);

}; // GeometryQuad3D

#endif // pylith_feassemble_geometryquad3d_hh


// End of file
