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
 * @file pylith/feassemble/GeometryQuad2D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 2-D
 * quadrilateral cell.
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

#if !defined(pylith_feassemble_geometryquad2d_hh)
#define pylith_feassemble_geometryquad2d_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace feassemble {
    class GeometryQuad2D;
  } // feassemble
} // pylith

class pylith::feassemble::GeometryQuad2D : public CellGeometry
{ // GeometryQuad2D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryQuad2D(void);

  /// Default destructor.
  ~GeometryQuad2D(void);

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

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  GeometryQuad2D(const GeometryQuad2D&); ///< Not implemented
  const GeometryQuad2D& operator=(const GeometryQuad2D&); ///< Not implemented

}; // GeometryQuad2D

#endif // pylith_feassemble_geometryquad2d_hh


// End of file
