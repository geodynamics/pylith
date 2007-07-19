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
 * @file pylith/feassemble/GeometryTri2D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 2-D
 * triangular cell.
 *
 * Reference cell:
 *
 * Vertex   x     y
 *    0   -1.0  -1.0
 *    1   +1.0  -1.0
 *    2   -1.0  +1.0
 */

#if !defined(pylith_feassemble_geometrytri2d_hh)
#define pylith_feassemble_geometrytri2d_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace feassemble {
    class GeometryTri2D;
  } // feassemble
} // pylith

class pylith::feassemble::GeometryTri2D : public CellGeometry
{ // GeometryTri2D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryTri2D(void);

  /// Default destructor.
  ~GeometryTri2D(void);

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

  GeometryTri2D(const GeometryTri2D&); ///< Not implemented
  const GeometryTri2D& operator=(const GeometryTri2D&); ///< Not implemented

}; // GeometryTri2D

#endif // pylith_feassemble_geometrytri2d_hh


// End of file
