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
 * @file pylith/feassemble/GeometryHex3D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 3-D
 * hexahedral cell.
 *
 * Reference cell:
 *
 * Vertex   x     y    z
 *    0   -1.0  -1.0  -1.0
 *    1   +1.0  -1.0  -1.0
 *    2   +1.0  +1.0  -1.0
 *    3   -1.0  +1.0  -1.0
 *    4   -1.0  -1.0  +1.0
 *    5   +1.0  -1.0  +1.0
 *    6   +1.0  +1.0  +1.0
 *    7   -1.0  +1.0  +1.0
 */

#if !defined(pylith_feassemble_geometryhex3d_hh)
#define pylith_feassemble_geometryhex3d_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace feassemble {
    class GeometryHex3D;
  } // feassemble
} // pylith

class pylith::feassemble::GeometryHex3D : public CellGeometry
{ // GeometryHex3D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryHex3D(void);

  /// Default destructor.
  ~GeometryHex3D(void);

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

}; // GeometryHex3D

#endif // pylith_feassemble_geometryhex3d_hh


// End of file
