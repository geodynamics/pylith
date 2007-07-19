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
 * @file pylith/feassemble/GeometryTet3D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 3-D
 * tetrahedral cell.
 *
 * Reference cell:
 *
 * Vertex   x     y     z
 *    0   -1.0  -1.0  -1.0
 *    1   +1.0  -1.0  -1.0
 *    2   -1.0  +1.0  -1.0
 *    3   -1.0  -1.0  +1.0
 */

#if !defined(pylith_feassemble_geometrytet3d_hh)
#define pylith_feassemble_geometrytet3d_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace feassemble {
    class GeometryTet3D;
  } // feassemble
} // pylith

class pylith::feassemble::GeometryTet3D : public CellGeometry
{ // GeometryTet3D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryTet3D(void);

  /// Default destructor.
  ~GeometryTet3D(void);

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

  GeometryTet3D(const GeometryTet3D&); ///< Not implemented
  const GeometryTet3D& operator=(const GeometryTet3D&); ///< Not implemented

}; // GeometryTet3D

#endif // pylith_feassemble_geometrytet3d_hh


// End of file
