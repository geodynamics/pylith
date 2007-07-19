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
 * @file pylith/feassemble/GeometryPoint2D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 0-D
 * point cell.
 */

#if !defined(pylith_feassemble_geometrypoint2d_hh)
#define pylith_feassemble_geometrypoint2d_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace feassemble {
    class GeometryPoint2D;
  } // feassemble
} // pylith

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

  GeometryPoint2D(const GeometryPoint2D&); ///< Not implemented
  const GeometryPoint2D& operator=(const GeometryPoint2D&); ///< Not implemented

}; // GeometryPoint2D

#endif // pylith_feassemble_geometrypoint2d_hh


// End of file
