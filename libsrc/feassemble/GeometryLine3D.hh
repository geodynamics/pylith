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
 * @file pylith/feassemble/GeometryLine3D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 1-D
 * line cell.
 */

#if !defined(pylith_feassemble_geometryline3d_hh)
#define pylith_feassemble_geometryline3d_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace feassemble {
    class GeometryLine3D;
  } // feassemble
} // pylith

class pylith::feassemble::GeometryLine3D : public CellGeometry
{ // GeometryLine3D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryLine3D(void);

  /// Default destructor.
  ~GeometryLine3D(void);

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
   * @param coordsGlobal Coordinates in global coordinate system.
   * @param coordsRef Coordinates in reference cell.
   * @param vertices Array of cell vertices in global coordinates.
   * @param dim Dimension of global coordinate system.
   */
  void coordsRefToGlobal(double* coordsGlobal,
			 const double* coordsRef,
			 const double* vertices,
			 const int dim) const;

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
   * @param location Location in reference cell at which to compute Jacobian.
   * @param dim Dimension of coordinate system.
   */
  void jacobian(double* jacobian,
		double* det,
		const double* vertices,
		const double* location,
		const int dim) const;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  GeometryLine3D(const GeometryLine3D&); ///< Not implemented
  const GeometryLine3D& operator=(const GeometryLine3D&); ///< Not implemented

}; // GeometryLine3D

#endif // pylith_feassemble_geometryline3d_hh


// End of file
