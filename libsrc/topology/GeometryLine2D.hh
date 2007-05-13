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
 * @file pylith/topology/GeometryLine2D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 1-D
 * line cell in 2-D space.
 */

#if !defined(pylith_topology_geometryline2d_hh)
#define pylith_topology_geometryline2d_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace topology {
    class GeometryLine2D;
  } // topology
} // pylith

class pylith::topology::GeometryLine2D : public CellGeometry
{ // GeometryLine2D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryLine2D(void);

  /// Default destructor.
  ~GeometryLine2D(void);

  /** Get cell geometry for lower dimension cell.
   *
   * @returns Pointer to cell geometry object corresponding to next
   * lower dimension, NULL if there is no lower dimension object.
   */
  CellGeometry* geometryLowerDim(void) const;

  /** Compute Jacobian at location in cell.
   *
   * @param jacobian Jacobian at location.
   * @param vertices Coordinates of vertices of cell.
   * @param location Location in reference cell at which to compute Jacobian.
   */
  void jacobian(double_array* jacobian,
		const double_array& vertices,
		const double_array& location) const;

}; // GeometryLine2D

#endif // pylith_topology_geometryline2d_hh


// End of file
