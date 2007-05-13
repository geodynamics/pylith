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
 * @file pylith/topology/GeometryQuad3D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 2-D
 * quadrilateral cell in 3-D.
 */

#if !defined(pylith_topology_geometryquad3d_hh)
#define pylith_topology_geometryquad3d_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace topology {
    class GeometryQuad3D;
  } // topology
} // pylith

class pylith::topology::GeometryQuad3D : public CellGeometry
{ // GeometryQuad3D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryQuad3D(void);

  /// Default destructor.
  ~GeometryQuad3D(void);

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

}; // GeometryQuad3D

#endif // pylith_topology_geometryquad3d_hh


// End of file
