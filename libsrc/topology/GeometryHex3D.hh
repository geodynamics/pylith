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
 * @file pylith/topology/GeometryHex3D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 3-D
 * hexahedral cell.
 */

#if !defined(pylith_topology_geometryhex3d_hh)
#define pylith_topology_geometryhex3d_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace topology {
    class GeometryHex3D;
  } // topology
} // pylith

class pylith::topology::GeometryHex3D : public CellGeometry
{ // GeometryHex3D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryHex3D(void);

  /// Default destructor.
  ~GeometryHex3D(void);

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

}; // GeometryHex3D

#endif // pylith_topology_geometryhex3d_hh


// End of file
