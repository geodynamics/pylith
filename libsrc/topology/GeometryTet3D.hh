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
 * @file pylith/topology/GeometryTet3D.hh
 *
 * @brief C++ implementation of cell geometry calculations for 3-D
 * tetrahedral cell.
 */

#if !defined(pylith_topology_geometrytet3d_hh)
#define pylith_topology_geometrytet3d_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace topology {
    class GeometryTet3D;
  } // topology
} // pylith

class pylith::topology::GeometryTet3D : public CellGeometry
{ // GeometryTet3D

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryTet3D(void);

  /// Default destructor.
  ~GeometryTet3D(void);

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

}; // GeometryTet3D

#endif // pylith_topology_geometrytet3d_hh


// End of file
