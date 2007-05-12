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
 * @file pylith/topology/GeometryTet.hh
 *
 * @brief C++ implementation of cell geometry calculations for 3-D
 * tetrahedral cell.
 */

#if !defined(pylith_topology_geometrytet_hh)
#define pylith_topology_geometrytet_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace topology {
    class GeometryTet;
  } // topology
} // pylith

class pylith::topology::GeometryTet
{ // GeometryTet

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryTet(void);

  /// Default destructor.
  ~GeometryTet(void);

  /** Compute Jacobian at location in cell.
   *
   * @param jacobian Jacobian at location.
   * @param vertices Coordinates of vertices of cell.
   * @param location Location in reference cell at which to compute Jacobian.
   */
  void jacobian(double_array* jacobian,
		const double_array& vertices,
		const double_array& location) const;

}; // GeometryTet

#endif // pylith_topology_geometrytet_hh


// End of file
