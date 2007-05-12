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
 * @file pylith/topology/GeometryHex.hh
 *
 * @brief C++ implementation of cell geometry calculations for 3-D
 * hexahedral cell.
 */

#if !defined(pylith_topology_geometryhex_hh)
#define pylith_topology_geometryhex_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace topology {
    class GeometryHex;
  } // topology
} // pylith

class pylith::topology::GeometryHex : public CellGeometry
{ // GeometryHex

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryHex(void);

  /// Default destructor.
  ~GeometryHex(void);

  /** Compute Jacobian at location in cell.
   *
   * @param jacobian Jacobian at location.
   * @param vertices Coordinates of vertices of cell.
   * @param location Location in reference cell at which to compute Jacobian.
   */
  void jacobian(double_array* jacobian,
		const double_array& vertices,
		const double_array& location) const;

}; // GeometryHex

#endif // pylith_topology_geometryhex_hh


// End of file
