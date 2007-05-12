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
 * @file pylith/topology/GeometryLine.hh
 *
 * @brief C++ implementation of cell geometry calculations for 1-D
 * line cell.
 */

#if !defined(pylith_topology_geometryline_hh)
#define pylith_topology_geometryline_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace topology {
    class GeometryLine;
  } // topology
} // pylith

class pylith::topology::GeometryLine
{ // GeometryLine

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryLine(void);

  /// Default destructor.
  ~GeometryLine(void);

  /** Compute Jacobian at location in cell.
   *
   * @param jacobian Jacobian at location.
   * @param vertices Coordinates of vertices of cell.
   * @param location Location in reference cell at which to compute Jacobian.
   */
  void jacobian(double_array* jacobian,
		const double_array& vertices,
		const double_array& location) const;

}; // GeometryLine

#endif // pylith_topology_geometryline_hh


// End of file
