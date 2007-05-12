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
 * @file pylith/topology/GeometryTri.hh
 *
 * @brief C++ implementation of cell geometry calculations for 2-D
 * triangular cell.
 */

#if !defined(pylith_topology_geometrytri_hh)
#define pylith_topology_geometrytri_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace topology {
    class GeometryTri;
  } // topology
} // pylith

class pylith::topology::GeometryTri : public CellGeometry
{ // GeometryTri

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryTri(void);

  /// Default destructor.
  ~GeometryTri(void);

  /** Compute Jacobian at location in cell.
   *
   * @param jacobian Jacobian at location.
   * @param vertices Coordinates of vertices of cell.
   * @param location Location in reference cell at which to compute Jacobian.
   */
  void jacobian(double_array* jacobian,
		const double_array& vertices,
		const double_array& location) const;

}; // GeometryTri

#endif // pylith_topology_geometrytri_hh


// End of file
