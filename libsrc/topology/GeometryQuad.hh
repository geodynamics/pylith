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
 * @file pylith/topology/GeometryQuad.hh
 *
 * @brief C++ implementation of cell geometry calculations for 2-D
 * quadrilateral cell.
 */

#if !defined(pylith_topology_geometryquad_hh)
#define pylith_topology_geometryquad_hh

#include "CellGeometry.hh" // ISA CellGeometry

namespace pylith {
  namespace topology {
    class GeometryQuad;
  } // topology
} // pylith

class pylith::topology::GeometryQuad : public CellGeometry
{ // GeometryQuad

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  GeometryQuad(void);

  /// Default destructor.
  ~GeometryQuad(void);

  /** Compute Jacobian at location in cell.
   *
   * @param jacobian Jacobian at location.
   * @param vertices Coordinates of vertices of cell.
   * @param location Location in reference cell at which to compute Jacobian.
   */
  void jacobian(double_array* jacobian,
		const double_array& vertices,
		const double_array& location) const;

}; // GeometryQuad

#endif // pylith_topology_geometryquad_hh


// End of file
