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
 * @file pylith/topology/CellGeometry.hh
 *
 * @brief C++ abstract base class for cell geometry calculations.
 */

#if !defined(pylith_topology_cellgeometry_hh)
#define pylith_topology_cellgeometry_hh

#include "pylith/utils/arrayfwd.hh" // USES double_array

namespace pylith {
  namespace topology {
    class CellGeometry;
  } // topology
} // pylith

class pylith::topology::CellGeometry
{ // CellGeometry

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Default constructor.
  CellGeometry(void);

  /// Default destructor.
  virtual
  ~CellGeometry(void);

  /** Compute Jacobian at location in cell.
   *
   * @param jacobian Jacobian at location.
   * @param vertices Coordinates of vertices of cell.
   * @param location Location in reference cell at which to compute Jacobian.
   */
  virtual
  void jacobian(double_array* jacobian,
		const double_array& vertices,
		const double_array& location) const = 0;

}; // CellGeometry

#endif // pylith_topology_cellgeometry_hh


// End of file
