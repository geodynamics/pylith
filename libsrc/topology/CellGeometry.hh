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
  
  /** Default constructor.
   *
   * @param cellDim Dimension of cell.
   * @param spaceDim Dimension of coordinate space.
   * @param numCorners Number of corners in cell.
   */
  CellGeometry(const int cellDim,
	       const int spaceDim,
	       const int numCorners);

  /// Default destructor.
  virtual
  ~CellGeometry(void);

  /** Get dimension of cell.
   *
   * @returns Spatial dimension of cell.
   */
  int cellDim(void) const;

  /** Get dimension of coordinate space.
   *
   * @returns Dimension of coordinate space.
   */
  int spaceDim(void) const;

  /** Get number of vertices in cell.
   *
   * @returns Number of vertices in cell.
   */
  int numCorners(void) const;

  /** Get cell geometry for lower dimension cell.
   *
   * @returns Pointer to cell geometry object corresponding to next
   * lower dimension, NULL if there is no lower dimension object.
   */
  virtual
  CellGeometry* geometryLowerDim(void) const = 0;

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

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  int _cellDim; ///< Dimension of cell.
  int _spaceDim; ///< Dimension of coordinate space.
  int _numCorners; ///< Number of corners in cell.
  
}; // CellGeometry

#include "CellGeometry.icc" // inline methods

#endif // pylith_topology_cellgeometry_hh


// End of file
