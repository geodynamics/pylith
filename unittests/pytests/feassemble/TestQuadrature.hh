// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/pytests/feassemble/TestQuadrature.hh
 *
 * @brief C++ TestQuadrature object
 *
 * Helper class for unit testing of Python Quadrature.
 */

#if !defined(pylith_feassemble_pytestquadrature_hh)
#define pylith_feassemble_pytestquadrature_hh

#include "pylith/feassemble/Quadrature.hh"

/// Namespace for spatialdata package
namespace pylith {
  namespace feassemble {
    class TestQuadrature;
  } // feassemble
} // pylith

/// Helper class for unit testing of Python Quadrature
class pylith::feassemble::TestQuadrature
{ // class TestQuadrature

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Get number of dimensions in reference cell.
   *
   * @returns Number of dimensions
   */
  static int cellDim(const Quadrature& q);

  /** Get number of vertices in cell.
   *
   * @returns Number of vertices
   */
  static int numBasis(const Quadrature& q);

  /** Get number of quadrature points.
   *
   * @returns Number of points
   */
  static int numQuadPts(const Quadrature& q);

  /** Get number of dimensions in coordinates of cell vertices.
   *
   * @returns Number of dimensions
   */
  static int spaceDim(const Quadrature& q);

  /** Get vertices in reference cell.
   *
   * @returns Array of coordinates of vertices in reference cell.
   */
  static const double* vertices(const Quadrature& q);

  /** Get basis functions evaluated at quadrature points.
   *
   * @returns Array of basis functions evaluated at quadrature points
   */
  static const double* basis(const Quadrature& q);

  /** Get derivatives of basis functions evaluated at quadrature points.
   *
   * @returns Array of derivatives of basis fns evaluated at quad pts
   */
  static const double* basisDeriv(const Quadrature& q);

  /** Get coordinates of quadrature points in reference cell.
   *
   * @returns Array of coordinates of quadrature points
   */
  static const double* quadPtsRef(const Quadrature& q);

  /** Get weights of quadrature points.
   *
   * @returns Array of weights of quadrature points
   */
  static const double* quadWts(const Quadrature& q);

}; // class TestQuadrature

#include "TestQuadrature.icc" // inline methods

#endif // pylith_feassemble_pytestquadrature_hh

// End of file 
