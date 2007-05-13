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
 * @file unittests/libtests/topology/TestCellGeometry.hh
 *
 * @brief C++ TestCellGeometry object
 *
 * C++ unit testing for CellGeometry.
 */

#if !defined(pylith_topology_testcellgeometry_hh)
#define pylith_topology_testcellgeometry_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace topology {
    class TestCellGeometry;
    class CellGeometry;

    class CellGeomData;
  } // topology
} // pylith

/// C++ unit testing for TestCellGeometry
class pylith::topology::TestCellGeometry : public CppUnit::TestFixture
{ // class TestCellGeometry

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Test cellDim().
   *
   * @param geometry CellGeometry object.
   * @param data Geometry data.
   */
  void _testCellDim(const CellGeometry& geometry,
		    const CellGeomData& data);

  /** Test spaceDim().
   *
   * @param geometry CellGeometry object.
   * @param data Geometry data.
   */
  void _testSpaceDim(const CellGeometry& geometry,
		     const CellGeomData& data);

  /** Test numCorners().
   *
   * @param geometry CellGeometry object.
   * @param data Geometry data.
   */
  void _testNumCorners(const CellGeometry& geometry,
		       const CellGeomData& data);

  /** Test jacobian().
   *
   * @param geometry CellGeometry object.
   * @param data Geometry data.
   */
  void _testJacobian(const CellGeometry* geometry,
		     const CellGeomData& data);

}; // class TestCellGeometry

#endif // pylith_topology_testcellgeometry_hh

// End of file 
