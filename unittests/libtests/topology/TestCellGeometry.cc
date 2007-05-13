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

#include <portinfo>

#include "TestCellGeometry.hh" // Implementation of class methods

#include "pylith/topology/CellGeometry.hh" // USES CellGeometry
#include "pylith/utils/array.hh" // USES double_array

#include "data/CellGeomData.hh" // USES CellGeomData

// ----------------------------------------------------------------------
// Test cellDim().
void
pylith::topology::TestCellGeometry::_testCellDim(const CellGeometry& geometry,
						 const CellGeomData& data)
{ // _testCellDim
  CPPUNIT_ASSERT_EQUAL(data.cellDim, geometry.cellDim());
} // _testCellDim

// ----------------------------------------------------------------------
// Test spaceDim().
void
pylith::topology::TestCellGeometry::_testSpaceDim(const CellGeometry& geometry,
						  const CellGeomData& data)
{ // _testSpaceDim
  CPPUNIT_ASSERT_EQUAL(data.spaceDim, geometry.spaceDim());
} // _testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners().
void
pylith::topology::TestCellGeometry::_testNumCorners(const CellGeometry& geometry,
						    const CellGeomData& data)
{ // _testNumCorners
  CPPUNIT_ASSERT_EQUAL(data.numCorners, geometry.numCorners());
} // _testNumConers

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::topology::TestCellGeometry::_testJacobian(const CellGeometry* geometry,
						  const CellGeomData& data)
{ // _testJacobian
  const int cellDim = data.cellDim;
  const int spaceDim = data.spaceDim;
  const int numCorners = data.numCorners;

  const int numLocs = data.numLocs;

  CPPUNIT_ASSERT_EQUAL(cellDim, geometry->cellDim());
  CPPUNIT_ASSERT_EQUAL(spaceDim, geometry->spaceDim());
  CPPUNIT_ASSERT_EQUAL(numCorners, geometry->numCorners());

  double_array jacobian(cellDim*spaceDim);
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    double_array vertices(data.vertices, numCorners*spaceDim);
    double_array location(&data.locations[iLoc*cellDim], cellDim);

    geometry->jacobian(&jacobian, vertices, location);

    const int size = jacobian.size();
    const int index = iLoc*cellDim*spaceDim;
    const double tolerance = 1.0e-06;
    for (int i=0; i < size; ++i)
      if (data.jacobian[index+i] < 1.0)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.jacobian[index+i], jacobian[i],
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, jacobian[i]/data.jacobian[index+i],
				     tolerance);
  } // for
} // _testJacobian


// End of file 
