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

#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/utils/array.hh" // USES double_array

#include "data/CellGeomData.hh" // USES CellGeomData

// ----------------------------------------------------------------------
// Setup data.
void 
pylith::feassemble::TestCellGeometry::setUp(void)
{ // setUp
  _object = 0;
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down test data.
void
pylith::feassemble::TestCellGeometry::tearDown(void)
{ // tearDown
  delete _object; _object = 0;
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test clone()
void
pylith::feassemble::TestCellGeometry::testClone(void)
{ // testClone
  CellGeometry* copy = _object->clone();
  CPPUNIT_ASSERT_EQUAL(_object->cellDim(), copy->cellDim());
  CPPUNIT_ASSERT_EQUAL(_object->spaceDim(), copy->spaceDim());
  CPPUNIT_ASSERT_EQUAL(_object->numCorners(), copy->numCorners());
  delete copy; copy = 0;
} // testClone

// ----------------------------------------------------------------------
// Test cellDim().
void
pylith::feassemble::TestCellGeometry::testCellDim(void)
{ // testCellDim
  CPPUNIT_ASSERT(0 != _object);
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(_data->cellDim, _object->cellDim());
} // testCellDim

// ----------------------------------------------------------------------
// Test spaceDim().
void
pylith::feassemble::TestCellGeometry::testSpaceDim(void)
{ // _testSpaceDim
  CPPUNIT_ASSERT(0 != _object);
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(_data->spaceDim, _object->spaceDim());
} // _testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners().
void
pylith::feassemble::TestCellGeometry::testNumCorners(void)
{ // testNumCorners
  CPPUNIT_ASSERT(0 != _object);
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(_data->numCorners, _object->numCorners());
} // testNumConers

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::feassemble::TestCellGeometry::testJacobian(void)
{ // testJacobian
  CPPUNIT_ASSERT(0 != _object);
  CPPUNIT_ASSERT(0 != _data);

  const int cellDim = _data->cellDim;
  const int spaceDim = _data->spaceDim;
  const int numCorners = _data->numCorners;

  const int numLocs = _data->numLocs;

  CPPUNIT_ASSERT_EQUAL(cellDim, _object->cellDim());
  CPPUNIT_ASSERT_EQUAL(spaceDim, _object->spaceDim());
  CPPUNIT_ASSERT_EQUAL(numCorners, _object->numCorners());

  double_array jacobian(cellDim*spaceDim);
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    double_array vertices(_data->vertices, numCorners*spaceDim);
    double_array location(&_data->locations[iLoc*cellDim], cellDim);

    _object->jacobian(&jacobian, vertices, location);

    const int size = jacobian.size();
    const int index = iLoc*cellDim*spaceDim;
    const double tolerance = 1.0e-06;
    for (int i=0; i < size; ++i)
      if (_data->jacobian[index+i] < 1.0)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->jacobian[index+i], jacobian[i],
				     tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, jacobian[i]/_data->jacobian[index+i],
				     tolerance);
  } // for
} // testJacobian


// End of file 
