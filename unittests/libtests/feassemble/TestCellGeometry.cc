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
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestCellGeometry );

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
// Test _orient1D().
void
pylith::feassemble::TestCellGeometry::testOrient1D(void)
{ // testOrient1D
  double_array jacobian;
  double jacobianDet;
  double_array upDir;
  double_array orientation(1);
  
  CellGeometry::_orient1D(&orientation, jacobian, jacobianDet, upDir);

  const int size = orientation.size();
  CPPUNIT_ASSERT_EQUAL(1, size);
  const double tolerance = 1.0e-6;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, orientation[i], tolerance);
} // testOrient1D

// ----------------------------------------------------------------------
// Test _orient2D().
void
pylith::feassemble::TestCellGeometry::testOrient2D(void)
{ // testOrient2D
  const int numLocs = 2;
  const int spaceDim = 2;
  const int orientSize = 4;

  const double jacobianVals[] = {
    -1.0, 2.0,
    -0.5, 1.0
  };
  const double orientationE[] = {
    -1.0,  2.0,  2.0, 1.0,
    -0.5,  1.0,  1.0, 0.5
  };

  const int jacobianSize = spaceDim*(spaceDim-1);
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    double_array jacobian(&jacobianVals[iLoc*jacobianSize], jacobianSize);
    double jacobianDet;
    double_array upDir;
    double_array orientation(orientSize);

    CellGeometry::_orient2D(&orientation, jacobian, jacobianDet, upDir);

    const int size = orientation.size();
    CPPUNIT_ASSERT_EQUAL(orientSize, size);
    const double tolerance = 1.0e-6;
    for (int i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(orientationE[iLoc*orientSize+i],
				   orientation[i], tolerance);
  } // for
} // testOrient2D

// ----------------------------------------------------------------------
// Test _orient3D().
void
pylith::feassemble::TestCellGeometry::testOrient3D(void)
{ // testOrient3D
  const int numLocs = 2;
  const int spaceDim = 3;
  const int orientSize = 9;

  const double jacobianVals[] = {
    2.0,  -0.5,
    1.0,  -0.2,
    0.5,   2.0,

    -1.0,  2.0,
    -3.0, -0.2,
    -0.3,  0.3,
  };
  const double jacobianDetVals[] = {
    1.3, 0.7
  };
  const double upDirVals[] = { 0.0, 0.0, 1.0 };
  const double orientationE[] = {
    1.1654847299258313, 0.57588657243394026, 0.0, 
    -0.012145479112634533, 0.024580136299379406, 1.2997108540889502, 
    0.57575848378190342, -1.1652255028919474, 0.027417070656281111,

    0.20879249519516274, -0.66813598462452073, 0.0,
    0.65951433797689429, 0.20609823061777949, 0.11209084413599232, 
    -0.10698846644884991, -0.033433895765265592, 0.69096717914882233
  };

  double_array upDir(upDirVals, 3);
  const int jacobianSize = spaceDim*(spaceDim-1);
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    double_array jacobian(&jacobianVals[iLoc*jacobianSize], jacobianSize);
    double jacobianDet = jacobianDetVals[iLoc];
    double_array orientation(orientSize);

    CellGeometry::_orient3D(&orientation, jacobian, jacobianDet, upDir);

    const int size = orientation.size();
    CPPUNIT_ASSERT_EQUAL(orientSize, size);
    const double tolerance = 1.0e-6;
    for (int i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(orientationE[iLoc*orientSize+i],
				   orientation[i], tolerance);
  } // for
} // testOrient3D

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
// Test orientFn.
void
pylith::feassemble::TestCellGeometry::testOrientFn(void)
{ // testOrientFn
  CPPUNIT_ASSERT(0 != _object);
  CPPUNIT_ASSERT(0 != _data);
  switch (_data->cellDim)
    { // switch
    case 0 :
      CPPUNIT_ASSERT(CellGeometry::_orient1D == _object->_orientFn);
      break;
    case 1 :
      CPPUNIT_ASSERT(CellGeometry::_orient2D == _object->_orientFn);
      break;
    case 2 :
      CPPUNIT_ASSERT(CellGeometry::_orient3D == _object->_orientFn);
      break;
    case 3 :
      CPPUNIT_ASSERT(0 == _object->_orientFn);
      break;
    default :
      assert(0);
    } // switch
} // testOrientFn

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
  double det = 0;
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    double_array vertices(_data->vertices, numCorners*spaceDim);
    double_array location(&_data->locations[iLoc*cellDim], cellDim);

    _object->jacobian(&jacobian, &det, vertices, location);

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
    CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->jacobianDet[iLoc], det, tolerance);
  } // for
} // testJacobian


// End of file 
