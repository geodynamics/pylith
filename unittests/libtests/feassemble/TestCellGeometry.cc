// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestCellGeometry.hh" // Implementation of class methods

#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include "data/CellGeomData.hh" // USES CellGeomData

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestCellGeometry );

// ----------------------------------------------------------------------
// Setup data.
void 
pylith::feassemble::TestCellGeometry::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _object = 0;
  _data = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down test data.
void
pylith::feassemble::TestCellGeometry::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _object; _object = 0;
  delete _data; _data = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test _orient1D().
void
pylith::feassemble::TestCellGeometry::testOrient1D(void)
{ // testOrient1D
  PYLITH_METHOD_BEGIN;

  const int numLocs = 2;
  const int spaceDim = 2;
  const int orientSize = 4;

  const PylithScalar jacobianVals[numLocs*spaceDim*(spaceDim-1)] = {
    -1.0, 2.0,
    -0.5, 1.0
  };
  const PylithScalar orientationE[numLocs*orientSize] = {
    -1.0,  2.0,  2.0, 1.0,
    -0.5,  1.0,  1.0, 0.5
  };

  const int jacobianSize = spaceDim*(spaceDim-1);
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    scalar_array jacobian(&jacobianVals[iLoc*jacobianSize], jacobianSize);
    PylithScalar jacobianDet;
    scalar_array upDir;
    scalar_array orientation(orientSize);

    CellGeometry::_orient1D(&orientation, jacobian, jacobianDet, upDir);

    const int size = orientation.size();
    CPPUNIT_ASSERT_EQUAL(orientSize, size);
    const PylithScalar tolerance = 1.0e-6;
    for (int i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(orientationE[iLoc*orientSize+i], orientation[i], tolerance);
  } // for

  PYLITH_METHOD_END;
} // testOrient1D

// ----------------------------------------------------------------------
// Test _orient2D().
void
pylith::feassemble::TestCellGeometry::testOrient2D(void)
{ // testOrient2D
  PYLITH_METHOD_BEGIN;

  const int numLocs = 2;
  const int spaceDim = 3;
  const int orientSize = 9;

  const PylithScalar jacobianVals[numLocs*spaceDim*(spaceDim-1)] = {
    2.0,  -0.5,
    1.0,  -0.2,
    0.5,   2.0,

    -1.0,  2.0,
    -3.0, -0.2,
    -0.3,  0.3,
  };
  const PylithScalar jacobianDetVals[numLocs] = {
    1.3, 0.7
  };
  const PylithScalar upDirVals[spaceDim] = { 0.0, 0.0, 1.0 };
  const PylithScalar orientationE[numLocs*orientSize] = {
    1.1654847299258313, 0.57588657243394026, 0.0, 
    -0.012145479112634533, 0.024580136299379406, 1.2997108540889502, 
    0.57575848378190342, -1.1652255028919474, 0.027417070656281111,

    0.20879249519516274, -0.66813598462452073, 0.0,
    0.65951433797689429, 0.20609823061777949, 0.11209084413599232, 
    -0.10698846644884991, -0.033433895765265592, 0.69096717914882233
  };

  scalar_array upDir(upDirVals, 3);
  const int jacobianSize = spaceDim*(spaceDim-1);
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    scalar_array jacobian(&jacobianVals[iLoc*jacobianSize], jacobianSize);
    PylithScalar jacobianDet = jacobianDetVals[iLoc];
    scalar_array orientation(orientSize);

    CellGeometry::_orient2D(&orientation, jacobian, jacobianDet, upDir);

    const int size = orientation.size();
    CPPUNIT_ASSERT_EQUAL(orientSize, size);
    const PylithScalar tolerance = 1.0e-6;
    for (int i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(orientationE[iLoc*orientSize+i], orientation[i], tolerance);
  } // for

  PYLITH_METHOD_END;
} // testOrient2D

// ----------------------------------------------------------------------
// Test clone()
void
pylith::feassemble::TestCellGeometry::testClone(void)
{ // testClone
  PYLITH_METHOD_BEGIN;

  CellGeometry* copy = _object->clone();
  CPPUNIT_ASSERT_EQUAL(_object->cellDim(), copy->cellDim());
  CPPUNIT_ASSERT_EQUAL(_object->spaceDim(), copy->spaceDim());
  CPPUNIT_ASSERT_EQUAL(_object->numCorners(), copy->numCorners());
  delete copy; copy = 0;

  PYLITH_METHOD_END;
} // testClone

// ----------------------------------------------------------------------
// Test cellDim().
void
pylith::feassemble::TestCellGeometry::testCellDim(void)
{ // testCellDim
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_object);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(_data->cellDim, _object->cellDim());

  PYLITH_METHOD_END;
} // testCellDim

// ----------------------------------------------------------------------
// Test spaceDim().
void
pylith::feassemble::TestCellGeometry::testSpaceDim(void)
{ // _testSpaceDim
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_object);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(_data->spaceDim, _object->spaceDim());

  PYLITH_METHOD_END;
} // _testSpaceDim

// ----------------------------------------------------------------------
// Test numCorners().
void
pylith::feassemble::TestCellGeometry::testNumCorners(void)
{ // testNumCorners
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_object);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(_data->numCorners, _object->numCorners());

  PYLITH_METHOD_END;
} // testNumConers

// ----------------------------------------------------------------------
// Test orientFn.
void
pylith::feassemble::TestCellGeometry::testOrientFn(void)
{ // testOrientFn
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_object);
  CPPUNIT_ASSERT(_data);

  // We can't compare function pointers, because they depend on how
  // the dynamic library is loaded (seems okay on Linux and Darwin but
  // fails with cygwin). Instead we settle for making sure the
  // function pointer is not 0 (NULL).

  switch (_data->cellDim)
    { // switch
    case 0 :
    case 1 :
    case 2 :
      CPPUNIT_ASSERT(_object->_orientFn);
      break;
    case 3 :
      CPPUNIT_ASSERT(0 == _object->_orientFn);
      break;
    default :
      CPPUNIT_ASSERT(0);
    } // switch

  PYLITH_METHOD_END;
} // testOrientFn

// ----------------------------------------------------------------------
// Test jacobian().
void
pylith::feassemble::TestCellGeometry::testJacobian(void)
{ // testJacobian
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_object);
  CPPUNIT_ASSERT(_data);

  const int cellDim = _data->cellDim;
  const int spaceDim = _data->spaceDim;
  const int numCorners = _data->numCorners;

  const int numLocs = _data->numLocs;

  CPPUNIT_ASSERT_EQUAL(cellDim, _object->cellDim());
  CPPUNIT_ASSERT_EQUAL(spaceDim, _object->spaceDim());
  CPPUNIT_ASSERT_EQUAL(numCorners, _object->numCorners());

  scalar_array jacobian(cellDim*spaceDim);
  PylithScalar det = 0;
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    _object->jacobian(&jacobian, &det, _data->vertices, numCorners, spaceDim, &_data->locations[iLoc*cellDim], cellDim);

    const int size = jacobian.size();
    const int index = iLoc*cellDim*spaceDim;
    const PylithScalar tolerance = 1.0e-06;

    for (int i=0; i < size; ++i) {
      if (_data->jacobian[index+i] < 1.0)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->jacobian[index+i], jacobian[i], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, jacobian[i]/_data->jacobian[index+i], tolerance);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->jacobianDet[iLoc], det, tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testJacobian


// End of file 
