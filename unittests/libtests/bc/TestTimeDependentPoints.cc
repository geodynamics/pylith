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

#include "TestTimeDependentPoints.hh" // Implementation of class methods

#include "pylith/bc/PointForce.hh" // USES PointForce

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitormesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestTimeDependentPoints );

// ----------------------------------------------------------------------
namespace pylith {
  namespace bc {
    namespace _TestTimeDependentPoints {
      const PylithScalar pressureScale = 4.0;
      const PylithScalar lengthScale = 1.5;
      const PylithScalar timeScale = 0.5;
      const int npointsIn = 2;
      const int pointsIn[npointsIn] = { 3, 5, };
      const int npointsOut = 2;
      const int pointsOut[npointsOut] = { 2, 4, };

      const int numBCDOF = 2;
      const int bcDOF[numBCDOF] = { 1, 0 };
      const PylithScalar initial[npointsIn*numBCDOF] = {
	0.3,  0.4,
	0.7,  0.6,
      };
      const PylithScalar rate[npointsIn*numBCDOF] = {
	-0.2,  -0.1,
	 0.4,   0.3,
      };
      const PylithScalar rateTime[npointsIn] = {
	0.5,
	0.8,
      };
      const PylithScalar change[npointsIn*numBCDOF] = {
	1.3,  1.4,
	1.7,  1.6,
      };
      const PylithScalar changeTime[npointsIn] = {
	2.0,
	2.4,
      };

      const PylithScalar tValue = 2.2;
      const PylithScalar tValue2 = 2.6;
      const PylithScalar valuesRate[npointsIn*numBCDOF] = {
	-0.34,  -0.17,
	 0.56,   0.42,
      };
      const PylithScalar valuesChange[npointsIn*numBCDOF] = {
	1.3,  1.4,
	0.0,  0.0,
      };
      const PylithScalar valuesChangeTH[npointsIn*numBCDOF] = {
	1.3*0.98,  1.4*0.98,
	0.0,  0.0,
      };
      const PylithScalar valuesIncrInitial[npointsIn*numBCDOF] = {
	0.0,  0.0,
	0.0,  0.0,
      };
      const PylithScalar valuesIncrRate[npointsIn*numBCDOF] = {
	-0.08,  -0.04,
	 0.16,   0.12,
      };
      const PylithScalar valuesIncrChange[npointsIn*numBCDOF] = {
	0.0,  0.0,
	1.7,  1.6,
      };
      const PylithScalar valuesIncrChangeTH[npointsIn*numBCDOF] = {
	1.3*-0.04,  1.4*-0.04,
	1.7*0.98,  1.6*0.98,
      };

      // Check values in section against expected values.
      static
      void _checkValues(const PylithScalar* valuesE,
			const int fiberDimE,
			const topology::Field& field,
			const PylithScalar scale);
    } // _TestTimeDependentPoints
  } // bc
} // pylith

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestTimeDependentPoints::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  const char* filename = "data/tri3.mesh";

  _mesh = new topology::Mesh();CPPUNIT_ASSERT(_mesh);
  meshio::MeshIOAscii iohandler;
  iohandler.filename(filename);
  iohandler.read(_mesh);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(_mesh->dimension());
  cs.initialize();
  _mesh->coordsys(&cs);

  spatialdata::units::Nondimensional normalizer;
  normalizer.pressureScale(_TestTimeDependentPoints::pressureScale);
  normalizer.lengthScale(_TestTimeDependentPoints::lengthScale);
  normalizer.timeScale(_TestTimeDependentPoints::timeScale);
  topology::MeshOps::nondimensionalize(_mesh, normalizer);

  _bc = new PointForce();CPPUNIT_ASSERT(_bc);
  _bc->label("bc");
  _bc->normalizer(normalizer);
  _bc->bcDOF(_TestTimeDependentPoints::bcDOF, _TestTimeDependentPoints::numBCDOF);
  _bc->_getPoints(*_mesh);

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestTimeDependentPoints::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _mesh; _mesh = 0;
  delete _bc; _bc = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestTimeDependentPoints::testBCDOF(void)
{ // testBCDOF
  PYLITH_METHOD_BEGIN;

  PointForce bc;

  const size_t numDOF = 4;
  const int fixedDOF[numDOF] = { 0, 2, 3, 5 };
  bc.bcDOF(fixedDOF, numDOF);

  CPPUNIT_ASSERT_EQUAL(numDOF, bc._bcDOF.size());
  for (int i=0; i < numDOF; ++i)
    CPPUNIT_ASSERT_EQUAL(fixedDOF[i], bc._bcDOF[i]);

  PYLITH_METHOD_END;
} // testBCDOF

// ----------------------------------------------------------------------
// Test _getLabel().
void
pylith::bc::TestTimeDependentPoints::testGetLabel(void)
{ // testGetLabel
  PYLITH_METHOD_BEGIN;

  PointForce bc;
  
  const std::string& label = "point force";
  bc.label(label.c_str());
  CPPUNIT_ASSERT_EQUAL(label, std::string(bc._getLabel()));

  PYLITH_METHOD_END;
} // testGetLabel

// ----------------------------------------------------------------------
// Test _queryDatabases().
void
pylith::bc::TestTimeDependentPoints::testQueryDatabases(void)
{ // testQueryDatabases
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_bc);

  spatialdata::spatialdb::SimpleDB dbInitial("TestTimeDependentPoints _queryDatabases initial");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/tri3_force.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbRate("TestTimeDependentPoints _queryDatabases rate");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/tri3_force_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbChange("TestTimeDependentPoints _queryDatabases change");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/tri3_force_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::TimeHistory th("TestTimeDependentPoints _queryDatabases time history");
  th.filename("data/tri3_force.timedb");

  _bc->dbInitial(&dbInitial);
  _bc->dbRate(&dbRate);
  _bc->dbChange(&dbChange);
  _bc->dbTimeHistory(&th);

  const PylithScalar pressureScale = _TestTimeDependentPoints::pressureScale;
  const PylithScalar lengthScale = _TestTimeDependentPoints::lengthScale;
  const PylithScalar timeScale = _TestTimeDependentPoints::timeScale;
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;
  const char* fieldName = "force";
  _bc->_queryDatabases(*_mesh, forceScale, fieldName);

  const PylithScalar tolerance = 1.0e-06;
  const int numBCDOF = _TestTimeDependentPoints::numBCDOF;
  CPPUNIT_ASSERT(_bc->_parameters);
  
  // Check initial values.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::initial, numBCDOF, _bc->_parameters->get("initial"), forceScale);

  // Check rate values.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::rate, numBCDOF, _bc->_parameters->get("rate"), forceScale/timeScale);

  // Check rate start time.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::rateTime, 1, _bc->_parameters->get("rate time"), timeScale);

  // Check change values.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::change, numBCDOF, _bc->_parameters->get("change"), forceScale);

  // Check change start time.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::changeTime, 1, _bc->_parameters->get("change time"), timeScale);
  th.close();

  PYLITH_METHOD_END;
} // testQueryDatabases

// ----------------------------------------------------------------------
// Test _calculateValue() with initial value.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueInitial(void)
{ // testCalculateValueInitial
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_bc);

  spatialdata::spatialdb::SimpleDB dbInitial("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/tri3_force.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  _bc->dbInitial(&dbInitial);

  const PylithScalar pressureScale = _TestTimeDependentPoints::pressureScale;
  const PylithScalar lengthScale = _TestTimeDependentPoints::lengthScale;
  const PylithScalar timeScale = _TestTimeDependentPoints::timeScale;
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;
  const char* fieldName = "force";
  _bc->_queryDatabases(*_mesh, forceScale, fieldName);
  _bc->_calculateValue(_TestTimeDependentPoints::tValue/timeScale);

  const PylithScalar tolerance = 1.0e-06;
  const int numBCDOF = _TestTimeDependentPoints::numBCDOF;
  CPPUNIT_ASSERT(_bc->_parameters);
  
  // Check values.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::initial, numBCDOF, _bc->_parameters->get("value"), forceScale);

  PYLITH_METHOD_END;
} // testCalculateValueInitial

// ----------------------------------------------------------------------
// Test _calculateValue() with rate.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueRate(void)
{ // testCalculateValueRate
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_bc);

  spatialdata::spatialdb::SimpleDB dbRate("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/tri3_force_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  _bc->dbRate(&dbRate);

  const PylithScalar pressureScale = _TestTimeDependentPoints::pressureScale;
  const PylithScalar lengthScale = _TestTimeDependentPoints::lengthScale;
  const PylithScalar timeScale = _TestTimeDependentPoints::timeScale;
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;
  const char* fieldName = "force";
  _bc->_queryDatabases(*_mesh, forceScale, fieldName);
  _bc->_calculateValue(_TestTimeDependentPoints::tValue/timeScale);

  const PylithScalar tolerance = 1.0e-06;
  const int numBCDOF = _TestTimeDependentPoints::numBCDOF;
  CPPUNIT_ASSERT(_bc->_parameters);
  
  // Check values.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesRate, numBCDOF, _bc->_parameters->get("value"), forceScale);

  PYLITH_METHOD_END;
} // testCalculateValueRate

// ----------------------------------------------------------------------
// Test _calculateValue() with temporal change.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueChange(void)
{ // testCalculateValueChange
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_bc);

  spatialdata::spatialdb::SimpleDB dbChange("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/tri3_force_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  _bc->dbChange(&dbChange);

  const PylithScalar pressureScale = _TestTimeDependentPoints::pressureScale;
  const PylithScalar lengthScale = _TestTimeDependentPoints::lengthScale;
  const PylithScalar timeScale = _TestTimeDependentPoints::timeScale;
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;
  const char* fieldName = "force";
  _bc->_queryDatabases(*_mesh, forceScale, fieldName);
  _bc->_calculateValue(_TestTimeDependentPoints::tValue/timeScale);

  const PylithScalar tolerance = 1.0e-06;
  const int numBCDOF = _TestTimeDependentPoints::numBCDOF;
  CPPUNIT_ASSERT(_bc->_parameters);
  
  // Check values.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesChange, numBCDOF, _bc->_parameters->get("value"), forceScale);

  PYLITH_METHOD_END;
} // testCalculateValueChange

// ----------------------------------------------------------------------
// Test _calculateValue() with temporal change w/time history.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueChangeTH(void)
{ // testCalculateValueChangeTH
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_bc);

  spatialdata::spatialdb::SimpleDB dbChange("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/tri3_force_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::TimeHistory th("TestTimeDependentPoints _queryDatabases");
  th.filename("data/tri3_force.timedb");

  _bc->dbChange(&dbChange);
  _bc->dbTimeHistory(&th);

  const PylithScalar pressureScale = _TestTimeDependentPoints::pressureScale;
  const PylithScalar lengthScale = _TestTimeDependentPoints::lengthScale;
  const PylithScalar timeScale = _TestTimeDependentPoints::timeScale;
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;
  const char* fieldName = "force";
  _bc->_queryDatabases(*_mesh, forceScale, fieldName);
  _bc->_calculateValue(_TestTimeDependentPoints::tValue/timeScale);

  const PylithScalar tolerance = 1.0e-06;
  const int numBCDOF = _TestTimeDependentPoints::numBCDOF;
  CPPUNIT_ASSERT(_bc->_parameters);
  
  // Check values.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesChangeTH, numBCDOF, _bc->_parameters->get("value"), forceScale);

  PYLITH_METHOD_END;
} // testCalculateValueChangeTH

// ----------------------------------------------------------------------
// Test _calculateValue() with initial, rate, and temporal change w/time history.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueAll(void)
{ // testCalculateValueAll
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_bc);

  spatialdata::spatialdb::SimpleDB dbInitial("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/tri3_force.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbRate("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/tri3_force_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbChange("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/tri3_force_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::TimeHistory th("TestTimeDependentPoints _queryDatabases");
  th.filename("data/tri3_force.timedb");

  _bc->dbInitial(&dbInitial);
  _bc->dbRate(&dbRate);
  _bc->dbChange(&dbChange);
  _bc->dbTimeHistory(&th);

  const PylithScalar pressureScale = _TestTimeDependentPoints::pressureScale;
  const PylithScalar lengthScale = _TestTimeDependentPoints::lengthScale;
  const PylithScalar timeScale = _TestTimeDependentPoints::timeScale;
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;
  const char* fieldName = "force";
  _bc->_queryDatabases(*_mesh, forceScale, fieldName);
  _bc->_calculateValue(_TestTimeDependentPoints::tValue/timeScale);

  const PylithScalar tolerance = 1.0e-06;
  const int numBCDOF = _TestTimeDependentPoints::numBCDOF;
  CPPUNIT_ASSERT(_bc->_parameters);
  
  // Check values.
  const int npoints = _TestTimeDependentPoints::npointsIn;
  scalar_array valuesE(npoints*numBCDOF);
  for (int i=0; i < valuesE.size(); ++i)
    valuesE[i] = 
      _TestTimeDependentPoints::initial[i] +
      _TestTimeDependentPoints::valuesRate[i] +
      _TestTimeDependentPoints::valuesChangeTH[i];

  _TestTimeDependentPoints::_checkValues(&valuesE[0], numBCDOF, _bc->_parameters->get("value"), forceScale);

  PYLITH_METHOD_END;
} // testCalculateValueAll

// ----------------------------------------------------------------------
// Test _calculateValueIncr() with initial value.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueIncrInitial(void)
{ // testCalculateValueIncrInitial
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_bc);

  spatialdata::spatialdb::SimpleDB dbInitial("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/tri3_force.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  _bc->dbInitial(&dbInitial);

  const PylithScalar pressureScale = _TestTimeDependentPoints::pressureScale;
  const PylithScalar lengthScale = _TestTimeDependentPoints::lengthScale;
  const PylithScalar timeScale = _TestTimeDependentPoints::timeScale;
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;
  const char* fieldName = "force";
  _bc->_queryDatabases(*_mesh, forceScale, fieldName);
  const PylithScalar t0 = _TestTimeDependentPoints::tValue / timeScale;
  const PylithScalar t1 = _TestTimeDependentPoints::tValue2 / timeScale;
  _bc->_calculateValueIncr(t0, t1);

  const PylithScalar tolerance = 1.0e-06;
  const int numBCDOF = _TestTimeDependentPoints::numBCDOF;
  CPPUNIT_ASSERT(_bc->_parameters);
  
  // Check values.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesIncrInitial, numBCDOF, _bc->_parameters->get("value"), forceScale);

  PYLITH_METHOD_END;
} // testCalculateValueIncrInitial

// ----------------------------------------------------------------------
// Test _calculateValueIncr() with rate.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueIncrRate(void)
{ // testCalculateValueIncrRate
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_bc);

  spatialdata::spatialdb::SimpleDB dbRate("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/tri3_force_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  _bc->dbRate(&dbRate);

  const PylithScalar pressureScale = _TestTimeDependentPoints::pressureScale;
  const PylithScalar lengthScale = _TestTimeDependentPoints::lengthScale;
  const PylithScalar timeScale = _TestTimeDependentPoints::timeScale;
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;
  const char* fieldName = "force";
  _bc->_queryDatabases(*_mesh, forceScale, fieldName);
  const PylithScalar t0 = _TestTimeDependentPoints::tValue / timeScale;
  const PylithScalar t1 = _TestTimeDependentPoints::tValue2 / timeScale;
  _bc->_calculateValueIncr(t0, t1);

  const PylithScalar tolerance = 1.0e-06;
  const int numBCDOF = _TestTimeDependentPoints::numBCDOF;
  CPPUNIT_ASSERT(_bc->_parameters);
  
  // Check values.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesIncrRate, numBCDOF, _bc->_parameters->get("value"), forceScale);

  PYLITH_METHOD_END;
} // testCalculateValueIncrRate

// ----------------------------------------------------------------------
// Test _calculateValueIncr() with temporal change.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueIncrChange(void)
{ // testCalculateValueIncrChange
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_bc);

  spatialdata::spatialdb::SimpleDB dbChange("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/tri3_force_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  _bc->dbChange(&dbChange);

  const PylithScalar pressureScale = _TestTimeDependentPoints::pressureScale;
  const PylithScalar lengthScale = _TestTimeDependentPoints::lengthScale;
  const PylithScalar timeScale = _TestTimeDependentPoints::timeScale;
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;
  const char* fieldName = "force";
  _bc->_queryDatabases(*_mesh, forceScale, fieldName);
  const PylithScalar t0 = _TestTimeDependentPoints::tValue / timeScale;
  const PylithScalar t1 = _TestTimeDependentPoints::tValue2 / timeScale;
  _bc->_calculateValueIncr(t0, t1);

  const PylithScalar tolerance = 1.0e-06;
  const int numBCDOF = _TestTimeDependentPoints::numBCDOF;
  CPPUNIT_ASSERT(_bc->_parameters);
  
  // Check values.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesIncrChange, numBCDOF, _bc->_parameters->get("value"), forceScale);

  PYLITH_METHOD_END;
} // testCalculateValueIncrChange

// ----------------------------------------------------------------------
// Test _calculateValueIncr() with temporal change w/time history.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueIncrChangeTH(void)
{ // testCalculateValueIncrChangeTH
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_bc);

  spatialdata::spatialdb::SimpleDB dbChange("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/tri3_force_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::TimeHistory th("TestTimeDependentPoints _queryDatabases");
  th.filename("data/tri3_force.timedb");

  _bc->dbChange(&dbChange);
  _bc->dbTimeHistory(&th);

  const PylithScalar pressureScale = _TestTimeDependentPoints::pressureScale;
  const PylithScalar lengthScale = _TestTimeDependentPoints::lengthScale;
  const PylithScalar timeScale = _TestTimeDependentPoints::timeScale;
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;
  const char* fieldName = "force";
  _bc->_queryDatabases(*_mesh, forceScale, fieldName);
  const PylithScalar t0 = _TestTimeDependentPoints::tValue / timeScale;
  const PylithScalar t1 = _TestTimeDependentPoints::tValue2 / timeScale;
  _bc->_calculateValueIncr(t0, t1);

  const PylithScalar tolerance = 1.0e-06;
  const int numBCDOF = _TestTimeDependentPoints::numBCDOF;
  CPPUNIT_ASSERT(_bc->_parameters);
  
  // Check values.
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesIncrChangeTH, numBCDOF, _bc->_parameters->get("value"), forceScale);

  PYLITH_METHOD_END;
} // testCalculateValueIncrChangeTH

// ----------------------------------------------------------------------
// Test _calculateValueIncr() with initial, rate, and temporal change w/time history.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueIncrAll(void)
{ // testCalculateValueIncrAll
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_bc);

  spatialdata::spatialdb::SimpleDB dbInitial("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbInitialIO;
  dbInitialIO.filename("data/tri3_force.spatialdb");
  dbInitial.ioHandler(&dbInitialIO);
  dbInitial.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbRate("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbRateIO;
  dbRateIO.filename("data/tri3_force_rate.spatialdb");
  dbRate.ioHandler(&dbRateIO);
  dbRate.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::SimpleDB dbChange("TestTimeDependentPoints _queryDatabases");
  spatialdata::spatialdb::SimpleIOAscii dbChangeIO;
  dbChangeIO.filename("data/tri3_force_change.spatialdb");
  dbChange.ioHandler(&dbChangeIO);
  dbChange.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

  spatialdata::spatialdb::TimeHistory th("TestTimeDependentPoints _queryDatabases");
  th.filename("data/tri3_force.timedb");

  _bc->dbInitial(&dbInitial);
  _bc->dbRate(&dbRate);
  _bc->dbChange(&dbChange);
  _bc->dbTimeHistory(&th);

  const PylithScalar pressureScale = _TestTimeDependentPoints::pressureScale;
  const PylithScalar lengthScale = _TestTimeDependentPoints::lengthScale;
  const PylithScalar timeScale = _TestTimeDependentPoints::timeScale;
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;
  const char* fieldName = "force";
  _bc->_queryDatabases(*_mesh, forceScale, fieldName);
  const PylithScalar t0 = _TestTimeDependentPoints::tValue / timeScale;
  const PylithScalar t1 = _TestTimeDependentPoints::tValue2 / timeScale;
  _bc->_calculateValueIncr(t0, t1);

  const PylithScalar tolerance = 1.0e-06;
  const int numBCDOF = _TestTimeDependentPoints::numBCDOF;
  CPPUNIT_ASSERT(_bc->_parameters);
  
  // Check values.
  const int npoints = _TestTimeDependentPoints::npointsIn;
  scalar_array valuesE(npoints*numBCDOF);
  for (int i=0; i < valuesE.size(); ++i)
    valuesE[i] = 
      _TestTimeDependentPoints::valuesIncrInitial[i] +
      _TestTimeDependentPoints::valuesIncrRate[i] +
      _TestTimeDependentPoints::valuesIncrChangeTH[i];

  _TestTimeDependentPoints::_checkValues(&valuesE[0], numBCDOF, _bc->_parameters->get("value"), forceScale);

  PYLITH_METHOD_END;
} // testCalculateValueIncrAll

// ----------------------------------------------------------------------
// Check values in section against expected values.
void
pylith::bc::_TestTimeDependentPoints::_checkValues(const PylithScalar* valuesE,
						   const int fiberDimE,
						   const topology::Field& field,
						   const PylithScalar scale)
{ // _checkValues
  PYLITH_METHOD_BEGIN;

  topology::VecVisitorMesh fieldVisitor(field);
  const PetscScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);
  
  const PylithScalar tolerance = 1.0e-06;

  // Check values at points associated with BC.
  const int npointsIn = _TestTimeDependentPoints::npointsIn;
  for (int i=0; i < npointsIn; ++i) {
    const int p_bc = _TestTimeDependentPoints::pointsIn[i];
    
    const PetscInt off = fieldVisitor.sectionOffset(p_bc);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fieldVisitor.sectionDof(p_bc));
    
    for (int iDim=0; iDim < fiberDimE; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i*fiberDimE+iDim]/scale, fieldArray[off+iDim], tolerance);
  } // for
  
  PYLITH_METHOD_END;
} // _checkValues


// End of file 
