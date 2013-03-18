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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestTimeDependentPoints.hh" // Implementation of class methods

#include "pylith/bc/PointForce.hh" // USES PointForce

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory

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
			PetscSection section,
			Vec vec,
			const PylithScalar scale);
    } // _TestTimeDependentPoints
  } // bc
} // pylith

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestTimeDependentPoints::setUp(void)
{ // setUp
  const char* filename = "data/tri3.mesh";

  _mesh = new topology::Mesh();
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
  _mesh->nondimensionalize(normalizer);

  _bc = new PointForce();
  _bc->label("bc");
  _bc->normalizer(normalizer);
  _bc->bcDOF(_TestTimeDependentPoints::bcDOF, _TestTimeDependentPoints::numBCDOF);
  _bc->_getPoints(*_mesh);
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestTimeDependentPoints::tearDown(void)
{ // tearDown
  delete _mesh; _mesh = 0;
  delete _bc; _bc = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestTimeDependentPoints::testBCDOF(void)
{ // testBCDOF
  PointForce bc;

  const size_t numDOF = 4;
  const int fixedDOF[numDOF] = { 0, 2, 3, 5 };
  bc.bcDOF(fixedDOF, numDOF);

  CPPUNIT_ASSERT_EQUAL(numDOF, bc._bcDOF.size());
  for (int i=0; i < numDOF; ++i)
    CPPUNIT_ASSERT_EQUAL(fixedDOF[i], bc._bcDOF[i]);
} // testBCDOF

// ----------------------------------------------------------------------
// Test _getLabel().
void
pylith::bc::TestTimeDependentPoints::testGetLabel(void)
{ // testGetLabel
  PointForce bc;
  
  const std::string& label = "point force";
  bc.label(label.c_str());
  CPPUNIT_ASSERT_EQUAL(label, std::string(bc._getLabel()));
} // testGetLabel

// ----------------------------------------------------------------------
// Test _queryDatabases().
void
pylith::bc::TestTimeDependentPoints::testQueryDatabases(void)
{ // testQueryDatabases
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
  PetscSection initialSection = _bc->_parameters->get("initial").petscSection();
  Vec          initialVec     = _bc->_parameters->get("initial").localVector();
  CPPUNIT_ASSERT(initialSection);CPPUNIT_ASSERT(initialVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::initial, numBCDOF, initialSection, initialVec, forceScale);

  // Check rate values.
  PetscSection rateSection = _bc->_parameters->get("rate").petscSection();
  Vec          rateVec     = _bc->_parameters->get("rate").localVector();
  CPPUNIT_ASSERT(rateSection);CPPUNIT_ASSERT(rateVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::rate, numBCDOF, rateSection, rateVec, forceScale/timeScale);

  // Check rate start time.
  PetscSection rateTimeSection = _bc->_parameters->get("rate time").petscSection();
  Vec          rateTimeVec     = _bc->_parameters->get("rate time").localVector();
  CPPUNIT_ASSERT(rateTimeSection);CPPUNIT_ASSERT(rateTimeVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::rateTime, 1, rateTimeSection, rateTimeVec, timeScale);

  // Check change values.
  PetscSection changeSection = _bc->_parameters->get("change").petscSection();
  Vec          changeVec     = _bc->_parameters->get("change").localVector();
  CPPUNIT_ASSERT(changeSection);CPPUNIT_ASSERT(changeVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::change, numBCDOF, changeSection, changeVec, forceScale);

  // Check change start time.
  PetscSection changeTimeSection = _bc->_parameters->get("change time").petscSection();
  Vec          changeTimeVec     = _bc->_parameters->get("change time").localVector();
  CPPUNIT_ASSERT(changeTimeSection);CPPUNIT_ASSERT(changeTimeVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::changeTime, 1, changeTimeSection, changeTimeVec, timeScale);
  th.close();
} // testQueryDatabases

// ----------------------------------------------------------------------
// Test _calculateValue() with initial value.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueInitial(void)
{ // testCalculateValueInitial
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
  PetscSection valueSection = _bc->_parameters->get("value").petscSection();
  Vec          valueVec     = _bc->_parameters->get("value").localVector();
  CPPUNIT_ASSERT(valueSection);CPPUNIT_ASSERT(valueVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::initial, numBCDOF, valueSection, valueVec, forceScale);
} // testCalculateValueInitial

// ----------------------------------------------------------------------
// Test _calculateValue() with rate.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueRate(void)
{ // testCalculateValueRate
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
  PetscSection valueSection = _bc->_parameters->get("value").petscSection();
  Vec          valueVec     = _bc->_parameters->get("value").localVector();
  CPPUNIT_ASSERT(valueSection);CPPUNIT_ASSERT(valueVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesRate, numBCDOF, valueSection, valueVec, forceScale);
} // testCalculateValueRate

// ----------------------------------------------------------------------
// Test _calculateValue() with temporal change.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueChange(void)
{ // testCalculateValueChange
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
  PetscSection valueSection = _bc->_parameters->get("value").petscSection();
  Vec          valueVec     = _bc->_parameters->get("value").localVector();
  CPPUNIT_ASSERT(valueSection);CPPUNIT_ASSERT(valueVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesChange, numBCDOF, valueSection, valueVec, forceScale);
} // testCalculateValueChange

// ----------------------------------------------------------------------
// Test _calculateValue() with temporal change w/time history.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueChangeTH(void)
{ // testCalculateValueChangeTH
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
  PetscSection valueSection = _bc->_parameters->get("value").petscSection();
  Vec          valueVec     = _bc->_parameters->get("value").localVector();
  CPPUNIT_ASSERT(valueSection);CPPUNIT_ASSERT(valueVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesChangeTH, numBCDOF, valueSection, valueVec, forceScale);
} // testCalculateValueChangeTH

// ----------------------------------------------------------------------
// Test _calculateValue() with initial, rate, and temporal change w/time history.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueAll(void)
{ // testCalculateValueAll
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

  PetscSection valueSection = _bc->_parameters->get("value").petscSection();
  Vec          valueVec     = _bc->_parameters->get("value").localVector();
  CPPUNIT_ASSERT(valueSection);CPPUNIT_ASSERT(valueVec);
  _TestTimeDependentPoints::_checkValues(&valuesE[0], numBCDOF, valueSection, valueVec, forceScale);
} // testCalculateValueAll

// ----------------------------------------------------------------------
// Test _calculateValueIncr() with initial value.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueIncrInitial(void)
{ // testCalculateValueIncrInitial
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
  PetscSection valueSection = _bc->_parameters->get("value").petscSection();
  Vec          valueVec     = _bc->_parameters->get("value").localVector();
  CPPUNIT_ASSERT(valueSection);CPPUNIT_ASSERT(valueVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesIncrInitial, numBCDOF, valueSection, valueVec, forceScale);
} // testCalculateValueIncrInitial

// ----------------------------------------------------------------------
// Test _calculateValueIncr() with rate.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueIncrRate(void)
{ // testCalculateValueIncrRate
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
  PetscSection valueSection = _bc->_parameters->get("value").petscSection();
  Vec          valueVec     = _bc->_parameters->get("value").localVector();
  CPPUNIT_ASSERT(valueSection);CPPUNIT_ASSERT(valueVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesIncrRate, numBCDOF, valueSection, valueVec, forceScale);
} // testCalculateValueIncrRate

// ----------------------------------------------------------------------
// Test _calculateValueIncr() with temporal change.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueIncrChange(void)
{ // testCalculateValueIncrChange
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
  PetscSection valueSection = _bc->_parameters->get("value").petscSection();
  Vec          valueVec     = _bc->_parameters->get("value").localVector();
  CPPUNIT_ASSERT(valueSection);CPPUNIT_ASSERT(valueVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesIncrChange, numBCDOF, valueSection, valueVec, forceScale);
} // testCalculateValueIncrChange

// ----------------------------------------------------------------------
// Test _calculateValueIncr() with temporal change w/time history.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueIncrChangeTH(void)
{ // testCalculateValueIncrChangeTH
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
  PetscSection valueSection = _bc->_parameters->get("value").petscSection();
  Vec          valueVec     = _bc->_parameters->get("value").localVector();
  CPPUNIT_ASSERT(valueSection);CPPUNIT_ASSERT(valueVec);
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesIncrChangeTH, numBCDOF, valueSection, valueVec, forceScale);
} // testCalculateValueIncrChangeTH

// ----------------------------------------------------------------------
// Test _calculateValueIncr() with initial, rate, and temporal change w/time history.
void
pylith::bc::TestTimeDependentPoints::testCalculateValueIncrAll(void)
{ // testCalculateValueIncrAll
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

  PetscSection valueSection = _bc->_parameters->get("value").petscSection();
  Vec          valueVec     = _bc->_parameters->get("value").localVector();
  CPPUNIT_ASSERT(valueSection);CPPUNIT_ASSERT(valueVec);
  _TestTimeDependentPoints::_checkValues(&valuesE[0], numBCDOF, valueSection, valueVec, forceScale);
} // testCalculateValueIncrAll

// ----------------------------------------------------------------------
// Check values in section against expected values.
void
pylith::bc::_TestTimeDependentPoints::_checkValues(const PylithScalar* valuesE,
						   const int fiberDimE,
						   PetscSection section,
                           Vec vec,
						   const PylithScalar scale)
{ // _checkValues
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
  PetscScalar   *array;
  PetscErrorCode err;
  
  const PylithScalar tolerance = 1.0e-06;

  // Check values at points associated with BC.
  const int npointsIn = _TestTimeDependentPoints::npointsIn;
  err = VecGetArray(vec, &array);CHECK_PETSC_ERROR(err);
  for (int i=0; i < npointsIn; ++i) {
    PetscInt dof, off;
    const int p_bc = _TestTimeDependentPoints::pointsIn[i];

    err = PetscSectionGetDof(section, p_bc, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(section, p_bc, &off);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, dof);
    for (int iDim=0; iDim < fiberDimE; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i*fiberDimE+iDim]/scale, array[off+iDim], tolerance);
  } // for
  err = VecRestoreArray(vec, &array);CHECK_PETSC_ERROR(err);
} // _checkValues


// End of file 
