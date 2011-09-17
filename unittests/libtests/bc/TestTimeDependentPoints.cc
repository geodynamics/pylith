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
// Copyright (c) 2010-2011 University of California, Davis
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
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestTimeDependentPoints );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

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
			const ALE::Obj<RealSection>& section,
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
  const ALE::Obj<RealSection>& initialSection = 
    _bc->_parameters->get("initial").section();
  CPPUNIT_ASSERT(!initialSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::initial,
					 numBCDOF, initialSection, forceScale);

  // Check rate values.
  const ALE::Obj<RealSection>& rateSection = 
    _bc->_parameters->get("rate").section();
  CPPUNIT_ASSERT(!rateSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::rate,
					 numBCDOF, rateSection, forceScale/timeScale);

  // Check rate start time.
  const ALE::Obj<RealSection>& rateTimeSection = 
    _bc->_parameters->get("rate time").section();
  CPPUNIT_ASSERT(!rateTimeSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::rateTime,
					 1, rateTimeSection, timeScale);

  // Check change values.
  const ALE::Obj<RealSection>& changeSection = 
    _bc->_parameters->get("change").section();
  CPPUNIT_ASSERT(!changeSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::change,
					 numBCDOF, changeSection, forceScale);

  // Check change start time.
  const ALE::Obj<RealSection>& changeTimeSection = 
    _bc->_parameters->get("change time").section();
  CPPUNIT_ASSERT(!changeTimeSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::changeTime,
					 1, changeTimeSection, timeScale);
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
  const ALE::Obj<RealSection>& valueSection = 
    _bc->_parameters->get("value").section();
  CPPUNIT_ASSERT(!valueSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::initial,
					 numBCDOF, valueSection, forceScale);
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
  const ALE::Obj<RealSection>& valueSection = 
    _bc->_parameters->get("value").section();
  CPPUNIT_ASSERT(!valueSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesRate,
					 numBCDOF, valueSection, forceScale);
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
  const ALE::Obj<RealSection>& valueSection = 
    _bc->_parameters->get("value").section();
  CPPUNIT_ASSERT(!valueSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesChange,
					 numBCDOF, valueSection, forceScale);
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
  const ALE::Obj<RealSection>& valueSection = 
    _bc->_parameters->get("value").section();
  CPPUNIT_ASSERT(!valueSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesChangeTH,
					 numBCDOF, valueSection, forceScale);
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

  const ALE::Obj<RealSection>& valueSection = 
    _bc->_parameters->get("value").section();
  CPPUNIT_ASSERT(!valueSection.isNull());
  _TestTimeDependentPoints::_checkValues(&valuesE[0],
					 numBCDOF, valueSection, forceScale);
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
  const ALE::Obj<RealSection>& valueSection = 
    _bc->_parameters->get("value").section();
  CPPUNIT_ASSERT(!valueSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesIncrInitial,
					 numBCDOF, valueSection, forceScale);
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
  const ALE::Obj<RealSection>& valueSection = 
    _bc->_parameters->get("value").section();
  CPPUNIT_ASSERT(!valueSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesIncrRate,
					 numBCDOF, valueSection, forceScale);
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
  const ALE::Obj<RealSection>& valueSection = 
    _bc->_parameters->get("value").section();
  CPPUNIT_ASSERT(!valueSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesIncrChange,
					 numBCDOF, valueSection, forceScale);
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
  const ALE::Obj<RealSection>& valueSection = 
    _bc->_parameters->get("value").section();
  CPPUNIT_ASSERT(!valueSection.isNull());
  _TestTimeDependentPoints::_checkValues(_TestTimeDependentPoints::valuesIncrChangeTH,
					 numBCDOF, valueSection, forceScale);
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

  const ALE::Obj<RealSection>& valueSection = 
    _bc->_parameters->get("value").section();
  CPPUNIT_ASSERT(!valueSection.isNull());
  _TestTimeDependentPoints::_checkValues(&valuesE[0],
					 numBCDOF, valueSection, forceScale);
} // testCalculateValueIncrAll

// ----------------------------------------------------------------------
// Check values in section against expected values.
void
pylith::bc::_TestTimeDependentPoints::_checkValues(const PylithScalar* valuesE,
						   const int fiberDimE,
						   const ALE::Obj<RealSection>& section,
						   const PylithScalar scale)
{ // _checkValues
  CPPUNIT_ASSERT(!section.isNull());
  
  const PylithScalar tolerance = 1.0e-06;

  // Check values at points associated with BC.
  const int npointsIn = _TestTimeDependentPoints::npointsIn;
  for (int i=0; i < npointsIn; ++i) {
    const int p_bc = _TestTimeDependentPoints::pointsIn[i];
    const int fiberDim = section->getFiberDimension(p_bc);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);

    const PylithScalar* values = section->restrictPoint(p_bc);
    CPPUNIT_ASSERT(values);
    for (int iDim=0; iDim < fiberDimE; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesE[i*fiberDimE+iDim]/scale,
				   values[iDim], tolerance);
  } // for
} // _checkValues


// End of file 
