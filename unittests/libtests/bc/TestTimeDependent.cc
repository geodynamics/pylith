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

#include "TestTimeDependent.hh" // Implementation of class methods

#include "pylith/bc/PointForce.hh" // USES PointForce
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "spatialdata/spatialdb/UniformDB.hh" // USES UniformDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestTimeDependent );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestTimeDependent::testBCDOF(void)
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
// Test dbInitial().
void
pylith::bc::TestTimeDependent::testDBInitial(void)
{ // testDBInitial
  PointForce bc;

  spatialdata::spatialdb::UniformDB db;
  bc.dbInitial(&db);

  CPPUNIT_ASSERT(0 != bc._dbInitial);
  CPPUNIT_ASSERT(0 == bc._dbRate);
  CPPUNIT_ASSERT(0 == bc._dbRateTime);
  CPPUNIT_ASSERT(0 == bc._dbChange);
  CPPUNIT_ASSERT(0 == bc._dbChangeTime);
  CPPUNIT_ASSERT(0 == bc._dbTimeHistory);
} // testDBInitial

// ----------------------------------------------------------------------
// Test dbRate().
void
pylith::bc::TestTimeDependent::testDBRate(void)
{ // testDBRate
  PointForce bc;

  spatialdata::spatialdb::UniformDB db;
  bc.dbRate(&db);

  CPPUNIT_ASSERT(0 == bc._dbInitial);
  CPPUNIT_ASSERT(0 != bc._dbRate);
  CPPUNIT_ASSERT(0 == bc._dbRateTime);
  CPPUNIT_ASSERT(0 == bc._dbChange);
  CPPUNIT_ASSERT(0 == bc._dbChangeTime);
  CPPUNIT_ASSERT(0 == bc._dbTimeHistory);
} // testDBRate

// ----------------------------------------------------------------------
// Test dbRateTime().
void
pylith::bc::TestTimeDependent::testDBRateTime(void)
{ // testDBRateTime
  PointForce bc;

  spatialdata::spatialdb::UniformDB db;
  bc.dbRateTime(&db);

  CPPUNIT_ASSERT(0 == bc._dbInitial);
  CPPUNIT_ASSERT(0 == bc._dbRate);
  CPPUNIT_ASSERT(0 != bc._dbRateTime);
  CPPUNIT_ASSERT(0 == bc._dbChange);
  CPPUNIT_ASSERT(0 == bc._dbChangeTime);
  CPPUNIT_ASSERT(0 == bc._dbTimeHistory);
} // testDBRateTime

// ----------------------------------------------------------------------
// Test dbChange().
void
pylith::bc::TestTimeDependent::testDBChange(void)
{ // testDBChange
  PointForce bc;

  spatialdata::spatialdb::UniformDB db;
  bc.dbChange(&db);

  CPPUNIT_ASSERT(0 == bc._dbInitial);
  CPPUNIT_ASSERT(0 == bc._dbRate);
  CPPUNIT_ASSERT(0 == bc._dbRateTime);
  CPPUNIT_ASSERT(0 != bc._dbChange);
  CPPUNIT_ASSERT(0 == bc._dbChangeTime);
  CPPUNIT_ASSERT(0 == bc._dbTimeHistory);
} // testDBChange

// ----------------------------------------------------------------------
// Test dbChangeTime().
void
pylith::bc::TestTimeDependent::testDBChangeTime(void)
{ // testDBChangeTime
  PointForce bc;

  spatialdata::spatialdb::UniformDB db;
  bc.dbChangeTime(&db);

  CPPUNIT_ASSERT(0 == bc._dbInitial);
  CPPUNIT_ASSERT(0 == bc._dbRate);
  CPPUNIT_ASSERT(0 == bc._dbRateTime);
  CPPUNIT_ASSERT(0 == bc._dbChange);
  CPPUNIT_ASSERT(0 != bc._dbChangeTime);
  CPPUNIT_ASSERT(0 == bc._dbTimeHistory);
} // testDBChangeTime

// ----------------------------------------------------------------------
// Test dbTimeHistory().
void
pylith::bc::TestTimeDependent::testDBTimeHistory(void)
{ // testDBTimeHistory
  PointForce bc;

  spatialdata::spatialdb::TimeHistory th;
  bc.dbTimeHistory(&th);

  CPPUNIT_ASSERT(0 == bc._dbInitial);
  CPPUNIT_ASSERT(0 == bc._dbRate);
  CPPUNIT_ASSERT(0 == bc._dbRateTime);
  CPPUNIT_ASSERT(0 == bc._dbChange);
  CPPUNIT_ASSERT(0 == bc._dbChangeTime);
  CPPUNIT_ASSERT(0 != bc._dbTimeHistory);
} // testDBTimeHistory

// ----------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::bc::TestTimeDependent::testVerifyConfiguration(void)
{ // testVerifyConfiguration
  topology::Mesh mesh;
  spatialdata::spatialdb::UniformDB db;
  spatialdata::spatialdb::TimeHistory th;

  { // initial
    PointForce bc;
    bc.dbInitial(&db);
    bc.TimeDependent::verifyConfiguration(mesh);
  } // initial

  { // rate
    PointForce bc;
    bc.dbRate(&db);
    bc.dbRateTime(&db);
    bc.TimeDependent::verifyConfiguration(mesh);
  } // rate

  { // change
    PointForce bc;
    bc.dbChange(&db);
    bc.dbChangeTime(&db);
    bc.dbTimeHistory(&th);
    bc.TimeDependent::verifyConfiguration(mesh);
  } // change

  { // change
    PointForce bc;
    bc.dbChange(&db);
    bc.dbChangeTime(&db);
    bc.TimeDependent::verifyConfiguration(mesh);
  } // change

  { // rate (missing start time)
    PointForce bc;
    bc.dbRate(&db);
    
    bool caught = false;
    try {
      bc.TimeDependent::verifyConfiguration(mesh);
    } catch ( const std::exception& err) {
      caught = true;
    } // catch
    CPPUNIT_ASSERT(caught);
  } // rate (missing start time)

  { // rate (missing rate)
    PointForce bc;
    bc.dbRateTime(&db);
    
    bool caught = false;
    try {
      bc.TimeDependent::verifyConfiguration(mesh);
    } catch ( const std::exception& err) {
      caught = true;
    } // catch
    CPPUNIT_ASSERT(caught);
  } // rate (missing rate)

  { // change (missing start time)
    PointForce bc;
    bc.dbChange(&db);
    
    bool caught = false;
    try {
      bc.TimeDependent::verifyConfiguration(mesh);
    } catch ( const std::exception& err) {
      caught = true;
    } // catch
    CPPUNIT_ASSERT(caught);
  } // change (missing start time)

  { // change (missing change)
    PointForce bc;
    bc.dbChangeTime(&db);
    
    bool caught = false;
    try {
      bc.TimeDependent::verifyConfiguration(mesh);
    } catch ( const std::exception& err) {
      caught = true;
    } // catch
    CPPUNIT_ASSERT(caught);
  } // change (missing change)

  { // change (missing start time)
    PointForce bc;
    bc.dbChange(&db);
    
    bool caught = false;
    try {
      bc.TimeDependent::verifyConfiguration(mesh);
    } catch ( const std::exception& err) {
      caught = true;
    } // catch
    CPPUNIT_ASSERT(caught);
  } // change (missing start time)

  { // change (missing change)
    PointForce bc;
    bc.dbTimeHistory(&th);
    bc.dbChangeTime(&db);

    bool caught = false;
    try {
      bc.TimeDependent::verifyConfiguration(mesh);
    } catch ( const std::exception& err) {
      caught = true;
    } // catch
    CPPUNIT_ASSERT(caught);
  } // change (missing change)

  { // change (missing start time)
    PointForce bc;
    bc.dbTimeHistory(&th);
    bc.dbChange(&db);
    
    bool caught = false;
    try {
      bc.TimeDependent::verifyConfiguration(mesh);
    } catch ( const std::exception& err) {
      caught = true;
    } // catch
    CPPUNIT_ASSERT(caught);
  } // change (missing start time)

} // testVerifyConfiguration


// End of file 
