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
  CPPUNIT_ASSERT(0 == bc._dbChange);
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
  CPPUNIT_ASSERT(0 == bc._dbChange);
  CPPUNIT_ASSERT(0 == bc._dbTimeHistory);
} // testDBRate

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
  CPPUNIT_ASSERT(0 != bc._dbChange);
  CPPUNIT_ASSERT(0 == bc._dbTimeHistory);
} // testDBChange

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
  CPPUNIT_ASSERT(0 == bc._dbChange);
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
    bc.TimeDependent::verifyConfiguration(mesh);
  } // rate

  { // change
    PointForce bc;
    bc.dbChange(&db);
    bc.dbTimeHistory(&th);
    bc.TimeDependent::verifyConfiguration(mesh);
  } // change

  { // change (missing change)
    PointForce bc;
    bc.dbTimeHistory(&th);

    bool caught = false;
    try {
      bc.TimeDependent::verifyConfiguration(mesh);
    } catch ( const std::exception& err) {
      caught = true;
    } // catch
    CPPUNIT_ASSERT(caught);
  } // change (missing change)

} // testVerifyConfiguration


// End of file 
