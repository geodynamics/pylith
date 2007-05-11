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

#include "TestBruneSlipFn.hh" // Implementation of class methods

#include "pylith/faults/BruneSlipFn.hh" // USES BruneSlipFn

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB

#include <stdexcept> // TEMPORARY
// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestBruneSlipFn );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestBruneSlipFn::testConstructor(void)
{ // testConstructor
  BruneSlipFn slipfn;
} // testConstructor

// ----------------------------------------------------------------------
// Test clone()
void
pylith::faults::TestBruneSlipFn::testClone(void)
{ // testClone
  BruneSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB dbFinalSlip("final slip");
  slipfn.dbFinalSlip(&dbFinalSlip);
  spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
  slipfn.dbSlipTime(&dbSlipTime);
  spatialdata::spatialdb::SimpleDB dbPeakRate("peak rate");
  slipfn.dbPeakRate(&dbPeakRate);

  SlipTimeFn* slipfnCopy = slipfn.clone();
  CPPUNIT_ASSERT(0 != slipfnCopy);

  BruneSlipFn* bruneCopy = dynamic_cast<BruneSlipFn*>(slipfnCopy);
  CPPUNIT_ASSERT(0 != bruneCopy);

  CPPUNIT_ASSERT(0 != bruneCopy->_dbFinalSlip);
  CPPUNIT_ASSERT_EQUAL(std::string("final slip"),
		       std::string(bruneCopy->_dbFinalSlip->label()));
  CPPUNIT_ASSERT(0 != bruneCopy->_dbSlipTime);
  CPPUNIT_ASSERT_EQUAL(std::string("slip time"),
		       std::string(bruneCopy->_dbSlipTime->label()));
  CPPUNIT_ASSERT(0 != bruneCopy->_dbPeakRate);
  CPPUNIT_ASSERT_EQUAL(std::string("peak rate"),
		       std::string(bruneCopy->_dbPeakRate->label()));

  delete slipfnCopy; slipfnCopy = 0;
} // testClone

// ----------------------------------------------------------------------
// Test dbFinalSlip().
void
pylith::faults::TestBruneSlipFn::testDbFinalSlip(void)
{ // testDbFinalSlip
  const char* label = "database ABC";
  BruneSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbFinalSlip(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbFinalSlip);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbFinalSlip->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbSlipTime);
  CPPUNIT_ASSERT(0 == slipfn._dbPeakRate);
} // testDbFinalSlip

// ----------------------------------------------------------------------
// Test dbSlipTime().
void
pylith::faults::TestBruneSlipFn::testDbSlipTime(void)
{ // testDbSlipTime
  const char* label = "database ABCD";
  BruneSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbSlipTime(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbSlipTime);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbSlipTime->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbFinalSlip);
  CPPUNIT_ASSERT(0 == slipfn._dbPeakRate);
} // testDbSlipTime

// ----------------------------------------------------------------------
// Test dbPeakRate().
void
pylith::faults::TestBruneSlipFn::testDbPeakRate(void)
{ // testDbPeakRate
  const char* label = "database ABCDE";
  BruneSlipFn slipfn;
  
  spatialdata::spatialdb::SimpleDB db(label);
  slipfn.dbPeakRate(&db);

  CPPUNIT_ASSERT(0 != slipfn._dbPeakRate);
  CPPUNIT_ASSERT_EQUAL(std::string(label),
		       std::string(slipfn._dbPeakRate->label()));
  CPPUNIT_ASSERT(0 == slipfn._dbFinalSlip);
  CPPUNIT_ASSERT(0 == slipfn._dbSlipTime);
} // testDbPeakRate

// ----------------------------------------------------------------------
// Test initialize() in 1-D.
void
pylith::faults::TestBruneSlipFn::testInitialize1D(void)
{ // testInitialize1D
  throw std::logic_error("Unit test not implemented.");
} // testInitialize1D

// ----------------------------------------------------------------------
// Test initialize() in 2-D.
void
pylith::faults::TestBruneSlipFn::testInitialize2D(void)
{ // testInitialize2D
  throw std::logic_error("Unit test not implemented.");
} // testInitialize2D

// ----------------------------------------------------------------------
// Test initialize() in 3-D.
void
pylith::faults::TestBruneSlipFn::testInitialize3D(void)
{ // testInitialize3D
  throw std::logic_error("Unit test not implemented.");
} // testInitialize3D

// ----------------------------------------------------------------------
// Test slip().
void
pylith::faults::TestBruneSlipFn::testSlip(void)
{ // testSlip
  throw std::logic_error("Unit test not implemented.");
} // testSlip

// ----------------------------------------------------------------------
// Test _slip().
void
pylith::faults::TestBruneSlipFn::testSlipTH(void)
{ // testSlipTH
  const double t = 0.734;
  const double finalSlip = 4.64;
  const double peakRate = 3.23;

  const double tau = finalSlip / (exp(1.0) * peakRate);
  const double slipE = finalSlip * (1.0 - exp(-t/tau) * (1.0 + t/tau));

  double slip = BruneSlipFn::_slip(t, finalSlip, peakRate);

  const double tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(slipE, slip, tolerance);

  slip = BruneSlipFn::_slip(-0.5, finalSlip, peakRate);
  CPPUNIT_ASSERT_EQUAL(0.0, slip);

  slip = BruneSlipFn::_slip(1.0e+10, finalSlip, peakRate);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(finalSlip, slip, tolerance);
} // testSlipTH



// End of file 
