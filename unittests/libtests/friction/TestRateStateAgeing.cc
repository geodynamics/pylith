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

#include "TestRateStateAgeing.hh" // Implementation of class methods

#include "data/RateStateAgeingData.hh" // USES RateStateAgeingData

#include "pylith/friction/RateStateAgeing.hh" // USES RateStateAgeing

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::friction::TestRateStateAgeing );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::friction::TestRateStateAgeing::setUp(void)
{ // setUp
  _friction = new RateStateAgeing();
  _data = new RateStateAgeingData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test hasProperty().
void
pylith::friction::TestRateStateAgeing::testHasProperty(void)
{ // testHasProperty
  RateStateAgeing material;

  CPPUNIT_ASSERT(material.hasProperty("reference_friction_coefficient"));
  CPPUNIT_ASSERT(material.hasProperty("reference_slip_rate"));
  CPPUNIT_ASSERT(material.hasProperty("characteristic_slip_distance"));
  CPPUNIT_ASSERT(material.hasProperty("constitutive_parameter_a"));
  CPPUNIT_ASSERT(material.hasProperty("constitutive_parameter_b"));
  CPPUNIT_ASSERT(!material.hasProperty("aaa"));
} // testHasProperty

// ----------------------------------------------------------------------
// Test hasStateVar().
void
pylith::friction::TestRateStateAgeing::testHasStateVar(void)
{ // testHasStateVar
  RateStateAgeing material;

  CPPUNIT_ASSERT(material.hasStateVar("state_variable"));
  CPPUNIT_ASSERT(!material.hasStateVar("aaa"));
} // testHasStateVar


// End of file 
