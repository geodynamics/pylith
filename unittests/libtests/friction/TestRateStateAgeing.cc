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
// Test properties metadata.
void
pylith::friction::TestRateStateAgeing::testPropertiesMetadata(void)
{ // testPropertiesMetadata
  RateStateAgeing model;

  CPPUNIT_ASSERT_EQUAL(6, model._metadata.numDBProperties());
  const char* const* names = model._metadata.dbProperties();
  CPPUNIT_ASSERT_EQUAL(std::string("reference-friction-coefficient"), 
		       std::string(names[0]));
  CPPUNIT_ASSERT_EQUAL(std::string("reference-slip-rate"), 
		       std::string(names[1]));
  CPPUNIT_ASSERT_EQUAL(std::string("characteristic-slip-distance"), 
		       std::string(names[2]));
  CPPUNIT_ASSERT_EQUAL(std::string("constitutive-parameter-a"), 
		       std::string(names[3]));
  CPPUNIT_ASSERT_EQUAL(std::string("constitutive-parameter-b"), 
		       std::string(names[4]));
  CPPUNIT_ASSERT_EQUAL(std::string("cohesion"),
		       std::string(names[5]));
} // testPropertiesMetadata

// ----------------------------------------------------------------------
// Test state variable metadata.
void
pylith::friction::TestRateStateAgeing::testStateVarsMetadata(void)
{ // testStateVarsMetadata
  RateStateAgeing model;

  CPPUNIT_ASSERT_EQUAL(1, model._metadata.numDBStateVars());
  const char* const* names = model._metadata.dbStateVars();
  CPPUNIT_ASSERT_EQUAL(std::string("state-variable"), 
		       std::string(names[0]));
} // testStateVarsMetadata

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
