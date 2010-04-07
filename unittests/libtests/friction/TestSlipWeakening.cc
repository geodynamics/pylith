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

#include "TestSlipWeakening.hh" // Implementation of class methods

#include "data/SlipWeakeningData.hh" // USES SlipWeakeningData

#include "pylith/friction/SlipWeakening.hh" // USES SlipWeakening

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::friction::TestSlipWeakening );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::friction::TestSlipWeakening::setUp(void)
{ // setUp
  _friction = new SlipWeakening();
  _data = new SlipWeakeningData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test properties metadata.
void
pylith::friction::TestSlipWeakening::testPropertiesMetadata(void)
{ // testPropertiesMetadata
  SlipWeakening model;

  CPPUNIT_ASSERT_EQUAL(4, model._metadata.numDBProperties());
  const char* const* names = model._metadata.dbProperties();
  CPPUNIT_ASSERT_EQUAL(std::string("static-coefficient"), 
		       std::string(names[0]));
  CPPUNIT_ASSERT_EQUAL(std::string("dynamic-coefficient"), 
		       std::string(names[1]));
  CPPUNIT_ASSERT_EQUAL(std::string("slip-weakening-parameter"), 
		       std::string(names[2]));
  CPPUNIT_ASSERT_EQUAL(std::string("cohesion"),
		       std::string(names[3]));
} // testPropertiesMetadata

// ----------------------------------------------------------------------
// Test state variable metadata.
void
pylith::friction::TestSlipWeakening::testStateVarsMetadata(void)
{ // testStateVarsMetadata
  SlipWeakening model;

  CPPUNIT_ASSERT_EQUAL(2, model._metadata.numDBStateVars());
  const char* const* names = model._metadata.dbStateVars();
  CPPUNIT_ASSERT_EQUAL(std::string("cumulative-slip"), 
		       std::string(names[0]));
  CPPUNIT_ASSERT_EQUAL(std::string("previous-slip"), 
		       std::string(names[1]));
} // testStateVarsMetadata

// ----------------------------------------------------------------------
// Test hasProperty().
void
pylith::friction::TestSlipWeakening::testHasProperty(void)
{ // testHasProperty
  SlipWeakening material;

  CPPUNIT_ASSERT(material.hasProperty("static_coefficient"));
  CPPUNIT_ASSERT(material.hasProperty("dynamic_coefficient"));
  CPPUNIT_ASSERT(material.hasProperty("slip_weakening_parameter"));
  CPPUNIT_ASSERT(!material.hasProperty("aaa"));
} // testHasProperty

// ----------------------------------------------------------------------
// Test hasStateVar().
void
pylith::friction::TestSlipWeakening::testHasStateVar(void)
{ // testHasStateVar
  SlipWeakening material;

  CPPUNIT_ASSERT(material.hasStateVar("cumulative_slip"));
  CPPUNIT_ASSERT(material.hasStateVar("previous_slip"));
  CPPUNIT_ASSERT(!material.hasStateVar("aaa"));
} // testHasStateVar


// End of file 
