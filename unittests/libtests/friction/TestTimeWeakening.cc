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

#include "TestTimeWeakening.hh" // Implementation of class methods

#include "data/TimeWeakeningData.hh" // USES TimeWeakeningData

#include "pylith/friction/TimeWeakening.hh" // USES TimeWeakening

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::friction::TestTimeWeakening );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::friction::TestTimeWeakening::setUp(void)
{ // setUp
  _friction = new TimeWeakening();
  _data = new TimeWeakeningData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test properties metadata.
void
pylith::friction::TestTimeWeakening::testPropertiesMetadata(void)
{ // testPropertiesMetadata
  TimeWeakening model;

  CPPUNIT_ASSERT_EQUAL(4, model._metadata.numDBProperties());
  const char* const* names = model._metadata.dbProperties();
  CPPUNIT_ASSERT_EQUAL(std::string("static-coefficient"), 
		       std::string(names[0]));
  CPPUNIT_ASSERT_EQUAL(std::string("dynamic-coefficient"), 
		       std::string(names[1]));
  CPPUNIT_ASSERT_EQUAL(std::string("time-weakening-parameter"), 
		       std::string(names[2]));
  CPPUNIT_ASSERT_EQUAL(std::string("cohesion"),
		       std::string(names[3]));
} // testPropertiesMetadata

// ----------------------------------------------------------------------
// Test state variable metadata.
void
pylith::friction::TestTimeWeakening::testStateVarsMetadata(void)
{ // testStateVarsMetadata
  TimeWeakening model;

  CPPUNIT_ASSERT_EQUAL(1, model._metadata.numDBStateVars());
  const char* const* names = model._metadata.dbStateVars();
  CPPUNIT_ASSERT_EQUAL(std::string("elapsed-time"), 
		       std::string(names[0]));
} // testStateVarsMetadata

// ----------------------------------------------------------------------
// Test hasProperty().
void
pylith::friction::TestTimeWeakening::testHasProperty(void)
{ // testHasProperty
  TimeWeakening material;

  CPPUNIT_ASSERT(material.hasProperty("static_coefficient"));
  CPPUNIT_ASSERT(material.hasProperty("dynamic_coefficient"));
  CPPUNIT_ASSERT(material.hasProperty("time_weakening_parameter"));
  CPPUNIT_ASSERT(!material.hasProperty("aaa"));
} // testHasProperty

// ----------------------------------------------------------------------
// Test hasStateVar().
void
pylith::friction::TestTimeWeakening::testHasStateVar(void)
{ // testHasStateVar
  TimeWeakening material;

  CPPUNIT_ASSERT(material.hasStateVar("elapsed_time"));
  CPPUNIT_ASSERT(!material.hasStateVar("aaa"));
} // testHasStateVar


// End of file 
