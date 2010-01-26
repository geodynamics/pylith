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
// Test hasProperty().
void
pylith::friction::TestSlipWeakening::testHasProperty(void)
{ // testHasProperty
  SlipWeakening material;

  CPPUNIT_ASSERT(material.hasProperty("static-coefficient"));
  CPPUNIT_ASSERT(material.hasProperty("dynamic-coefficient"));
  CPPUNIT_ASSERT(material.hasProperty("slip-weakening-parameter"));
  CPPUNIT_ASSERT(!material.hasProperty("aaa"));
} // testHasProperty

// ----------------------------------------------------------------------
// Test hasStateVar().
void
pylith::friction::TestSlipWeakening::testHasStateVar(void)
{ // testHasStateVar
  SlipWeakening material;

  CPPUNIT_ASSERT(material.hasStateVar("cumulative-slip"));
  CPPUNIT_ASSERT(material.hasStateVar("previous-slip"));
  CPPUNIT_ASSERT(!material.hasStateVar("aaa"));
} // testHasStateVar


// End of file 
