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

#include "TestStaticFriction.hh" // Implementation of class methods

#include "data/StaticFrictionData.hh" // USES StaticFrictionData

#include "pylith/friction/StaticFriction.hh" // USES StaticFriction

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::friction::TestStaticFriction );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::friction::TestStaticFriction::setUp(void)
{ // setUp
  _friction = new StaticFriction();
  _data = new StaticFrictionData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test hasProperty().
void
pylith::friction::TestStaticFriction::testHasProperty(void)
{ // testHasProperty
  StaticFriction material;

  CPPUNIT_ASSERT(material.hasProperty("friction_coefficient"));
  CPPUNIT_ASSERT(!material.hasProperty("aaa"));
} // testHasProperty

// ----------------------------------------------------------------------
// Test hasStateVar().
void
pylith::friction::TestStaticFriction::testHasStateVar(void)
{ // testHasStateVar
  StaticFriction material;

  CPPUNIT_ASSERT(!material.hasStateVar("aaa"));
} // testHasStateVar


// End of file 
