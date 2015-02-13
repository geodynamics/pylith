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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
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
// Test properties metadata.
void
pylith::friction::TestStaticFriction::testPropertiesMetadata(void)
{ // testPropertiesMetadata
  StaticFriction model;

  CPPUNIT_ASSERT_EQUAL(2, model._metadata.numDBProperties());
  const char* const* names = model._metadata.dbProperties();
  CPPUNIT_ASSERT_EQUAL(std::string("friction-coefficient"), 
		       std::string(names[0]));
  CPPUNIT_ASSERT_EQUAL(std::string("cohesion"),
		       std::string(names[1]));
} // testPropertiesMetadata

// ----------------------------------------------------------------------
// Test state variable metadata.
void
pylith::friction::TestStaticFriction::testStateVarsMetadata(void)
{ // testStateVarsMetadata
  StaticFriction model;

  CPPUNIT_ASSERT_EQUAL(0, model._metadata.numDBStateVars());
} // testStateVarsMetadata

// ----------------------------------------------------------------------
// Test hasProperty().
void
pylith::friction::TestStaticFriction::testHasPropStateVar(void)
{ // testHasPropStateVar
  StaticFriction model;

  CPPUNIT_ASSERT(model.hasPropStateVar("friction_coefficient"));
  CPPUNIT_ASSERT(!model.hasPropStateVar("aaa"));
} // testHasPropStateVar


// End of file 
