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

#include "TestSlipWeakeningTime.hh" // Implementation of class methods

#include "data/SlipWeakeningTimeData.hh" // USES SlipWeakeningTimeData

#include "pylith/friction/SlipWeakeningTime.hh" // USES SlipWeakeningTime

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::friction::TestSlipWeakeningTime );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::friction::TestSlipWeakeningTime::setUp(void)
{ // setUp
  _friction = new SlipWeakeningTime();
  _data = new SlipWeakeningTimeData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test properties metadata.
void
pylith::friction::TestSlipWeakeningTime::testPropertiesMetadata(void)
{ // testPropertiesMetadata
  SlipWeakeningTime model;

  CPPUNIT_ASSERT_EQUAL(5, model._metadata.numDBProperties());
  const char* const* names = model._metadata.dbProperties();
  CPPUNIT_ASSERT_EQUAL(std::string("static-coefficient"), 
		       std::string(names[0]));
  CPPUNIT_ASSERT_EQUAL(std::string("dynamic-coefficient"), 
		       std::string(names[1]));
  CPPUNIT_ASSERT_EQUAL(std::string("slip-weakening-parameter"), 
		       std::string(names[2]));
  CPPUNIT_ASSERT_EQUAL(std::string("cohesion"),
		       std::string(names[3]));
  CPPUNIT_ASSERT_EQUAL(std::string("weakening-time"),
		       std::string(names[4]));
} // testPropertiesMetadata

// ----------------------------------------------------------------------
// Test state variable metadata.
void
pylith::friction::TestSlipWeakeningTime::testStateVarsMetadata(void)
{ // testStateVarsMetadata
  SlipWeakeningTime model;

  CPPUNIT_ASSERT_EQUAL(2, model._metadata.numDBStateVars());
  const char* const* names = model._metadata.dbStateVars();
  CPPUNIT_ASSERT_EQUAL(std::string("cumulative-slip"), 
		       std::string(names[0]));
  CPPUNIT_ASSERT_EQUAL(std::string("previous-slip"), 
		       std::string(names[1]));
} // testStateVarsMetadata

// ----------------------------------------------------------------------
// Test hasPropStateVar().
void
pylith::friction::TestSlipWeakeningTime::testHasPropStateVar(void)
{ // testHasPropStateVar
  SlipWeakeningTime material;

  CPPUNIT_ASSERT(material.hasPropStateVar("static_coefficient"));
  CPPUNIT_ASSERT(material.hasPropStateVar("dynamic_coefficient"));
  CPPUNIT_ASSERT(material.hasPropStateVar("slip_weakening_parameter"));
  CPPUNIT_ASSERT(material.hasPropStateVar("cohesion"));
  CPPUNIT_ASSERT(material.hasPropStateVar("weakening_time"));
  CPPUNIT_ASSERT(!material.hasPropStateVar("aaa"));
  CPPUNIT_ASSERT(material.hasPropStateVar("cumulative_slip"));
  CPPUNIT_ASSERT(material.hasPropStateVar("previous_slip"));
} // testHasPropStateVar


// End of file 
