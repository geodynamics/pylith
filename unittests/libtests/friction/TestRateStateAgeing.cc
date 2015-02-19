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
// Test cutoff for linear slip rate.
void
pylith::friction::TestRateStateAgeing::testLinearSlipRate(void)
{ // testLinearSlipRate
  RateStateAgeing model;

  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0e-12), model._linearSlipRate); // default

  const PylithScalar value = 1.0e-20;
  model.linearSlipRate(value);
  CPPUNIT_ASSERT_EQUAL(value, model._linearSlipRate);
} // testLinearSlipRate
  
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
// Test hasPropStateVar().
void
pylith::friction::TestRateStateAgeing::testHasPropStateVar(void)
{ // testHasProperty
  RateStateAgeing material;

  CPPUNIT_ASSERT(material.hasPropStateVar("reference_friction_coefficient"));
  CPPUNIT_ASSERT(material.hasPropStateVar("reference_slip_rate"));
  CPPUNIT_ASSERT(material.hasPropStateVar("characteristic_slip_distance"));
  CPPUNIT_ASSERT(material.hasPropStateVar("constitutive_parameter_a"));
  CPPUNIT_ASSERT(material.hasPropStateVar("constitutive_parameter_b"));
  CPPUNIT_ASSERT(!material.hasPropStateVar("aaa"));
  CPPUNIT_ASSERT(material.hasPropStateVar("state_variable"));
} // testHasPropStateVar


// End of file 
