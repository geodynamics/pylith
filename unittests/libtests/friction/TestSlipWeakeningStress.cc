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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestSlipWeakeningStress.hh" // Implementation of class methods

#include "data/SlipWeakeningStressData.hh" // USES SlipWeakeningStressData

#include "pylith/friction/SlipWeakeningStress.hh" // USES SlipWeakeningStress

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::friction::TestSlipWeakeningStress );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::friction::TestSlipWeakeningStress::setUp(void)
{ // setUp
  _friction = new SlipWeakeningStress();
  _data = new SlipWeakeningStressData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test properties metadata.
void
pylith::friction::TestSlipWeakeningStress::testPropertiesMetadata(void)
{ // testPropertiesMetadata
  SlipWeakeningStress model;

  CPPUNIT_ASSERT_EQUAL(4, model._metadata.numDBProperties());
  const char* const* names = model._metadata.dbProperties();
  CPPUNIT_ASSERT_EQUAL(std::string("static-stress"), std::string(names[0]));
  CPPUNIT_ASSERT_EQUAL(std::string("dynamic-stress"), std::string(names[1]));
  CPPUNIT_ASSERT_EQUAL(std::string("slip-weakening-parameter"), std::string(names[2]));
  CPPUNIT_ASSERT_EQUAL(std::string("weakening-time"), std::string(names[3]));
} // testPropertiesMetadata

// ----------------------------------------------------------------------
// Test state variable metadata.
void
pylith::friction::TestSlipWeakeningStress::testStateVarsMetadata(void)
{ // testStateVarsMetadata
  SlipWeakeningStress model;

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
pylith::friction::TestSlipWeakeningStress::testHasPropStateVar(void)
{ // testHasPropStateVar
  SlipWeakeningStress material;

  CPPUNIT_ASSERT(material.hasPropStateVar("static_stress"));
  CPPUNIT_ASSERT(material.hasPropStateVar("dynamic_stress"));
  CPPUNIT_ASSERT(material.hasPropStateVar("slip_weakening_parameter"));
  CPPUNIT_ASSERT(material.hasPropStateVar("weakening_time"));
  CPPUNIT_ASSERT(!material.hasPropStateVar("aaa"));
  CPPUNIT_ASSERT(material.hasPropStateVar("cumulative_slip"));
  CPPUNIT_ASSERT(material.hasPropStateVar("previous_slip"));
} // testHasPropStateVar


// End of file 
