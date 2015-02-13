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

#include "TestMetadata.hh" // Implementation of class methods

#include "pylith/materials/Metadata.hh" // USES Metadata
#include "pylith/utils/array.hh" // USES scalar_array

#include <cstring> // USES strcmp()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestMetadata );

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _TestMetadata {
      const Metadata::ParamDescription properties[] = {
	{ "one", 1, pylith::topology::FieldBase::SCALAR },
	{ "two", 2, pylith::topology::FieldBase::VECTOR },
	{ "three", 3, pylith::topology::FieldBase::TENSOR },
      };
      const int numProperties = 3;
      const char* dbProperties[] = {
	"prop_one", "prop_two", "prop_three", "prop_four",
      };
      const int numDBProperties = 4;
      const Metadata::ParamDescription stateVars[] = {
	{ "var1", 1, pylith::topology::FieldBase::TENSOR },
	{ "var2", 3, pylith::topology::FieldBase::VECTOR },
      };
      const int numStateVars = 2;
      const char* dbStateVars[] = {
	"var_one", "var_two", "var_three",
      };
      const int numDBStateVars = 3;
    } // _TestMetadata
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Setup test data.
void
pylith::materials::TestMetadata::setUp(void)
{ // setUp
  _metadata = new Metadata(_TestMetadata::properties,
			   _TestMetadata::numProperties,
			   _TestMetadata::dbProperties,
			   _TestMetadata::numDBProperties,
			   _TestMetadata::stateVars,
			   _TestMetadata::numStateVars,
			   _TestMetadata::dbStateVars,
			   _TestMetadata::numDBStateVars);
  CPPUNIT_ASSERT(0 != _metadata);
} // setUp

// ----------------------------------------------------------------------
// Tear down test data.
void
pylith::materials::TestMetadata::tearDown(void)
{ // tearDown
  delete _metadata; _metadata = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::materials::TestMetadata::testConstructor(void)
{ // testConstructor
  CPPUNIT_ASSERT(0 != _metadata);
} // testConstructor

// ----------------------------------------------------------------------
// Test copy constructor.
void
pylith::materials::TestMetadata::testCopyConstructor(void)
{ // testCopyConstructor
  Metadata m(*_metadata);

  delete _metadata; _metadata = &m;
  testProperties();
  testStateVars();
  testDBProperties();
  testDBStateVars();

  _metadata = 0;
} // testCopyConstructor

// ----------------------------------------------------------------------
// Test properties().
void
pylith::materials::TestMetadata::testProperties(void)
{ // testProperties
  CPPUNIT_ASSERT(_metadata);

  const int numProperties = _TestMetadata::numProperties;
  CPPUNIT_ASSERT_EQUAL(numProperties, _metadata->numProperties());
  for (int i=0; i < numProperties; ++i) {
    const Metadata::ParamDescription& property = _metadata->getProperty(i);
    CPPUNIT_ASSERT_EQUAL(std::string(_TestMetadata::properties[i].name),
			 property.name);
    CPPUNIT_ASSERT_EQUAL(_TestMetadata::properties[i].fiberDim, 
			 property.fiberDim);
    CPPUNIT_ASSERT_EQUAL(_TestMetadata::properties[i].fieldType, 
			 property.fieldType);
  } // for
} // testProperties

// ----------------------------------------------------------------------
// Test stateVars().
void
pylith::materials::TestMetadata::testStateVars(void)
{ // testStateVars
  CPPUNIT_ASSERT(_metadata);

  const int numStateVars = _TestMetadata::numStateVars;
  CPPUNIT_ASSERT_EQUAL(numStateVars, _metadata->numStateVars());
  for (int i=0; i < numStateVars; ++i) {
    const Metadata::ParamDescription& stateVar = _metadata->getStateVar(i);
    CPPUNIT_ASSERT_EQUAL(std::string(_TestMetadata::stateVars[i].name),
			 stateVar.name);
    CPPUNIT_ASSERT_EQUAL(_TestMetadata::stateVars[i].fiberDim, 
			 stateVar.fiberDim);
    CPPUNIT_ASSERT_EQUAL(_TestMetadata::stateVars[i].fieldType, 
			 stateVar.fieldType);
  } // for
} // testStateVars

// ----------------------------------------------------------------------
// Test dbProperties().
void
pylith::materials::TestMetadata::testDBProperties(void)
{ // testDBProperties
  CPPUNIT_ASSERT(0 != _metadata);

  const int numDBProperties = _TestMetadata::numDBProperties;
  CPPUNIT_ASSERT_EQUAL(numDBProperties, _metadata->numDBProperties());

  const char* const* dbPropertiesE = _TestMetadata::dbProperties;
  const char* const* dbProperties = _metadata->dbProperties();
 
  for (int i=0; i < numDBProperties; ++i)
    CPPUNIT_ASSERT(0 == strcmp(dbPropertiesE[i], dbProperties[i]));
} // testDBProperties

// ----------------------------------------------------------------------
// Test dbStateVars().
void
pylith::materials::TestMetadata::testDBStateVars(void)
{ // testDBStateVars
  CPPUNIT_ASSERT(0 != _metadata);

  const int numDBStateVars = _TestMetadata::numDBStateVars;
  CPPUNIT_ASSERT_EQUAL(numDBStateVars, _metadata->numDBStateVars());

  const char* const* dbStateVarsE = _TestMetadata::dbStateVars;
  const char* const* dbStateVars = _metadata->dbStateVars();
 
  for (int i=0; i < numDBStateVars; ++i)
    CPPUNIT_ASSERT(0 == strcmp(dbStateVarsE[i], dbStateVars[i]));
} // testDBStateVars


// End of file 
