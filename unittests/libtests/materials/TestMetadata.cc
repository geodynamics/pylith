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

#include "TestMetadata.hh" // Implementation of class methods

#include "pylith/materials/Metadata.hh" // USES Metadata
#include "pylith/utils/array.hh" // USES double_array

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
  testFiberDim();
  testFieldType();
  testDBProperties();
  testDBStateVars();

  _metadata = 0;
} // testCopyConstructor

// ----------------------------------------------------------------------
// Test properties().
void
pylith::materials::TestMetadata::testProperties(void)
{ // testProperties
  CPPUNIT_ASSERT(0 != _metadata);

  const string_vector& properties = _metadata->properties();
  const size_t numProperties = _TestMetadata::numProperties;
  CPPUNIT_ASSERT_EQUAL(numProperties, properties.size());
  for (size_t i=0; i < numProperties; ++i)
    CPPUNIT_ASSERT_EQUAL(std::string(_TestMetadata::properties[i].name),
			 properties[i]);
} // testProperties

// ----------------------------------------------------------------------
// Test stateVars().
void
pylith::materials::TestMetadata::testStateVars(void)
{ // testStateVars
  CPPUNIT_ASSERT(0 != _metadata);

  const string_vector& stateVars = _metadata->stateVars();
  const size_t numStateVars = _TestMetadata::numStateVars;
  CPPUNIT_ASSERT_EQUAL(numStateVars, stateVars.size());
  for (size_t i=0; i < numStateVars; ++i)
    CPPUNIT_ASSERT_EQUAL(std::string(_TestMetadata::stateVars[i].name),
			 stateVars[i]);
} // testStateVars

// ----------------------------------------------------------------------
// Test fiberDim().
void
pylith::materials::TestMetadata::testFiberDim(void)
{ // testFiberDim
  CPPUNIT_ASSERT(0 != _metadata);

  { // check property
  const int index = 1;
  const char* property = _TestMetadata::properties[index].name;
  const int fiberDimE = _TestMetadata::properties[index].fiberDim;
  const int fiberDim = _metadata->fiberDim(property, Metadata::PROPERTY);
  CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
  } // check property

  { // check state variable
  const int index = 1;
  const char* stateVar = _TestMetadata::stateVars[index].name;
  const int fiberDimE = _TestMetadata::stateVars[index].fiberDim;
  const int fiberDim = _metadata->fiberDim(stateVar, Metadata::STATEVAR);
  CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
  } // check state variable
} // testFiberDim

// ----------------------------------------------------------------------
// Test fieldType().
void
pylith::materials::TestMetadata::testFieldType(void)
{ // testFieldType
  CPPUNIT_ASSERT(0 != _metadata);

  { // check property
  const int index = 2;
  const char* property = _TestMetadata::properties[index].name;
  const topology::FieldBase::VectorFieldEnum fieldTypeE = 
    _TestMetadata::properties[index].fieldType;
  const topology::FieldBase::VectorFieldEnum fieldType =
    _metadata->fieldType(property, Metadata::PROPERTY);
  CPPUNIT_ASSERT_EQUAL(fieldTypeE, fieldType);
  } // check property

  { // check state variable
  const int index = 0;
  const char* stateVar = _TestMetadata::stateVars[index].name;
  const topology::FieldBase::VectorFieldEnum fieldTypeE = 
    _TestMetadata::stateVars[index].fieldType;
  const topology::FieldBase::VectorFieldEnum fieldType =
    _metadata->fieldType(stateVar, Metadata::STATEVAR);
  CPPUNIT_ASSERT_EQUAL(fieldTypeE, fieldType);
  } // check state variable
} // testFieldType

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
