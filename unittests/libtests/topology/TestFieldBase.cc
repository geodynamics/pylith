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

#include "TestFieldBase.hh" // Implementation of class methods
#include "pylith/topology/FieldBase.hh" // USES Field


#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldBase );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldBase::testConstructor(void)
{ // testConstructor
  FieldBase field;
} // testConstructor

// ----------------------------------------------------------------------
// Test name().
void 
pylith::topology::TestFieldBase::testName(void)
{ // testName
  FieldBase field;

  CPPUNIT_ASSERT_EQUAL(std::string("unknown"), std::string(field.name()));

  const std::string& name = "field A";
  field.name(name.c_str());
  CPPUNIT_ASSERT_EQUAL(name, std::string(field.name()));
} // testName

// ----------------------------------------------------------------------
// Test vectorFieldType().
void
pylith::topology::TestFieldBase::testVectorFieldType(void)
{ // testVectorFieldType
  FieldBase field;

  CPPUNIT_ASSERT_EQUAL(FieldBase::OTHER, field.vectorFieldType());

  const FieldBase::VectorFieldEnum ftype = FieldBase::TENSOR;
  field.vectorFieldType(ftype);
  CPPUNIT_ASSERT_EQUAL(ftype, field.vectorFieldType());
} // testVectorFieldType

// ----------------------------------------------------------------------
// Test scale().
void
pylith::topology::TestFieldBase::testScale(void)
{ // testScale
  FieldBase field;

  CPPUNIT_ASSERT_EQUAL(1.0, field.scale());

  const double scale = 4.0;
  field.scale(scale);
  CPPUNIT_ASSERT_EQUAL(scale, field.scale());
} // testScale

// ----------------------------------------------------------------------
// Test addDimensionOkay().
void
pylith::topology::TestFieldBase::testAddDimensionOkay(void)
{ // testAddDimensionOkay
  FieldBase field;

  CPPUNIT_ASSERT_EQUAL(false, field.addDimensionOkay());

  field.addDimensionOkay(true);
  CPPUNIT_ASSERT_EQUAL(true, field.addDimensionOkay());
} // testAddDimensionOkay


// End of file 
