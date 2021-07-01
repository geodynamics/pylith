// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFieldBase.hh" // Implementation of class methods

#include "pylith/topology/FieldBase.hh" // USES Field

#include <string> // USES std::string

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestFieldBase );

// ----------------------------------------------------------------------
// Test vectorFieldString()
void
pylith::topology::TestFieldBase::testVectorFieldString(void)
{ // testVectorFieldString
  CPPUNIT_ASSERT_EQUAL(std::string("scalar"),
		       std::string(FieldBase::vectorFieldString(FieldBase::SCALAR)));
  CPPUNIT_ASSERT_EQUAL(std::string("vector"),
		       std::string(FieldBase::vectorFieldString(FieldBase::VECTOR)));
  CPPUNIT_ASSERT_EQUAL(std::string("tensor"),
		       std::string(FieldBase::vectorFieldString(FieldBase::TENSOR)));
  CPPUNIT_ASSERT_EQUAL(std::string("other"),
		       std::string(FieldBase::vectorFieldString(FieldBase::OTHER)));
  CPPUNIT_ASSERT_EQUAL(std::string("multi_scalar"),
		       std::string(FieldBase::vectorFieldString(FieldBase::MULTI_SCALAR)));
  CPPUNIT_ASSERT_EQUAL(std::string("multi_vector"),
		       std::string(FieldBase::vectorFieldString(FieldBase::MULTI_VECTOR)));
  CPPUNIT_ASSERT_EQUAL(std::string("multi_tensor"),
		       std::string(FieldBase::vectorFieldString(FieldBase::MULTI_TENSOR)));
  CPPUNIT_ASSERT_EQUAL(std::string("multi_other"),
		       std::string(FieldBase::vectorFieldString(FieldBase::MULTI_OTHER)));
} // testVectorFieldString

// ----------------------------------------------------------------------
// Test parseVectorFieldString()
void
pylith::topology::TestFieldBase::testParseVectorFieldString(void)
{ // testParseVectorFieldString
  CPPUNIT_ASSERT_EQUAL(FieldBase::SCALAR,
		       FieldBase::parseVectorFieldString("scalar"));
  CPPUNIT_ASSERT_EQUAL(FieldBase::VECTOR,
		       FieldBase::parseVectorFieldString("vector"));
  CPPUNIT_ASSERT_EQUAL(FieldBase::TENSOR,
		       FieldBase::parseVectorFieldString("tensor"));
  CPPUNIT_ASSERT_EQUAL(FieldBase::OTHER,
		       FieldBase::parseVectorFieldString("other"));
  CPPUNIT_ASSERT_EQUAL(FieldBase::MULTI_SCALAR,
		       FieldBase::parseVectorFieldString("multi_scalar"));
  CPPUNIT_ASSERT_EQUAL(FieldBase::MULTI_VECTOR,
		       FieldBase::parseVectorFieldString("multi_vector"));
  CPPUNIT_ASSERT_EQUAL(FieldBase::MULTI_TENSOR,
		       FieldBase::parseVectorFieldString("multi_tensor"));
  CPPUNIT_ASSERT_EQUAL(FieldBase::MULTI_OTHER,
		       FieldBase::parseVectorFieldString("multi_other"));
} // testParseVectorFieldString


// End of file 
