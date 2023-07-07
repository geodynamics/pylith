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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/FieldBase.hh" // USES Field

#include <string> // USES std::string

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace topology {
        class TestFieldBase;
    } // topology
} // pylith

class pylith::topology::TestFieldBase : public pylith::utils::GenericComponent {
    // PUBLIC METHODS ///////////////////////////////////////////////////////////////////////////////
public:

    /// Test vectorFieldString().
    static
    void testVectorFieldString(void);

    /// Test parseVectorFieldString().
    static
    void testParseVectorFieldString(void);

}; // class TestFieldBase

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestFieldBase::testVectorFieldString", "[TestFieldBase]") {
    pylith::topology::TestFieldBase::testVectorFieldString();
}
TEST_CASE("TestFieldBase::testParVectorFieldString", "[TestFieldBase]") {
    pylith::topology::TestFieldBase::testParseVectorFieldString();
}

// ------------------------------------------------------------------------------------------------
// Test vectorFieldString()
void
pylith::topology::TestFieldBase::testVectorFieldString(void) {
    CHECK(std::string("scalar") == std::string(FieldBase::vectorFieldString(FieldBase::SCALAR)));
    CHECK(std::string("vector") == std::string(FieldBase::vectorFieldString(FieldBase::VECTOR)));
    CHECK(std::string("tensor") == std::string(FieldBase::vectorFieldString(FieldBase::TENSOR)));
    CHECK(std::string("other") == std::string(FieldBase::vectorFieldString(FieldBase::OTHER)));
    CHECK(std::string("multi_scalar") == std::string(FieldBase::vectorFieldString(FieldBase::MULTI_SCALAR)));
    CHECK(std::string("multi_vector") == std::string(FieldBase::vectorFieldString(FieldBase::MULTI_VECTOR)));
    CHECK(std::string("multi_tensor") == std::string(FieldBase::vectorFieldString(FieldBase::MULTI_TENSOR)));
    CHECK(std::string("multi_other") == std::string(FieldBase::vectorFieldString(FieldBase::MULTI_OTHER)));
} // testVectorFieldString


// ----------------------------------------------------------------------
// Test parseVectorFieldString()
void
pylith::topology::TestFieldBase::testParseVectorFieldString(void) {
    CHECK(FieldBase::SCALAR == FieldBase::parseVectorFieldString("scalar"));
    CHECK(FieldBase::VECTOR == FieldBase::parseVectorFieldString("vector"));
    CHECK(FieldBase::TENSOR == FieldBase::parseVectorFieldString("tensor"));
    CHECK(FieldBase::OTHER == FieldBase::parseVectorFieldString("other"));
    CHECK(FieldBase::MULTI_SCALAR == FieldBase::parseVectorFieldString("multi_scalar"));
    CHECK(FieldBase::MULTI_VECTOR == FieldBase::parseVectorFieldString("multi_vector"));
    CHECK(FieldBase::MULTI_TENSOR == FieldBase::parseVectorFieldString("multi_tensor"));
    CHECK(FieldBase::MULTI_OTHER == FieldBase::parseVectorFieldString("multi_other"));
} // testParseVectorFieldString


// End of file
