// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file tests/libtests/bc/TestAbsorbingDampers.hh
 *
 * @brief C++ TestAbsorbingDampers object.
 *
 * C++ unit testing for AbsorbingDampers.
 */
#pragma once

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/bc/bcfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/topology/Field.hh" // HOLDSA Discretization

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA UserFunctionDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA CoordSys
#include "pylith/scales/scalesfwd.hh" // HOLDSA Scales

/// Namespace for pylith package
namespace pylith {
    namespace bc {
        class TestAbsorbingDampers;

        class TestAbsorbingDampers_Data;
    } // bc
} // pylith

/// C++ unit testing for DirichletBC.
class pylith::bc::TestAbsorbingDampers : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestAbsorbingDampers);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testAccessors);
    CPPUNIT_TEST(testAuxFieldDiscretization);
    CPPUNIT_TEST(testAuxFieldDB);
    CPPUNIT_TEST(testScales);
    CPPUNIT_TEST(testVerifyConfiguration);
    CPPUNIT_TEST(testInitialize);
    CPPUNIT_TEST(testComputeRHSResidual);
    CPPUNIT_TEST(testAuxFieldSetup);

    CPPUNIT_TEST_SUITE_END_ABSTRACT();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test constructor.
    void testConstructor(void);

    /// Test accessors (field).
    void testAccessors(void);

    /// Test auxFieldDiscretization().
    void testAuxFieldDiscretization(void);

    /// Test auxFieldDB().
    void testAuxFieldDB(void);

    /// Test scales().
    void testScales(void);

    /// Test verifyConfiguration().
    void testVerifyConfiguration(void);

    /// Test initialize().
    void testInitialize(void);

    /// Test computeRHSResidual().
    void testComputeRHSResidual(void);

    /// Test _auxiliaryFieldsSetup().
    void testAuxFieldSetup(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    TestAbsorbingDampers_Data* _data; ///< Data for testing

    pylith::bc::AbsorbingDampers* _bc; /// Test subject.
    pylith::topology::Mesh* _mesh; /// Mesh used in testing.
    pylith::topology::Field* _solution; ///< Solution field used in testing.

    static const double FILL_VALUE; ///< Fill value for unconstrained values.

    // PRIVATE METHODS ////////////////////////////////////////////////////
private:

    /// Initializer boundary condition for testing.
    void _initialize(void);

    /// Setup solution field.
    void _setupSolutionField(void);

}; // class TestAbsorbingDampers

// ======================================================================
class pylith::bc::TestAbsorbingDampers_Data {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    TestAbsorbingDampers_Data(void);

    /// Destructor
    ~TestAbsorbingDampers_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    const char* meshFilename; ///< Name of file with ASCII mesh.
    const char* bcLabel; ///< Label defining cells associated with material.

    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    pylith::scales::Scales* scales; ///< Scales for nondimensionalization.

    const char* field; ///< Name of solution field constrained.
    pylith::topology::FieldBase::VectorFieldEnum vectorFieldType; ///< Vector field type for constrained field.

    int numAuxSubfields; ///< Number of auxiliary subfields.
    const char** auxSubfields; ///< Names of auxiliary subfields.
    pylith::topology::Field::Discretization* auxDiscretizations; ///< Discretizations for auxiliary fields.
    spatialdata::spatialdb::UserFunctionDB* auxDB; ///< Spatial database with auxiliary field.

    PylithReal t; ///< Time associated with setting solution.
    int solnNumSubfields; ///< Number of solution subfields.
    pylith::topology::FieldBase::Discretization* solnDiscretizations; ///< Discretizations for solution fields.
    spatialdata::spatialdb::UserFunctionDB* solnDB; ///< Spatial database with solution.

}; // class TestAbsorbingDampers_Data

// End of file
