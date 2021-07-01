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

/**
 * @file tests/libtests/problems/TestSolutionFactory.hh
 *
 * @brief C++ TestSolutionFactory object.
 *
 * C++ unit testing for SolutionFactory.
 */

#if !defined(pylith_problems_testsolutionfactory_hh)
#define pylith_problems_testsolutionfactory_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/problems/problemsfwd.hh" // HOLDSA SolutionFactory
#include "pylith/topology/Field.hh" // HOLDSA Field::SubfieldInfo
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA Coordsys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional

#include <map> // USES std::map

/// Namespace for pylith package
namespace pylith {
    namespace problems {
        class TestSolutionFactory;
        class TestSolutionFactory_Data;
    } // problems
} // pylith

class pylith::problems::TestSolutionFactory : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestSolutionFactory);

    CPPUNIT_TEST(testDispVel);
    CPPUNIT_TEST(testDispLagrangeFault);
    CPPUNIT_TEST(testPressure);
    CPPUNIT_TEST(testDispTemp);
    CPPUNIT_TEST(testSetValues);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test adding displacement and velocity subfields.
    void testDispVel(void);

    /// Test adding displacement and fault Lagrange multiplier subfields.
    void testDispLagrangeFault(void);

    /// Test adding pressure and fluid pressure subfields.
    void testPressure(void);

    /// Test adding displacement and temperature subfields.
    void testDispTemp(void);

    /// Test setValues().
    void testSetValues(void);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Initialze mesh, coordinate system, solution, and factory.
    void _initialize(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    SolutionFactory* _factory; ///< Test subject.
    TestSolutionFactory_Data* _data; ///< Test data.

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.
    pylith::topology::Field* _solution; ///< Solution field for test subject.

}; // class TestSolutionFactory

// =====================================================================================================================
class pylith::problems::TestSolutionFactory_Data {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestSolutionFactory_Data(void);

    /// Destructor
    ~TestSolutionFactory_Data(void);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    size_t dimension; ///< Spatial dimension.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    spatialdata::units::Nondimensional* normalizer; ///< Scales for nondimensionalization.

    std::map<std::string, pylith::topology::Field::SubfieldInfo> subfields;
    spatialdata::spatialdb::UserFunctionDB* solutionDB; ///< Spatial database with values for solution.

}; // class TestSolutionFactory_Data

#endif // pylith_problems_testsolutionfactory_hh

// End of file
