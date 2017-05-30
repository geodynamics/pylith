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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/bc/TestDirichletTimeDependent.hh
 *
 * @brief C++ TestDirichletTimeDependent object.
 *
 * C++ unit testing for DirichletBC.
 */

#if !defined(pylith_bc_testdirichlettimedependent_hh)
#define pylith_bc_testdirichlettimedependent_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/bc/bcfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/topology/Field.hh" // HOLDSA Discretization

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SimpleGridDB
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional


/// Namespace for pylith package
namespace pylith {
    namespace bc {
        class TestDirichletTimeDependent;

        class TestDirichletTimeDependent_Data;
    } // bc
} // pylith

/// C++ unit testing for DirichletBC.
class pylith::bc::TestDirichletTimeDependent : public CppUnit::TestFixture {

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestDirichletTimeDependent);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testAccessors);
    CPPUNIT_TEST(testAuxFieldDiscretization);
    CPPUNIT_TEST(testAuxFieldsDB);
    CPPUNIT_TEST(testNormalizer);
    CPPUNIT_TEST(testVerifyConfiguration);
    CPPUNIT_TEST(testInitialize);
    CPPUNIT_TEST(testPrestep);
    CPPUNIT_TEST(testSetSolution);
    CPPUNIT_TEST(testAuxFieldsSetup);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test constructor.
    void testConstructor(void);

    /// Test accessors (field, constrainedDOF, dbTimeHistory, useInitial, useRate, useTimeHistory).
    void testAccessors(void);

    /// Test auxFieldDiscretization().
    void testAuxFieldDiscretization(void);

    /// Test auxFieldsDB().
    void testAuxFieldsDB(void);

    /// Test normalizer().
    void testNormalizer(void);

    /// Test verifyConfiguration().
    void testVerifyConfiguration(void);

    /// Test initialize().
    void testInitialize(void);

    /// Test prestep().
    void testPrestep(void);

    /// Test setSolution().
    void testSetSolution(void);

    /// Test _auxFieldsSetup().
    void testAuxFieldsSetup(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    TestDirichletTimeDependent_Data* _data; ///< Data for testing

    pylith::bc::DirichletTimeDependent* _bc; /// Test subject.
    pylith::topology::Mesh* _mesh; /// Mesh used in testing.
    pylith::topology::Field* _solution; ///< Solution field used in testing.

    // PRIVATE METHODS ////////////////////////////////////////////////////
private:

    /// Initializer boundary condition for testing.
    void _initialize(void);

    /// Setup solution field.
    void _setupSolutionField(void);

}; // class TestDirichletTimeDependent


// ======================================================================
class pylith::bc::TestDirichletTimeDependent_Data {

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    TestDirichletTimeDependent_Data(void);

    /// Destructor
    ~TestDirichletTimeDependent_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    const char* meshFilename; ///< Name of file with ASCII mesh.
    const char* bcLabel; ///< Label defining cells associated with material.

    spatialdata::units::Nondimensional* normalizer; ///< Scales for nondimensionalization.

    const char* field; ///< Name of solution field constrained.
    pylith::topology::FieldBase::VectorFieldEnum vectorFieldType; ///< Vector field type for constrained field.
    PylithReal scale; ///< Scale of constrained field.
    int numConstrainedDOF; ///< Number of constrained DOF;
    const int* constrainedDOF; ///< Array of constrained DOF.

    bool useInitial; ///< True if using initial value.
    bool useRate; ///< True if using initial rate.
    bool useTimeHistory; ///< True if using time history;
    const char* thFilename; ///< Name of file with time history.

    PylithReal t; ///< Time associated with setting solution.
    int solnNumFields; ///< Number of solution fields.
    pylith::topology::FieldBase::Discretization* solnDiscretizations; ///< Discretizations for solution fields.
    const char* solnDBFilename; ///< Name of file with data for solution.

    int numAuxFields; ///< Number of auxiliary fields.
    const char** auxFields; ///< Names of auxiliary fields.
    pylith::topology::Field::Discretization* auxDiscretizations; ///< Discretizations for auxiliary fields.
    const char* auxDBFilename; ///< Name of file with data for auxFieldsDB.

}; // class TestDirichletTimeDependent_Data


#endif // pylith_bc_dirichlettimedependent_hh


// End of file
