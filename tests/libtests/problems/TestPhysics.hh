// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/problems/problemsfwd.hh" // HOLDSA Physics
#include "pylith/topology/topologyfwd.hh" // HOLDSA Mesh, Field

/// Namespace for pylith package
namespace pylith {
    namespace problems {
        class TestPhysics;
    } // problems
} // pylith

class pylith::problems::TestPhysics : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestPhysics);

    CPPUNIT_TEST(testSetNormalizer);
    CPPUNIT_TEST(testSetAuxiliaryFieldDB);
    CPPUNIT_TEST(testSetAuxiliarySubfieldDiscretization);
    CPPUNIT_TEST(testObservers);
    CPPUNIT_TEST(testGetKernelConstants);
    CPPUNIT_TEST(testVerifyConfiguration);
    CPPUNIT_TEST(testCreateIntegrator);
    CPPUNIT_TEST(testCreateConstraint);
    CPPUNIT_TEST(testCreateAuxiliaryField);
    CPPUNIT_TEST(testCreateDerivedField);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test setNormalizer().
    void testSetNormalizer(void);

    /// Test setAuxiliaryFieldDB().
    void testSetAuxiliaryFieldDB(void);

    /// Test setAuxiliarySubfieldDiscretization();
    void testSetAuxiliarySubfieldDiscretization(void);

    /// Test registerObserver(), removeObserver(), getObservers().
    void testObservers(void);

    /// Test getKernelConstants().
    void testGetKernelConstants(void);

    /// Test verifyConfiguration().
    void testVerifyConfiguration(void);

    /// Test createIntegrator().
    void testCreateIntegrator(void);

    /// Test createConstraint().
    void testCreateConstraint(void);

    /// Test createAuxiliaryField().
    void testCreateAuxiliaryField(void);

    /// Test createDerivedField().
    void testCreateDerivedField(void);

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::problems::Physics* _physics; ///< Test subject.
    pylith::topology::Mesh* _mesh; // Mesh for test subject.
    pylith::topology::Field* _solution; ///< Solution field for test subject.

}; // class TestPhysics

// End of file
