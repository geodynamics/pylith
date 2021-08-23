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
 * @file tests/libtests/faults/TestInterfacePatches.hh
 *
 * C++ unit testing for InterfacePatches.
 */

#if !defined(pylith_feassemble_testinterfacepatches_hh)
#define pylith_feassemble_testinterfacepatches_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // HASA Mesh
#include "pylith/testing/testingfwd.hh" // HASA FaultCohesiveStub

namespace pylith {
    namespace feassemble {
        class TestInterfacePatches;
        class TestInterfacePatches_Data;
    } // feassemble
} // pylith

class pylith::feassemble::TestInterfacePatches : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestInterfacePatches);

    CPPUNIT_TEST(testAccessors);
    CPPUNIT_TEST(testCreateMaterialPairs);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test getLabelName().
    void testAccessors(void);

    // Test createMaterialPairs().
    void testCreateMaterialPairs(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////
protected:

    TestInterfacePatches_Data* _data; ///< Data for testing.
    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.
    pylith::faults::FaultCohesiveStub* _fault; ///< Fault to create interface.

    // PRIVATE METHODS //////////////////////////////////////////////////////
private:

    /// Setup mesh.
    void _initialize();

}; // class TestInterfacePatches

class pylith::feassemble::TestInterfacePatches_Data {
    // PUBLIC METHODS //////////////////////////////////////////////////////////
public:

    /// Constructor
    TestInterfacePatches_Data(void);

    /// Destructor
    ~TestInterfacePatches_Data(void);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////
public:

    const char* filename; ///< Name of mesh file.
    const char* faultLabel; ///< Label for fault.
    const char* edgeLabel; ///< Label for buried fault edges.

    // Information for integration patches for material pairs.
    struct KeyValues {
        PetscInt negative_value;
        PetscInt positive_value;
    }; // KeyValues

    PylithInt numPatches; ///< Number of integration patches.
    KeyValues* patchKeys; ///< Weak form keys for integration patches.
    PylithInt* patchNumCells; ///< Number of cohesive cells in each integration patch.
    PylithInt** patchCells; ///< List of cohesive cells in each integration patch.

}; // TestInterfacePatches_Data

#endif // pylith_feassemble_testinterfacepatches_hh

// End of file
