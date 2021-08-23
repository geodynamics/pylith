// -*- C++ -*-
//
// -----------------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// -----------------------------------------------------------------------------
//

#include <portinfo>

#include "TestInterfacePatches.hh" // Implementation of class methods

// -----------------------------------------------------------------------------
namespace pylith {
    namespace feassemble {
        // ---------------------------------------------------------------------
        class TestInterfacePatches_Quad : public TestInterfacePatches {
            CPPUNIT_TEST_SUB_SUITE(TestInterfacePatches_Quad, TestInterfacePatches);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestInterfacePatches::setUp();

                _data->filename = "data/quad_patches_a.mesh";
                _data->faultLabel = "fault";
                _data->edgeLabel = NULL;

                static const size_t numPatches = 6;
                _data->numPatches = numPatches;

                static const TestInterfacePatches_Data::KeyValues patchKeys[numPatches] = {
                    {1, 1},
                    {1, 2},
                    {3, 2},
                    {3, 4},
                    {4, 4},
                    {4, 3},
                };
                _data->patchKeys = const_cast<TestInterfacePatches_Data::KeyValues*>(patchKeys);

                static const PylithInt patchNumCells[numPatches] = {
                    2, 1, 1, 1, 1, 1,
                };
                _data->patchNumCells = const_cast<PylithInt*>(patchNumCells);

                static const PylithInt patchCells11[2] = { 14, 15 };
                static const PylithInt patchCells12[1] = { 16 };
                static const PylithInt patchCells32[1] = { 17 };
                static const PylithInt patchCells34[1] = { 18 };
                static const PylithInt patchCells44[1] = { 19 };
                static const PylithInt patchCells43[1] = { 20 };
                static const PylithInt* patchCells[numPatches] = {
                    patchCells11,
                    patchCells12,
                    patchCells32,
                    patchCells34,
                    patchCells44,
                    patchCells43,
                };
                _data->patchCells = const_cast<PylithInt**>(patchCells);
            } // setUp

        }; // TestInterfacePatches_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestInterfacePatches_Quad);

    } // faults
} // pylith

// End of file
