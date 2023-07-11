// -*- C++ -*-
//
// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
//

#include <portinfo>

#include "TestInterfacePatches.hh" // Implementation of class methods

namespace pylith {
    namespace feassemble {
        class TestInterfacePatches_Quad;
    }
}

class pylith::feassemble::TestInterfacePatches_Quad {
public:

    // Data factory methods
    static TestInterfacePatches_Data* caseA(void);

}; // TestInterfacePatches_Quad

// ------------------------------------------------------------------------------------------------
pylith::feassemble::TestInterfacePatches_Data*
pylith::feassemble::TestInterfacePatches_Quad::caseA(void) {
    pylith::feassemble::TestInterfacePatches_Data* data = new pylith::feassemble::TestInterfacePatches_Data();

    data->filename = "data/quad_patches_a.mesh";
    data->faultLabel = "fault";
    data->edgeLabel = NULL;

    static const size_t numPatches = 6;
    data->numPatches = numPatches;

    static const TestInterfacePatches_Data::KeyValues patchKeys[numPatches] = {
        {1, 1},
        {1, 2},
        {3, 2},
        {3, 4},
        {4, 4},
        {4, 3},
    };
    data->patchKeys = const_cast<TestInterfacePatches_Data::KeyValues*>(patchKeys);

    static const PylithInt patchNumCells[numPatches] = {
        2, 1, 1, 1, 1, 1,
    };
    data->patchNumCells = const_cast<PylithInt*>(patchNumCells);

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
    data->patchCells = const_cast<PylithInt**>(patchCells);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
#include "catch2/catch_test_macros.hpp"

TEST_CASE("TestInterfacePatches::testAccessors", "[TestInterfacePatches]") {
    pylith::feassemble::TestInterfacePatches::testAccessors();
}
TEST_CASE("TestInterfacePatches_QuadA", "[TestInterfacePatches][Quad]") {
    pylith::feassemble::TestInterfacePatches(pylith::feassemble::TestInterfacePatches_Quad::caseA()).testCreateMaterialPairs();
}

// End of file
