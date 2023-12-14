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

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

namespace pylith {
    namespace meshio {
        class TestExodusII;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestExodusII : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Test constructor
    void testConstructor(void);

    /// Test filename()
    void testFilename(void);

    /// Test open() and close().
    void testOpenClose(void);

    /// Test hasDim()
    void testHasDim(void);

    /// Test hasAtt()
    void testHasAtt(void);

    /// Test hasVar()
    void testHasVar(void);

    /// Test getVar(PylithScalar*)
    void testGetVarDouble(void);

    /// Test getVar(int*)
    void testGetVarInt(void);

    /// Test getVar(string_vector)
    void testGetVarString(void);

}; // class TestExodusII

// End of file
