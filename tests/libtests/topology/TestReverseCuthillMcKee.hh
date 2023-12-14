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

#include "pylith/topology/topologyfwd.hh" // USES Mesh

namespace pylith {
    namespace topology {
        class TestReverseCuthillMcKee;
        class TestReverseCuthillMcKee_Data;
    } // topology
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestReverseCuthillMcKee : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestReverseCuthillMcKee(TestReverseCuthillMcKee_Data* data);

    /// Destructor.
    ~TestReverseCuthillMcKee(void);

    /// Test reorder().
    void testReorder(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    TestReverseCuthillMcKee_Data* _data; ///< Data for testing.
    Mesh* _mesh; ///< Finite-element mesh.

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    /// Setup mesh.
    void _initialize();

}; // class TestReverseCuthillMcKee

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestReverseCuthillMcKee_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestReverseCuthillMcKee_Data(void);

    /// Destructor
    ~TestReverseCuthillMcKee_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    const char* filename; ///< Name of mesh file.
    const char* faultLabel; ///< Label for fault (use NULL for no fault).

};  // TestReverseCuthillMcKee_Data

// End of file
