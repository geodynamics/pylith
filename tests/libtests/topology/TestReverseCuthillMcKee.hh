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

/**
 * @file tests/libtests/topology/TestReverseCuthillMcKee.hh
 *
 * @brief C++ TestReverseCuthillMcKee object
 *
 * C++ unit testing for ReverseCuthillMcKee.
 */

#if !defined(pylith_topology_testreversecuthillmckee_hh)
#define pylith_topology_testreversecuthillmckee_hh

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

#endif // pylith_topology_testreversecuthillmckee_hh

// End of file
