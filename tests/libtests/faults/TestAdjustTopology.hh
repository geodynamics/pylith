// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/topology/topologyfwd.hh" // HOLDSA Mesh

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

/// Namespace for pylith package
namespace pylith {
    namespace faults {
        class TestAdjustTopology;
        class TestAdjustTopology_Data;
    } // faults
} // pylith

/// C++ unit testing of FaultCohesive::adjustTopology()
class pylith::faults::TestAdjustTopology : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestAdjustTopology(TestAdjustTopology_Data* data);

    /// Destructor.
    ~TestAdjustTopology(void);

    /// Test adjustTopology().
    void run(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    TestAdjustTopology_Data* _data; ///< Data for testing.
    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    /// Setup mesh.
    void _initialize();

}; // class TestAdjustTopology

class pylith::faults::TestAdjustTopology_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestAdjustTopology_Data(void);

    /// Destructor
    ~TestAdjustTopology_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    const char* filename; ///< Name of mesh file.

    size_t numFaults; ///< Number of faults
    const char** faultSurfaceLabels; ///< Labels marking fault surfaces.
    const char** faultEdgeLabels; ///< Labels for buried edges.
    const int* interfaceIds; ///< Label values for interfaces.

    size_t cellDim; ///< Number of dimensions associated with cells.
    size_t spaceDim; ///< Spatial dimension for vertex coordinates.
    size_t numVertices; // Number of vertices in updated mesh.
    size_t numCells; ///< Number of cells in updated mesh.

    const int* numCorners; ///< Number of vertices in each cell in updated mesh.
    const int* materialIds; ///< Material identifief for each cell in updated mesh.

    size_t numGroups; ///< Number of groups.
    const int* groupSizes; ///< Array of sizes in each group.
    char** groupNames; ///< Names of groups.
    char** groupTypes; ///< Type of group.

    bool failureExpected; ///< Flag indicating adjust topology should fail.

}; // TestAdjustTopology_Data

// End of file
