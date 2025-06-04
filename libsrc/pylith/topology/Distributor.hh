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

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/meshio/meshiofwd.hh" // USES DataWriter
#include "pylith/faults/faultsfwd.hh" // USES FaultCohesive

// Distributor ----------------------------------------------------------
/// Distribute mesh among processors.
class pylith::topology::Distributor : public pylith::utils::GenericComponent {
    friend class TestDistributor; // unit testing

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /// Constructor
    Distributor(void);

    /// Destructor
    ~Distributor(void);

    /** Distribute mesh among processors.
     *
     * @param[in] newMesh Distributed mesh (result).
     * @param[in] origMesh Mesh to distribute.
     * @param[in] faults Array of fault interfaces.
     * @param[in] numFaults Number of fault interfaces.
     * @param[in] partitionerName Name of PETSc partitioner to use in distributing mesh.
     * @param[in] useEdgeWeighting Use edge weighting when partitioning (parmetis only).
     */
    static
    void distribute(pylith::topology::Mesh* const newMesh,
                    const pylith::topology::Mesh& origMesh,
                    pylith::faults::FaultCohesive* faults[],
                    const int numFaults,
                    const char* partitionerName,
                    const bool useEdgeWeighting);

    /** Write partitioning info for distributed mesh.
     *
     * @param writer Data writer for partition information.
     * @param mesh Distributed mesh.
     * @param cs Coordinate system for mesh.
     */
    static
    void write(meshio::DataWriter* const writer,
               const topology::Mesh& mesh);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    Distributor(const Distributor&); ///< Not implemented
    const Distributor& operator=(const Distributor&); ///< Not implemented

}; // Distributor

// End of file
