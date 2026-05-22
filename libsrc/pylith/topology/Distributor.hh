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

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/meshio/meshiofwd.hh" // USES DataWriter
#include "pylith/faults/faultsfwd.hh" // USES FaultCohesive
#include "pylith/utils/petscfwd.h" // USES PetscDM

#include <vector> // USES std::vector

class pylith::topology::Distributor : public pylith::utils::PyreComponent {
    friend class TestDistributor; // unit testing

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    Distributor(void);

    /// Destructor
    ~Distributor(void);

    /** Set data writer.
     *
     * @param[in] writer Data writer.
     */
    void setDataWriter(pylith::meshio::DataWriter* writer);

    /** Set edge weighting.
     *
     * @param[in] useEdgeWeighting If true, weight edges in distribution.
     */
    void setUseEdgeWeighting(const bool flag);

    /** Set partitioner.
     *
     * @param[in] partitioner Name of mesh partitioner.
     */
    void setPartitioner(const char* partitioner);

    /** Distribute mesh.
     *
     * @param[in] mesh Original mesh.
     * @param[in] faults Fault interfaces.
     */
    pylith::topology::Mesh* distribute(const pylith::topology::Mesh& mesh,
                                       const std::vector<pylith::faults::FaultCohesive*>& faults) const;

    /** Distribute custom overlap based on PETSc labels.
     *
     * The overlap excludes cohesive cells but includes cells adjacent to faults.
     * This is a custom version of DMPlexDistributeOverlap()
     *
     * @param[out] dmOverlap PETSc DM for the overlap.
     * @param[in] dmMesh PETSc DM for the current mesh.
     * @param[in] faults Array of fault interfaces.
     *
     * @returns PETSc error code (0==success).
     */
    static
    PetscErrorCode distributeOverlap(PetscDM* dmOverlap,
                                     PetscDM dmMesh,
                                     const std::vector<pylith::faults::FaultCohesive*>& faults);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::meshio::DataWriter* _writer; ///< Data writer for partition.
    std::string _partitioner; ///< Name of mesh partitioner
    bool _useEdgeWeighting; ///< Use edge weighting in partitioning.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    Distributor(const Distributor&); ///< Not implemented
    const Distributor& operator=(const Distributor&); ///< Not implemented

}; // Distributor

// End of file
