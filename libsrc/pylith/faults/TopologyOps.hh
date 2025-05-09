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

#include "pylith/faults/faultsfwd.hh" // forward declarations

#include "pylith/topology/Mesh.hh" // USES Mesh
#include <set> // USES std::set

// TopologyOps ----------------------------------------------------------
/// Helper object for creation of cohesive cells.
class pylith::faults::TopologyOps {
    // PUBLIC TYPEDEFS ////////////////////////////////////////////////////
public:

    typedef std::set < PetscInt > PointSet;

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /** Create the fault mesh.
     *
     * @param faultMesh Finite-element mesh of fault (output).
     * @param mesh Finite-element mesh of domain.
     * @param surfaceLabel Label for points on fault surface.
     * @param surfaceLabelValue Value for label for points on fault surface.
     */
    static
    void createFault(topology::Mesh* faultMesh,
                     const topology::Mesh& mesh,
                     PetscDMLabel surfaceLabel,
                     const int surfaceLabelValue);

    /** Create cohesive cells in an interpolated mesh.
     *
     * @param fault Finite-element mesh of fault (output)
     * @param mesh Finite-element mesh
     * @param materialId Material id for cohesive elements.
     */
    static
    void create(topology::Mesh* mesh,
                const topology::Mesh& faultMesh,
                PetscDMLabel faultBdLabel,
                const int faultBdLabelValue,
                const int cohesiveLabelValue);

    /** Create (distributed) fault mesh from cohesive cells.
     *
     * @param faultMesh Finite-element mesh of fault (output).
     * @param mesh Finite-element mesh.
     * @param labelValue Value of label associated with integration domain.
     * @param labelName Name of label associated with integration domain.
     * @param surfaceLabel Name of label for interface surface.
     */
    static
    void createFaultParallel(topology::Mesh* faultMesh,
                             const topology::Mesh& mesh,
                             const int labelValue,
                             const char* labelName,
                             const char* surfaceLabel);

    /** Classify cells adjacent to the fault as to the side of the fault each cell is on.
     */
    static
    void classifyCellsDM(PetscDM dmMesh,
                         PetscInt vertex,
                         const int depth,
                         const int faceSize,
                         PetscInt firstCohesiveCell,
                         PointSet& replaceCells,
                         PointSet& noReplaceCells,
                         const int debug);

    /** Get name of PETSc DM label for interfaces.
     *
     * @returns PETSc Label name.
     */
    static
    const char* getInterfacesLabelName(void);

    /** Get PETSc DM label for interfaces, creating if necessary.
     *
     * @param[inout] dm PETSc DM holding interfaces label.
     * @returns PETSc DM label for interfaces.
     */
    static
    PetscDMLabel getInterfacesLabel(PetscDM dm);

    /** Get cells adjacent to cohesive cell on negative and positive sides of the fault.
     *
     * @param[out] adjacentCellNegative Adjacent cell on negative side of the fault.
     * @param[out] adjacentCellPositive Adjacent cell on positive side of the fault.
     * @param[in] dmMesh DM for finite-element mesh.
     * @param[in] cohesiveCell Cohesive cell.
     */
    static
    void getAdjacentCells(PylithInt* adjacentCellNegative,
                          PylithInt* adjacentCellPositive,
                          PetscDM dmMesh,
                          const PylithInt cohesiveCell);

}; // class TopologyOps

// End of file
