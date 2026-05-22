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

    /** Complete labels and change value to match dimension of point.
     */
    static
    void updateCohesiveLabel(const pylith::topology::Mesh* mesh,
                             const char* labelName,
                             const int labelValue);

    /** Set label and label value for newly created cohesive cells.
     *
     * @param[inout] dmMeshNew PETSc DM with cohesive cells.
     * @param[in] dmMesh PETSc DM without cohesive cells.
     * @param[in] labelName Name of label for cohesive cells.
     * @param[in] labelValue Label value for cohesive cells.
     */
    static
    void setCohesiveCellLabel(PetscDM dmMeshNew,
                              PetscDM dmMesh,
                              const char* labelName,
                              const int labelValue);

    /** Remove cohesive points from labels.
     *
     * @param[inout] dmMeshNew PETSc DM with cohesive cells.
     */
    static
    void labelsRemoveCohesivePoints(PetscDM dmMeshNew);

    /** Create buried edge label.
     *
     * @param[inout] dmMeshNew PETSc DM with cohesive cells.
     * @param[in] dmMesh PETSc DM without cohesive cells.
     * @param[in] buriedEdgeLabelName Name of label for buried edge.
     * @param[in] buriedEdgeLabelValue Label value for buried edge.
     * @param[in] surfaceLabel Label identifying fault surface.
     * @param[in] transform Transform creating cohesive cells.
     */
    static
    void createBuriedEdgeLabel(PetscDM dmMeshNew,
                               PetscDM dmMesh,
                               const char* buriedEdgeLabelName,
                               const PylithInt buriedEdgeLabelValue,
                               PetscDMLabel surfaceLabel,
                               DMPlexTransform transform);

    /** Create (distributed) fault mesh from cohesive cells.
     *
     * @param faultMesh Finite-element mesh of fault (output).
     * @param mesh Finite-element mesh.
     * @param labelValue Value of label associated with integration domain.
     * @param labelName Name of label associated with integration domain.
     * @param surfaceLabel Name of label for interface surface.
     */
    static
    void createFaultFromCohesiveCells(pylith::topology::Mesh* faultMesh,
                                      const topology::Mesh& mesh,
                                      const char* labelName,
                                      const int labelValue,
                                      const char* surfaceLabel);

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
