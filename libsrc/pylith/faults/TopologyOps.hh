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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/TopologyOps.hh
 *
 * @brief C++ helper object for creating cohesive cells.
 */

#if !defined(pylith_faults_topologyops_hh)
#define pylith_faults_topologyops_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

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
     * @param groupdField Group of vertices assocated with faces of
     *   cells defining fault surface
     */
    static
    void createFault(topology::Mesh* faultMesh,
                     const topology::Mesh& mesh,
                     DMLabel groupField);

    /** Create cohesive cells in an interpolated mesh.
     *
     * If firstFaultVertex == 0, then firstFaultVertex is set to the first point
     * not currently used in the mesh, and firstFaultCell is incremented with this
     * point. These values are updated as new fault vertices and cells are added.
     *
     * @param fault Finite-element mesh of fault (output)
     * @param mesh Finite-element mesh
     * @param materialId Material id for cohesive elements.
     */
    static
    void create(topology::Mesh* mesh,
                const topology::Mesh& faultMesh,
                PetscDMLabel faultBdLabel,
                const int materialId);

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

#endif // pylith_faults_topologyops_hh

// End of file
