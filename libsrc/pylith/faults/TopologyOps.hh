// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
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

}; // class TopologyOps

#endif // pylith_faults_topologyops_hh

// End of file
