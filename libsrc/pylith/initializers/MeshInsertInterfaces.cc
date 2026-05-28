// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/initializers/MeshInsertInterfaces.hh" // implementation of class methods

#include "pylith/problems/Problem.hh" // USES Problem
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Distributor.hh" // USES Distributor
#include "pylith/materials/Material.hh" // USES Material
#include "pylith/faults/FaultCohesive.hh" // USES FaultCohesive
#include "pylith/utils/error.hh" // USES PylithCallPetsc()
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ------------------------------------------------------------------------------------------------
// Default constructor
pylith::initializers::MeshInsertInterfaces::MeshInsertInterfaces(void) {}


// ------------------------------------------------------------------------------------------------
// Default destructor
pylith::initializers::MeshInsertInterfaces::~MeshInsertInterfaces(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::initializers::MeshInsertInterfaces::deallocate(void) {
    InitializePhase::deallocate();
} // deallocate


// ------------------------------------------------------------------------------------------------
// Run initialization phase.
pylith::topology::Mesh*
pylith::initializers::MeshInsertInterfaces::run(pylith::topology::Mesh* mesh,
                                                const pylith::problems::Problem& problem) {
    PYLITH_METHOD_BEGIN;
    if (problem.getInterfaces().size() > 0) {
        PYLITH_INFO_ROOT(pylith::journal::application_flow, "Inserting cohesive cells.");
    } // if
    assert(mesh);

    if (!problem.getInterfaces().size()) {
        PetscDM dmOrig = mesh->getDM();assert(dmOrig);
        PylithCallPetsc(PetscObjectReference((PetscObject) dmOrig));
        pylith::topology::Mesh* meshNew = new pylith::topology::Mesh(dmOrig, *mesh);
        PYLITH_METHOD_RETURN(meshNew);
    } // if
    PylithCallPetsc(DMPlexCheckGeometry(mesh->getDM()));

    pythia::journal::debug_t debug("initialize_mesh");
    if (debug.state()) {
        mesh->view(":mesh_domain_before_faults.txt:ascii_info_detail");
        pylith::topology::Mesh* meshExploded = pylith::topology::MeshOps::explode(*mesh);
        meshExploded->view(":mesh_domain_before_faults.tex:ascii_latex");
        delete meshExploded;meshExploded = nullptr;
    } // if

    // Determine starting label value for cohesive cells.
    PylithInt cohesiveLabelValue = 100;
    for (auto material : problem.getMaterials()) {
        const PylithInt materialLabelValue = material->getLabelValue();
        cohesiveLabelValue = std::max(cohesiveLabelValue, materialLabelValue);
    } // for

    pylith::topology::Mesh* meshTmp = mesh;
    pylith::topology::Mesh* meshNew = nullptr;
    for (auto interface : problem.getInterfaces()) {
        interface->setCohesiveLabelValue(cohesiveLabelValue);
        meshNew = interface->transformTopology(meshTmp);
        if (meshTmp != mesh) {
            delete meshTmp;meshTmp = nullptr;
        } // if
        meshTmp = meshNew;
        cohesiveLabelValue += 1;
    } // for

    if (debug.state()) {
        DMPlexCheckTransform(meshNew->getDM());
        meshNew->view(":mesh_domain_before_overlap.txt:ascii_info_detail");
        pylith::topology::Mesh* meshExploded = pylith::topology::MeshOps::explode(*meshNew);
        meshExploded->view(":mesh_domain_before_overlap.tex:ascii_latex");
        delete meshExploded;meshExploded = nullptr;
    } // if

    PylithCallPetsc(DMPlexCheckSkeleton(meshNew->getDM(), 0));
    PylithCallPetsc(DMPlexCheckGeometry(meshNew->getDM()));

    PetscDM dmNew = nullptr;
    // Set overlap since cohesive cells can be put in the SF
    PylithCallPetsc(DMPlexSetOverlap(meshNew->getDM(), nullptr, 1));
    pylith::topology::Distributor::distributeOverlap(&dmNew, meshNew->getDM(), problem.getInterfaces());
    meshNew->setDM(dmNew);

    /* Need to reorder supports of cohesive cells after migration */
    DMPlexTransform tr;
    PylithCallPetsc(DMPlexTransformCreate(meshNew->getComm(), &tr));
    PylithCallPetsc(DMPlexTransformSetType(tr, DMPLEXCOHESIVEEXTRUDE));
    // PylithCallPetsc(DMPlexTransformSetUp(tr));
    PylithCallPetsc(DMPlexTransformOrderSupports(tr, dmNew, dmNew));
    PylithCallPetsc(DMPlexTransformDestroy(&tr));

    PylithCallPetsc(DMPlexCheckGeometry(meshNew->getDM()));

    PYLITH_METHOD_RETURN(meshNew);
} // run


// End of file
