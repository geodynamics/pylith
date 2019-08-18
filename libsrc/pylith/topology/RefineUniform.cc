// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "RefineUniform.hh" // implementation of class methods

#include "Mesh.hh" // USES Mesh
#include "MeshOps.hh" // USES MeshOps

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::topology::RefineUniform::RefineUniform(void) {}


// ----------------------------------------------------------------------
// Destructor
pylith::topology::RefineUniform::~RefineUniform(void) {
    deallocate();
}


// ----------------------------------------------------------------------
// Deallocate data structures.
void
pylith::topology::RefineUniform::deallocate(void) {}


// ----------------------------------------------------------------------
// Refine mesh.
void
pylith::topology::RefineUniform::refine(Mesh* const newMesh,
                                        const Mesh& mesh,
                                        const int levels) {
    PYLITH_METHOD_BEGIN;

    if (levels < 1) {
        PYLITH_METHOD_END;
    } // if

    assert(newMesh);

    PetscErrorCode err;
    PetscDM dmOrig = mesh.dmMesh();assert(dmOrig);

    PetscInt meshDepth = 0;
    err = DMPlexGetDepth(dmOrig, &meshDepth);

    const int meshDim = mesh.dimension();
    if (( meshDim > 0) && ( meshDepth != meshDim) ) {
        std::ostringstream msg;
        msg << "Mesh refinement for uninterpolated meshes not supported.\n"
            << "Turn on interpolated meshes using 'interpolate' mesh generator property.";
        throw std::runtime_error(msg.str());
    } // if

    // Refine, keeping original mesh intact.
    PetscDM dmNew = NULL;
    err = DMPlexSetRefinementUniform(dmOrig, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
    err = DMRefine(dmOrig, mesh.comm(), &dmNew);PYLITH_CHECK_ERROR(err);

    for (int i = 1; i < levels; ++i) {
        PetscDM dmCur = dmNew;dmNew = NULL;
        err = DMPlexSetRefinementUniform(dmCur, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
        err = DMRefine(dmCur, mesh.comm(), &dmNew);PYLITH_CHECK_ERROR(err);

        err = DMDestroy(&dmCur);PYLITH_CHECK_ERROR(err);
    } // for

    newMesh->dmMesh(dmNew);

    // Remove all non-cells from material-id label
    DMLabel mlabel;
    PetscIS vIS;
    const PetscInt *values;
    PetscInt cStart, cEnd, nv;
    err = DMPlexGetHeightStratum(dmNew, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmNew, "material-id", &mlabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetNumValues(mlabel, &nv);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetValueIS(mlabel, &vIS);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(vIS, &values);PYLITH_CHECK_ERROR(err);
    for (PetscInt v = 0; v < nv; ++v) {
      PetscIS sIS;
      const PetscInt *points;
      const PetscInt value = values[v];
      PetscInt np;
      err = DMLabelGetStratumSize(mlabel, value, &np);PYLITH_CHECK_ERROR(err);
      err = DMLabelGetStratumIS(mlabel, value, &sIS);PYLITH_CHECK_ERROR(err);
      err = ISGetIndices(sIS, &points);PYLITH_CHECK_ERROR(err);
      for (PetscInt p = 0; p < np; ++p) {
        const PetscInt point = points[p];
        if (point < cStart || point >= cEnd) {
          err = DMLabelClearValue(mlabel, point, value);PYLITH_CHECK_ERROR(err);
        }
      }
      err = ISRestoreIndices(sIS, &points);PYLITH_CHECK_ERROR(err);
      err = ISDestroy(&sIS);PYLITH_CHECK_ERROR(err);
    }
    err = ISRestoreIndices(vIS, &values);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&vIS);PYLITH_CHECK_ERROR(err);

    // Check consistency
    topology::MeshOps::checkTopology(*newMesh);

    // newMesh->view("REFINED_MESH", "::ascii_info_detail");

    PYLITH_METHOD_END;
} // refine


// End of file
