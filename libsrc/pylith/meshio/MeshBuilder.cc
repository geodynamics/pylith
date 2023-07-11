// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "MeshBuilder.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Set vertices and cells in mesh.
void
pylith::meshio::MeshBuilder::buildMesh(topology::Mesh* mesh,
                                       scalar_array* coordinates,
                                       const int numVertices,
                                       int spaceDim,
                                       const int_array& cells,
                                       const int numCells,
                                       const int numCorners,
                                       const int meshDim,
                                       const bool isParallel) { // buildMesh
    PYLITH_METHOD_BEGIN;

    assert(mesh);
    assert(coordinates);
    MPI_Comm comm = mesh->getComm();
    PetscInt dim = meshDim;
    PetscErrorCode err;

    { // Check to make sure every vertex is in at least one cell.
      // This is required by PETSc
        std::vector<bool> vertexInCell(numVertices, false);
        const int size = cells.size();
        for (int i = 0; i < size; ++i) {
            vertexInCell[cells[i]] = true;
        }
        int count = 0;
        for (int i = 0; i < numVertices; ++i) {
            if (!vertexInCell[i]) {
                ++count;
            }
        }
        if (count > 0) {
            std::ostringstream msg;
            msg << "Mesh contains " << count << " vertices that are not in any cells.";
            throw std::runtime_error(msg.str());
        } // if
    } // check

    /* DMPlex */
    PetscDM dmMesh = NULL;
    PetscBool interpolate = PETSC_TRUE;

    err = MPI_Bcast(&dim, 1, MPIU_INT, 0, comm);PYLITH_CHECK_ERROR(err);
    err = MPI_Bcast(&spaceDim, 1, MPIU_INT, 0, comm);PYLITH_CHECK_ERROR(err);
    const PetscInt bound = numCells*numCorners;
    for (PetscInt coff = 0; coff < bound; coff += numCorners) {
        DMPolytopeType ct;

        if (dim < 3) { continue;}
        switch (numCorners) {
        case 4: ct = DM_POLYTOPE_TETRAHEDRON;break;
        case 6: ct = DM_POLYTOPE_TRI_PRISM;break;
        case 8: ct = DM_POLYTOPE_HEXAHEDRON;break;
        default: continue;
        }
        err = DMPlexInvertCell(ct, (int *) &cells[coff]);PYLITH_CHECK_ERROR(err);
    }
    err = DMPlexCreateFromCellListPetsc(comm, dim, numCells, numVertices, numCorners, interpolate, &cells[0], spaceDim, &(*coordinates)[0], &dmMesh);PYLITH_CHECK_ERROR(err);
    mesh->setDM(dmMesh);

    PYLITH_METHOD_END;
} // buildMesh


// ----------------------------------------------------------------------
// Build a point group as an int section.
void
pylith::meshio::MeshBuilder::setGroup(pylith::topology::Mesh* mesh,
                                      const char* name,
                                      const GroupPtType groupType,
                                      const int_array& points) {
    PYLITH_METHOD_BEGIN;
    assert(mesh);

    PetscDM dmMesh = mesh->getDM();assert(dmMesh);
    const PetscInt numPoints = points.size();
    DMLabel label;
    PetscErrorCode err;

    err = DMCreateLabel(dmMesh, name);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmMesh, name, &label);PYLITH_CHECK_ERROR(err);
    if (CELL == groupType) {
        for (PetscInt p = 0; p < numPoints; ++p) {
            err = DMLabelSetValue(label, points[p], 1);PYLITH_CHECK_ERROR(err);
        } // for
    } else if (VERTEX == groupType) {
        PetscInt cStart, cEnd, vStart, vEnd, numCells;

        err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
        numCells = cEnd - cStart;
        for (PetscInt p = 0; p < numPoints; ++p) {
            err = DMLabelSetValue(label, numCells+points[p], 1);PYLITH_CHECK_ERROR(err);
        } // for
          // Also add any non-cells which have all vertices marked
        for (PetscInt p = 0; p < numPoints; ++p) {
            const PetscInt vertex = numCells+points[p];
            PetscInt      *star = NULL, starSize, s;

            err = DMPlexGetTransitiveClosure(dmMesh, vertex, PETSC_FALSE, &starSize, &star);PYLITH_CHECK_ERROR(err);
            for (s = 0; s < starSize*2; s += 2) {
                const PetscInt point = star[s];
                PetscInt      *closure = NULL, closureSize, c, value;
                PetscBool marked = PETSC_TRUE;

                if ((point >= cStart) && (point < cEnd)) { continue;}
                err = DMPlexGetTransitiveClosure(dmMesh, point, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
                for (c = 0; c < closureSize*2; c += 2) {
                    if ((closure[c] >= vStart) && (closure[c] < vEnd)) {
                        err = DMLabelGetValue(label, closure[c], &value);PYLITH_CHECK_ERROR(err);
                        if (value != 1) {marked = PETSC_FALSE;break;}
                    }
                }
                err = DMPlexRestoreTransitiveClosure(dmMesh, point, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
                if (marked) {err = DMLabelSetValue(label, point, 1);PYLITH_CHECK_ERROR(err);}
            }
            err = DMPlexRestoreTransitiveClosure(dmMesh, vertex, PETSC_FALSE, &starSize, &star);PYLITH_CHECK_ERROR(err);
        }
    } // if/else

    PYLITH_METHOD_END;
} // setGroup


// End of file
