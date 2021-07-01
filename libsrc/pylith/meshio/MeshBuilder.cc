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
// Copyright (c) 2010-2021 University of California, Davis
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
    PetscBool interpolate = PETSC_TRUE; /* interpolate = interpolate ? PETSC_TRUE : PETSC_FALSE; */

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


// End of file
