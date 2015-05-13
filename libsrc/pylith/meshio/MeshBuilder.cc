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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "MeshBuilder.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

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
				       const bool interpolate,
				       const bool isParallel)
{ // buildMesh
  PYLITH_METHOD_BEGIN;

  assert(mesh);
  assert(coordinates);
  MPI_Comm comm  = mesh->comm();
  PetscInt dim  = meshDim;
  PetscErrorCode err;

  { // Check to make sure every vertex is in at least one cell.
    // This is required by PETSc
    std::vector<bool> vertexInCell(numVertices, false);
    const int size = cells.size();
    for (int i=0; i < size; ++i)
      vertexInCell[cells[i]] = true;
    int count = 0;
    for (int i=0; i < numVertices; ++i)
      if (!vertexInCell[i])
        ++count;
    if (count > 0) {
      std::ostringstream msg;
      msg << "Mesh contains " << count << " vertices that are not in any cells.";
      throw std::runtime_error(msg.str());
    } // if
  } // check

  /* DMPlex */
  PetscDM   dmMesh;
  PetscBool pInterpolate = PETSC_TRUE; /* pInterpolate = interpolate ? PETSC_TRUE : PETSC_FALSE; */
  PetscInt  bound        = numCells*numCorners, coff;

  err = MPI_Bcast(&dim, 1, MPIU_INT, 0, comm);PYLITH_CHECK_ERROR(err);
  err = MPI_Bcast(&spaceDim, 1, MPIU_INT, 0, comm);PYLITH_CHECK_ERROR(err);
  for (coff = 0; coff < bound; coff += numCorners) {
    err = DMPlexInvertCell(dim, numCorners, (int *) &cells[coff]);PYLITH_CHECK_ERROR(err);
  }
  err = DMPlexCreateFromCellList(comm, dim, numCells, numVertices, numCorners, pInterpolate, &cells[0], spaceDim, &(*coordinates)[0], &dmMesh);PYLITH_CHECK_ERROR(err);
  mesh->dmMesh(dmMesh);

  PYLITH_METHOD_END;
} // buildMesh

// End of file 
