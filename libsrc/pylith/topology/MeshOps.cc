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

#include "MeshOps.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/utils/array.hh" // USES int_array

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

#include <algorithm> // USES std::sort, std::find
#include <map> // USES std::map


// ----------------------------------------------------------------------
// Create PETSc DM mesh.
void
pylith::topology::MeshOps::createDMMesh(Mesh* const mesh,
					const int dim,
					const MPI_Comm& comm,
					const char* label)
{ // createDMMesh
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err;
  PetscDM dmMesh = NULL;
  err = DMCreate(comm, &dmMesh);PYLITH_CHECK_ERROR(err);
  err = DMSetType(dmMesh, DMPLEX);PYLITH_CHECK_ERROR(err);
  err = DMSetDimension(dmMesh, dim);PYLITH_CHECK_ERROR(err);
  mesh->dmMesh(dmMesh, label);

  PYLITH_METHOD_END;
} // createDMMesh

// ----------------------------------------------------------------------
// Nondimensionalize the finite-element mesh.
void 
pylith::topology::MeshOps::nondimensionalize(Mesh* const mesh,
					     const spatialdata::units::Nondimensional& normalizer)
{ // nondimensionalize
  PYLITH_METHOD_BEGIN;

  PetscVec coordVec;
  const PylithScalar lengthScale = normalizer.lengthScale();
  PetscErrorCode err;

  PetscDM dmMesh = mesh->dmMesh();assert(dmMesh);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);PYLITH_CHECK_ERROR(err);assert(coordVec);
  // There does not seem to be an advantage to calling nondimensionalize()
  err = VecScale(coordVec, 1.0/lengthScale);PYLITH_CHECK_ERROR(err);
  err = DMPlexSetScale(dmMesh, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);
  err = DMViewFromOptions(dmMesh, NULL, "-pylith_nondim_dm_view");PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // nondimensionalize


// ----------------------------------------------------------------------
// Check topology of mesh.
void
pylith::topology::MeshOps::checkTopology(const Mesh& mesh)
{ // checkTopology
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);

  const int cellDim = mesh.dimension();
  const int numCorners = mesh.numCorners();
  PetscBool isSimplexMesh = PETSC_TRUE;
  if ((cellDim == 2 && numCorners == 4) ||
      (cellDim == 3 && numCorners == 8)) {
    isSimplexMesh = PETSC_FALSE;
  } // if
  DMLabel subpointMap;
  PetscErrorCode ierr = DMPlexGetSubpointMap(dmMesh, &subpointMap);PYLITH_CHECK_ERROR(ierr);
  PetscInt cellHeight = subpointMap ? 1 : 0;

  PetscErrorCode err;
  err = DMViewFromOptions(dmMesh, NULL, "-pylith_checktopo_dm_view");PYLITH_CHECK_ERROR(err);
  err = DMPlexCheckSymmetry(dmMesh);PYLITH_CHECK_ERROR_MSG(err, "Error in topology of mesh associated with symmetry of adjacency information.");

  err = DMPlexCheckSkeleton(dmMesh, isSimplexMesh, cellHeight);PYLITH_CHECK_ERROR_MSG(err, "Error in topology of mesh cells.");
} // checkTopology


// ----------------------------------------------------------------------
void
pylith::topology::MeshOps::checkMaterialIds(const Mesh& mesh,
					    int* const materialIds,
					    const int numMaterials)
{ // checkMaterialIds
  PYLITH_METHOD_BEGIN;

  assert((!numMaterials && !materialIds) || (numMaterials && materialIds));
  PetscErrorCode err;

  // Create map with indices for each material
  std::map<int, int> materialIndex;
  for (int i=0; i < numMaterials; ++i) {
    materialIndex[materialIds[i]] = i;
  } // for

  int_array matCellCounts(numMaterials);
  matCellCounts = 0;

  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  Stratum cellsStratum(dmMesh, Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  PetscDMLabel materialsLabel = NULL;
  err = DMPlexGetLabel(dmMesh, "material-id", &materialsLabel);PYLITH_CHECK_ERROR(err);assert(materialsLabel);

  int *matBegin = materialIds;
  int *matEnd = materialIds + numMaterials;
  std::sort(matBegin, matEnd);

  for (PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt matId;

    err = DMLabelGetValue(materialsLabel, c, &matId);PYLITH_CHECK_ERROR(err);
    if (matId < 0) {
      // :KLUDGE: Skip cells that are probably hybrid cells in halo
      // around fault that we currently ignore when looping over
      // materials (including cohesive cells).
      continue;
    } // if
    const int *result = std::find(matBegin, matEnd, matId);
    if (result == matEnd) {
      std::ostringstream msg;
      msg << "Material id '" << matId << "' for cell '" << c
          << "' does not match the id of any available materials or interfaces.";
      throw std::runtime_error(msg.str());
    } // if

    const int matIndex = materialIndex[matId];
    assert(0 <= matIndex && matIndex < numMaterials);
    ++matCellCounts[matIndex];
  } // for

  // Make sure each material has cells.
  int_array matCellCountsAll(matCellCounts.size());
  err = MPI_Allreduce(&matCellCounts[0], &matCellCountsAll[0],
                      matCellCounts.size(), MPI_INT, MPI_SUM, mesh.comm());PYLITH_CHECK_ERROR(err);
  for (int i=0; i < numMaterials; ++i) {
    const int matId = materialIds[i];
    const int matIndex = materialIndex[matId];
    assert(0 <= matIndex && matIndex < numMaterials);
    if (matCellCountsAll[matIndex] <= 0) {
      std::ostringstream msg;
      msg << "No cells associated with material with id '" << matId << "'.";
      throw std::runtime_error(msg.str());
    } // if
  } // for
  
  PYLITH_METHOD_END;
} // checkMaterialIds


// ----------------------------------------------------------------------
int
pylith::topology::MeshOps::numMaterialCells(const Mesh& mesh,
					    int materialId)
{ // numMaterialCells
  PYLITH_METHOD_BEGIN;

  PetscInt ncells = 0;

  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  PetscErrorCode err = DMPlexGetStratumSize(dmMesh, "material-id", materialId, &ncells);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_RETURN(ncells);
} // numMaterialCells


// End of file 
