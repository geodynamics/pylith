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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "MeshOps.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/array.hh" // USES int_array

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

#include <algorithm> // USES std::sort, std::find
#include <map> // USES std::map


// ----------------------------------------------------------------------
int
pylith::topology::MeshOps::numMaterialCells(const Mesh& mesh,
					    int materialId)
{ // numMaterialCells
  PetscInt ncells = 0;

  DM dmMesh = mesh.dmMesh();
  assert(dmMesh);
  PetscErrorCode err = DMPlexGetStratumSize(dmMesh, "material-id", materialId, &ncells);CHECK_PETSC_ERROR(err);
  return ncells;
} // numMaterialCells


// ----------------------------------------------------------------------
void
pylith::topology::MeshOps::checkMaterialIds(const Mesh& mesh,
					    int* const materialIds,
					    const int numMaterials)
{ // checkMaterialIds
  assert((!numMaterials && !materialIds) || (numMaterials && materialIds));
  PetscErrorCode err;

  // Create map with indices for each material
  std::map<int, int> materialIndex;
  for (int i=0; i < numMaterials; ++i)
    materialIndex[materialIds[i]] = i;

  int_array matCellCounts(numMaterials);
  matCellCounts = 0;

  DM       dmMesh = mesh.dmMesh();
  PetscInt cStart, cEnd;

  assert(dmMesh);
  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  DMLabel materialsLabel;
  err = DMPlexGetLabel(dmMesh, "material-id", &materialsLabel);CHECK_PETSC_ERROR(err);
  assert(materialsLabel);

  int *matBegin = materialIds;
  int *matEnd   = materialIds + numMaterials;
  std::sort(matBegin, matEnd);

  for (PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt matId;

    err = DMLabelGetValue(materialsLabel, c, &matId);CHECK_PETSC_ERROR(err);
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

  // Make sure each material has 
  int_array matCellCountsAll(matCellCounts.size());
  err = MPI_Allreduce(&matCellCounts[0], &matCellCountsAll[0],
                      matCellCounts.size(), MPI_INT, MPI_SUM, mesh.comm());CHECK_PETSC_ERROR(err);
  for (int i=0; i < numMaterials; ++i) {
    const int matId    = materialIds[i];
    const int matIndex = materialIndex[matId];
    assert(0 <= matIndex && matIndex < numMaterials);
    if (matCellCountsAll[matIndex] <= 0) {
      std::ostringstream msg;
      msg << "No cells associated with material with id '" << matId
	  << "'.";
      throw std::runtime_error(msg.str());
    } // if
  } // for
  
} // checkMaterialIds


// End of file 
