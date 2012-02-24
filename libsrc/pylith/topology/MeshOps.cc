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
  int ncells = 0;

  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  if (!sieveMesh.isNull()) {
    const ALE::Obj<Mesh::SieveMesh::label_sequence>& cells = 
      sieveMesh->getLabelStratum("material-id", materialId);
    if (!cells.isNull())
      ncells = cells->size();
  } // if
} // numMaterialCells


// ----------------------------------------------------------------------
void
pylith::topology::MeshOps::checkMaterialIds(const Mesh& mesh,
					    int* const materialIds,
					    const int numMaterials)
{ // checkMaterialIds
  typedef Mesh::SieveMesh SieveMesh;

  assert( (0 == numMaterials && 0 == materialIds) ||
	  (0 < numMaterials && 0 != materialIds) );

  // Create map with indices for each material
  std::map<int, int> materialIndex;
  for (int i=0; i < numMaterials; ++i)
    materialIndex[materialIds[i]] = i;

  int_array matCellCounts(numMaterials);
  matCellCounts = 0;

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<SieveMesh::label_type>& materialsLabel = 
    sieveMesh->getLabel("material-id");

  //mesh.view("===== MESH BEFORE CHECKING MATERIALS =====");
  //materialsLabel->view("material-id");

  int* matBegin = materialIds;
  int* matEnd = materialIds + numMaterials;
  std::sort(matBegin, matEnd);

  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellId = sieveMesh->getValue(materialsLabel, *c_iter);
    const int* result = std::find(matBegin, matEnd, cellId);
    if (result == matEnd) {
      std::ostringstream msg;
      msg << "Material id '" << cellId << "' for cell '" << *c_iter
	  << "' does not match the id of any available materials or interfaces.";
      throw std::runtime_error(msg.str());
    } // if

    const int matIndex = materialIndex[cellId];
    assert(0 <= matIndex && matIndex < numMaterials);
    ++matCellCounts[matIndex];
  } // for

  // Make sure each material has 
  int_array matCellCountsAll(matCellCounts.size());
  MPI_Allreduce(&matCellCounts[0], &matCellCountsAll[0],
		matCellCounts.size(), MPI_INT, MPI_SUM, mesh.comm());
  for (int i=0; i < numMaterials; ++i) {
    const int matId = materialIds[i];
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
