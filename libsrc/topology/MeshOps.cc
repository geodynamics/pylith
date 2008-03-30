// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "MeshOps.hh" // implementation of class methods

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <assert.h> // USES assert()

#include <algorithm> // USES std::sort, std::find

// ----------------------------------------------------------------------
void
pylith::topology::MeshOps::checkMaterialIds(const ALE::Obj<Mesh>& mesh,
					    int* materialIds,
					    const int numMaterials)
{ // checkMaterialIds
  assert( (0 == numMaterials && 0 == materialIds) ||
	  (0 < numMaterials && 0 != materialIds) );

  const ALE::Obj<ALE::Mesh::label_sequence>& cells = mesh->heightStratum(0);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<Mesh::label_type>& materialsLabel = mesh->getLabel("material-id");

  int* matBegin = materialIds;
  int* matEnd = materialIds + numMaterials;
  std::sort(matBegin, matEnd);

  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellId = mesh->getValue(materialsLabel, *c_iter);
    const int* result = std::find(matBegin, matEnd, cellId);
    if (result == matEnd) {
      std::ostringstream msg;
      msg << "Material id '" << cellId << "' for cell '" << *c_iter
	  << "' does not match the id of any available materials or interfaces.";
      throw std::runtime_error(msg.str());
    } // if
  } // for
} // checkMaterialIds


// End of file 
