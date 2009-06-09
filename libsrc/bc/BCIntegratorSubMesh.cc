// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "BCIntegratorSubMesh.hh" // implementation of object methods

#include "pylith/topology/SubMesh.hh" // USES SubMesh

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::BCIntegratorSubMesh::BCIntegratorSubMesh(void) :
  _boundaryMesh(0),
  _parameters(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::BCIntegratorSubMesh::~BCIntegratorSubMesh(void)
{ // destructor
  delete _boundaryMesh; _boundaryMesh = 0;
  delete _parameters; _parameters = 0;
} // destructor

// ----------------------------------------------------------------------
// Get submesh associated with boundary condition.
void
pylith::bc::BCIntegratorSubMesh::createSubMesh(const topology::Mesh& mesh)
{ // _createSubMesh
  delete _boundaryMesh; _boundaryMesh = 0;
  delete _parameters; _parameters = 0;

  _boundaryMesh = new topology::SubMesh(mesh, _label.c_str());
  assert(0 != _boundaryMesh);
} // _createSubMesh

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::BCIntegratorSubMesh::verifyConfiguration(const topology::Mesh& mesh) const 
{ // verifyConfiguration
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);

  BoundaryCondition::verifyConfiguration(mesh);

  // check compatibility of quadrature and boundary mesh
  if (_quadrature->cellDim() != _boundaryMesh->dimension()) {
    std::ostringstream msg;
    msg << "Quadrature is incompatible with cells for boundary condition '"
	<< _label << "'.\n"
	<< "Dimension of boundary mesh: " << _boundaryMesh->dimension()
	<< ", dimension of quadrature: " << _quadrature->cellDim()
	<< ".";
    throw std::runtime_error(msg.str());
  } // if
  const int numCorners = _quadrature->numBasis();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh>& sieveSubMesh = _boundaryMesh->sieveMesh();
  assert(!sieveSubMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    sieveSubMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Make sure surface cells are compatible with quadrature.
  const int boundaryDepth = sieveSubMesh->depth()-1; // depth of bndry cells
  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellNumCorners = 
      sieveSubMesh->getNumCellCorners(*c_iter, boundaryDepth);
    if (numCorners != cellNumCorners) {
      std::ostringstream msg;
      msg << "Quadrature is incompatible with cell for boundary condition '"
	  << _label << "'.\n"
	  << "Cell " << *c_iter << " has " << cellNumCorners
	  << " vertices but quadrature reference cell has "
	  << numCorners << " vertices.";
      throw std::runtime_error(msg.str());
    } // if
  } // for

} // verifyConfiguration


// End of file 
