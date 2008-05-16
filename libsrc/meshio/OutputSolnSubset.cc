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

#include "OutputSolnSubset.hh" // implementation of class methods

#include <Selection.hh> // USES submesh algorithms

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnSubset::OutputSolnSubset(void) :
  _label("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnSubset::~OutputSolnSubset(void)
{ // destructor
} // destructor  

// ----------------------------------------------------------------------
// Set label identifier for subdomain.
void
pylith::meshio::OutputSolnSubset::label(const char* value)
{ // label
  _label = value;
} // label

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputSolnSubset::verifyConfiguration(const ALE::Obj<Mesh>& mesh) const
{ // verifyConfiguration
  assert(!mesh.isNull());

  if (!mesh->hasIntSection(_label)) {
    std::ostringstream msg;
    msg << "Mesh missing group of vertices '" << _label
	<< " for subdomain output.";
    throw std::runtime_error(msg.str());
  } // if
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get mesh associated with subdomain.
const ALE::Obj<pylith::Mesh>&
pylith::meshio::OutputSolnSubset::subdomainMesh(const ALE::Obj<Mesh>& mesh)
{ // subdomainMesh
  _mesh =
    ALE::Selection<Mesh>::submeshV(mesh, mesh->getIntSection(_label));
  if (_mesh.isNull()) {
    std::ostringstream msg;
    msg << "Could not construct mesh of subdomain " << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  _mesh->setRealSection("coordinates", 
			mesh->getRealSection("coordinates"));

  return _mesh;
} // subdomainMesh


// End of file 
