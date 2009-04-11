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

#include "pylith/topology/Mesh.hh" // USES Mesh

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnSubset::OutputSolnSubset(void) :
  _label(""),
  _submesh(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnSubset::~OutputSolnSubset(void)
{ // destructor
  delete _submesh; _submesh = 0;
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
pylith::meshio::OutputSolnSubset::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  const ALE::Obj<topology::Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  if (!sieveMesh->hasIntSection(_label)) {
    std::ostringstream msg;
    msg << "Mesh missing group of vertices '" << _label
	<< " for subdomain output.";
    throw std::runtime_error(msg.str());
  } // if
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get mesh associated with subdomain.
const pylith::topology::SubMesh&
pylith::meshio::OutputSolnSubset::subdomainMesh(const topology::Mesh& mesh)
{ // subdomainMesh
  delete _submesh; _submesh = new topology::SubMesh(mesh, _label.c_str());
  assert(0 != _submesh);
  return *_submesh;
} // subdomainMesh

// ----------------------------------------------------------------------
// Append finite-element vertex field to file.
void
pylith::meshio::OutputSolnSubset::appendVertexField(const double t,
						    const topology::Field<topology::Mesh>& field)
{ // appendVertexField
  // How do we call DataWriter<SubMesh>::appendVertexField() with a Field<Mesh>?
} // appendVertexField


// End of file 
