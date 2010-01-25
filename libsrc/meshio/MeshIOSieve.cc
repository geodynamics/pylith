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

#include "MeshIOSieve.hh" // implementation of class methods

#include "petscmesh.hh"

#include "MeshBuilder.hh" // USES MeshBuilder
#include "pylith/topology/Mesh.hh" // USES Mesh

#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOSieve::MeshIOSieve(void) :
  _filename("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOSieve::~MeshIOSieve(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::MeshIOSieve::deallocate(void)
{ // deallocate
  MeshIO::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Read mesh.
void
pylith::meshio::MeshIOSieve::_read(void)
{ // _read
  MPI_Comm comm = _mesh->comm();
  int rank = 0;
  int meshDim = 0;
  int spaceDim = 0;
  int numVertices = 0;
  int numCells = 0;
  int numCorners = 0;


  // :TODO: STUFF GOES HERE
  assert(false);
} // read

// ----------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIOSieve::_write(void) const
{ // write

  ALE::MeshSerializer::writeMesh(_filename, *_mesh->sieveMesh());

} // write

  
// End of file 
