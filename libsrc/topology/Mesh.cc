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

#include "Mesh.hh" // implementation of class methods

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

// ----------------------------------------------------------------------
// Default constructor
pylith::topology::Mesh::Mesh(void) :
  _coordsys(0),
  _comm(PETSC_COMM_WORLD),
  _debug(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default constructor
pylith::topology::Mesh::Mesh(const int dim,
			     const MPI_Comm& comm) :
  _mesh(new SieveMesh(comm, dim)),
  _coordsys(0),
  _comm(comm),
  _debug(false)
{ // constructor
  assert(!_mesh->getFactory().isNull());
  _mesh->getFactory()->clear();
} // constructor

// ----------------------------------------------------------------------
// Default destructor
pylith::topology::Mesh::~Mesh(void)
{ // destructor
  delete _coordsys; _coordsys = 0;
} // destructor

// ----------------------------------------------------------------------
// Create Sieve mesh.
void
pylith::topology::Mesh::createSieveMesh(const int dim)
{ // createSieveMesh
  _mesh.destroy();
  _mesh = new SieveMesh(_comm, dim);
  _mesh->setDebug(_debug);
  assert(!_mesh->getFactory().isNull());
  _mesh->getFactory()->clear();
} // createSieveMesh

// ----------------------------------------------------------------------
// Set coordinate system.
void
pylith::topology::Mesh::coordsys(const spatialdata::geocoords::CoordSys* cs)
{ // coordsys
  delete _coordsys; _coordsys = (0 != cs) ? cs->clone() : 0;
  if (0 != _coordsys)
    _coordsys->initialize();
} // coordsys

// ----------------------------------------------------------------------
// Initialize the finite-element mesh.
void 
pylith::topology::Mesh::initialize(void)
{ // initialize
} // initialize


// End of file 
