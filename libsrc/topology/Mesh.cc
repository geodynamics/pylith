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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>
#include <stdexcept>

#include "Mesh.hh" // implementation of class methods

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "pylith/utils/array.hh" // USES double_array

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

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
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Mesh::deallocate(void)
{ // deallocate
  delete _coordsys; _coordsys = 0;
} // deallocate
  
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
// Nondimensionalizer the finite-element mesh.
void 
pylith::topology::Mesh::nondimensionalize(const spatialdata::units::Nondimensional& normalizer)
{ // initialize
  // Get coordinates (currently dimensioned).
  assert(!_mesh.isNull());
  const ALE::Obj<RealSection>& coordsSection =
    _mesh->getRealSection("coordinates");
  assert(!coordsSection.isNull());

  // Get field for dimensioned coordinates.
  const ALE::Obj<RealSection>& coordsDimSection =
    _mesh->getRealSection("coordinates_dimensioned");
  assert(!coordsDimSection.isNull());
  coordsDimSection->setAtlas(coordsSection->getAtlas());
  coordsDimSection->allocateStorage();
  coordsDimSection->setBC(coordsSection->getBC());

  const double lengthScale = normalizer.lengthScale();
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    _mesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = 
    vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = 
    vertices->end();

  double coordsVertex[3];
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
      v_iter != verticesEnd;
      ++v_iter) {
    const int spaceDim = coordsSection->getFiberDimension(*v_iter);
    assert(spaceDim <= 3);
    const double* coordsDimVertex = coordsSection->restrictPoint(*v_iter);
    
    // Update section with dimensioned coordinates
    assert(spaceDim == 
	   coordsDimSection->getFiberDimension(*v_iter));
    coordsDimSection->updatePoint(*v_iter, coordsDimVertex);

    // Copy coordinates to array for nondimensionalization.
    for (int i=0; i < spaceDim; ++i)
      coordsVertex[i] = coordsDimVertex[i];

    // Nondimensionalize original coordinates.
    normalizer.nondimensionalize(&coordsVertex[0], spaceDim, lengthScale);
    
    // Update section with nondimensional coordinates
    assert(spaceDim == 
	   coordsSection->getFiberDimension(*v_iter));
    coordsSection->updatePoint(*v_iter, coordsVertex);
  } // for
} // nondimensionalize

// ----------------------------------------------------------------------
// Return the names of all vertex groups.
void
pylith::topology::Mesh::groups(int *numNames, char ***outNames)
{ // groups
  const ALE::Obj<std::set<std::string> >& sectionNames =  _mesh->getIntSections();
  
  *numNames = sectionNames->size();
  PetscErrorCode ierr = PetscMalloc((*numNames) * sizeof(char *), outNames);
  const std::set<std::string>::const_iterator namesEnd = sectionNames->end();
  int i = 0;
  for (std::set<std::string>::const_iterator name = sectionNames->begin(); name != namesEnd; ++name) {
    char *newName;

    ierr = PetscStrallocpy(name->c_str(), &newName);
    (*outNames)[i++] = newName;
  }
} // groups

// ----------------------------------------------------------------------
// Return the names of all vertex groups.
int
pylith::topology::Mesh::groupSize(const char *name)
{ // groupSize
  if (!_mesh->hasIntSection(name)) {
    std::ostringstream msg;
    msg << "Cannot get size of group '" << name
	<< "'. Group missing from mesh.";
    throw std::runtime_error(msg.str());
  } // if

  const ALE::Obj<IntSection>&            group    = _mesh->getIntSection(name);
  const IntSection::chart_type&          chart    = group->getChart();
  IntSection::chart_type::const_iterator chartEnd = chart.end();
  int                                    size     = 0;

  for(IntSection::chart_type::const_iterator c_iter = chart.begin();
      c_iter != chartEnd;
      ++c_iter) {
    if (group->getFiberDimension(*c_iter))
      size++;
  } // for

  return size;
} // groupSize


// End of file 
