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
  const ALE::Obj<IntSection>&            group    = _mesh->getIntSection(name);
  const IntSection::chart_type&          chart    = group->getChart();
  IntSection::chart_type::const_iterator chartEnd = chart.end();
  int                                    size     = 0;

  for(IntSection::chart_type::const_iterator c_iter = chart.begin(); c_iter != chartEnd; ++c_iter) {
    if (group->getFiberDimension(*c_iter)) {
      size++;
    }
  }
  return size;
} // groupSize


// End of file 
