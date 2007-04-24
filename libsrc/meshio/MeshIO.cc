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

#include "MeshIO.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES double_array, int_array

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIO::MeshIO(void) :
  _useIndexZero(true),
  _debug(false),
  _interpolate(false),
  _mesh(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIO::~MeshIO(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Get spatial dimension of mesh.
int
pylith::meshio::MeshIO::getMeshDim(void) const
{ // getMeshDim
  assert(0 != _mesh);
  assert(!_mesh->isNull());
  return (*_mesh)->getDimension();
} // getMeshDim

// ----------------------------------------------------------------------
// Read mesh from file.
void 
pylith::meshio::MeshIO::read(ALE::Obj<Mesh>* mesh)
{ // read
  assert(0 == _mesh);

  _mesh = mesh;
  _read();
  _mesh = 0;
} // read

// ----------------------------------------------------------------------
// Write mesh to file.
void 
pylith::meshio::MeshIO::write(ALE::Obj<Mesh>* mesh)
{ // write
  assert(0 == _mesh);

  _mesh = mesh;
  _write();
  _mesh = 0;
} // write

// ----------------------------------------------------------------------
// Create cube boundary mesh.
void 
pylith::meshio::MeshIO::createCubeBoundary(ALE::Obj<Mesh>* mesh)
{ // read
  double lower[3] = {0.0, 0.0, 0.0};
  double upper[3] = {1.0, 1.0, 1.0};
  int    faces[3] = {1, 1, 1};

  *mesh = ALE::MeshBuilder::createCubeBoundary(PETSC_COMM_WORLD, lower, upper, faces, _debug);
  mesh->addRef();
} // createCubeBoundary

// ----------------------------------------------------------------------
// Set vertices in mesh.
void
pylith::meshio::MeshIO::_buildMesh(const double_array& coordinates,
				   const int numVertices,
				   const int spaceDim,
				   const int_array& cells,
				   const int numCells,
				   const int numCorners,
				   const int meshDim)
{ // _buildMesh
  assert(0 != _mesh);
  assert(coordinates.size() == numVertices*spaceDim);
  assert(cells.size() == numCells*numCorners);

  *_mesh = new Mesh(PETSC_COMM_WORLD, meshDim);
  _mesh->addRef();
  assert(!_mesh->isNull());
  (*_mesh)->setDebug(_debug);

  ALE::Obj<sieve_type> sieve = new sieve_type((*_mesh)->comm());

  ALE::SieveBuilder<Mesh>::buildTopology(sieve, meshDim, 
                                         numCells, 
                                         const_cast<int*>(&cells[0]), 
                                         numVertices, 
                                         _interpolate, numCorners, -1,
                                         (*_mesh)->getArrowSection("orientation"));
  (*_mesh)->setSieve(sieve);
  (*_mesh)->stratify();
  ALE::SieveBuilder<Mesh>::buildCoordinates(*_mesh, spaceDim, &coordinates[0]);
} // _buildMesh

// ----------------------------------------------------------------------
// Get coordinates of vertices in mesh.
void
pylith::meshio::MeshIO::_getVertices(double_array* coordinates,
				     int* numVertices,
				     int* spaceDim) const
{ // _getVertices
  assert(0 != coordinates);
  assert(0 != numVertices);
  assert(0 != spaceDim);
  assert(0 != _mesh);
  assert(!_mesh->isNull());

  const ALE::Obj<Mesh::label_sequence>& vertices = (*_mesh)->depthStratum(0);
  assert(!vertices.isNull());
  const ALE::Obj<real_section_type>& coordsField =
    (*_mesh)->getRealSection("coordinates");
  assert(!coordsField.isNull());

  *numVertices = vertices->size();
  *spaceDim = coordsField->getFiberDimension(*vertices->begin());

  const int size = (*numVertices) * (*spaceDim);
  coordinates->resize(size);

  int i = 0;
  for(Mesh::label_sequence::iterator v_iter = vertices->begin();
      v_iter != vertices->end();
      ++v_iter) {
    const real_section_type::value_type *vertexCoords = 
      coordsField->restrictPoint(*v_iter);
    for (int iDim=0; iDim < *spaceDim; ++iDim)
      (*coordinates)[i++] = vertexCoords[iDim];
  } // for
} // _getVertices

// ----------------------------------------------------------------------
// Get cells in mesh.
void
pylith::meshio::MeshIO::_getCells(int_array* cells,
				  int* numCells,
				  int* numCorners,
				  int* meshDim) const
{ // _getCells
  assert(0 != cells);
  assert(0 != numCells);
  assert(0 != meshDim);
  assert(0 != _mesh);
  assert(!_mesh->isNull());

  const ALE::Obj<sieve_type>& sieve = (*_mesh)->getSieve();
  assert(!sieve.isNull());
  const ALE::Obj<Mesh::label_sequence>& meshCells = (*_mesh)->heightStratum(0);
  assert(!meshCells.isNull());

  *meshDim = (*_mesh)->getDimension();
  *numCells = meshCells->size();
  *numCorners = sieve->nCone(*meshCells->begin(), (*_mesh)->depth())->size();
  
  const ALE::Obj<Mesh::numbering_type>& vNumbering = 
    (*_mesh)->getFactory()->getLocalNumbering(*_mesh, 0);

  const int size = (*numCells) * (*numCorners);
  cells->resize(size);
    
  const int offset = (useIndexZero()) ? 0 : 1;
  int i = 0;
  for(Mesh::label_sequence::iterator e_iter = meshCells->begin();
      e_iter != meshCells->end();
      ++e_iter) {
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = 
      sieve->cone(*e_iter);
    for(sieve_type::traits::coneSequence::iterator c_iter = cone->begin();
	c_iter != cone->end();
	++c_iter)
      (*cells)[i++] = vNumbering->getIndex(*c_iter) + offset;
  } // for
} // _getCells

// ----------------------------------------------------------------------
// Tag cells in mesh with material identifiers.
void
pylith::meshio::MeshIO::_setMaterials(const int_array& materialIds)
{ // _setMaterials
  assert(0 != _mesh);
  assert(!_mesh->isNull());
  
  const ALE::Obj<Mesh::label_sequence>& cells = (*_mesh)->heightStratum(0);
  assert(!cells.isNull());

  const int numCells = materialIds.size();
  if (cells->size() != numCells) {
    std::ostringstream msg;
    msg << "Mismatch in size of materials identifier array ("
	<< numCells << ") and number of cells in mesh ("
	<< cells->size() << ").";
    throw std::runtime_error(msg.str());
  } // if

  const ALE::Obj<Mesh::label_type>& labelMaterials = 
    (*_mesh)->createLabel("material-id");
  
  int i = 0;
  for(Mesh::label_sequence::iterator e_iter = cells->begin();
      e_iter != cells->end();
      ++e_iter)
    (*_mesh)->setValue(labelMaterials, *e_iter, materialIds[i++]);
} // _setMaterials

// ----------------------------------------------------------------------
// Get material identifiers for cells.
void
pylith::meshio::MeshIO::_getMaterials(int_array* materialIds) const
{ // _getMaterials
  assert(0 != materialIds);
  assert(0 != _mesh);
  assert(!_mesh->isNull());

  const ALE::Obj<Mesh::label_sequence>& cells = (*_mesh)->heightStratum(0);
  assert(!cells.isNull());
  const int numCells = cells->size();

  const int size = numCells;
  materialIds->resize(size);
  
  const ALE::Obj<Mesh::label_type>& labelMaterials = 
    (*_mesh)->getLabel("material-id");
  const int idDefault = 0;

  int i = 0;
  for(Mesh::label_sequence::iterator e_iter = cells->begin();
      e_iter != cells->end();
      ++e_iter)
    (*materialIds)[i++] = 
      (*_mesh)->getValue(labelMaterials, *e_iter, idDefault);
} // _getMaterials

// ----------------------------------------------------------------------
// Build a point group as an int sectio
void
pylith::meshio::MeshIO::_setGroup(const std::string& name,
				  const GroupPtType type,
				  const int_array& points)
{ // _setGroup
  assert(0 != _mesh);
  assert(!_mesh->isNull());

  const ALE::Obj<int_section_type>& groupField = (*_mesh)->getIntSection(name);
  assert(!groupField.isNull());

  const int numPoints = points.size();
  if (CELL == type)
    for(int i=0; i < numPoints; ++i)
      groupField->setFiberDimension(points[i], 1);
  else if (VERTEX == type) {
    const int numCells = (*_mesh)->heightStratum(0)->size();
    for(int i=0; i < numPoints; ++i)
      groupField->setFiberDimension(points[i]+numCells, 1);
  } // if/else
  (*_mesh)->allocate(groupField);
} // _setGroup

// ----------------------------------------------------------------------
// Get names of all groups in mesh.
void
pylith::meshio::MeshIO::_getGroupNames(string_vector* names) const
{ // _getGroups
  assert(0 != names);
  assert(0 != _mesh);
  assert(!_mesh->isNull());
  
  const ALE::Obj<std::set<std::string> >& sectionNames = 
    (*_mesh)->getIntSections();
  
  const int numGroups = sectionNames->size();
  names->resize(numGroups);
  int i=0;
  for (std::set<std::string>::const_iterator name=sectionNames->begin();
       name != sectionNames->end();
       ++name)
    (*names)[i++] = *name;
} // _getGroups

// ----------------------------------------------------------------------
// Get group entities
void
pylith::meshio::MeshIO::_getGroup(int_array* points,
				  GroupPtType* type,
				  const char *name) const
{ // _getMaterials
  assert(0 != points);
  assert(0 != type);
  assert(0 != _mesh);
  assert(!_mesh->isNull());

  const ALE::Obj<int_section_type>& groupField = (*_mesh)->getIntSection(name);
  assert(!groupField.isNull());
  const int_section_type::chart_type& chart = groupField->getChart();
  const Mesh::point_type firstPoint = *chart.begin();
  ALE::Obj<Mesh::numbering_type> numbering;

  if ((*_mesh)->height(firstPoint) == 0) {
    *type = CELL;
    numbering = (*_mesh)->getFactory()->getNumbering(*_mesh, 
						     (*_mesh)->depth());
  } else {
    *type = VERTEX;
    numbering = (*_mesh)->getFactory()->getNumbering(*_mesh, 0);
  } // if/else
  const int numPoints = chart.size();
  points->resize(numPoints);
  int i = 0;

  for(int_section_type::chart_type::iterator c_iter = chart.begin();
      c_iter != chart.end();
      ++c_iter) {
    assert(!numbering.isNull());
    (*points)[i++] = numbering->getIndex(*c_iter);
  } // for
} // _getMaterials


// End of file 
