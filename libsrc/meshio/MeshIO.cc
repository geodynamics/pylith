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

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIO::MeshIO(void) :
  _useIndexZero(true),
  _interpolate(false),
  _mesh(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIO::~MeshIO(void)
{ // destructor
  assert(0 == _mesh);
} // destructor
  
// ----------------------------------------------------------------------
// Get spatial dimension of mesh.
int
pylith::meshio::MeshIO::getMeshDim(void) const
{ // getMeshDim
  assert(0 != _mesh);
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
// Set vertices in mesh.
void
pylith::meshio::MeshIO::_buildMesh(const double* coordinates,
				   const int numVertices,
				   const int spaceDim,
				   const int* cells,
				   const int numCells,
				   const int numCorners,
				   const int meshDim)
{ // _buildMesh
  assert(0 != _mesh);

  *_mesh = new Mesh(PETSC_COMM_WORLD, meshDim);
  ALE::Obj<Mesh>& mesh = *_mesh;
  mesh.addRef();
  
  ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());

  ALE::New::SieveBuilder<Mesh>::buildTopology(sieve, meshDim, 
						    numCells, 
						    const_cast<int*>(cells), 
						    numVertices, 
						    _interpolate, numCorners);
  mesh->setSieve(sieve);
  mesh->stratify();
  ALE::New::SieveBuilder<Mesh>::buildCoordinatesNew(mesh, spaceDim, coordinates);
} // _buildMesh

// ----------------------------------------------------------------------
// Get coordinates of vertices in mesh.
void
pylith::meshio::MeshIO::_getVertices(double** pCoordinates,
				     int* pNumVertices,
				     int* pSpaceDim) const
{ // _getVertices
  assert(0 != _mesh);
  ALE::Obj<Mesh>& mesh = *_mesh;
  mesh.addRef();

  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const ALE::Obj<Mesh::real_section_type>& coordsField =
    mesh->getRealSection("coordinates");

  const int numVertices = vertices->size();
  const int spaceDim = 
    coordsField->getFiberDimension(*vertices->begin());

  double* coordinates = 0;
  const int size = numVertices * spaceDim;
  if (0 != pCoordinates && size > 0) {
    coordinates = new double[size];

    int i = 0;
    for(Mesh::label_sequence::iterator v_iter = vertices->begin();
        v_iter != vertices->end();
	    ++v_iter) {
      const Mesh::real_section_type::value_type *vertexCoords = 
	coordsField->restrictPoint(*v_iter);
      for (int iDim=0; iDim < spaceDim; ++iDim)
	coordinates[i++] = vertexCoords[iDim];
    } // for
  } // if

  if (0 != pCoordinates)
    *pCoordinates = coordinates;
  if (0 != pNumVertices)
    *pNumVertices = numVertices;
  if (0 != pSpaceDim)
    *pSpaceDim = spaceDim;
} // _getVertices

// ----------------------------------------------------------------------
// Get cells in mesh.
void
pylith::meshio::MeshIO::_getCells(int** pCells,
				  int* pNumCells,
				  int* pNumCorners,
				  int* pMeshDim) const
{ // _getCells
  assert(0 != _mesh);
  ALE::Obj<Mesh>& mesh = *_mesh;
  mesh.addRef();

  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

  const int meshDim = mesh->getDimension();
  const int numCells = cells->size();
  const int numCorners = sieve->nCone(*cells->begin(), 
				      mesh->depth())->size();

  const ALE::Obj<Mesh::numbering_type>& vNumbering = 
    mesh->getFactory()->getLocalNumbering(mesh, 0);

  int* cellsArray = 0;
  const int size = numCells * numCorners;
  if (0 != pCells && size > 0) {
    cellsArray = new int[size];
    
    const int offset = (useIndexZero()) ? 0 : 1;
    int i = 0;
    for(Mesh::label_sequence::iterator e_iter = cells->begin();
	e_iter != cells->end();
	++e_iter) {
      const ALE::Obj<sieve_type::traits::coneSequence>& cone = 
	sieve->cone(*e_iter);
      for(sieve_type::traits::coneSequence::iterator c_iter = cone->begin();
	  c_iter != cone->end();
	  ++c_iter)
	cellsArray[i++] = vNumbering->getIndex(*c_iter) + offset;
    } // for
  } // if  

  if (0 != pCells)
    *pCells = cellsArray;
  if (0 != pNumCells)
    *pNumCells = numCells;
  if (0 != pNumCorners)
    *pNumCorners = numCorners;
  if (0 != pMeshDim)
    *pMeshDim = meshDim;
} // _getCells

// ----------------------------------------------------------------------
// Tag cells in mesh with material identifiers.
void
pylith::meshio::MeshIO::_setMaterials(const int* materialIds,
				      const int numCells)
{ // _setMaterials
  assert(0 != _mesh);
  ALE::Obj<Mesh>& mesh = *_mesh;
  mesh.addRef();
  
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

  if (cells->size() != numCells) {
    std::ostringstream msg;
    msg << "Mismatch in size of materials identifier array ("
	<< numCells << ") and number of cells in mesh ("
	<< cells->size() << ").";
    throw std::runtime_error(msg.str());
  } // if

  const ALE::Obj<label_type>& labelMaterials = mesh->createLabel("material-id");
  
  int i = 0;
  for(Mesh::label_sequence::iterator e_iter = cells->begin();
      e_iter != cells->end();
      ++e_iter)
    mesh->setValue(labelMaterials, *e_iter, materialIds[i++]);
} // _setMaterials

// ----------------------------------------------------------------------
// Get material identifiers for cells.
void
pylith::meshio::MeshIO::_getMaterials(int** pMaterialIds,
				      int* pNumCells) const
{ // _getMaterials
  assert(0 != _mesh);
  ALE::Obj<Mesh>& mesh = *_mesh;
  mesh.addRef();

  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);
  const int numCells = cells->size();

  int* materialsArray = 0;
  const int size = numCells;
  if (0 != pMaterialIds && size > 0) {
    materialsArray = new int[size];
  
    const ALE::Obj<label_type>& labelMaterials = mesh->getLabel("material-id");
    const int idDefault = 0;

    int i = 0;
    for(Mesh::label_sequence::iterator e_iter = cells->begin();
	e_iter != cells->end();
	++e_iter)
      materialsArray[i++] = 
	mesh->getValue(labelMaterials, *e_iter, idDefault);
  } // if  

  if (0 != pMaterialIds)
    *pMaterialIds = materialsArray;
  if (0 != pNumCells)
    *pNumCells = numCells;
} // _getMaterials

// ----------------------------------------------------------------------
// Build a point group as an int sectio
void
pylith::meshio::MeshIO::_buildGroup(const std::string& name,
                   const PointType type,
                   const int numPoints,
				   const int* points)
{ // _buildMesh
  assert(0 != _mesh);

  ALE::Obj<Mesh>& mesh = *_mesh;
  mesh.addRef();

  const ALE::Obj<Mesh::int_section_type>& groupField = mesh->getIntSection(name);

  if (CELL == type) {
    for(int i = 0; i < numPoints; ++i) {
      groupField->setFiberDimension(points[i], 1);
    }
  } else if (VERTEX == type) {
    const int numCells = mesh->heightStratum(0)->size();
    for(int i = 0; i < numPoints; ++i) {
      groupField->setFiberDimension(points[i]+numCells, 1);
    }
  } // if
  mesh->allocate(groupField);
  groupField->view("Group Field");
} // _buildMesh

// ----------------------------------------------------------------------
// Get group names
ALE::Obj<std::set<std::string> >
pylith::meshio::MeshIO::_getGroups() const
{
  return (*_mesh)->getIntSections();
}

// ----------------------------------------------------------------------
// Get group entities
void
pylith::meshio::MeshIO::_getGroup(const char *name,
                      PointType& type,
                      int& numPoints,
				      int *points[]) const
{ // _getMaterials
  assert(0 != _mesh);
  ALE::Obj<Mesh>& mesh = *_mesh;
  mesh.addRef();

  const ALE::Obj<Mesh::int_section_type>&   groupField = mesh->getIntSection(name);
  const Mesh::int_section_type::chart_type& chart      = groupField->getChart();
  const Mesh::point_type                    firstPoint = *chart.begin();
  ALE::Obj<Mesh::numbering_type>            numbering;

  if (mesh->height(firstPoint) == 0) {
    type      = CELL;
    numbering = mesh->getFactory()->getNumbering(mesh, mesh->depth());
  } else {
    type      = VERTEX;
    numbering = mesh->getFactory()->getNumbering(mesh, 0);
  }
  numPoints = chart.size();
  int *indices = new int[numPoints];
  int  i       = 0;

  for(Mesh::int_section_type::chart_type::iterator c_iter = chart.begin(); c_iter != chart.end(); ++c_iter) {
    indices[i++] = numbering->getIndex(*c_iter);
  }
  *points = indices;
} // _getMaterials

// End of file 
