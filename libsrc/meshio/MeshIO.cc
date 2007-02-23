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
  
  ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());
  ALE::Obj<topology_type> topology = new topology_type(mesh->comm());

  ALE::New::SieveBuilder<sieve_type>::buildTopology(sieve, meshDim, 
						    numCells, 
						    const_cast<int*>(cells), 
						    numVertices, 
						    _interpolate, numCorners);
  sieve->stratify();
  topology->setPatch(0, sieve);
  topology->stratify();
  mesh->setTopology(topology);
  ALE::New::SieveBuilder<sieve_type>::buildCoordinates(
		      mesh->getRealSection("coordinates"), 
		      spaceDim, coordinates);
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

  const Mesh::real_section_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = mesh->getTopology();

  const ALE::Obj<Mesh::topology_type::label_sequence>& vertices = 
    topology->depthStratum(patch, 0);
  const ALE::Obj<Mesh::real_section_type>& coordsField =
    mesh->getRealSection("coordinates");

  const int numVertices = vertices->size();
  const int spaceDim = 
    coordsField->getFiberDimension(patch, *vertices->begin());

  double* coordinates = 0;
  const int size = numVertices * spaceDim;
  if (0 != pCoordinates && size > 0) {
    coordinates = new double[size];

    int i = 0;
    for(Mesh::topology_type::label_sequence::iterator v_iter = 
	  vertices->begin();
	v_iter != vertices->end();
	++v_iter) {
      const Mesh::real_section_type::value_type *vertexCoords = 
	coordsField->restrict(patch, *v_iter);
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

  const topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = mesh->getTopology();

  const ALE::Obj<sieve_type>& sieve = topology->getPatch(patch);
  const ALE::Obj<Mesh::topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);

  const int meshDim = mesh->getDimension();
  const int numCells = cells->size();
  const int numCorners = sieve->nCone(*cells->begin(), 
				      topology->depth())->size();

  const ALE::Obj<Mesh::numbering_type>& vNumbering = 
    mesh->getFactory()->getLocalNumbering(topology, patch, 0);

  int* cellsArray = 0;
  const int size = numCells * numCorners;
  if (0 != pCells && size > 0) {
    cellsArray = new int[size];
    
    const int offset = (useIndexZero()) ? 0 : 1;
    int i = 0;
    for(Mesh::topology_type::label_sequence::iterator e_iter = cells->begin();
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
  
  const topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = mesh->getTopology();
  const ALE::Obj<Mesh::topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);

  if (cells->size() != numCells) {
    std::ostringstream msg;
    msg << "Mismatch in size of materials identifier array ("
	<< numCells << ") and number of cells in mesh ("
	<< cells->size() << ").";
    throw std::runtime_error(msg.str());
  } // if

  const ALE::Obj<patch_label_type>& labelMaterials = 
    topology->createLabel(patch, "material-id");
  
  int i = 0;
  for(Mesh::topology_type::label_sequence::iterator e_iter = cells->begin();
      e_iter != cells->end();
      ++e_iter)
    topology->setValue(labelMaterials, *e_iter, materialIds[i++]);
} // _setMaterials

// ----------------------------------------------------------------------
// Get material identifiers for cells.
void
pylith::meshio::MeshIO::_getMaterials(int** pMaterialIds,
				      int* pNumCells)
{ // _getMaterials
  assert(0 != _mesh);
  ALE::Obj<Mesh>& mesh = *_mesh;

  const topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type>& topology = mesh->getTopology();
  const ALE::Obj<Mesh::topology_type::label_sequence>& cells = 
    topology->heightStratum(patch, 0);
  const int numCells = cells->size();

  int* materialsArray = 0;
  const int size = numCells;
  if (0 != pMaterialIds && size > 0) {
    materialsArray = new int[size];
  
    const ALE::Obj<patch_label_type>& labelMaterials = 
      topology->getLabel(patch, "material-id");
    const int idDefault = 0;

    int i = 0;
    for(Mesh::topology_type::label_sequence::iterator e_iter = cells->begin();
	e_iter != cells->end();
	++e_iter)
      materialsArray[i++] = 
	topology->getValue(labelMaterials, *e_iter, idDefault);
  } // if  

  if (0 != pMaterialIds)
    *pMaterialIds = materialsArray;
  if (0 != pNumCells)
    *pNumCells = numCells;
} // _getMaterials


// End of file 
