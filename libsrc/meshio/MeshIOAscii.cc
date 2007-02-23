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

#include "MeshIOAscii.hh" // implementation of class methods

#include <fstream> // USES std::ifstream, std::ofstream
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <assert.h> // USES assert()
#include <iomanip> // USES setw(), setiosflags(), resetiosflags()

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOAscii::MeshIOAscii(void) :
  _filename("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOAscii::~MeshIOAscii(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Unpickle mesh
void
pylith::meshio::MeshIOAscii::_read(void)
{ // _read
  int meshDim = 0;
  int spaceDim = 0;
  int numVertices = 0;
  int numCells = 0;
  int numCorners = 0;
  double* coordinates = 0;
  int* cells = 0;
  int* materialIds = 0;
  
  std::ifstream filein(_filename.c_str());
  if (!filein.is_open() || !filein.good()) {
    std::ostringstream msg;
    msg << "Could not open mesh file '" << _filename
	<< "' for reading.\n";
    throw std::runtime_error(msg.str());
  } // if

  std::string token;
  const int maxIgnore = 1024;
  filein >> token;
  if (0 != strcasecmp(token.c_str(), "mesh")) {
    std::ostringstream msg;
    msg << "Expected 'mesh' token but encountered '" << token << "'\n";
    throw std::runtime_error(msg.str());
  }

  bool readDim = false;
  bool readCells = false;
  bool readVertices = false;
  bool builtMesh = false;

  filein.ignore(maxIgnore, '{');
  filein >> token;
  while (filein.good() && token != "}") {
    if (0 == strcasecmp(token.c_str(), "dimension")) {
      filein.ignore(maxIgnore, '=');
      filein >> meshDim;
      readDim = true;
    } else if (0 == strcasecmp(token.c_str(), "use-index-zero")) {
      filein.ignore(maxIgnore, '=');
      std::string flag = "";
      filein >> flag;
      if (0 == strcasecmp(flag.c_str(), "true"))
        useIndexZero(true);
      else
        useIndexZero(false);
    } else if (0 == strcasecmp(token.c_str(), "vertices")) {
      filein.ignore(maxIgnore, '{');
      _readVertices(filein, &coordinates, &numVertices, &spaceDim);
      readVertices = true;
    } else if (0 == strcasecmp(token.c_str(), "cells")) {
      filein.ignore(maxIgnore, '{');
      _readCells(filein, &cells, &materialIds, &numCells, &numCorners);
      readCells = true;
    } else if (0 == strcasecmp(token.c_str(), "group")) {
      if (!builtMesh)
	throw std::runtime_error("Both 'vertices' and 'cells' must "
				 "precede any groups in mesh file.");
      filein.ignore(maxIgnore, '{');
      //_readGroup(filein, mesh);
    } else {
      std::ostringstream msg;
      msg << "Could not parse '" << token << "' into a mesh setting.";
      throw std::runtime_error(msg.str());  
    } // else

    if (readDim && readCells && readVertices && !builtMesh) {
      // Can now build mesh
      _buildMesh(coordinates, numVertices, spaceDim,
		 cells, numCells, numCorners, meshDim);
      _setMaterials(materialIds, numCells);
      delete[] coordinates; coordinates = 0;
      delete[] cells; cells = 0;
      delete[] materialIds; materialIds = 0;
      builtMesh = true;
    } // if

    filein >> token;
  } // while
  filein.close();
} // read

// ----------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIOAscii::_write(void) const
{ // write
  std::ofstream fileout(_filename.c_str());
  if (!fileout.is_open() || !fileout.good()) {
    std::ostringstream msg;
    msg << "Could not open mesh file '" << _filename
	<< "' for writing.\n";
    throw std::runtime_error(msg.str());
  } // if

  fileout
    << "mesh = {\n"
    << "  dimension = " << getMeshDim() << "\n"
    << "  use-index-zero = " << (useIndexZero() ? "true" : "false") << "\n";

  _writeVertices(fileout);
  _writeCells(fileout);

  // LOOP OVER GROUPS
  // _writeGroup(fileout, mesh, nameIter->c_str());

  fileout << "}\n";
  fileout.close();
} // write

// ----------------------------------------------------------------------
// Read mesh vertices.
void
pylith::meshio::MeshIOAscii::_readVertices(std::istream& filein,
					   double** pCoordinates,
					   int* pNumVertices, 
					   int* pNumDims) const
{ // _readVertices
  double* coordinates = 0;
  int numDims = 0;
  int numVertices = 0;

  std::string token;
  const int maxIgnore = 1024;
  filein >> token;
  while (filein.good() && token != "}") {
    if (0 == strcasecmp(token.c_str(), "dimension")) {
      filein.ignore(maxIgnore, '=');
      filein >> numDims;
    } else if (0 == strcasecmp(token.c_str(), "count")) {
      filein.ignore(maxIgnore, '=');
      filein >> numVertices;
    } else if (0 == strcasecmp(token.c_str(), "coordinates")) {
      const int size = numVertices*numDims;
      if (0 == size) {
	const char* msg = 
	  "Tokens 'dimension' and 'count' must precede 'coordinates'.";
	throw std::runtime_error(msg);
      } // if
      
      filein.ignore(maxIgnore, '{');
      delete[] coordinates; coordinates = new double[size];
      assert(0 != coordinates);
      for (int i=0; i < size; ++i)
	filein >> coordinates[i];
      filein.ignore(maxIgnore, '}');
    } else {
      std::ostringstream msg;
      msg << "Could not parse '" << token << "' into a vertices setting.";
      throw std::runtime_error(msg.str());
    } // else
    filein >> token;
  } // while
  if (!filein.good())
    throw std::runtime_error("I/O error while parsing vertices settings.");

  if (0 != pCoordinates)
    *pCoordinates = coordinates;
  if (0 != pNumVertices)
    *pNumVertices = numVertices;
  if (0 != pNumDims)
    *pNumDims = numDims;
} // _readVertices

// ----------------------------------------------------------------------
// Write mesh vertices.
void
pylith::meshio::MeshIOAscii::_writeVertices(std::ostream& fileout) const
{ // _writeVertices
  int spaceDim = 0;
  int numVertices = 0;
  double* coordinates = 0;
  _getVertices(&coordinates, &numVertices, &spaceDim);

  fileout
    << "  vertices = {\n"
    << "    dimension = " << spaceDim << "\n"
    << "    count = " << numVertices << "\n"
    << "    coordinates = {\n"
    << std::resetiosflags(std::ios::fixed)
    << std::setiosflags(std::ios::scientific)
    << std::setprecision(6);
  for(int iVertex=0, i=0; iVertex < numVertices; ++iVertex) {
    fileout << "      ";
    for(int iDim=0; iDim < spaceDim; ++iDim)
      fileout << std::setw(18) << coordinates[i++];
    fileout << "\n";
  } // for
  fileout
    << "    }\n"
    << "  }\n";

  delete[] coordinates; coordinates = 0;
} // _writeVertices
  
// ----------------------------------------------------------------------
// Read mesh cells.
void
pylith::meshio::MeshIOAscii::_readCells(std::istream& filein,
					int** pCells,
					int** pMaterialIds,
					int* pNumCells, 
					int* pNumCorners) const
{ // _readCells
  int* cells = 0;
  int* materialIds = 0;
  int numCells = 0;
  int numCorners = 0;
  int dimension = 0;

  std::string token;
  const int maxIgnore = 1024;
  filein >> token;
  while (filein.good() && token != "}") {
    if (0 == strcasecmp(token.c_str(), "num-corners")) {
      filein.ignore(maxIgnore, '=');
      filein >> numCorners;
    } else if (0 == strcasecmp(token.c_str(), "count")) {
      filein.ignore(maxIgnore, '=');
      filein >> numCells;
    } else if (0 == strcasecmp(token.c_str(), "simplices")) {
      const int size = numCells*numCorners;
      if (0 == size) {
	const char* msg = 
	  "Tokens 'num-corners' and 'count' must precede 'cells'.";
	throw std::runtime_error(msg);
      } // if
      
      filein.ignore(maxIgnore, '{');
      delete[] cells; cells = new int[size];
      assert(0 != cells);
      for (int i=0; i < size; ++i)
	filein >> cells[i];
      if (!useIndexZero()) {
	// if files begins with index 1, then decrement to index 0
	// for compatibility with Sieve
	for (int i=0; i < size; ++i)
	  --cells[i];
      } // if
      
      filein.ignore(maxIgnore, '}');
    } else if (0 == strcasecmp(token.c_str(), "material-ids")) {
      if (0 == numCells) {
	const char* msg =
	  "Token 'count' must precede 'material-ids'.";
	throw std::runtime_error(msg);
      } // if
      const int size = numCells;
      filein.ignore(maxIgnore, '{');
      delete[] materialIds; materialIds = new int[size];
      assert(0 != materialIds);
      for (int i=0; i < size; ++i)
	filein >> materialIds[i];      
      filein.ignore(maxIgnore, '}');
    } else {
      std::ostringstream msg;
      msg << "Could not parse '" << token << "' into an cells setting.";
      throw std::runtime_error(msg.str());
    } // else
    filein >> token;
  } // while
  if (!filein.good())
    throw std::runtime_error("I/O error while parsing cells settings.");

  // If no materials given, assign each cell material identifier of 0
  if (0 == materialIds && numCells > 0) {
    const int size = numCells;
    materialIds = new int[size];
    for (int i=0; i < size; ++i)
      materialIds[i] = 0;
  } // if

  if (0 != pCells)
    *pCells = cells;
  if (0 != pNumCells)
    *pNumCells = numCells;
  if (0 != pNumCorners)
    *pNumCorners = numCorners;
  if (0 != pMaterialIds)
    *pMaterialIds = materialIds;
} // _readCells

// ----------------------------------------------------------------------
// Write mesh cells.
void
pylith::meshio::MeshIOAscii::_writeCells(std::ostream& fileout) const
{ // _writeCells
  int meshDim = 0;
  int numCells = 0;
  int numCorners = 0;
  int* cells = 0;
  _getCells(&cells, &numCells, &numCorners, &meshDim);
  
  fileout
    << "  cells = {\n"
    << "    count = " << numCells << "\n"
    << "    num-corners = " << numCorners << "\n"
    << "    simplices = {\n";

  for(int iCell=0, i=0; iCell < numCells; ++iCell) {
    fileout << "      ";
    for (int iCorner=0; iCorner < numCorners; ++iCorner)
      fileout << std::setw(8)
	      << cells[i++];
    fileout << "\n";
  } // for
  fileout
    << "    }\n";
  delete[] cells; cells = 0;

  // Write material identifiers
  int* materialIds = 0;
  _getMaterials(&materialIds, &numCells);
  assert( (0 != materialIds && 0 < numCells) ||
	  (0 == materialIds && 0 == numCells) );
  fileout
    << "    material-ids = {\n";
  for(int iCell=0, i=0; iCell < numCells; ++iCell)
    fileout << "      " << materialIds[iCell] << "\n";
  fileout
    << "    }\n";  
  delete[] materialIds; materialIds = 0;

  fileout
    << "  }\n";
} // _writeCells

#if 0
// ----------------------------------------------------------------------
// Read mesh group.
void
pylith::meshio::MeshIOAscii::_readGroup(std::istream& filein,
					const ALE::Obj<Mesh>& mesh) const
{ // _readGroup
  std::string name = ""; // Name of group
  int dimension = 0; // Topology dimension associated with group
  int count = 0; // Number of entities in group
  int* indices = 0; // Indices of entities in group

  std::string token;
  const int maxIgnore = 1024;
  filein >> token;
  while (filein.good() && token != "}") {
    if (0 == strcasecmp(token.c_str(), "name")) {
      filein.ignore(maxIgnore, '=');
      filein >> name;
    } else if (0 == strcasecmp(token.c_str(), "dimension")) {
      filein.ignore(maxIgnore, '=');
      filein >> dimension;
    } else if (0 == strcasecmp(token.c_str(), "count")) {
      filein.ignore(maxIgnore, '=');
      filein >> count;
    } else if (0 == strcasecmp(token.c_str(), "indices")) {
      if (0 == count) {
	std::ostringstream msg;
	msg << "Tokens 'count' must precede 'indices'.";
	throw std::runtime_error(msg.str());
      } // if
      
      filein.ignore(maxIgnore, '{');
      delete[] indices; indices = new int[count];
      assert(0 != indices);
      for (int i = 0; i < count; ++i)
	filein >> indices[i];
      filein.ignore(maxIgnore, '}');
    } else {
      std::ostringstream msg;
      msg << "Could not parse '" << token << "' into a group setting.";
      throw std::runtime_error(msg.str());
    } // else
    filein >> token;
  } // while
  if (!filein.good())
    throw std::runtime_error("I/O error while parsing group settings.");

  assert(!mesh.isNull());
  ALE::Obj<Mesh::field_type> groupField = mesh->getField(name);
  const int meshDim = mesh->getDimension();
  ALE::Obj<std::list<Mesh::point_type> > patchPoints = 
    std::list<Mesh::point_type>();
  Mesh::field_type::patch_type patch;

  patchPoints->clear();
  if (meshDim == dimension) {
    for (int i=0; i < count; ++i)
      patchPoints->push_back(Mesh::point_type(0, indices[i]));
    groupField->setPatch(patchPoints, patch);
  } else if (0 == dimension) {
  } // if
  groupField->setFiberDimensionByHeight(patch, 0, 1);
  groupField->orderPatches();
  const double zero = 0;
  for (int i=0; i < count; ++i)
    groupField->update(patch, Mesh::point_type(0, i), &zero);
} // _readGroup

// ----------------------------------------------------------------------
// Write mesh group.
void
pylith::meshio::MeshIOAscii::_writeGroup(std::ostream& fileout,
					 const ALE::Obj<Mesh>& mesh,
					 const char* name) const
{ // _writeGroup
  //_writeGroup(fileout, mesh);
  // ADD STUFF HERE
  int count = 0; // TEMPORARY
  int dimension = 0; // TEMPORARY

  fileout
    << "  group = {\n"
    << "    name = " << name << "\n"
    << "    dimension = " << dimension << "\n"
    << "    count = " << count << "\n"
    << "    indices = {\n";
    
  fileout
    << "    }\n"
    << "  }\n";
} // _writeGroup
#endif
  
// End of file 
