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

const char *pylith::meshio::MeshIOAscii::groupTypeNames[] = {"vertices", "cells"};

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
      std::string name;
      pylith::meshio::MeshIO::PointType type;
      int  numPoints;
      int *points;

      if (!builtMesh)
        throw std::runtime_error("Both 'vertices' and 'cells' must "
				 "precede any groups in mesh file.");
      filein.ignore(maxIgnore, '{');
      _readGroup(filein, name, type, numPoints, &points);
      _buildGroup(name, type, numPoints, points);
      delete [] points;
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
  ALE::Obj<std::set<std::string> > groups = _getGroups();

  for(std::set<std::string>::const_iterator name = groups->begin(); name != groups->end(); ++name) {
    _writeGroup(fileout, name->c_str());
  }

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
      int label;
      for (int iVertex=0, i=0; iVertex < numVertices; ++iVertex) {
	filein >> label;
	for (int iDim=0; iDim < numDims; ++iDim)
	  filein >> coordinates[i++];
      } // for
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
    fileout << std::setw(8) << iVertex;
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
      int label;
      for (int iCell=0, i=0; iCell < numCells; ++iCell) {
	filein >> label;
	for (int iCorner=0; iCorner < numCorners; ++iCorner)
	  filein >> cells[i++];
      } // for
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
      int label = 0;
      for (int iCell=0; iCell < numCells; ++iCell) {
	filein >> label;
	filein >> materialIds[iCell];
      } // for
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
    fileout << "      " << std::setw(8) << iCell;
    for (int iCorner=0; iCorner < numCorners; ++iCorner)
      fileout << std::setw(8) << cells[i++];
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
  for(int iCell=0, i=0; iCell < numCells; ++iCell) {
    fileout << "      " << std::setw(8) << iCell;
    fileout << std::setw(4) << materialIds[iCell] << "\n";
  } // for
  fileout
    << "    }\n";  
  delete[] materialIds; materialIds = 0;

  fileout
    << "  }\n";
} // _writeCells

// ----------------------------------------------------------------------
// Read mesh group.
void
pylith::meshio::MeshIOAscii::_readGroup(std::istream& filein,
             std::string& name,
             pylith::meshio::MeshIO::PointType& type,
             int& numPoints,
             int *points[]) const
{ // _readGroup
  int* indices = 0; // Indices of entities in group

  std::string token;
  const int maxIgnore = 1024;
  numPoints = -1;
  filein >> token;
  while (filein.good() && token != "}") {
    if (0 == strcasecmp(token.c_str(), "name")) {
      filein.ignore(maxIgnore, '=');
      filein >> name;
    } else if (0 == strcasecmp(token.c_str(), "type")) {
      std::string typeName;
      filein.ignore(maxIgnore, '=');
      filein >> typeName;
      if (typeName == groupTypeNames[VERTEX]) {
        type = VERTEX;
      } else if (typeName == groupTypeNames[CELL]) {
        type = CELL;
      } else {
        std::ostringstream msg;
        msg << "Invalid point type " << typeName << ".";
        throw std::runtime_error(msg.str());
      }
    } else if (0 == strcasecmp(token.c_str(), "count")) {
      filein.ignore(maxIgnore, '=');
      filein >> numPoints;
    } else if (0 == strcasecmp(token.c_str(), "indices")) {
      if (-1 == numPoints) {
        std::ostringstream msg;
        msg << "Tokens 'count' must precede 'indices'.";
        throw std::runtime_error(msg.str());
      } // if
      
      filein.ignore(maxIgnore, '{');
      delete[] indices; indices = new int[numPoints];
      assert(0 != indices);
      for (int i = 0; i < numPoints; ++i)
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
  *points = indices;
} // _readGroup

// ----------------------------------------------------------------------
// Write mesh group.
void
pylith::meshio::MeshIOAscii::_writeGroup(std::ostream& fileout,
					 const char* name) const
{ // _writeGroup
  pylith::meshio::MeshIO::PointType         type;
  int numPoints;
  int *points;

  _getGroup(name, type, numPoints, &points);
  fileout
    << "  group = {\n"
    << "    name = " << name << "\n"
    << "    type = " << groupTypeNames[type] << "\n"
    << "    count = " << numPoints << "\n"
    << "    indices = {\n";
  for(int i = 0; i < numPoints; ++i) {
    fileout << "      " << points[i] << "\n";
  }
  fileout
    << "    }\n"
    << "  }\n";
} // _writeGroup
  
// End of file 
