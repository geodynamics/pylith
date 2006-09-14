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

#include "MeshIO.hh" // MeshIOAscii ISA MeshIO
#include "MeshIOAscii.hh" // implementation of class methods

#include "PetscMesh.hh"

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
pylith::meshio::MeshIOAscii::read(ALE::Obj<ALE::PetscMesh>* pMesh)
{ // read
  assert(0 != pMesh);

  int meshDim = 0;
  int numDims = 0;
  int numVertices = 0;
  int numElements = 0;
  int numCorners = 0;
  double* coordinates = 0;
  int* elements = 0;

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
  bool readElements = false;
  bool readVertices = false;
  bool builtTopology = false;

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
      _readVertices(filein, &coordinates, &numVertices, &numDims);
      readVertices = true;
    } else if (0 == strcasecmp(token.c_str(), "elements")) {
      filein.ignore(maxIgnore, '{');
      _readElements(filein, &elements, &numElements, &numCorners);
      readElements = true;
    } else if (0 == strcasecmp(token.c_str(), "chart")) {
      if (!builtTopology)
	throw std::runtime_error("Both 'vertices' and 'elements' must "
				 "precede any charts in mesh file.");
      filein.ignore(maxIgnore, '{');
      _readChart(filein, pMesh);
    } else {
      std::ostringstream msg;
      msg << "Could not parse '" << token << "' into a mesh setting.";
      throw std::runtime_error(msg.str());  
    } // else

    if (readDim && readElements && readVertices && !builtTopology) {
      // Can now build topology
      *pMesh = ALE::PetscMesh(PETSC_COMM_WORLD, meshDim);
      (*pMesh)->debug = true;
      bool interpolate = false;

#if 1
// allow mesh to have different dimension than coordinates
      ALE::Obj<ALE::PetscMesh::sieve_type> topology = (*pMesh)->getTopology();
      topology->setStratification(false);
      (*pMesh)->buildTopology(numElements, elements, numVertices, 
			      interpolate, numCorners);
      topology->stratify();
      topology->setStratification(true);
      (*pMesh)->createVertexBundle(numElements, elements, 0, numCorners);
      (*pMesh)->createSerialCoordinates(numDims, numElements, coordinates);
#else
// require mesh to have same dimension as coordinates
      (*pMesh)->populate(numElements, elements, numVertices, coordinates, 
			 interpolate, numCorners);
#endif
      delete[] coordinates; coordinates = 0;
      delete[] elements; elements = 0;
      builtTopology = true;
    } // if

    filein >> token;
  } // while
  filein.close();
} // read

// ----------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIOAscii::write(const ALE::Obj<ALE::PetscMesh>& mesh) const
{ // write
  std::ofstream fileout(_filename.c_str());
  if (!fileout.is_open() || !fileout.good()) {
    std::ostringstream msg;
    msg << "Could not open mesh file '" << _filename
	<< "' for writing.\n";
    throw std::runtime_error(msg.str());
  } // if

  const int dimension = mesh->getDimension();

  fileout
    << "mesh = {\n"
    << "  dimension = " << dimension << "\n"
    << "  use-index-zero = " << (useIndexZero() ? "true" : "false") << "\n";

  _writeVertices(fileout, mesh);
  _writeElements(fileout, mesh);

  // LOOP OVER CHARTS
  // _writeChart(fileout, mesh, nameIter->c_str());

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
	std::ostringstream msg;
	msg << "Tokens 'dimension' and 'count' must precede 'coordinates'.";
	throw std::runtime_error(msg.str());
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
pylith::meshio::MeshIOAscii::_writeVertices(std::ostream& fileout,
			       const ALE::Obj<ALE::PetscMesh>& mesh) const
{ // _writeVertices
  ALE::Obj<ALE::PetscMesh::field_type> coords_field = mesh->getCoordinates();
  ALE::Obj<ALE::PetscMesh::bundle_type> vertexBundle = mesh->getBundle(0);
  ALE::PetscMesh::field_type::patch_type patch;
  const double* coordinates = coords_field->restrict(patch);
  const int numVertices = (vertexBundle->getGlobalOffsets()) ?
    vertexBundle->getGlobalOffsets()[mesh->commSize()] :
    mesh->getTopology()->depthStratum(0)->size();
  const int numDims = coords_field->getFiberDimension(patch, 
			      *mesh->getTopology()->depthStratum(0)->begin());

  fileout
    << "  vertices = {\n"
    << "    dimension = " << numDims << "\n"
    << "    count = " << numVertices << "\n"
    << "    coordinates = {\n";
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    fileout << "      ";
    
    fileout
      << std::resetiosflags(std::ios::fixed)
      << std::setiosflags(std::ios::scientific)
      << std::setprecision(6);
    for (int iDim=0; iDim < numDims; ++iDim) 
      fileout << std::setw(18) << coordinates[numDims*iVertex+iDim];
    fileout << "\n";
  } // for
  fileout
    << "    }\n"
    << "  }\n";
} // _writeVertices
  
// ----------------------------------------------------------------------
// Read mesh elements.
void
pylith::meshio::MeshIOAscii::_readElements(std::istream& filein,
					   int** pElements,
					   int* pNumElements, 
					   int* pNumCorners) const
{ // _readElements
  int* elements = 0;
  int numElements = 0;
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
      filein >> numElements;
    } else if (0 == strcasecmp(token.c_str(), "simplices")) {
      const int size = numElements*numCorners;
      if (0 == size) {
	std::ostringstream msg;
	msg << "Tokens 'num-corners' and 'count' must precede 'elements'.";
	throw std::runtime_error(msg.str());
      } // if
      
      filein.ignore(maxIgnore, '{');
      delete[] elements; elements = new int[size];
      assert(0 != elements);
      for (int i=0; i < size; ++i)
	filein >> elements[i];
      if (!useIndexZero()) {
	// if files begins with index 1, then decrement to index 0
	// for compatibility with Sieve
	for (int i=0; i < size; ++i)
	  --elements[i];
      } // if
      
      filein.ignore(maxIgnore, '}');
    } else {
      std::ostringstream msg;
      msg << "Could not parse '" << token << "' into an elements setting.";
      throw std::runtime_error(msg.str());
    } // else
    filein >> token;
  } // while
  if (!filein.good())
    throw std::runtime_error("I/O error while parsing elements settings.");

  if (0 != pElements)
    *pElements = elements;
  if (0 != pNumElements)
    *pNumElements = numElements;
  if (0 != pNumCorners)
    *pNumCorners = numCorners;
} // _readElements

// ----------------------------------------------------------------------
// Write mesh elements.
void
pylith::meshio::MeshIOAscii::_writeElements(std::ostream& fileout,
				  const ALE::Obj<ALE::PetscMesh>& mesh) const
{ // _writeElements
  ALE::Obj<ALE::PetscMesh::sieve_type> topology = mesh->getTopology();
  ALE::Obj<ALE::PetscMesh::sieve_type::traits::heightSequence> elements = 
    topology->heightStratum(0);
  ALE::Obj<ALE::PetscMesh::bundle_type> vertexBundle =  mesh->getBundle(0);
  ALE::PetscMesh::bundle_type::patch_type patch;
  std::string orderName("element");

  assert(0 != topology);
  assert(0 != elements);
  assert(0 != vertexBundle);

  int numCorners = 
    topology->nCone(*elements->begin(), topology->depth())->size();
  const int numElements = mesh->getTopology()->heightStratum(0)->size();

  fileout
    << "  elements = {\n"
    << "    count = " << numElements << "\n"
    << "    num-corners = " << numCorners << "\n"
    << "    simplices = {\n";

  const int offset = (useIndexZero()) ? 0 : 1;
  for(ALE::PetscMesh::sieve_type::traits::heightSequence::iterator e_itor = 
	elements->begin(); 
      e_itor != elements->end();
      ++e_itor) {
    fileout << "      ";
    ALE::Obj<ALE::PetscMesh::bundle_type::order_type::coneSequence> cone = 
      vertexBundle->getPatch(orderName, *e_itor);
    assert(0 != cone);
    for(ALE::PetscMesh::bundle_type::order_type::coneSequence::iterator c_itor = 
	  cone->begin(); 
	c_itor != cone->end(); 
	++c_itor)
      fileout
	<< std::setw(8)
	<< offset + vertexBundle->getIndex(patch, *c_itor).prefix;
    fileout << "\n";
  } // for
  fileout
    << "    }\n"
    << "  }\n";
} // _writeElements

// ----------------------------------------------------------------------
// Read mesh charts.
void
pylith::meshio::MeshIOAscii::_readChart(std::istream& filein,
					ALE::Obj<ALE::PetscMesh>* pMesh) const
{ // _readChart
  std::string name = ""; // Name of chart
  int dimension = 0; // Topology dimension associated with chart
  int count = 0; // Number of entities in chart
  int* indices = 0; // Indices of entities in chart

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
      for (int i=0; i < count; ++i)
	filein >> indices[i];
      filein.ignore(maxIgnore, '}');
    } else {
      std::ostringstream msg;
      msg << "Could not parse '" << token << "' into a chart setting.";
      throw std::runtime_error(msg.str());
    } // else
    filein >> token;
  } // while
  if (!filein.good())
    throw std::runtime_error("I/O error while parsing chart settings.");

#if 0
  assert(0 != pMesh);
  ALE::Obj<ALE::PetscMesh::field_type> chartField = (*pMesh)->getField(name);
  const int meshDim = (*pMesh)->getDimension();
  ALE::Obj<std::list<ALE::PetscMesh::point_type> > patchPoints = 
    std::list<ALE::PetscMesh::point_type>();
  ALE::PetscMesh::field_type::patch_type patch;

  patchPoints->clear();
  if (meshDim == dimension) {
    for (int i=0; i < count; ++i)
      patchPoints->push_back(ALE::PetscMesh::point_type(0, indices[i]));
    chartField->setPatch(patchPoints, patch);
  } else if (0 == dimension) {
  } // if
  chartField->setFiberDimensionByHeight(patch, 0, 1);
  chartField->orderPatches();
  const double zero = 0;
  for (int i=0; i < count; ++i)
    chartField->update(patch, ALE::PetscMesh::point_type(0, i), &zero);
#endif
} // _readChart

// ----------------------------------------------------------------------
// Write mesh chart.
void
pylith::meshio::MeshIOAscii::_writeChart(std::ostream& fileout,
					 const ALE::Obj<ALE::PetscMesh>& mesh,
					 const char* name) const
{ // _writeChart
  //_writeChart(fileout, mesh);
  // ADD STUFF HERE
  int count = 0; // TEMPORARY
  int dimension = 0; // TEMPORARY

  fileout
    << "  chart = {\n"
    << "    name = " << name << "\n"
    << "    dimension = " << dimension << "\n"
    << "    count = " << count << "\n"
    << "    indices = {\n";
    
  fileout
    << "    }\n"
    << "  }\n";
} // _writeChart
  
// End of file 
