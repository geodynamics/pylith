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

#include <fstream> // USES std::ifstream, std::ofstream
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <assert.h> // USES assert()
#include <iomanip> // USES setw(), setiosflags(), resetiosflags()

// ----------------------------------------------------------------------
// Constructor
pylith::meshIO::MeshIOAscii::MeshIOAscii(void) :
  _filename("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshIO::MeshIOAscii::~MeshIOAscii(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Unpickle mesh
void
pylith::meshIO::MeshIOAscii::read(Obj<Mesh>& mesh)
{ // read
  int meshDim = 0;
  int numDims = 0;
  int numVertices = 0;
  int numCells = 0;
  int numCorners = 0;
  double* coordinates = 0;
  int* cells = 0;

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
    } else if (0 == strcasecmp(token.c_str(), "cells")) {
      filein.ignore(maxIgnore, '{');
      _readCells(filein, &cells, &numCells, &numCorners);
      readCells = true;
    } else if (0 == strcasecmp(token.c_str(), "chart")) {
      if (!builtTopology)
	throw std::runtime_error("Both 'vertices' and 'cells' must "
				 "precede any charts in mesh file.");
      filein.ignore(maxIgnore, '{');
      _readChart(filein, mesh);
    } else {
      std::ostringstream msg;
      msg << "Could not parse '" << token << "' into a mesh setting.";
      throw std::runtime_error(msg.str());  
    } // else

    if (readDim && readCells && readVertices && !builtTopology) {
      // Can now build topology
      mesh = Mesh(PETSC_COMM_WORLD, meshDim);
      mesh->debug = true;
      bool interpolate = false;

      // allow mesh to have different dimension than coordinates
      Obj<sieve_type>    sieve    = new sieve_type(mesh->comm(), mesh->debug);
      Obj<topology_type> topology = new topology_type(mesh->comm(), mesh->debug);

      ALE::New::SieveBuilder<sieve_type>::buildTopology(sieve, meshDim, numCells, cells, numVertices, interpolate, numCorners);
      sieve->stratify();
      topology->setPatch(0, sieve);
      topology->stratify();
      mesh->setTopologyNew(topology);
      ALE::New::SieveBuilder<sieve_type>::buildCoordinates(mesh->getSection("coordinates"), meshDim, coordinates);
      delete[] coordinates; coordinates = NULL;
      delete[] cells; cells = NULL;
      builtTopology = true;
    } // if

    filein >> token;
  } // while
  filein.close();
} // read

// ----------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshIO::MeshIOAscii::write(const Obj<Mesh>& mesh) const
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
  _writeCells(fileout, mesh);

  // LOOP OVER CHARTS
  // _writeChart(fileout, mesh, nameIter->c_str());

  fileout << "}\n";
  fileout.close();
} // write

// ----------------------------------------------------------------------
// Read mesh vertices.
void
pylith::meshIO::MeshIOAscii::_readVertices(std::istream& filein,
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
pylith::meshIO::MeshIOAscii::_writeVertices(std::ostream& fileout,
			       const Obj<Mesh>& mesh) const
{ // _writeVertices
  const Obj<Mesh::section_type>&       coords_field = mesh->getSection("coordinates");
  const Obj<Mesh::topology_type>&      topology     = mesh->getTopologyNew();
  const Mesh::section_type::patch_type patch        = 0;
  const Obj<Mesh::topology_type::label_sequence>& vertices = topology->depthStratum(patch, 0);
  const int                            embedDim     = coords_field->getFiberDimension(patch, *vertices->begin());

  fileout
    << "  vertices = {\n"
    << "    dimension = " << embedDim << "\n"
    << "    count = " << vertices->size() << "\n"
    << "    coordinates = {\n"
    << std::resetiosflags(std::ios::fixed)
    << std::setiosflags(std::ios::scientific)
    << std::setprecision(6);
  for(Mesh::topology_type::label_sequence::iterator v_iter = vertices->begin(); v_iter != vertices->end(); ++v_iter) {
    const Mesh::section_type::value_type *coordinates = coords_field->restrict(patch, *v_iter);

    fileout << "      ";
    for(int d = 0; d < embedDim; ++d) {
      fileout << std::setw(18) << coordinates[d];
    }
    fileout << "\n";
  } // for
  fileout
    << "    }\n"
    << "  }\n";
} // _writeVertices
  
// ----------------------------------------------------------------------
// Read mesh cells.
void
pylith::meshIO::MeshIOAscii::_readCells(std::istream& filein,
					   int** pCells,
					   int* pNumCells, 
					   int* pNumCorners) const
{ // _readCells
  int* cells = 0;
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
	std::ostringstream msg;
	msg << "Tokens 'num-corners' and 'count' must precede 'cells'.";
	throw std::runtime_error(msg.str());
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
    } else {
      std::ostringstream msg;
      msg << "Could not parse '" << token << "' into an cells setting.";
      throw std::runtime_error(msg.str());
    } // else
    filein >> token;
  } // while
  if (!filein.good())
    throw std::runtime_error("I/O error while parsing cells settings.");

  if (0 != pCells)
    *pCells = cells;
  if (0 != pNumCells)
    *pNumCells = numCells;
  if (0 != pNumCorners)
    *pNumCorners = numCorners;
} // _readCells

// ----------------------------------------------------------------------
// Write mesh cells.
void
pylith::meshIO::MeshIOAscii::_writeCells(std::ostream& fileout,
				  const Obj<Mesh>& mesh) const
{ // _writeCells
  const Obj<topology_type>&       topology   = mesh->getTopologyNew();
  const topology_type::patch_type patch      = 0;
  const Obj<sieve_type>&          sieve      = topology->getPatch(patch);
  const Obj<Mesh::topology_type::label_sequence>& cells = topology->heightStratum(patch, 0);
  const int                       numCorners = sieve->nCone(*cells->begin(), topology->depth())->size();

  fileout
    << "  cells = {\n"
    << "    count = " << cells->size() << "\n"
    << "    num-corners = " << numCorners << "\n"
    << "    simplices = {\n";

  const int offset = (useIndexZero()) ? 0 : 1;
  for(Mesh::topology_type::label_sequence::iterator e_iter = cells->begin(); e_iter != cells->end(); ++e_iter) {
    fileout << "      ";
    const Obj<sieve_type::traits::coneSequence>& cone = sieve->cone(*e_iter);

    for(sieve_type::traits::coneSequence::iterator c_iter = cone->begin(); c_iter != cone->end(); ++c_iter) {
      fileout << std::setw(8) << *c_iter + offset;
    }
    fileout << "\n";
  } // for
  fileout
    << "    }\n"
    << "  }\n";
} // _writeCells

// ----------------------------------------------------------------------
// Read mesh charts.
void
pylith::meshIO::MeshIOAscii::_readChart(std::istream& filein,
					const Obj<Mesh>& mesh) const
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
      for (int i = 0; i < count; ++i)
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
  assert(!mesh.isNull());
  Obj<Mesh::field_type> chartField = mesh->getField(name);
  const int meshDim = mesh->getDimension();
  Obj<std::list<Mesh::point_type> > patchPoints = 
    std::list<Mesh::point_type>();
  Mesh::field_type::patch_type patch;

  patchPoints->clear();
  if (meshDim == dimension) {
    for (int i=0; i < count; ++i)
      patchPoints->push_back(Mesh::point_type(0, indices[i]));
    chartField->setPatch(patchPoints, patch);
  } else if (0 == dimension) {
  } // if
  chartField->setFiberDimensionByHeight(patch, 0, 1);
  chartField->orderPatches();
  const double zero = 0;
  for (int i=0; i < count; ++i)
    chartField->update(patch, Mesh::point_type(0, i), &zero);
#endif
} // _readChart

// ----------------------------------------------------------------------
// Write mesh chart.
void
pylith::meshIO::MeshIOAscii::_writeChart(std::ostream& fileout,
					 const Obj<Mesh>& mesh,
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
