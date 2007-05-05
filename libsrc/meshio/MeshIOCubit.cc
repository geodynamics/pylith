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

#include "MeshIOCubit.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES double_array, int_array, string_vector

#include "journal/info.h" // USES journal::info_t

#include <netcdfcpp.h> // USES netcdf

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOCubit::MeshIOCubit(void) :
  _filename("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOCubit::~MeshIOCubit(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Unpickle mesh
void
pylith::meshio::MeshIOCubit::_read(void)
{ // _read
  try {
    NcFile ncfile(_filename.c_str());
    if (!ncfile.is_valid()) {
      std::ostringstream msg;
      msg << "Could not open Cubit Exodus file '" << _filename
	  << "' for reading.\n";
      throw std::runtime_error(msg.str());
    } // if

    int meshDim = 0;
    int spaceDim = 0;
    int numVertices = 0;
    int numCells = 0;
    int numCorners = 0;
    double_array coordinates;
    int_array cells;
    int_array materialIds;

    NcDim* num_dim = ncfile.get_dim("num_dim");
    if (0 == num_dim)
      throw std::runtime_error("Could not get dimension 'num_dim'.");
    meshDim = num_dim->size();

    _readVertices(ncfile, &coordinates, &numVertices, &spaceDim);
    _readCells(ncfile, &cells, &materialIds, &numCells, &numCorners);
    _buildMesh(coordinates, numVertices, spaceDim,
	       cells, numCells, numCorners, meshDim);
    _setMaterials(materialIds);

    _readGroups(ncfile);
  } catch (std::exception& err) {
    std::ostringstream msg;
    msg << "Error while reading Cubit Exodus file '" << _filename << "'.\n"
	<< err.what();
    throw std::runtime_error(msg.str());
  } catch (...) {
    std::ostringstream msg;
    msg << "Unknown error while reading Cubit Exodus file '" << _filename
	<< "'.";
    throw std::runtime_error(msg.str());
  } // try/catch
} // read

// ----------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIOCubit::_write(void) const
{ // write
  NcFile ncfile(_filename.c_str());
  if (!ncfile.is_valid()) {
    std::ostringstream msg;
    msg << "Could not open Cubit Exodus file '" << _filename
	<< "' for writing.\n";
    throw std::runtime_error(msg.str());
  } // if

  _writeDimensions(ncfile);
  _writeVariables(ncfile);
  _writeAttributes(ncfile);
} // write

// ----------------------------------------------------------------------
// Read mesh vertices.
void
pylith::meshio::MeshIOCubit::_readVertices(NcFile& ncfile,
					   double_array* coordinates,
					   int* numVertices, 
					   int* numDims) const
{ // _readVertices
  assert(0 != coordinates);
  assert(0 != numVertices);
  assert(0 != numDims);

  journal::info_t info("meshiocubit");
    
  NcDim* num_nodes = ncfile.get_dim("num_nodes");
  if (0 == num_nodes)
    throw std::runtime_error("Could not get dimension 'num_nodes'.");
  *numVertices = num_nodes->size();
  info << "Reading " << *numVertices << " vertices." << journal::endl;

  NcVar* coord = ncfile.get_var("coord");
  if (0 == coord)
    throw std::runtime_error("Could not get variable 'coord'.");
  if (2 != coord->num_dims())
    throw std::runtime_error("Number of dimensions of variable 'coord' "
			     "must be 2.");
  const int size = coord->num_vals();

  NcDim* space_dim = coord->get_dim(0);
  if (0 == space_dim)
    throw std::runtime_error("Could not get dimensions of coordinates.");
  *numDims = space_dim->size();

  assert((*numVertices)*(*numDims) == size);
  coordinates->resize(size);
  long* counts = coord->edges();
  bool ok = coord->get(&(*coordinates)[0], counts);
  delete[] counts; counts = 0;
  if (!ok)
    throw std::runtime_error("Could not get coordinate values.");
} // _readVertices

// ----------------------------------------------------------------------
// Read mesh cells.
void
pylith::meshio::MeshIOCubit::_readCells(NcFile& ncfile,
					int_array* cells,
					int_array* materialIds,
					int* numCells, 
					int* numCorners) const
{ // _readCells
  assert(0 != cells);
  assert(0 != materialIds);
  assert(0 != numCells);
  assert(0 != numCorners);

  journal::info_t info("meshiocubit");

  NcDim* num_elem = ncfile.get_dim("num_elem");
  if (0 == num_elem)
    throw std::runtime_error("Could not get dimension 'num_elem'.");
  *numCells = num_elem->size();
  NcDim* num_el_blk = ncfile.get_dim("num_el_blk");
  if (0 == num_el_blk)
    throw std::runtime_error("Could not get dimension 'num_el_blk'.");
  const int numMaterials = num_el_blk->size();

  info << "Reading " << numCells << " cells in " << numMaterials 
       << " blocks." << journal::endl;

  NcVar* eb_prop1 = ncfile.get_var("eb_prop1");
  if (0 == eb_prop1) 
    throw std::runtime_error("Could not get variable 'eb_prop1'.");
  std::valarray<int> blockIds(numMaterials);
  long* counts = eb_prop1->edges();
  bool ok = eb_prop1->get(&blockIds[0], counts);
  delete[] counts; counts = 0;
  materialIds->resize(*numCells);

  *numCorners = 0;
  for (int iMaterial=0, index=0; iMaterial < numMaterials; ++iMaterial) {
    std::ostringstream varname;
    varname << "num_nod_per_el" << iMaterial+1;
    NcDim* num_nod_per_el = ncfile.get_dim(varname.str().c_str());
    if (0 == num_nod_per_el)
      throw std::runtime_error("Could not get dimension 'num_nod_per_el'.");
    if (0 == *numCorners) {
      *numCorners = num_nod_per_el->size();
      const int size = (*numCells) * (*numCorners);
      cells->resize(size);
    } else
      assert(num_nod_per_el->size() == *numCorners);

    varname.str("");
    varname << "num_el_in_blk" << iMaterial+1;
    NcDim* num_el_in_blk = ncfile.get_dim(varname.str().c_str());
    if (0 == num_el_in_blk)
      throw std::runtime_error("Could not get dimension 'num_el_in_blk'.");
    const int blockSize = num_el_in_blk->size();
	
    varname.str("");
    varname << "connect" << iMaterial+1;
    NcVar* connect = ncfile.get_var(varname.str().c_str());
    if (0 == connect)
      throw std::runtime_error("Could not get variable 'connect'.");
    if (2 != connect->num_dims())
      throw std::runtime_error("Number of dimensions of variable "
			       "'connect' must be 2.");
    const int size = connect->num_vals();
    assert(blockSize * (*numCorners) == size);
    long* counts = connect->edges();
    bool ok = connect->get(&(*cells)[index * (*numCorners)], counts);
    delete[] counts; counts = 0;
    if (!ok)
      throw std::runtime_error("Could not get cell values.");
	
    for (int i=0; i < blockSize; ++i)
      materialIds[index+i] = blockIds[iMaterial];
    
    index += blockSize;
  } // for
  cells -= 1; // use zero index
} // _readCells

// ----------------------------------------------------------------------
// Read mesh groups.
void
pylith::meshio::MeshIOCubit::_readGroups(NcFile& ncfile)
{ // _readGroups
  journal::info_t info("meshiocubit");

  NcDim* num_node_sets = ncfile.get_dim("num_node_sets");
  if (0 == num_node_sets)
    throw std::runtime_error("Could not get dimension 'num_node_sets'.");
  const int numGroups = num_node_sets->size();
  info << "Found " << numGroups << " node sets." << journal::endl;
      
  NcVar* ns_prop1 = ncfile.get_var("ns_prop1");
  if (0 == ns_prop1) 
    throw std::runtime_error("Could not get variable 'ns_prop1'.");
  std::valarray<int> ids(numGroups);
  long* counts = ns_prop1->edges();
  bool ok = ns_prop1->get(&ids[0], counts);
  delete[] counts; counts = 0;
      
  for (int iGroup=0; iGroup < numGroups; ++iGroup) {
    std::valarray<int> points;
	
    std::ostringstream varname;
    varname << "node_ns" << iGroup+1;
    NcVar* node_ns = ncfile.get_var(varname.str().c_str());
    if (0 == node_ns)
      throw std::runtime_error("Could not get node set.");
    const int size = node_ns->num_vals();
    info << "Reading node set " << ids[iGroup] << " with "
	 << size << " nodes." << journal::endl;

    points.resize(size);
    long* counts = node_ns->edges();
    bool ok = node_ns->get(&points[0], counts);
    delete[] counts; counts = 0;
    if (!ok)
      throw std::runtime_error("Could not get node set.");
    std::sort(&points[0], &points[size]);

    GroupPtType type = VERTEX;
    std::string name = "";
    _setGroup(name, type, points);
  } // for  
} // _readGroups

// ----------------------------------------------------------------------
// Write mesh dimensions.
void
pylith::meshio::MeshIOCubit::_writeDimensions(NcFile& ncfile) const
{ // _writeDimensions
  throw std::logic_error("MeshIOCubit::_writeDimensions() not implemented.");
} // _writeDimensions
  
// ----------------------------------------------------------------------
// Write mesh variables.
void
pylith::meshio::MeshIOCubit::_writeVariables(NcFile& ncfile) const
{ // _writeVariables
  throw std::logic_error("MeshIOCubit::_writeVariables() not implemented.");
} // _writeVariables
  
// ----------------------------------------------------------------------
// Write mesh attributes.
void
pylith::meshio::MeshIOCubit::_writeAttributes(NcFile& ncfile) const
{ // _writeAttributes
  throw std::logic_error("MeshIOCubit::_writeAttributes() not implemented.");
} // _writeAttributes

  
// End of file 
