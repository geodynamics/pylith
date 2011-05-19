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

#include "MeshIOCubit.hh" // implementation of class methods

#include "MeshBuilder.hh" // USES MeshBuilder

#include "pylith/utils/array.hh" // USES double_array, int_array, string_vector

#include "journal/info.h" // USES journal::info_t

#include <petsc.h> // USES MPI_Comm

#include <netcdfcpp.h> // USES netcdf

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOCubit::MeshIOCubit(void) :
  _filename(""),
  _useNodesetNames(true)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOCubit::~MeshIOCubit(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::MeshIOCubit::deallocate(void)
{ // deallocate
  MeshIO::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Unpickle mesh
void
pylith::meshio::MeshIOCubit::_read(void)
{ // _read
  assert(_mesh);

  int rank = 0;
  int meshDim = 0;
  int spaceDim = 0;
  int numVertices = 0;
  int numCells = 0;
  int numCorners = 0;
  double_array coordinates;
  int_array cells;
  int_array materialIds;

  MPI_Comm_rank(_mesh->comm(), &rank);
  if (0 == rank) {
    try {
      NcFile ncfile(_filename.c_str());
      if (!ncfile.is_valid()) {
        std::ostringstream msg;
        msg << "Could not open Cubit Exodus file '" << _filename
            << "' for reading.\n";
        throw std::runtime_error(msg.str());
      } // if


      NcDim* num_dim = ncfile.get_dim("num_dim");
      if (0 == num_dim)
        throw std::runtime_error("Could not get dimension 'num_dim'.");
      meshDim = num_dim->size();

      _readVertices(ncfile, &coordinates, &numVertices, &spaceDim);
      _readCells(ncfile, &cells, &materialIds, &numCells, &numCorners);
      _orientCells(&cells, numCells, numCorners, meshDim);
      MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim,
			     cells, numCells, numCorners, meshDim,
			     _interpolate);
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
  } else {
    MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim,
			   cells, numCells, numCorners, meshDim,
			   _interpolate);
    _setMaterials(materialIds);
  }
  _distributeGroups();
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
  info << journal::at(__HERE__)
       << "Reading " << *numVertices << " vertices." << journal::endl;

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
  long* counts = coord->edges();
  double_array buffer(size);
  bool ok = coord->get(&buffer[0], counts);
  delete[] counts; counts = 0;
  if (!ok)
    throw std::runtime_error("Could not get coordinate values.");

  coordinates->resize(size);
  for (int iVertex=0; iVertex < *numVertices; ++iVertex)
    for (int iDim=0; iDim < *numDims; ++iDim)
      (*coordinates)[iVertex*(*numDims)+iDim] = 
	buffer[iDim*(*numVertices)+iVertex];
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

  info << journal::at(__HERE__)
       << "Reading " << *numCells << " cells in " << numMaterials 
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
    } else if (num_nod_per_el->size() != *numCorners) {
      std::ostringstream msg;
      msg << "All materials must have the same number of vertices per cell.\n"
	  << "Expected " << *numCorners << " vertices per cell, but block "
	  << blockIds[iMaterial] << " has " << num_nod_per_el->size()
	  << " vertices.";
      throw std::runtime_error(msg.str());
    } // if

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
      (*materialIds)[index+i] = blockIds[iMaterial];
    
    index += blockSize;
  } // for

  *cells -= 1; // use zero index
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
  info << journal::at(__HERE__)
       << "Found " << numGroups << " node sets." << journal::endl;
      
  NcVar* ns_prop1 = ncfile.get_var("ns_prop1");
  if (0 == ns_prop1) 
    throw std::runtime_error("Could not get variable 'ns_prop1'.");
  std::valarray<int> ids(numGroups);
  long* counts = ns_prop1->edges();
  bool ok = ns_prop1->get(&ids[0], counts);
  delete[] counts; counts = 0;
      
  string_vector groupNames(numGroups);
  if (_useNodesetNames) {
    NcVar* ns_names = ncfile.get_var("ns_names");
    if (0 == ns_names) 
      throw std::runtime_error("Could not get variable 'ns_names'.");
    long* counts = ns_names->edges();
    if (ns_names->num_dims() != 2)
      throw std::runtime_error("Expected variable 'ns_names' to have "
			       "2 dimensions.");
    const int bufferSize = counts[1];
    char* buffer = (bufferSize > 0) ? new char[bufferSize] : 0;
    for (int i=0; i < numGroups; ++i) {
      ns_names->set_cur(i);
      bool ok = ns_names->get(buffer, 1, bufferSize, 0, 0, 0);
      if (!ok) {
	std::ostringstream msg;
	msg << "Could not read name of nodeset " << ids[i] << ".";
	throw std::runtime_error(msg.str());
      } // if
      groupNames[i] = buffer;
      //std::cout << "GROUP: '" << groupNames[i] << "'." << std::endl;
    } // for
    delete[] buffer; buffer = 0;
    delete[] counts; counts = 0;
  } // if

  for (int iGroup=0; iGroup < numGroups; ++iGroup) {
    std::valarray<int> points;
	
    std::ostringstream varname;
    varname << "node_ns" << iGroup+1;
    NcVar* node_ns = ncfile.get_var(varname.str().c_str());
    if (0 == node_ns)
      throw std::runtime_error("Could not get node set.");
    const int size = node_ns->num_vals();
    info << journal::at(__HERE__)
	 << "Reading node set '" << groupNames[iGroup] << "' with id "
	 << ids[iGroup] << " containing "
	 << size << " nodes." << journal::endl;

    points.resize(size);
    long* counts = node_ns->edges();
    bool ok = node_ns->get(&points[0], counts);
    delete[] counts; counts = 0;
    if (!ok)
      throw std::runtime_error("Could not get node set.");
    std::sort(&points[0], &points[size]);
    points -= 1; // use zero index

    GroupPtType type = VERTEX;
    if (_useNodesetNames)
      _setGroup(groupNames[iGroup], type, points);
    else {
      std::ostringstream name;
      name << ids[iGroup];
      _setGroup(name.str().c_str(), type, points);
    } // if/else
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

// ----------------------------------------------------------------------
// Reorder vertices in cells to match PyLith conventions.
void
pylith::meshio::MeshIOCubit::_orientCells(int_array* const cells,
					  const int numCells,
					  const int numCorners,
					  const int meshDim)
{ // _orientCells
  assert(0 != cells);
  assert(cells->size() == numCells*numCorners);

  if (2 == meshDim && 4 == numCorners) { // QUAD4
    ; // do nothing

  } else if (3 == meshDim && 8 == numCorners) { // HEX8
    ; // do nothing

  } else if (2 == meshDim && 6 == numCorners) { // TRI6
    // CUBIT
    // corners, 
    // bottom edges, middle edges, top edges

    // PyLith (Sieve)
    // bottom edge, right edge, left edge, corners

    // Permutation: 3, 4, 5, 0, 1, 2
    int tmp = 0;
    for (int iCell=0; iCell < numCells; ++iCell) {
      const int ii = iCell*numCorners;
      tmp = (*cells)[ii+0];
      (*cells)[ii+0] = (*cells)[ii+3];
      (*cells)[ii+3] = tmp;

      tmp = (*cells)[ii+1];
      (*cells)[ii+1] = (*cells)[ii+4];
      (*cells)[ii+4] = tmp;

      tmp = (*cells)[ii+2];
      (*cells)[ii+2] = (*cells)[ii+5];
      (*cells)[ii+5] = tmp;
    } // for

  } else if (3 == meshDim && 27 == numCorners) { // HEX27
    // CUBIT
    // corners, 
    // bottom edges, middle edges, top edges
    // interior
    // bottom/top, left/right, front/back

    // PyLith
    // corners, 
    // bottom edges, top edges, middle edges
    // left/right, front/back, bottom/top
    // interior
    int tmp = 0;
    for (int iCell=0; iCell < numCells; ++iCell) {
      const int i12 = iCell*numCorners+12;
      const int i13 = iCell*numCorners+13;
      const int i14 = iCell*numCorners+14; 
      const int i15 = iCell*numCorners+15; 
      const int i16 = iCell*numCorners+16; 
      const int i17 = iCell*numCorners+17; 
      const int i18 = iCell*numCorners+18; 
      const int i19 = iCell*numCorners+19; 
      const int i20 = iCell*numCorners+20; 
      const int i21 = iCell*numCorners+21; 
      const int i22 = iCell*numCorners+22; 
      const int i23 = iCell*numCorners+23; 
      const int i24 = iCell*numCorners+24; 
      const int i25 = iCell*numCorners+25; 
      const int i26 = iCell*numCorners+26; 

      tmp = (*cells)[i12];
      (*cells)[i12] = (*cells)[i16];
      (*cells)[i16] = tmp;

      tmp = (*cells)[i13];
      (*cells)[i13] = (*cells)[i17];
      (*cells)[i17] = tmp;

      tmp = (*cells)[i14];
      (*cells)[i14] = (*cells)[i18];
      (*cells)[i18] = tmp;

      tmp = (*cells)[i15];
      (*cells)[i15] = (*cells)[i19];
      (*cells)[i19] = tmp;

      tmp = (*cells)[i20];
      (*cells)[i20] = (*cells)[i23];
      (*cells)[i23] = (*cells)[i26];
      (*cells)[i26] = tmp;

      tmp = (*cells)[i21];
      (*cells)[i21] = (*cells)[i24];
      (*cells)[i24] = tmp;

      tmp = (*cells)[i22];
      (*cells)[i22] = (*cells)[i25];
      (*cells)[i25] = tmp;
    } // for
  } // if/else

} // _orientCells
  

// End of file 
