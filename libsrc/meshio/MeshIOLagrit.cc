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

#include "MeshIOLagrit.hh" // implementation of class methods

#include "GMVFileAscii.hh" // USES GMVFileAscii
#include "GMVFileBinary.hh" // USES GMVFileBinary
#include "PsetFileAscii.hh" // USES PsetFileAscii
#include "PsetFileBinary.hh" // USES PsetFileBinary
#include "MeshBuilder.hh" // USES MeshBuilder

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/array.hh" // USES double_array, int_array

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error()

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOLagrit::MeshIOLagrit(void) :
  _filenameGmv(""),
  _filenamePset(""),
  _flipEndian(false),
  _ioInt32(true),
  _isRecordHeader32Bit(true)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOLagrit::~MeshIOLagrit(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Unpickle mesh
void
pylith::meshio::MeshIOLagrit::_read(void)
{ // _read
  MPI_Comm comm = _mesh->comm();
  int rank;
  int meshDim = 0;
  int spaceDim = 0;
  int numVertices = 0;
  int numCells = 0;
  int numCorners = 0;
  double_array coordinates;
  int_array cells;
  int_array materialIds;

  MPI_Comm_rank(comm, &rank);
  if (!rank) {
    if (GMVFile::isAscii(_filenameGmv.c_str())) {
      GMVFileAscii filein(_filenameGmv.c_str());
      filein.read(&coordinates, &cells, &materialIds, 
                  &meshDim, &spaceDim, &numVertices, &numCells, &numCorners);
      _orientCellsAscii(&cells, numCells, numCorners, meshDim);
    } else {
      GMVFileBinary filein(_filenameGmv.c_str(), _flipEndian);
      filein.read(&coordinates, &cells, &materialIds, 
                  &meshDim, &spaceDim, &numVertices, &numCells, &numCorners);
      _orientCellsBinary(&cells, numCells, numCorners, meshDim);
    } // if/else
  }
  MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim,
			 cells, numCells, numCorners, meshDim,
			 _interpolate, *_normalizer);
  _setMaterials(materialIds);

  if (0 == rank) {
    std::vector<PsetFile::Pset> groups;
    if (PsetFile::isAscii(_filenamePset.c_str())) {
      PsetFileAscii filein(_filenamePset.c_str());
      filein.read(&groups);
    } else {
      PsetFileBinary filein(_filenamePset.c_str(), 
			    _flipEndian,
			    _ioInt32,
			    _isRecordHeader32Bit);
      filein.read(&groups);
    } // if/else
    GroupPtType type = VERTEX;
    const int numGroups = groups.size();
    for (int iGroup=0; iGroup < numGroups; ++iGroup)
      _setGroup(groups[iGroup].name, type, groups[iGroup].points);
  }
  _distributeGroups();
} // _read

// ----------------------------------------------------------------------
// Pickle mesh
void
pylith::meshio::MeshIOLagrit::_write(void) const
{ // _write
  throw std::logic_error("MeshIOLagrit::_write not implemented.");
} // _write

// ----------------------------------------------------------------------
// Reorder vertices in cells from ASCII GMV file to match PyLith
// conventions.
void
pylith::meshio::MeshIOLagrit::_orientCellsAscii(int_array* const cells,
						const int numCells,
						const int numCorners,
						const int meshDim)
{ // _orientCellsAscii
  assert(0 != cells);
  assert(cells->size() == numCells*numCorners);

  if (3 == meshDim && 4 == numCorners) // TET
    for (int iCell=0; iCell < numCells; ++iCell) {
      const int i1 = iCell*numCorners+1;
      const int i2 = iCell*numCorners+2;
      const int tmp = (*cells)[i1];
      (*cells)[i1] = (*cells)[i2];
      (*cells)[i2] = tmp;
    } // for
} // _orientCellsAscii
  
// ----------------------------------------------------------------------
// Reorder vertices in cells from binary GMV file to match PyLith
// conventions.
void
pylith::meshio::MeshIOLagrit::_orientCellsBinary(int_array* const cells,
						 const int numCells,
						 const int numCorners,
						 const int meshDim)
{ // _orientCellsBinary
  assert(0 != cells);
  assert(cells->size() == numCells*numCorners);

  if (3 == meshDim && 4 == numCorners)  // TET
    ; // do nothing
} // _orientCellsBinary
  

// End of file 
