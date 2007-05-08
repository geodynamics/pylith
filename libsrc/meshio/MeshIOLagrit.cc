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

#include "pylith/utils/array.hh" // USES double_array, int_array

#include <stdexcept> // TEMPORARY

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOLagrit::MeshIOLagrit(void) :
  _filenameGmv(""),
  _filenamePset(""),
  _flipEndian(false)
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
  int meshDim = 0;
  int spaceDim = 0;
  int numVertices = 0;
  int numCells = 0;
  int numCorners = 0;
  double_array coordinates;
  int_array cells;
  int_array materialIds;

  if (GMVFile::isAscii(_filenameGmv.c_str())) {
    GMVFileAscii filein(_filenameGmv.c_str());
    filein.read(&coordinates, &cells, &materialIds, 
		&meshDim, &spaceDim, &numVertices, &numCells, &numCorners);
  } else {
    GMVFileBinary filein(_filenameGmv.c_str(), _flipEndian);
    filein.read(&coordinates, &cells, &materialIds, 
		&meshDim, &spaceDim, &numVertices, &numCells, &numCorners);
  } // if/else
  _buildMesh(coordinates, numVertices, spaceDim,
	     cells, numCells, numCorners, meshDim);
  _setMaterials(materialIds);

  std::vector<PsetFile::Pset> groups;
  if (PsetFile::isAscii(_filenamePset.c_str())) {
    PsetFileAscii filein(_filenamePset.c_str());
    filein.read(&groups);
  } else {
    PsetFileBinary filein(_filenamePset.c_str(), _flipEndian);
    filein.read(&groups);
  } // if/else
  GroupPtType type = VERTEX;
  const int numGroups = groups.size();
  for (int iGroup=0; iGroup < numGroups; ++iGroup)
    _setGroup(groups[iGroup].name, type, groups[iGroup].points);
} // _read

// ----------------------------------------------------------------------
// Pickle mesh
void
pylith::meshio::MeshIOLagrit::_write(void) const
{ // _write
  throw std::logic_error("MeshIOLagrit::_write not implemented.");
} // _write

  
// End of file 
