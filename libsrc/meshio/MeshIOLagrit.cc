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

#include "pylith/utils/array.hh" // USES double_array, int_array

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

#if 0
  if (PsetFile::_isPsetAscii(_filenamePset)) {
    PsetFileAscii filein(_filenamePset);
    filein.read();
  } else {
    PsetFileBinary filein(_filenamePset);
    filein.read();
  } // if/else
#endif
} // _read

  
// End of file 
