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

#include "GMVFileAscii.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES double_array, int_array

#include "journal/info.h" // USES journal::info_t

#include <fstream> // USES std::ifstream
#include <iomanip> // USES std::setw()
#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::exception

// ----------------------------------------------------------------------
const char* pylith::meshio::GMVFileAscii::_HEADER = "gmvinput ascii";

// ----------------------------------------------------------------------
// Constructor with name of GMV file.
pylith::meshio::GMVFileAscii::GMVFileAscii(const char* filename) :
  GMVFile(filename)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor 
pylith::meshio::GMVFileAscii::~GMVFileAscii(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Read ASCII GMV file.
void
pylith::meshio::GMVFileAscii::read(double_array* coordinates,
				   int_array* cells,
				   int_array* materialIds,
				   int* meshDim,
				   int* spaceDim,
				   int* numVertices,
				   int* numCells,
				   int* numCorners)
{ // read
  assert(0 != coordinates);
  assert(0 != cells);
  assert(0 != materialIds);
  assert(0 != meshDim);
  assert(0 != spaceDim);
  assert(0 != numVertices);
  assert(0 != numCells);
  assert(0 != numCorners);

  *meshDim = 3;

  journal::info_t info("gmvfile");

  std::ifstream fin(_filename.c_str(), std::ios::in);
  if (!(fin.is_open() && fin.good())) {
    std::ostringstream msg;
    msg
      << "Could not open ASCII GMV file '" << _filename
      << "' for reading.";
    throw std::runtime_error(msg.str());
  } // if
    
  info << journal::at(__HERE__)
       << "Reading ASCII GMV file '" << _filename << "'." << journal::endl;

  _readHeader(fin);

  std::string token;
  while (fin >> token) {
    if (token == "nodes")
      _readVertices(fin, coordinates, numVertices, spaceDim);
    else if (token == "cells")
      _readCells(fin, cells, numCells, numCorners);
    else if (token == "variable")
      _readVariables(fin, *numVertices, *numCells);
    else if (token == "flags")
      _readFlags(fin, *numVertices, *numCells);
    else if (token == "material")
      _readMaterials(fin, materialIds, *numVertices, *numCells);
  } // while

  assert(coordinates->size() == (*numVertices) * (*spaceDim));
  assert(cells->size() == (*numCells) * (*numCorners));
  assert(materialIds->size() == *numCells);
} // read

// ----------------------------------------------------------------------
// Write ASCII GMV file.
void
pylith::meshio::GMVFileAscii::write(const double_array& coordinates,
				    const int_array& cells,
				    const int_array& materialIds,
				    const int meshDim,
				    const int spaceDim,
				    const int numVertices,
				    const int numCells,
				    const int numCorners)
{ // write
  assert(coordinates.size() == numVertices * spaceDim);
  assert(cells.size() == numCells * numCorners);
  assert(materialIds.size() == numCells);

#if 0
  _writeHeader();
  _writeVertices(coordinates);
  _writeCells(cells);
  _writeMaterials(materialIds);
#endif
} // write

// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileAscii::_readHeader(std::ifstream& fin)
{ // _readHeader
  const int headerLen = strlen(_HEADER)+1;
  char buffer[headerLen];
  fin.get(buffer, headerLen, '\n');
  if (0 != strcmp(_HEADER, buffer)) {
    std::ostringstream msg;
    msg
      << "Header in ASCII GMV file '" << buffer
      << "' does not match anticipated header '"
      << _HEADER << "'";
    throw std::runtime_error(msg.str());
  } // if
} // _readHeader

// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileAscii::_readVertices(std::ifstream& fin,
					    double_array* coordinates,
					    int* numVertices,
					    int* spaceDim)
{ // _readVertices
  assert(0 != coordinates);
  assert(0 != numVertices);
  assert(0 != spaceDim);

  *spaceDim = 3;

  journal::info_t info("gmvfile");

  fin >> *numVertices;
  info << journal::at(__HERE__)
       << "Reading " << *numVertices << " nodes." << journal::endl;

  coordinates->resize(*numVertices * (*spaceDim));
  // NOTE: Order of loops is different than what we usually have
  for (int iDim=0; iDim < *spaceDim; ++iDim)
    for (int iVertex=0; iVertex < *numVertices; ++iVertex)
      fin >> (*coordinates)[iVertex*(*spaceDim)+iDim];

  info << journal::at(__HERE__)
       << "Done." << journal::endl;
} // readVertices

#include <iostream>
// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileAscii::_readCells(std::ifstream& fin,
					 int_array* cells,
					 int* numCells,
					 int* numCorners)
{ // readCells
  assert(0 != cells);
  assert(0 != numCells);
  assert(0 != numCorners);

  journal::info_t info("gmvfile");

  *numCorners = 0;

  fin >> *numCells;
  std::string cellString = "";
  info << journal::at(__HERE__)
       << "Reading " << numCells << " cells." << journal::endl;
  for (int iCell=0; iCell < *numCells; ++iCell) {
    std::string cellStringCur;
    int numCornersCur = 0;
    fin >> cellStringCur;
    fin >> numCornersCur;
    if (0 != *numCorners) {
      if (cellStringCur != cellString) {
	std::ostringstream msg;
	msg 
	  << "Mutiple element types not supported. Found element types '"
	  << cellString << "' and '" << cellStringCur << "' in GMV file '"
	  << _filename << "'.";
	throw std::runtime_error(msg.str());
      } // if
      assert(*numCorners == numCornersCur);
    } else {
      cellString = cellStringCur;
      *numCorners = numCornersCur;
      cells->resize((*numCells) * (*numCorners));
    } // if/else
    for (int iCorner=0; iCorner < *numCorners; ++iCorner)
      fin >> (*cells)[iCell*(*numCorners)+iCorner];
    fin >> std::ws;
  } // for

  *cells -= 1; // use zero base

  info << journal::at(__HERE__)
       << "Done." << journal::endl;
} // readCells

// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileAscii::_readVariables(std::ifstream& fin,
					     const int numVertices,
					     const int numCells)
{ // _readVariables
  journal::info_t info("gmvfile");

  info << journal::at(__HERE__)
       << "Reading variables..." << journal::endl;

  std::string varName;
  fin >> varName;
  while("endvars" != varName && !fin.eof() && fin.good()) {
    int varType = 0;
    fin >> varType;
    if (1 == varType) { // variables/attributes associated with vertices
      const int numVars = 1;
      double_array vals(numVertices*numVars);
      for (int iVertex=0; iVertex < numVertices; ++iVertex)
	fin >> vals[iVertex];
    } else { // variables/attributes associated with cells
      const int numVars = 1;
      double_array vals(numCells*numVars);
      for (int iCell=0; iCell < numCells; ++iCell)
	fin >> vals[iCell];
    } // else
    fin >> varName;
  } // while

  info << journal::at(__HERE__)
       << "Done." << journal::endl;
} // _readVariables

// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileAscii::_readFlags(std::ifstream& fin,
					 const int numVertices,
					 const int numCells)
{ // _readFlags
  journal::info_t info("gmvfile");

  info << journal::at(__HERE__)
       << "Reading flags..." << journal::endl;

  std::string varName;
  fin >> varName;
  while("endflag" != varName && !fin.eof() && fin.good()) {
    int varType = 0;
    int numFlags = 0;
    fin >> numFlags >> varType;
    for (int iFlag=0; iFlag < numFlags; ++iFlag) {
      std::string flagName;
      fin >> flagName;
    } // for
    if (1 == varType) { // flag associated with vertices
      int flag;
      for (int iVertex=0; iVertex < numVertices; ++iVertex)
	fin >> flag;
    } else { // flag associated with cells
      int flag;
      for (int iCell=0; iCell < numCells; ++iCell)
	fin >> flag;
    } // else
    fin >> varName;
  } // while

  info << journal::at(__HERE__)
       << "Done." << journal::endl;
} // _readFlags

// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileAscii::_readMaterials(std::ifstream& fin,
					     int_array* materialIds,
					     const int numVertices,
					     const int numCells)
{ // _readMaterials
  assert(0 != materialIds);

  journal::info_t info("gmvfile");
  info << journal::at(__HERE__)
       << "Reading materials..." << journal::endl;

  int numMaterials = 0;
  int dataType = 0;
  fin >> numMaterials >> dataType;

  std::string name;
  for (int iMat=0; iMat < numMaterials; ++iMat)
    fin >> name;

  if (0 == dataType) { // material associated with elements
    materialIds->resize(numCells);
    for (int iCell=0; iCell < numCells; ++iCell)
      fin >> (*materialIds)[iCell];
  } else { // material associated with nodes
    int_array materials(numVertices);
    for (int iVertex=0; iVertex < numVertices; ++iVertex)
      fin >> materials[iVertex];
  } // else

  info << journal::at(__HERE__)
       << "Done." << journal::endl;
} // _readMaterials


// End of file 
