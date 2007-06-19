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

#include "PsetFileAscii.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES double_array, int_array

#include "journal/info.h" // USES journal::info_t

#include <fstream> // USES std::ifstream
#include <iomanip> // USES std::setw()
#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::exception

// ----------------------------------------------------------------------
const char* pylith::meshio::PsetFileAscii::_HEADER = "pset ascii";

// ----------------------------------------------------------------------
// Constructor with name of Pset file.
pylith::meshio::PsetFileAscii::PsetFileAscii(const char* filename) :
  PsetFile(filename)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor 
pylith::meshio::PsetFileAscii::~PsetFileAscii(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Read ASCII Pset file.
void
pylith::meshio::PsetFileAscii::read(std::vector<Pset>* groups)
{ // read
  assert(0 != groups);

  journal::info_t info("psetfile");

  std::ifstream fin(_filename.c_str(), std::ios::in);
  if (!(fin.is_open() && fin.good())) {
    std::ostringstream msg;
    msg
      << "Could not open ASCII Pset file '" << _filename
      << "' for reading.";
    throw std::runtime_error(msg.str());
  } // if
    
  info << journal::at(__HERE__)
       << "Reading ASCII Pset file '" << _filename << "'." << journal::endl;

  _readHeader(fin);

  // Read number of psets
  int numGroups = 0;
  fin >> numGroups;
  groups->resize(numGroups);

  // Read groups
  info << journal::at(__HERE__)
       << "Reading " << numGroups << " point sets from file." << journal::endl;
  for (int iGroup=0; iGroup < numGroups; ++iGroup)
    _readPset(fin, &(*groups)[iGroup]);
} // read

// ----------------------------------------------------------------------
// Write ASCII Pset file.
void
pylith::meshio::PsetFileAscii::write(const std::vector<Pset>& groups)
{ // write
  journal::info_t info("psetfile");

  std::ofstream fout(_filename.c_str(), std::ios::out);
  if (!(fout.is_open() && fout.good())) {
    std::ostringstream msg;
    msg
      << "Could not open ASCII Pset file '" << _filename
      << "' for writing.";
    throw std::runtime_error(msg.str());
  } // if
    
  info << journal::at(__HERE__)
       << "Writing ASCII Pset file '" << _filename << "'." << journal::endl;

  _writeHeader(fout);

  // Write number of groups
  const int numGroups = groups.size();
  fout << std::setw(4) << numGroups << std::endl;

  // Write groups
  info << journal::at(__HERE__)
       << "Writing " << numGroups << " point sets to file." << journal::endl;
  for (int iGroup=0; iGroup < numGroups; ++iGroup)
    _writePset(fout, groups[iGroup]);
} // write

// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileAscii::_readHeader(std::ifstream& fin)
{ // _readHeader
  const int headerLen = strlen(_HEADER)+1;
  char buffer[headerLen];
  fin.get(buffer, headerLen, '\n');
  if (0 != strcmp(_HEADER, buffer)) {
    std::ostringstream msg;
    msg
      << "Header in ASCII Pset file '" << buffer
      << "' does not match anticipated header '"
      << _HEADER << "'";
    throw std::runtime_error(msg.str());
  } // if
} // _readHeader

// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileAscii::_writeHeader(std::ofstream& fout)
{ // _writeHeader
  fout << _HEADER << std::endl;
} // _writeHeader

// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileAscii::_readPset(std::ifstream& fin,
					 Pset* group)
{ // _readPset
  assert(0 != group);

  journal::info_t info("psetfile");

  int size = 0;
  fin >> group->name >> group->id >> size;
  info << journal::at(__HERE__)
       << "Reading point set '" << group->name << "' with " << size
       << " points." << journal::endl;

  group->points.resize(size);
  for (int i=0; i < size; ++i)
    fin >> group->points[i];

  group->points -= 1; // use zero base

  info << journal::at(__HERE__)
       << "Done." << journal::endl;
} // _readPset

// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileAscii::_writePset(std::ofstream& fout,
					  const Pset& group)
{ // _writePset
  journal::info_t info("psetfile");
  const int size = group.points.size();

  info << journal::at(__HERE__)
       << "Writing point set '" << group.name << "' with " << size
       << " points." << journal::endl;

  fout << group.name << " " << group.id << "  " << size << std::endl;

  const int numCols = 10;
  for (int i=0, iCol=0; i < size; ++i) {
    fout << std::setw(8) << 1+group.points[i];
    if (++iCol == numCols) {
      fout << std::endl;
      iCol = 0;
    } // if
  } // for

  info << journal::at(__HERE__)
       << "Done." << journal::endl;
} // _writePset


// End of file 
