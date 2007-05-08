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

#include "PsetFileBinary.hh" // implementation of class methods

#include "BinaryIO.hh" // USES readString()

#include "pylith/utils/array.hh" // USES double_array, int_array

#include "journal/info.h" // USES journal::info_t

#include <fstream> // USES std::ifstream
#include <iomanip> // USES std::setw()
#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::exception

// ----------------------------------------------------------------------
const char* pylith::meshio::PsetFileBinary::_HEADER = "pset binary";

// ----------------------------------------------------------------------
// Constructor with name of Pset file.
pylith::meshio::PsetFileBinary::PsetFileBinary(const char* filename,
					       const bool flipEndian) :
  PsetFile(filename),
  _flipEndian(flipEndian)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor 
pylith::meshio::PsetFileBinary::~PsetFileBinary(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Read binary Pset file.
void
pylith::meshio::PsetFileBinary::read(std::vector<Pset>* groups)
{ // read
  assert(0 != groups);

  journal::info_t info("psetfile");

  std::ifstream fin(_filename.c_str(), std::ios::in);
  if (!(fin.is_open() && fin.good())) {
    std::ostringstream msg;
    msg
      << "Could not open binary Pset file '" << _filename
      << "' for reading.";
    throw std::runtime_error(msg.str());
  } // if
    
  info << "Reading binary Pset file '" << _filename << "'." << journal::endl;

  _readHeader(fin);

  // Read number of psets
  int numGroups = 0;
  fin.read((char*) numGroups, sizeof(int));
  if (_flipEndian)
    BinaryIO::swapByteOrder((char*) &numGroups, 1, sizeof(int));
  assert(numGroups >= 0);
  groups->resize(numGroups);

  // Read groups
  info << "Reading " << numGroups << " point sets from file." << journal::endl;
  for (int iGroup=0; iGroup < numGroups; ++iGroup)
    _readPset(fin, &(*groups)[iGroup]);
} // read

// ----------------------------------------------------------------------
// Write binary Pset file.
void
pylith::meshio::PsetFileBinary::write(const std::vector<Pset>& groups)
{ // write
  journal::info_t info("psetfile");

  std::ofstream fout(_filename.c_str(), std::ios::out);
  if (!(fout.is_open() && fout.good())) {
    std::ostringstream msg;
    msg
      << "Could not open binary Pset file '" << _filename
      << "' for writing.";
    throw std::runtime_error(msg.str());
  } // if
    
  info << "Writing binary Pset file '" << _filename << "'." << journal::endl;

  _writeHeader(fout);

  // Write number of groups
  int numGroups = groups.size();
  if (_flipEndian)
    BinaryIO::swapByteOrder((char*) &numGroups, 1, sizeof(numGroups));
  fout.write((char*) &numGroups, sizeof(numGroups));

  // Write groups
  info << "Writing " << numGroups << " point sets to file." << journal::endl;
  for (int iGroup=0; iGroup < numGroups; ++iGroup)
    _writePset(fout, groups[iGroup]);
} // write

// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileBinary::_readHeader(std::ifstream& fin)
{ // _readHeader
  std::string header = BinaryIO::readString(fin, strlen(_HEADER));
  std::string headerE = _HEADER;
  headerE = headerE.substr(0, headerE.find_first_of(" "));
  if (headerE != header) {
    std::ostringstream msg;
    msg
      << "Header in binary Pset file '" << header
      << "' does not match anticipated header '" << headerE << "'.";
    throw std::runtime_error(msg.str());
  } // if
} // _readHeader

// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileBinary::_writeHeader(std::ofstream& fout)
{ // _writeHeader
  fout.write((char*) _HEADER, strlen(_HEADER));
} // _writeHeader

// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileBinary::_readPset(std::ifstream& fin,
					 Pset* group)
{ // _readPset
  assert(0 != group);

  journal::info_t info("psetfile");

  group->name = BinaryIO::readString(fin, 32);

  int id = 0;
  fin.read((char*) &id, sizeof(int));
  if (_flipEndian)
    BinaryIO::swapByteOrder((char*) &id, 1, sizeof(id));

  int size = 0;
  fin.read((char*) &size, sizeof(int));
  if (_flipEndian)
    BinaryIO::swapByteOrder((char*) &size, 1, sizeof(size));
  assert(size >= 0);

  info << "Reading point set '" << group->name << "' with " << size
       << " points." << journal::endl;

  group->points.resize(size);
  fin.read((char*) &group->points[0], size*sizeof(int));
  if (_flipEndian)
    BinaryIO::swapByteOrder((char*) &group->points[0], size, sizeof(int));

  info << "Done." << journal::endl;
} // _readPset

// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileBinary::_writePset(std::ofstream& fout,
					   const Pset& group)
{ // _writePset
  journal::info_t info("psetfile");
  const int size = group.points.size();
  info << "Writing point set '" << group.name << "' with " << size
       << " points." << journal::endl;

  fout.write((char*) group.name.c_str(), 32);

  int id = group.id;
  if (_flipEndian)
    BinaryIO::swapByteOrder((char*) &id, 1, sizeof(id));
  fout.write((char*) &id, sizeof(int));

  int sizeIO = size;
  if (_flipEndian)
    BinaryIO::swapByteOrder((char*) &sizeIO, 1, sizeof(sizeIO));
  fout.write((char*) &sizeIO, sizeof(int));

  int_array pointsIO(group.points);
  if (_flipEndian)
    BinaryIO::swapByteOrder((char*) &pointsIO[0], size, sizeof(int));
  fout.write((char*) &pointsIO[0], size*sizeof(int));

  info << "Done." << journal::endl;
} // _writePset


// End of file 
