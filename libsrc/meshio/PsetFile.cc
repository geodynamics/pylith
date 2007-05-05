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

#include "PsetFile.hh" // implementation of class methods

#include "PsetFileAscii.hh"

#include <fstream> // uses std::fstream
#include <sstream> // uses std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------
pylith::meshio::PsetFile::PsetFile(const char* filename) :
  _filename(filename)
{ // constructor
} // constructor

// ----------------------------------------------------------------
pylith::meshio::PsetFile::~PsetFile(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------
bool
pylith::meshio::PsetFile::isAscii(const char* filename)
{ // isAscii
  std::ifstream fin(filename);
  if (!(fin.is_open() && fin.good())) {
    std::ostringstream msg;
    msg << "Could not open Pset file '" << filename << "' for reading.";
    throw std::runtime_error(msg.str());
  } // if
  const int headerLen = strlen(PsetFileAscii::header())+1;
  char buffer[headerLen];
  fin.get(buffer, headerLen, '\n');
  fin.close();
  return (0 == strcmp(PsetFileAscii::header(), buffer)) ? true : false;
} // isAscii


// End of file 
