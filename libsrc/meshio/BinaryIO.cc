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

#include "BinaryIO.hh" // implementation of class methods

#include <fstream> // USES std::ifstream
#include <assert.h> // USES assert()

#if defined(WORDS_BIGENDIAN)
#define NATIVE_BIG_ENDIAN
#else
#define NATIVE_LITTLE_ENDIAN
#endif

// ----------------------------------------------------------------------
// Read fixed length string from file.
std::string
pylith::meshio::BinaryIO::readString(std::ifstream& fin,
				     const int numChars)
{ // readString
  char* buffer = (numChars > 0) ? new char[numChars+1] : 0;
  fin.read(buffer, sizeof(char)*numChars);
  buffer[numChars] = '\0';

  // get string from buffer
  std::string bufstring(buffer);
  delete[] buffer; buffer = 0;

  // remove whitespace
  const int iLast = bufstring.find_first_of(" ");

  return bufstring.substr(0, iLast);
} // readString

// ----------------------------------------------------------------------
// Change endian type by swapping byte order.
void
pylith::meshio::BinaryIO::swapByteOrder(char* vals,
					const int numVals,
					const int typesize)
{ // swapByteOrder
  assert(0 != vals);
  const int numSwaps = sizeof(typesize) / 2;
  for (int iVal=0; iVal < numVals; ++iVal) {
    char* buf = (char*) (vals + iVal*typesize);
    for (int iSwap=0, jSwap=typesize-1; iSwap < numSwaps; ++iSwap, --jSwap) {
      char tmp = buf[iSwap];
      buf[iSwap] = buf[jSwap];
      buf[jSwap] = tmp;
    } // for
  } // for
} // swapByteOrder


// End of file 
