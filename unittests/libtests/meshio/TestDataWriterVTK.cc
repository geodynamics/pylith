// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterVTK.hh" // Implementation of class methods

#include "data/DataWriterVTKData.hh" // USES DataWriterVTKData

#include <string.h> // USES strcmp()

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTK::setUp(void)
{ // setUp
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTK::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Check VTK file against archived file.
void
pylith::meshio::TestDataWriterVTK::checkFile(const char* filenameRoot,
					     const double t,
					     const char* timeFormat)
{ // checkFile

  const std::string& fileroot = filenameRoot;

  std::ostringstream buffer;
  const int indexExt = fileroot.find(".vtk");
  // Add time stamp to filename
  char sbuffer[256];
  sprintf(sbuffer, timeFormat, t);
  std::string timestamp(sbuffer);
  const int pos = timestamp.find(".");
  if (pos != timestamp.length())
    timestamp.erase(pos, 1);
  buffer
    << std::string(fileroot, 0, indexExt) << "_t" << timestamp << ".vtk";
  
  const std::string& filename = buffer.str();
  const std::string filenameE = "data/" + filename;

  std::ifstream fileInE(filenameE.c_str());
  CPPUNIT_ASSERT(fileInE.is_open());

  std::ifstream fileIn(filename.c_str());
  CPPUNIT_ASSERT(fileIn.is_open());

  const int maxLen = 256;
  char line[maxLen];
  char lineE[maxLen];

  int i = 1;
  while(!fileInE.eof()) {
    fileInE.getline(lineE, maxLen);
    fileIn.getline(line, maxLen);
    if (0 != strcmp(line, lineE)) {
      std::cerr << "Line " << i << " of file '" << filename << "' is incorrect."
		<< std::endl;
      CPPUNIT_ASSERT(false);
    } // if
    ++i;
  } // while

  fileInE.close();
  fileIn.close();
} // checkFile


// End of file 
