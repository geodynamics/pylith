// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestXdmf.hh" // Implementation of class methods

#include "pylith/meshio/Xdmf.hh" // USES Xdmf

#include <cstring> // USES strcmp()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestXdmf );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::meshio::TestXdmf::testConstructor(void)
{ // testConstructor
  Xdmf one;
} // testConstructor

// ----------------------------------------------------------------------
// Test write() with 2D mesh and vertex data.
void
pylith::meshio::TestXdmf::testWrite2DVertex(void)
{ // testWrite2DVertex

  const char* filenameHDF5 = "data/tri3_vertex.h5";
  const char* filenameXdmf = "tri3_vertex.xmf";

  Xdmf metafile;
  metafile.write(filenameXdmf, filenameHDF5);

  _checkFile(filenameXdmf);
} // testWrite2DVertex

// ----------------------------------------------------------------------
// Test write() with 2D mesh and cell data.
void
pylith::meshio::TestXdmf::testWrite2DCell(void)
{ // testWrite2DCell
  const char* filenameHDF5 = "data/tri3_cell.h5";
  const char* filenameXdmf = "tri3_cell.xmf";

  Xdmf metafile;
  metafile.write(filenameXdmf, filenameHDF5);

  _checkFile(filenameXdmf);
} // testWrite2DCell

// ----------------------------------------------------------------------
// Test write() with 3D mesh and vertex data.
void
pylith::meshio::TestXdmf::testWrite3DVertex(void)
{ // testWrite3DVertex
  const char* filenameHDF5 = "data/tet4_vertex.h5";
  const char* filenameXdmf = "tet4_vertex.xmf";

  Xdmf metafile;
  metafile.write(filenameXdmf, filenameHDF5);

  _checkFile(filenameXdmf);
} // testWrite3DVertex

// ----------------------------------------------------------------------
// Test write() with 3D mesh and cell data.
void
pylith::meshio::TestXdmf::testWrite3DCell(void)
{ // testWrite3DCell
  const char* filenameHDF5 = "data/tet4_cell.h5";
  const char* filenameXdmf = "tet4_cell.xmf";

  Xdmf metafile;
  metafile.write(filenameXdmf, filenameHDF5);

  _checkFile(filenameXdmf);
} // testWrite3DCell

// ----------------------------------------------------------------------
// Check VTK file against archived file.
void
pylith::meshio::TestXdmf::_checkFile(const char* filename)
{ // _checkFile
  const std::string filenameE = std::string("data/") + std::string(filename);

  std::ifstream fileInE(filenameE.c_str());
  CPPUNIT_ASSERT(fileInE.is_open());

  std::ifstream fileIn(filename);
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
} // _checkFile


// End of file 
