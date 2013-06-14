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
// Copyright (c) 2010-2013 University of California, Davis
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
// Test write() with tri3 mesh and vertex data.
void
pylith::meshio::TestXdmf::testWriteTri3Vertex(void)
{ // testWriteTri3Vertex

  const char* filenameHDF5 = "data/tri3_vertex.h5";
  const char* filenameXdmf = "tri3_vertex.xmf";

  Xdmf metafile;
  metafile.write(filenameXdmf, filenameHDF5);

  _checkFile(filenameXdmf);
} // testWriteTri3Vertex

// ----------------------------------------------------------------------
// Test write() with Tri3 mesh and cell data.
void
pylith::meshio::TestXdmf::testWriteTri3Cell(void)
{ // testWriteTri3Cell
  const char* filenameHDF5 = "data/tri3_cell.h5";
  const char* filenameXdmf = "tri3_cell.xmf";

  Xdmf metafile;
  metafile.write(filenameXdmf, filenameHDF5);

  _checkFile(filenameXdmf);
} // testWriteTri3Cell

// ----------------------------------------------------------------------
// Test write() with quad4 mesh and vertex data.
void
pylith::meshio::TestXdmf::testWriteQuad4Vertex(void)
{ // testWriteQuad4Vertex

  const char* filenameHDF5 = "data/quad4_vertex.h5";
  const char* filenameXdmf = "quad4_vertex.xmf";

  Xdmf metafile;
  metafile.write(filenameXdmf, filenameHDF5);

  _checkFile(filenameXdmf);
} // testWriteQuad4Vertex

// ----------------------------------------------------------------------
// Test write() with Quad4 mesh and cell data.
void
pylith::meshio::TestXdmf::testWriteQuad4Cell(void)
{ // testWriteQuad4Cell
  const char* filenameHDF5 = "data/quad4_cell.h5";
  const char* filenameXdmf = "quad4_cell.xmf";

  Xdmf metafile;
  metafile.write(filenameXdmf, filenameHDF5);

  _checkFile(filenameXdmf);
} // testWriteQuad4Cell

// ----------------------------------------------------------------------
// Test write() with tet4 mesh and vertex data.
void
pylith::meshio::TestXdmf::testWriteTet4Vertex(void)
{ // testWriteTet4Vertex
  const char* filenameHDF5 = "data/tet4_vertex.h5";
  const char* filenameXdmf = "tet4_vertex.xmf";

  Xdmf metafile;
  metafile.write(filenameXdmf, filenameHDF5);

  _checkFile(filenameXdmf);
} // testWriteTet4Vertex

// ----------------------------------------------------------------------
// Test write() with tet4 mesh and cell data.
void
pylith::meshio::TestXdmf::testWriteTet4Cell(void)
{ // testWriteTet4Cell
  const char* filenameHDF5 = "data/tet4_cell.h5";
  const char* filenameXdmf = "tet4_cell.xmf";

  Xdmf metafile;
  metafile.write(filenameXdmf, filenameHDF5);

  _checkFile(filenameXdmf);
} // testWriteTet4Cell

// ----------------------------------------------------------------------
// Test write() with hex8 mesh and vertex data.
void
pylith::meshio::TestXdmf::testWriteHex8Vertex(void)
{ // testWriteHex8Vertex
  const char* filenameHDF5 = "data/hex8_vertex.h5";
  const char* filenameXdmf = "hex8_vertex.xmf";

  Xdmf metafile;
  metafile.write(filenameXdmf, filenameHDF5);

  _checkFile(filenameXdmf);
} // testWriteHex8Vertex

// ----------------------------------------------------------------------
// Test write() with hex8 mesh and cell data.
void
pylith::meshio::TestXdmf::testWriteHex8Cell(void)
{ // testWriteHex8Cell
  const char* filenameHDF5 = "data/hex8_cell.h5";
  const char* filenameXdmf = "hex8_cell.xmf";

  Xdmf metafile;
  metafile.write(filenameXdmf, filenameHDF5);

  _checkFile(filenameXdmf);
} // testWriteHex8Cell

// ----------------------------------------------------------------------
// Check Xdmf file against archived file.
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
