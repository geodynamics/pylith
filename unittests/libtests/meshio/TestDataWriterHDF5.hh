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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/meshio/TestDataWriterHDF5.hh
 *
 * @brief C++ TestDataWriterHDF5 object
 *
 * C++ unit testing for DataWriterHDF5.
 */

#if !defined(pylith_meshio_testdatawriterhdf5_hh)
#define pylith_meshio_testdatawriterhdf5_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterHDF5;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5
{ // class TestDataWriterHDF5

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Check HDF5 file against archived file.
   *
   * @param filename Name of file to check.
   */
  static
  void checkFile(const char* filename);

}; // class TestDataWriterHDF5

#endif // pylith_meshio_testdatawriterhdf5_hh


// End of file 
