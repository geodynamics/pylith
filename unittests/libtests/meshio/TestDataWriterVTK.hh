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
 * @file unittests/libtests/meshio/TestDataWriterVTK.hh
 *
 * @brief C++ TestDataWriterVTK object
 *
 * C++ unit testing for DataWriterVTK.
 */

#if !defined(pylith_meshio_testdatawritervtk_hh)
#define pylith_meshio_testdatawritervtk_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTK;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTK
{ // class TestDataWriterVTK

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Check VTK file against archived file.
   *
   * @param filename Name of file to check.
   * @param t Time for file.
   * @param timeFormat Format of timestamp in filename.
   */
  static
  void checkFile(const char* filename,
		 const PylithScalar t,
		 const char* timeFormat);
  
}; // class TestDataWriterVTK

#endif // pylith_meshio_testdatawritervtk_hh


// End of file 
