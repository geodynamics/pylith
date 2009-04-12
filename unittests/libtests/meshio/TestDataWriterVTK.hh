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

/**
 * @file unittests/libtests/meshio/TestDataWriterVTK.hh
 *
 * @brief C++ TestDataWriterVTK object
 *
 * C++ unit testing for DataWriterVTK.
 */

#if !defined(pylith_meshio_testdatawritervtk_hh)
#define pylith_meshio_testdatawritervtk_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterVTK;

    class DataWriterVTKData;
  } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTK : public CppUnit::TestFixture
{ // class TestDataWriterVTK

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

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
		 const double t,
		 const char* timeFormat);
  
  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  DataWriterVTKData* _data; ///< Data for testing

}; // class TestDataWriterVTK

#endif // pylith_meshio_testdatawritervtk_hh


// End of file 
