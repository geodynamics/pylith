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
 * @file unittests/libtests/utils/TestEventLogger.hh
 *
 * @brief C++ TestEventLogger object
 *
 * C++ unit testing for EventLogger.
 */

#if !defined(pylith_utils_testeventlogger_hh)
#define pylith_utils_testeventlogger_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace utils {
    class TestEventLogger;
  } // utils
} // pylith

/// C++ unit testing for TestEventLogger
class pylith::utils::TestEventLogger : public CppUnit::TestFixture
{ // class TestEventLogger

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestEventLogger );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testClassName );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testRegisterEvent );
  CPPUNIT_TEST( testEventId );
  CPPUNIT_TEST( testEventLogging );

  CPPUNIT_TEST_SUITE_END();

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test className().
  void testClassName(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test registerEvent().
  void testRegisterEvent(void);

  /// Test eventId().
  void testEventId(void);

  /// Test eventBegin() and eventEnd().
  void testEventLogging(void);

}; // class TestEventLogging

#endif // pylith_utils_testeventlogger_hh


// End of file 
