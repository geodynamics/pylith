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

#include <cppunit/extensions/TestFactoryRegistry.h>

#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/TextOutputter.h>

int
main(int argc,
     char* argv[])
{ // main
  // Create event manager and test controller
  CppUnit::TestResult controller;

  // Add listener to collect test results
  CppUnit::TestResultCollector result;
  controller.addListener(&result);

  // Add listener to show progress as tests run
  CppUnit::BriefTestProgressListener progress;
  controller.addListener(&progress);

  // Add top suite to test runner
  CppUnit::TestRunner runner;
  runner.addTest(CppUnit::TestFactoryRegistry::getRegistry().makeTest());
  runner.run(controller);

  // Print tests
  CppUnit::TextOutputter outputter(&result, std::cerr);
  outputter.write();

  return (result.wasSuccessful() ? 0 : 1);
} // main

// End of file
