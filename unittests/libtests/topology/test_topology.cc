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

#include <petsc.h>
#include <Python.h>

#include <cppunit/extensions/TestFactoryRegistry.h>

#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/TextOutputter.h>

#include <stdlib.h> // USES abort()

#include "journal/info.h"

int
main(int argc,
     char* argv[])
{ // main
  CppUnit::TestResultCollector result;

  try {
    // Initialize PETSc
    PetscErrorCode err = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    CHKERRQ(err);

    // Initialize Python
    Py_Initialize();

    journal::info_t info("gmvfile");
    //info.activate();

    // Create event manager and test controller
    CppUnit::TestResult controller;

    // Add listener to collect test results
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

    // Finalize Python
    Py_Finalize();

    // Finalize PETSc
    err = PetscFinalize();
    CHKERRQ(err);
  } catch (...) {
    abort();
  } // catch

  return (result.wasSuccessful() ? 0 : 1);
} // main

// End of file
