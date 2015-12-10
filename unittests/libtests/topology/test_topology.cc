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

#define MALLOC_DUMP

int
main(int argc,
     char* argv[])
{ // main
  CppUnit::TestResultCollector result;

  try {
    // Initialize PETSc
    PetscErrorCode err = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(err);
#if defined(MALLOC_DUMP)
    err = PetscOptionsSetValue(NULL, "-malloc_dump", "");CHKERRQ(err);
#endif

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

    // Finalize PETSc
    err = PetscFinalize();
    CHKERRQ(err);
  } catch (...) {
    abort();
  } // catch

#if !defined(MALLOC_DUMP)
  std::cout << "WARNING -malloc dump is OFF\n" << std::endl;
#endif

  return (result.wasSuccessful() ? 0 : 1);
} // main

// End of file
