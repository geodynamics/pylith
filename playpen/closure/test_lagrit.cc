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

#include "petsc.h"
#include <Python.h>

#include "TestClosure.hh"

#include "pylith/meshio/MeshIOLagrit.hh"
#include "pylith/topology/Mesh.hh"

#include "spatialdata/geocoords/CSCart.hh"

#include <stdlib.h> // USES abort()

int
main(int argc,
     char* argv[])
{ // main
  if (3 != argc) {
    std::cerr << "Usage: test_lagrit filenameGMV filenamePset" << std::endl;
    return 1;
  } // if

  const char* filenameGmv = argv[1];
  const char* filenamePset = argv[2];

  try {
    // Initialize PETSc
    PetscErrorCode err = PetscInitialize(&argc, &argv,
					 PETSC_NULL, PETSC_NULL);
    CHKERRQ(err);
    err = PetscLogBegin(); CHKERRQ(err);

    // Initialize Python.
    Py_Initialize();

    pylith::meshio::MeshIOLagrit reader;
    reader.filenameGmv(filenameGmv);
    reader.filenamePset(filenamePset);
    
    pylith::topology::Mesh mesh;
    spatialdata::geocoords::CSCart cs;
    mesh.coordsys(&cs);
    reader.read(&mesh);

    pylith::playpen::TestClosure test;
    test.testRestrictClosure(mesh);

    // Finalize Python
    Py_Finalize();

    // Finalize PETSc
    err = PetscFinalize(); CHKERRQ(err);
  } catch (...) {
    abort();
  } // catch

  return 0;
} // main

// End of file
