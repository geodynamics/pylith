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

#include "pylith/meshio/MeshIOLagrit.hh"
#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/utils/array.hh" // USES int_array

int
main(int argc,
     char** argv)
{ // main
  try {
    // Initialize PETSc
    PetscErrorCode err = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    CHKERRQ(err);
    
    // Initialize Python
    Py_Initialize();

    const char* filenameGmv = "strikeslip_tet4_1000m.gmv";
    const char* filenamePset = "strikeslip_tet4_1000m.pset";

    pylith::meshio::MeshIOLagrit iohandler;
    iohandler.filenameGmv(filenameGmv);
    iohandler.filenamePset(filenamePset);
    ALE::Obj<Mesh> mesh;
    iohandler.read(&mesh);

    // Finalize Python
    Py_Finalize();

    // Finalize PETSc
    err = PetscFinalize();
    CHKERRQ(err);
  } catch (...) {
    abort();
  } // catch
  
  return 0;
} // main


// End of file
