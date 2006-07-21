// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include "MeshIOAscii.hh"

#include "petsc.h" // USES PetscInitialize(), PetscFinalize()
#include "PetscMesh.hh" // USES PetscMesh

#include <iostream> // USES std::cerr

// ----------------------------------------------------------------------
int
main(int argc,
     char** argv)
{ // main
  PetscErrorCode err;

  PetscInitialize(&argc, &argv, 0, 0);
  
  if (3 != argc) {
    std::cerr << "usage: testascii MESHIN MESHOUT" << std::endl;
    return -1;
  } // if

  ALE::Obj<ALE::PetscMesh> mesh;

  MeshIOAscii iohandler;
  iohandler.filename(argv[1]);
  iohandler.read(&mesh);

  iohandler.filename(argv[2]);
  iohandler.write(mesh);

  err = PetscFinalize(); CHKERRQ(err);
  
  return err;
} // main

// version
// $Id$

// End of file 
