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

#include "MeshIO.hh"
#include "MeshIOAscii.hh"

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

  ALE::Obj<ALE::Mesh> mesh;

  pylith::meshIO::MeshIOAscii iohandler;
  iohandler.filename(argv[1]);
  iohandler.read(mesh);

  iohandler.filename(argv[2]);
  iohandler.write(mesh);

  err = PetscFinalize(); CHKERRQ(err);
  
  return err;
} // main

// version
// $Id$

// End of file 
