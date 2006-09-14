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

  if (argc < 3) {
    std::cerr << "usage: testascii MESHIN MESHOUT [options]" << std::endl;
    return -1;
  } // if

  try {
    ALE::Obj<ALE::Mesh> mesh;

    pylith::meshio::MeshIOAscii iohandler;
    iohandler.filename(argv[1]);
    iohandler.read(mesh);

    iohandler.filename(argv[2]);
    iohandler.write(mesh);
  } catch(const ALE::Exception& err) {
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    std::cout << "[" << rank << "]: " << err << std::endl;
  }
  err = PetscFinalize(); CHKERRQ(err);
  
  return err;
} // main

// version
// $Id$

// End of file 
