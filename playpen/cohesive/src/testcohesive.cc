// -*- C++ -*-
#include "pylith/meshio/MeshIOAscii.hh"
#include <iostream> // USES std::cerr

#include "src/dm/mesh/meshpylith.h"

int
main(int argc,
     char** argv)
{ // main
  PetscErrorCode err;

  PetscInitialize(&argc, &argv, 0, 0);

  if (argc < 2) {
    std::cerr << "usage: testcohesive MESHIN [options]" << std::endl;
    return -1;
  } // if

  try {
    ALE::Obj<ALE::Mesh> mesh;

    pylith::meshio::MeshIOAscii iohandler;
    iohandler.filename(argv[1]);
    iohandler.read(&mesh);

    const ALE::Mesh::topology_type::patch_type patch = 0;
    const ALE::Obj<ALE::Mesh::real_section_type>& coords = mesh->getRealSection("coordinates");

    mesh->view("Original Mesh");
    // For the tractest mesh, we will split a face on the midplane
    // Elem 2: 17-22-19-(18) Elem 23: 22-19-(25)-17
    std::set<ALE::Mesh::point_type> faultVertices;

    faultVertices.insert(17+41-1);
    faultVertices.insert(19+41-1);
    faultVertices.insert(22+41-1);
    ALE::PyLith::Builder::createCohesiveElements(mesh, faultVertices);
    mesh->view("New Mesh with Cohesive Elements");
  } catch(ALE::Exception e) {
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    std::cout <<"["<<rank<<"]: " << e << std::endl;
  }
  err = PetscFinalize(); CHKERRQ(err);
  
  return err;
} // main

// version
// $Id$

// End of file
