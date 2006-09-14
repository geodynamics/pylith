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
#include "elemVector.hh"
#include <iostream> // USES std::cerr

#include "Integration.hh"

// ----------------------------------------------------------------------
ALE::Mesh::section_type::value_type f(ALE::Mesh::section_type::value_type coords[])
{
  return 1.0;
}

int
main(int argc,
     char** argv)
{ // main
  PetscErrorCode err;

  PetscInitialize(&argc, &argv, 0, 0);

  if (argc < 2) {
    std::cerr << "usage: testintegrate MESHIN [options]" << std::endl;
    return -1;
  } // if

  try {
    ALE::Obj<ALE::Mesh> mesh;

    pylith::meshIO::MeshIOAscii iohandler;
    iohandler.filename(argv[1]);
    iohandler.read(mesh, false);

    const ALE::Mesh::topology_type::patch_type patch = 0;
    const Obj<ALE::Mesh::section_type>& field = mesh->getSection("field");
    pylith::feassemble::Integrator      integrator(mesh->getDimension(),
                                                   NUM_QUADRATURE_POINTS, points, weights,
                                                   NUM_BASIS_FUNCTIONS, Basis);

    field->setFiberDimensionByDepth(patch, 0, 1);
    field->allocate();
    integrator.integrateFunction(field, mesh->getSection("coordinates"), f);
    field->view("Weak form of f");
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
