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

#include "pylith/meshio/MeshIOAscii.hh"
#include "elemVector.hh"
#include <iostream> // USES std::cerr

#include "Integration.hh"

// ----------------------------------------------------------------------
ALE::Mesh::section_type::value_type zeroF(const ALE::Mesh::section_type::value_type coords[])
{
  return 0.0;
}

ALE::Mesh::section_type::value_type constantF(const ALE::Mesh::section_type::value_type coords[])
{
  return 1.0;
}

ALE::Mesh::section_type::value_type linearF(const ALE::Mesh::section_type::value_type coords[])
{
  return coords[0];
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

    pylith::meshio::MeshIOAscii iohandler;
    iohandler.filename(argv[1]);
    iohandler.read(mesh, false);

    const ALE::Mesh::topology_type::patch_type patch = 0;
    const Obj<ALE::Mesh::section_type>& coords = mesh->getSection("coordinates");
    const Obj<ALE::Mesh::section_type>& X = mesh->getSection("X");
    const Obj<ALE::Mesh::section_type>& F = mesh->getSection("F");
    pylith::feassemble::Integrator      integrator(mesh->getDimension(),
                                                   NUM_QUADRATURE_POINTS, points, weights,
                                                   NUM_BASIS_FUNCTIONS, Basis, BasisDerivatives);

    coords->view("Coordinates");
    X->setFiberDimensionByDepth(patch, 0, 1);
    X->allocate();
    F->setFiberDimensionByDepth(patch, 0, 1);
    F->allocate();
    integrator.integrateFunction(X, coords, constantF);
    X->view("Weak form of a constant");
    integrator.fillSection(X, coords, constantF);
    integrator.integrateLaplacianAction(X, F, coords);
    F->view("Weak form of \Delta constant");
    integrator.fillSection(X, coords, zeroF);
    integrator.integrateFunction(X, coords, linearF);
    X->view("Weak form of a linear");
    integrator.fillSection(X, coords, linearF);
    integrator.integrateLaplacianAction(X, F, coords);
    F->view("Weak form of \Delta linear");
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
