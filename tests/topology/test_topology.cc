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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

// Application for using valgrind to test for memory leaks in mesh
// related operations, including reading, reordering, distribution,
// and refinement.

#include <portinfo>

#include "pylith/topology/Mesh.hh"
#include "pylith/meshio/MeshIOCubit.hh"
#include "pylith/topology/ReverseCuthillMcKee.hh"
#include "pylith/topology/Distributor.hh"
#include "pylith/topology/RefineUniform.hh"

#include <petsc.h>
#include <Python.h>

#include <stdlib.h> // USES abort()

int
main(int argc,
     char* argv[])
{ // main
  std::string meshFilename = "data/tet4.exo";
    // Initialize PETSc
    PetscErrorCode err = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    CHKERRQ(err);

    // Initialize Python
    Py_Initialize();

    if (argc == 2)
      meshFilename = argv[1];

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    pylith::meshio::MeshIOCubit reader;
    reader.filename(meshFilename.c_str());
    reader.useNodesetNames(false);

    // Read mesh
    std::cout << "["<<rank<<"] Reading mesh" << std::endl;
    pylith::topology::Mesh mesh;
    reader.read(&mesh);

    // Reorder mesh
    std::cout << "["<<rank<<"] Reordering mesh" << std::endl;
    pylith::topology::ReverseCuthillMcKee order;
    order.reorder(&mesh);

    int nprocs = 0;
    MPI_Comm_size(mesh.comm(), &nprocs);
    if (nprocs > 1) {
      // Distribute mesh
      std::cout << "["<<rank<<"] Distributing mesh" << std::endl;
      pylith::topology::Mesh dmesh(mesh.dimension(), mesh.comm());
      pylith::topology::Distributor::distribute(&dmesh, dmesh, "chaco");
      mesh.deallocate();

      // Refine mesh
      std::cout << "["<<rank<<"] Refining mesh" << std::endl;
      pylith::topology::RefineUniform refiner;
      pylith::topology::Mesh rmesh(dmesh.dimension(), dmesh.comm());
      refiner.refine(&rmesh, dmesh);
      rmesh.deallocate();
    } else {
      // Refine mesh
      std::cout << "["<<rank<<"] Refining mesh" << std::endl;
      pylith::topology::RefineUniform refiner;
      pylith::topology::Mesh rmesh(mesh.dimension(), mesh.comm());
      refiner.refine(&rmesh, mesh);
      rmesh.deallocate();
    } // if/else

    // Finalize Python
    Py_Finalize();

    // Finalize PETSc
    err = PetscFinalize();
    CHKERRQ(err);

    return 0;
} // main


// End of file
