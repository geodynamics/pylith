// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "MeshBuilder.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "Selection.hh" // USES boundary()

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
// Set vertices and cells in mesh.
void
pylith::meshio::MeshBuilder::buildMesh(topology::Mesh* mesh,
				       scalar_array* coordinates,
				       const int numVertices,
				       const int spaceDim,
				       const int_array& cells,
				       const int numCells,
				       const int numCorners,
				       const int meshDim,
				       const bool interpolate,
				       const bool isParallel)
{ // buildMesh
  PYLITH_METHOD_BEGIN;

  assert(mesh);
  assert(coordinates);
  MPI_Comm comm  = mesh->comm();
  PetscInt dim  = meshDim;
  PetscMPIInt commRank = mesh->commRank();
  PetscErrorCode err;
  const bool useSieve = true; // :TODO: Remove sieve

  { // Check to make sure every vertex is in at least one cell.
    // This is required by Sieve
    std::vector<bool> vertexInCell(numVertices, false);
    const int size = cells.size();
    for (int i=0; i < size; ++i)
      vertexInCell[cells[i]] = true;
    int count = 0;
    for (int i=0; i < numVertices; ++i)
      if (!vertexInCell[i])
        ++count;
    if (count > 0) {
      std::ostringstream msg;
      msg << "Mesh contains " << count
          << " vertices that are not in any cells.";
      throw std::runtime_error(msg.str());
    } // if
  } // check

  /* Sieve */
  if (useSieve) {
    MPI_Bcast(&dim, 1, MPI_INT, 0, comm);
    mesh->createSieveMesh(dim);
    const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
    assert(!sieveMesh.isNull());

    ALE::Obj<SieveMesh::sieve_type> sieve = new SieveMesh::sieve_type(mesh->comm());
    sieveMesh->setSieve(sieve);

    if (0 == commRank || isParallel) {
      assert(coordinates->size() == numVertices*spaceDim);
      assert(cells.size() == numCells*numCorners);
      if (!interpolate) {
	// Create the ISieve
	sieve->setChart(SieveMesh::sieve_type::chart_type(0, numCells+numVertices));
	// Set cone and support sizes
	for (int c = 0; c < numCells; ++c)
	  sieve->setConeSize(c, numCorners);
	sieve->symmetrizeSizes(numCells, numCorners, const_cast<int*>(&cells[0]), numCells);
	// Allocate point storage
	sieve->allocate();
	// Fill up cones
	int *cone  = new int[numCorners];
	int *coneO = new int[numCorners];
	for (int v = 0; v < numCorners; ++v)
	  coneO[v] = 1;
	for (int c = 0; c < numCells; ++c) {
	  for (int v = 0; v < numCorners; ++v)
	    cone[v] = cells[c*numCorners+v]+numCells;
	  sieve->setCone(cone, c);
	  sieve->setConeOrientation(coneO, c);
	} // for
	delete[] cone; cone = 0;
	delete[] coneO; coneO = 0;
	// Symmetrize to fill up supports
	sieve->symmetrize();
      } else {
	// Same old thing
	ALE::Obj<SieveFlexMesh::sieve_type> s = new SieveFlexMesh::sieve_type(sieve->comm(), sieve->debug());
	ALE::Obj<SieveFlexMesh::arrow_section_type> orientation = new SieveFlexMesh::arrow_section_type(sieve->comm(), sieve->debug());

	s->setDebug(2);
	ALE::SieveBuilder<SieveFlexMesh>::buildTopology(s, meshDim, numCells, const_cast<int*>(&cells[0]), numVertices, interpolate,
							numCorners, -1, orientation);
	std::map<SieveFlexMesh::point_type,SieveFlexMesh::point_type> renumbering;
	ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
      } // if/else
      if (!interpolate) {
	// Optimized stratification
	const ALE::Obj<SieveMesh::label_type>& height = sieveMesh->createLabel("height");
	const ALE::Obj<SieveMesh::label_type>& depth = sieveMesh->createLabel("depth");
	
	for(int c = 0; c < numCells; ++c) {
	  height->setCone(0, c);
	  depth->setCone(1, c);
	} // for
	for(int v = numCells; v < numCells+numVertices; ++v) {
	  height->setCone(1, v);
	  depth->setCone(0, v);
	} // for
	sieveMesh->setHeight(1);
	sieveMesh->setDepth(1);
      } else {
	sieveMesh->stratify();
      } // if/else
    } else {
      sieveMesh->getSieve()->setChart(SieveMesh::sieve_type::chart_type());
      sieveMesh->getSieve()->allocate();
      sieveMesh->stratify();
    } // if/else

    ALE::SieveBuilder<SieveMesh>::buildCoordinates(sieveMesh, spaceDim, &(*coordinates)[0]);
    sieveMesh->getFactory()->clear();
  } // if

  /* DMPlex */
  PetscDM dmMesh;
  PetscBool pInterpolate = interpolate ? PETSC_TRUE : PETSC_FALSE;

  err = DMPlexCreateFromCellList(comm, dim, numCells, numVertices, numCorners, pInterpolate, &cells[0], spaceDim, &(*coordinates)[0], &dmMesh);CHECK_PETSC_ERROR(err);
  mesh->setDMMesh(dmMesh);

  PYLITH_METHOD_END;
} // buildMesh

// End of file 
