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
  assert(mesh);

  assert(coordinates);
  MPI_Comm comm = mesh->comm();
  int dim = meshDim;
  const int commRank = mesh->commRank();
  PetscErrorCode err;

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

  MPI_Bcast(&dim, 1, MPI_INT, 0, comm);
  mesh->createSieveMesh(dim);
  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
  assert(!sieveMesh.isNull());
  /* DMComplex */
  mesh->createDMMesh(dim);
  DM complexMesh = mesh->dmMesh();
  assert(complexMesh);

  ALE::Obj<SieveMesh::sieve_type> sieve = 
    new SieveMesh::sieve_type(mesh->comm());
  sieveMesh->setSieve(sieve);

  // Memory debugging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //logger.setDebug(1);

  logger.stagePush("MeshCreation");
  if (0 == commRank || isParallel) {
    assert(coordinates->size() == numVertices*spaceDim);
    assert(cells.size() == numCells*numCorners);
    if (!interpolate) {
      // Create the ISieve
      sieve->setChart(SieveMesh::sieve_type::chart_type(0, 
							numCells+numVertices));
      // Set cone and support sizes
      for (int c = 0; c < numCells; ++c)
	sieve->setConeSize(c, numCorners);
      sieve->symmetrizeSizes(numCells, numCorners, 
			     const_cast<int*>(&cells[0]), numCells);
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
      /* DMComplex */
      err = DMComplexSetChart(complexMesh, 0, numCells+numVertices);CHECK_PETSC_ERROR(err);
      for(PetscInt c = 0; c < numCells; ++c) {
        err = DMComplexSetConeSize(complexMesh, c, numCorners);CHECK_PETSC_ERROR(err);
      }
      err = DMSetUp(complexMesh);CHECK_PETSC_ERROR(err);
      PetscInt *cone2 = new PetscInt[numCorners];
      for(PetscInt c = 0; c < numCells; ++c) {
        for(PetscInt v = 0; v < numCorners; ++v) {
          cone2[v] = cells[c*numCorners+v]+numCells;
        }
        err = DMComplexSetCone(complexMesh, c, cone2);CHECK_PETSC_ERROR(err);
      } // for
      delete [] cone2; cone2 = 0;
      err = DMComplexSymmetrize(complexMesh);CHECK_PETSC_ERROR(err);
      err = DMComplexStratify(complexMesh);CHECK_PETSC_ERROR(err);
    } else {
      // Same old thing
      ALE::Obj<SieveFlexMesh::sieve_type> s =
	new SieveFlexMesh::sieve_type(sieve->comm(), sieve->debug());
      ALE::Obj<SieveFlexMesh::arrow_section_type> orientation = new SieveFlexMesh::arrow_section_type(sieve->comm(), sieve->debug());

      s->setDebug(2);
      ALE::SieveBuilder<SieveFlexMesh>::buildTopology(s, meshDim, 
                                                      numCells, 
                                                      const_cast<int*>(&cells[0]), 
                                                      numVertices, 
                                                      interpolate,
                                                      numCorners,
                                                      -1,
                                                      orientation);
      std::map<SieveFlexMesh::point_type,SieveFlexMesh::point_type> renumbering;
      ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
    } // if/else
    logger.stagePop();
    logger.stagePush("MeshStratification");
    if (!interpolate) {
      // Optimized stratification
      const ALE::Obj<SieveMesh::label_type>& height = 
	sieveMesh->createLabel("height");
      const ALE::Obj<SieveMesh::label_type>& depth =
	sieveMesh->createLabel("depth");

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
    logger.stagePop();
  } else {
    logger.stagePop();
    logger.stagePush("MeshStratification");
    sieveMesh->getSieve()->setChart(SieveMesh::sieve_type::chart_type());
    sieveMesh->getSieve()->allocate();
    sieveMesh->stratify();
    logger.stagePop();
  } // if/else

  logger.stagePush("MeshCoordinates");
  ALE::SieveBuilder<SieveMesh>::buildCoordinates(sieveMesh, spaceDim, 
						 &(*coordinates)[0]);
  /* DMComplex */
  PetscSection coordSection;
  Vec          coordVec;
  PetscScalar *coords;
  PetscInt     coordSize;

  err = DMComplexGetCoordinateSection(complexMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMComplexGetCoordinateVec(complexMesh, &coordVec);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(coordSection, numCells, numCells+numVertices);CHECK_PETSC_ERROR(err);
  for(PetscInt v = numCells; v < numCells+numVertices; ++v) {
    err = PetscSectionSetDof(coordSection, v, spaceDim);CHECK_PETSC_ERROR(err);
  }
  err = PetscSectionSetUp(coordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetStorageSize(coordSection, &coordSize);CHECK_PETSC_ERROR(err);
  err = VecSetSizes(coordVec, coordSize, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(coordVec);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  for(PetscInt v = 0; v < numVertices; ++v) {
    PetscInt off;

    err = PetscSectionGetOffset(coordSection, v+numCells, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < spaceDim; ++d) {
      coords[off+d] = (*coordinates)[v*spaceDim+d];
    }
  }
  err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  logger.stagePop(); // Coordinates

  sieveMesh->getFactory()->clear();
} // buildMesh

// ----------------------------------------------------------------------
// Set vertices and cells for fault mesh.
void
pylith::meshio::MeshBuilder::buildFaultMesh(const ALE::Obj<SieveMesh>& fault,
					    ALE::Obj<SieveFlexMesh>& faultBd,
					    const scalar_array& coordinates,
					    const int numVertices,
					    const int spaceDim,
					    const int_array& cells,
					    const int numCells,
					    const int numCorners,
					    const int firstCell,
					    const int_array& faceCells,
					    const int meshDim)
{ // buildFaultMesh
  int dim  = meshDim;

  assert(!fault.isNull());

  ALE::Obj<SieveMesh::sieve_type> sieve = 
    new SieveMesh::sieve_type(fault->comm());
  fault->setDebug(fault->debug());
  fault->setSieve(sieve);
  
  const int commRank = fault->commRank();

  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //logger.setDebug(fault->debug()/2);
  logger.stagePush("Creation");

  if (0 == commRank) {
    assert(coordinates.size() == numVertices*spaceDim);
    assert(cells.size() == numCells*numCorners);
    ALE::Obj<SieveFlexMesh::sieve_type> s = 
      new SieveFlexMesh::sieve_type(sieve->comm(), sieve->debug());
    
    ALE::SieveBuilder<SieveFlexMesh>::buildTopology(s, meshDim, 
						numCells, 
						const_cast<int*>(&cells[0]), 
						numVertices, 
						true,
						numCorners,
						0,
						fault->getArrowSection("orientation"),
						firstCell);

    // Add in cells
    for(int c = 0; c < numCells; ++c) {
      s->addArrow(c+firstCell, faceCells[c*2+0]);
      s->addArrow(c+firstCell, faceCells[c*2+1]);
    } // for
    
    SieveFlexMesh::renumbering_type& renumbering = fault->getRenumbering();
    ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering, false);
    ALE::ISieveConverter::convertOrientation(*s, *sieve, renumbering,
				fault->getArrowSection("orientation").ptr());
    
    Obj<SieveFlexMesh> tmpMesh = 
      new SieveFlexMesh(fault->comm(), dim, fault->debug());
    faultBd = ALE::Selection<SieveFlexMesh>::boundary(tmpMesh);

    logger.stagePop();
    logger.stagePush("Stratification");
    fault->stratify();
    logger.stagePop();
  } else {
    Obj<SieveFlexMesh> tmpMesh = 
      new SieveFlexMesh(fault->comm(), dim, fault->debug());
    faultBd = ALE::Selection<SieveFlexMesh>::boundary(tmpMesh);

    logger.stagePop();
    logger.stagePush("Stratification");
    fault->getSieve()->setChart(SieveMesh::sieve_type::chart_type());
    fault->getSieve()->allocate();
    fault->stratify();
    logger.stagePop();
  } // if/else

#if defined(ALE_MEM_LOGGING)
  std::cout
    << std::endl
    << "FaultCreation " << logger.getNumAllocations("Creation")
    << " allocations " << logger.getAllocationTotal("Creation")
    << " bytes" << std::endl
    
    << "FaultCreation " << logger.getNumDeallocations("Creation")
    << " deallocations " << logger.getDeallocationTotal("Creation")
    << " bytes" << std::endl
    
    << "FaultStratification " << logger.getNumAllocations("Stratification")
    << " allocations " << logger.getAllocationTotal("Stratification")
    << " bytes" << std::endl
    
    << "FaultStratification " << logger.getNumDeallocations("Stratification")
    << " deallocations " << logger.getDeallocationTotal("Stratification")
    << " bytes" << std::endl << std::endl;
#endif

} // buildFaultMesh


// End of file 
