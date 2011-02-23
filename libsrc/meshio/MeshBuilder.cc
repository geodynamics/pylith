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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "MeshBuilder.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/array.hh" // USES double_array, int_array

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
				       double_array* coordinates,
				       const int numVertices,
				       const int spaceDim,
				       const int_array& cells,
				       const int numCells,
				       const int numCorners,
				       const int meshDim,
				       const bool interpolate)
{ // buildMesh
  assert(0 != mesh);

  assert(0 != coordinates);
  MPI_Comm comm = mesh->comm();
  int dim = meshDim;
  int rank = 0;

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

  ALE::Obj<SieveMesh::sieve_type> sieve = 
    new SieveMesh::sieve_type(mesh->comm());
  sieveMesh->setSieve(sieve);

  MPI_Comm_rank(comm, &rank);
  // Memory debugging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.setDebug(mesh->debug()/2);

  logger.stagePush("Mesh");
  logger.stagePush("MeshCreation");
  if (0 == rank) {
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
    } else {
      // Same old thing
      ALE::Obj<FlexMesh::sieve_type> s =
	new FlexMesh::sieve_type(sieve->comm(), sieve->debug());

      ALE::SieveBuilder<FlexMesh>::buildTopology(s, meshDim, 
                                                  numCells, 
                                                  const_cast<int*>(&cells[0]), 
                                                  numVertices, 
                                                  interpolate,
                                                  numCorners);
      std::map<FlexMesh::point_type,FlexMesh::point_type> renumbering;
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

#ifdef IMESH_NEW_LABELS
      height->setChart(sieveMesh->getSieve()->getChart());
      depth->setChart(sieveMesh->getSieve()->getChart());
      for(int c = 0; c < numCells+numVertices; ++c) {
        height->setConeSize(c, 1);
        depth->setConeSize(c, 1);
      }
      if (numCells+numVertices)
	height->setSupportSize(0, numCells+numVertices);
      if (numCells+numVertices)
	depth->setSupportSize(0, numCells+numVertices);
      height->allocate();
      depth->allocate();
#endif
      for(int c = 0; c < numCells; ++c) {
        height->setCone(0, c);
        depth->setCone(1, c);
      } // for
      for(int v = numCells; v < numCells+numVertices; ++v) {
        height->setCone(1, v);
        depth->setCone(0, v);
      } // for
#ifdef IMESH_NEW_LABELS
      height->recalculateLabel();
      depth->recalculateLabel();
#endif
      sieveMesh->setHeight(1);
      sieveMesh->setDepth(1);
    } else {
      sieveMesh->stratify();
    } // if/else
    logger.stagePop();
  } else {
    logger.stagePush("MeshStratification");
    sieveMesh->getSieve()->setChart(SieveMesh::sieve_type::chart_type());
    sieveMesh->getSieve()->allocate();
    sieveMesh->stratify();
    logger.stagePop();
  } // if/else

  logger.stagePush("MeshCoordinates");
  ALE::SieveBuilder<SieveMesh>::buildCoordinates(sieveMesh, spaceDim, 
						 &(*coordinates)[0]);
  logger.stagePop(); // MeshCoordinates

  logger.stagePop(); // Mesh

  sieveMesh->getFactory()->clear();
} // buildMesh

// ----------------------------------------------------------------------
// Set vertices and cells for fault mesh.
void
pylith::meshio::MeshBuilder::buildFaultMesh(const ALE::Obj<SieveMesh>& fault,
					    ALE::Obj<FlexMesh>& faultBd,
					    const double_array& coordinates,
					    const int numVertices,
					    const int spaceDim,
					    const int_array& cells,
					    const int numCells,
					    const int numCorners,
					    const int firstCell,
					    const int_array& faceCells,
					    const int meshDim)
{ // buildFaultMesh
  MPI_Comm comm = PETSC_COMM_WORLD;
  int dim  = meshDim;
  int rank = 0;

  assert(!fault.isNull());

  ALE::Obj<SieveMesh::sieve_type> sieve = 
    new SieveMesh::sieve_type(fault->comm());
  fault->setDebug(fault->debug());
  fault->setSieve(sieve);

  MPI_Comm_rank(comm, &rank);

  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.setDebug(fault->debug()/2);
  logger.stagePush("FaultCreation");

  if (0 == rank) {
    assert(coordinates.size() == numVertices*spaceDim);
    assert(cells.size() == numCells*numCorners);
    ALE::Obj<FlexMesh::sieve_type> s = 
      new FlexMesh::sieve_type(sieve->comm(), sieve->debug());
    
    ALE::SieveBuilder<FlexMesh>::buildTopology(s, meshDim, 
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
    
    FlexMesh::renumbering_type& renumbering = fault->getRenumbering();
    ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering, false);
    ALE::ISieveConverter::convertOrientation(*s, *sieve, renumbering,
				fault->getArrowSection("orientation").ptr());
    
    Obj<FlexMesh> tmpMesh = 
      new FlexMesh(fault->comm(), dim, fault->debug());
    faultBd = ALE::Selection<FlexMesh>::boundary(tmpMesh);

    logger.stagePop();
    logger.stagePush("FaultStratification");
    fault->stratify();
    logger.stagePop();
  } else {
    Obj<FlexMesh> tmpMesh = 
      new FlexMesh(fault->comm(), dim, fault->debug());
    faultBd = ALE::Selection<FlexMesh>::boundary(tmpMesh);

    logger.stagePop();
    logger.stagePush("FaultStratification");
    fault->getSieve()->setChart(SieveMesh::sieve_type::chart_type());
    fault->getSieve()->allocate();
    fault->stratify();
    logger.stagePop();
  } // if/else

#if defined(ALE_MEM_LOGGING)
  std::cout
    << std::endl
    << "FaultCreation " << logger.getNumAllocations("FaultCreation")
    << " allocations " << logger.getAllocationTotal("FaultCreation")
    << " bytes" << std::endl
    
    << "FaultCreation " << logger.getNumDeallocations("FaultCreation")
    << " deallocations " << logger.getDeallocationTotal("FaultCreation")
    << " bytes" << std::endl
    
    << "FaultStratification " << logger.getNumAllocations("FaultStratification")
    << " allocations " << logger.getAllocationTotal("FaultStratification")
    << " bytes" << std::endl
    
    << "FaultStratification " << logger.getNumDeallocations("FaultStratification")
    << " deallocations " << logger.getDeallocationTotal("FaultStratification")
    << " bytes" << std::endl << std::endl;
#endif

} // buildFaultMesh


// End of file 
