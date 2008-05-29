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

#include <portinfo>

#include "Distributor.hh" // implementation of class methods

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/utils/vectorfields.hh" // USES SCALAR_FIELD
#include "pylith/meshio/OutputManager.hh" // USES OutputManager

#include <string.h> // USES strlen()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::topology::Distributor::Distributor(void)
{ // constructor
} // constructor
 
// ----------------------------------------------------------------------
// Destructor
pylith::topology::Distributor::~Distributor(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Distribute mesh among processors.
void
pylith::topology::Distributor::distribute(ALE::Obj<Mesh>* const newMesh,
					  const ALE::Obj<Mesh>& origMesh,
					  const char* partitioner)
{ // distribute
  typedef Mesh::point_type                  point_type;
  typedef ALE::DistributionNew<Mesh>        distribution_type;
  typedef distribution_type::partition_type partition_type;

  const Obj<Mesh::sieve_type>        newSieve        = new Mesh::sieve_type(origMesh->comm(), origMesh->debug());
  const Obj<Mesh::send_overlap_type> sendMeshOverlap = new Mesh::send_overlap_type(origMesh->comm(), origMesh->debug());
  const Obj<Mesh::recv_overlap_type> recvMeshOverlap = new Mesh::recv_overlap_type(origMesh->comm(), origMesh->debug());

  *newMesh = new Mesh(origMesh->comm(), origMesh->getDimension(), origMesh->debug());
  (*newMesh)->setSieve(newSieve);
  // IMESH_TODO
  //   This might be unnecessary, since the overlap for submeshes is just the restriction of the overlaps
  // std::map<point_type,point_type>    renumbering;
  Mesh::renumbering_type&            renumbering     = (*newMesh)->getRenumbering();
  // Distribute the mesh
  if (strlen(partitioner) != 0) {
    std::cout << "ERROR: Using default partitioner instead of " << partitioner << std::endl;
  }
  Obj<partition_type> partition = distribution_type::distributeMeshV(origMesh, (*newMesh), renumbering, sendMeshOverlap, recvMeshOverlap);
  // Distribute the coordinates
  const Obj<real_section_type>& coordinates         = origMesh->getRealSection("coordinates");
  const Obj<real_section_type>& parallelCoordinates = (*newMesh)->getRealSection("coordinates");

  (*newMesh)->setupCoordinates(parallelCoordinates);
  distribution_type::distributeSection(coordinates, partition, renumbering, sendMeshOverlap, recvMeshOverlap, parallelCoordinates);
  // Distribute other sections
  if (origMesh->getRealSections()->size() > 1) {
    throw ALE::Exception("Need to distribute more real sections");
  }
  if (origMesh->getIntSections()->size() > 0) {
    Obj<std::set<std::string> > names = origMesh->getIntSections();

    for(std::set<std::string>::const_iterator n_iter = names->begin(); n_iter != names->end(); ++n_iter) {
      const Obj<Mesh::int_section_type>& origSection = origMesh->getIntSection(*n_iter);
      const Obj<Mesh::int_section_type>& newSection  = (*newMesh)->getIntSection(*n_iter);

      // We assume all integer sections are complete sections
      newSection->setChart((*newMesh)->getSieve()->getChart());
      distribution_type::distributeSection(origSection, partition, renumbering, sendMeshOverlap, recvMeshOverlap, newSection);
#if 0
      std::string serialName("Serial ");
      std::string parallelName("Parallel ");
      serialName   += *n_iter;
      parallelName += *n_iter;
      origSection->view(serialName.c_str());
      newSection->view(parallelName.c_str());
#endif
    }
  }
  if (origMesh->getArrowSections()->size() > 1) {
    throw ALE::Exception("Need to distribute more arrow sections");
  }
  // Distribute labels
  const Mesh::labels_type& labels = origMesh->getLabels();

  for(Mesh::labels_type::const_iterator l_iter = labels.begin(); l_iter != labels.end(); ++l_iter) {
    if ((*newMesh)->hasLabel(l_iter->first)) continue;
    const Obj<Mesh::label_type>& origLabel = l_iter->second;
    const Obj<Mesh::label_type>& newLabel  = (*newMesh)->createLabel(l_iter->first);
    // Get remote labels
    ALE::New::Completion<Mesh,Mesh::point_type>::scatterCones(origLabel, newLabel, sendMeshOverlap, recvMeshOverlap, renumbering);
    // Create local label
    newLabel->add(origLabel, (*newMesh)->getSieve(), renumbering);
#if 0
    std::string serialName("Serial ");
    std::string parallelName("Parallel ");
    serialName   += l_iter->first;
    parallelName += l_iter->first;
    origLabel->view(serialName.c_str());
    newLabel->view(parallelName.c_str());
#endif
  }
  // Create the parallel overlap
  Obj<Mesh::send_overlap_type> sendParallelMeshOverlap = (*newMesh)->getSendOverlap();
  Obj<Mesh::recv_overlap_type> recvParallelMeshOverlap = (*newMesh)->getRecvOverlap();
  //   Can I figure this out in a nicer way?
  ALE::SetFromMap<std::map<point_type,point_type> > globalPoints(renumbering);

  ALE::OverlapBuilder<>::constructOverlap(globalPoints, renumbering, sendParallelMeshOverlap, recvParallelMeshOverlap);
  (*newMesh)->setCalculatedOverlap(true);
} // distribute

// ----------------------------------------------------------------------
// Write partitioning info for distributed mesh.
void
pylith::topology::Distributor::write(meshio::OutputManager* const output,
				     const ALE::Obj<Mesh>& mesh,
				     const spatialdata::geocoords::CoordSys* cs)
{ // write
  
  // Setup and allocate field
  const int fiberDim = 1;
  ALE::Obj<real_section_type> partition = 
    new real_section_type(mesh->comm(), mesh->debug());
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);
  assert(!cells.isNull());
  partition->setFiberDimension(cells, fiberDim);
  mesh->allocate(partition);

  const int rank  = mesh->commRank();
  double rankReal = double(rank);
  const Mesh::label_sequence::iterator  cellsEnd = cells->end();
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    partition->updatePoint(*c_iter, &rankReal);
  } // for


  partition->view("PARTITION");
  const double t = 0.0;
  const int numTimeSteps = 0;
  output->open(mesh, cs, numTimeSteps);
  output->openTimeStep(t, mesh, cs);
  output->appendCellField(t, "partition", partition, SCALAR_FIELD, mesh);
  output->closeTimeStep();
  output->close();
} // write


// End of file 
