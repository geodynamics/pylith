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

#include "Distributor.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field<Mesh>
#include "pylith/meshio/DataWriter.hh" // USES DataWriter

#include "journal/info.h" // USES journal::info_t

#include <cstring> // USES strlen()
#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::IntSection IntSection;

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
pylith::topology::Distributor::distribute(topology::Mesh* const newMesh,
					  const topology::Mesh& origMesh,
					  const char* partitioner)
{ // distribute
  assert(0 != newMesh);

  journal::info_t info("distributor");
    
  newMesh->coordsys(origMesh.coordsys());
  
  if (0 == strcasecmp(partitioner, "")) {
    info << journal::at(__HERE__)
	 << "Distributing mesh using dumb partitioner." << journal::endl;
    _distribute<ALE::DistributionNew<SieveMesh> >(newMesh, origMesh);
#if defined(PETSC_HAVE_CHACO)
  } else if (0 == strcasecmp(partitioner, "chaco")) {
    info << journal::at(__HERE__)
	 << "Distributing mesh using 'chaco' partitioner." << journal::endl;
    _distribute<ALE::DistributionNew<SieveMesh, ALE::Partitioner<ALE::Chaco::Partitioner<> > > >(newMesh, origMesh);
#endif
#if defined(PETSC_HAVE_PARMETIS)
  } else if (0 == strcasecmp(partitioner, "parmetis")) {
    info << journal::at(__HERE__)
	 << "Distributing mesh using 'parmetis' partitioner." << journal::endl;
   _distribute<ALE::DistributionNew<SieveMesh, ALE::Partitioner<ALE::ParMetis::Partitioner<> > > >(newMesh, origMesh);
#endif
  } else {
    info << journal::at(__HERE__)
	 << "Unknown partitioner '" << partitioner
	 << "', distribution mesh using dumb partitioner." << journal::endl;
    _distribute<ALE::DistributionNew<SieveMesh> >(newMesh, origMesh);
  } // else

} // distribute

// ----------------------------------------------------------------------
// Write partitioning info for distributed mesh.
void
pylith::topology::Distributor::write(meshio::DataWriter<topology::Mesh, topology::Field<topology::Mesh> >* const writer,
				     const topology::Mesh& mesh)
{ // write
  
  journal::info_t info("distributor");
    
  info << journal::at(__HERE__)
       << "Writing partition." << journal::endl;

  // Setup and allocate field
  const int fiberDim = 1;
  topology::Field<topology::Mesh> partition(mesh);
  partition.scale(1.0);
  partition.label("partition");
  partition.vectorFieldType(topology::FieldBase::SCALAR);
  partition.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
  partition.allocate();
  const ALE::Obj<RealSection>& partitionSection = partition.section();
  assert(!partitionSection.isNull());

  const ALE::Obj<SieveMesh> sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  double rankReal = double(sieveMesh->commRank());
  assert(sieveMesh->height() > 0);
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    partitionSection->updatePoint(*c_iter, &rankReal);
  } // for

  //partition->view("PARTITION");
  const double t = 0.0;
  const int numTimeSteps = 0;
  writer->open(mesh, numTimeSteps);
  writer->openTimeStep(t, mesh);
  writer->writeCellField(t, partition);
  writer->closeTimeStep();
  writer->close();
} // write

// ----------------------------------------------------------------------
template<typename DistributionType>
void
pylith::topology::Distributor::_distribute(topology::Mesh* const newMesh,
					   const topology::Mesh& origMesh)
{ // distribute
  typedef typename SieveMesh::point_type point_type;
  typedef typename DistributionType::partitioner_type partitioner_type;
  typedef typename DistributionType::partition_type   partition_type;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //logger.setDebug(1);
  logger.stagePush("Mesh");
  logger.stagePush("MeshCreation");

  journal::info_t info("distributor");
    
  info << journal::at(__HERE__)
       << "Partitioning and distributing the mesh." << journal::endl;

  ALE::Obj<SieveMesh>& newSieveMesh = newMesh->sieveMesh();
  assert(!newSieveMesh.isNull());
  const ALE::Obj<SieveMesh>& origSieveMesh = origMesh.sieveMesh();
  assert(!origSieveMesh.isNull());

  const ALE::Obj<SieveMesh::sieve_type> newSieve =
    new SieveMesh::sieve_type(origSieveMesh->comm(), origSieveMesh->debug());
  assert(!newSieve.isNull());
  const ALE::Obj<SieveMesh::send_overlap_type> sendMeshOverlap = 
    new SieveMesh::send_overlap_type(origSieveMesh->comm(), 
				     origSieveMesh->debug());
  assert(!sendMeshOverlap.isNull());
  const ALE::Obj<SieveMesh::recv_overlap_type> recvMeshOverlap = 
    new SieveMesh::recv_overlap_type(origSieveMesh->comm(), 
				     origSieveMesh->debug());
  assert(!recvMeshOverlap.isNull());

  newSieveMesh = new SieveMesh(origSieveMesh->comm(), 
			       origSieveMesh->getDimension(), 
			       origSieveMesh->debug());
  assert(!newSieveMesh.isNull());
  newSieveMesh->setSieve(newSieve);
  // IMESH_TODO This might be unnecessary, since the overlap for
  //   submeshes is just the restriction of the overlaps
  //   std::map<point_type,point_type> renumbering;
  SieveMesh::renumbering_type& renumbering = newSieveMesh->getRenumbering();
  // Distribute the mesh
  ALE::Obj<partition_type> partition = 
    DistributionType::distributeMeshV(origSieveMesh, newSieveMesh, 
				      renumbering, 
				      sendMeshOverlap, recvMeshOverlap);
  if (origSieveMesh->debug()) {
    std::cout << "["<<origSieveMesh->commRank()<<"]: Mesh Renumbering:"
	      << std::endl;
    for (SieveMesh::renumbering_type::const_iterator r_iter = renumbering.begin();
	 r_iter != renumbering.end();
	 ++r_iter) {
      std::cout << "["<<origSieveMesh->commRank()<<"]:   global point " 
		<< r_iter->first << " --> " << " local point " 
		<< r_iter->second << std::endl;
    } // for
  } // if

  info << journal::at(__HERE__)
       << "Checking the overlap." << journal::endl;

  // Check overlap
  int localSendOverlapSize = 0, sendOverlapSize;
  int localRecvOverlapSize = 0, recvOverlapSize;
  const int commSize = sendMeshOverlap->commSize();
  for (int p = 0; p < commSize; ++p) {
    localSendOverlapSize += sendMeshOverlap->cone(p)->size();
    localRecvOverlapSize += recvMeshOverlap->support(p)->size();
  } // for
  MPI_Allreduce(&localSendOverlapSize, &sendOverlapSize, 1, MPI_INT, MPI_SUM,
		sendMeshOverlap->comm());
  MPI_Allreduce(&localRecvOverlapSize, &recvOverlapSize, 1, MPI_INT, MPI_SUM,
		recvMeshOverlap->comm());
  if (sendOverlapSize != recvOverlapSize) {
    std::cout <<"["<<sendMeshOverlap->commRank()<<"]: Size mismatch " << 
      sendOverlapSize << " != " << recvOverlapSize << std::endl;
    sendMeshOverlap->view("Send Overlap");
    recvMeshOverlap->view("Recv Overlap");
    throw ALE::Exception("Invalid Overlap");
  } // if

  logger.stagePop();
  logger.stagePush("MeshCoordinates");

  // Distribute the coordinates
  info << journal::at(__HERE__)
       << "Distribution the vertex coordinates." << journal::endl;

  const ALE::Obj<RealSection>& coordinates = 
    origSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const ALE::Obj<RealSection>& parallelCoordinates = 
    newSieveMesh->getRealSection("coordinates");
  assert(!parallelCoordinates.isNull());

  newSieveMesh->setupCoordinates(parallelCoordinates);
  DistributionType::distributeSection(coordinates, partition, renumbering, 
				      sendMeshOverlap, recvMeshOverlap, 
				      parallelCoordinates);

  logger.stagePop();
  logger.stagePush("MeshRealSections");

  // Distribute other sections
  info << journal::at(__HERE__)
       << "Distribution other sections." << journal::endl;

  if (origSieveMesh->getRealSections()->size() > 1) {
    ALE::Obj<std::set<std::string> > names = origSieveMesh->getRealSections();
    assert(!names.isNull());
    int n = 0;

    const std::set<std::string>::const_iterator namesBegin = names->begin();
    const std::set<std::string>::const_iterator namesEnd = names->end();
    for (std::set<std::string>::const_iterator n_iter = namesBegin; 
	 n_iter != namesEnd;
	 ++n_iter) {
      if (*n_iter == "coordinates") continue; // already copied
      if (*n_iter == "replaced_cells") continue; // ignore
      std::cerr << "ERROR: Did not distribute real section '" << *n_iter
		<< "'." << std::endl;
      ++n;
    } // if
    if (n)
      throw std::logic_error("Need to distribute more real sections");
  }

  logger.stagePop();
  logger.stagePush("MeshIntSections");

  if (origSieveMesh->getIntSections()->size() > 0) {
    ALE::Obj<std::set<std::string> > names = origSieveMesh->getIntSections();
    assert(!names.isNull());

    std::set<std::string>::const_iterator namesBegin = names->begin();
    std::set<std::string>::const_iterator namesEnd = names->end();
    for (std::set<std::string>::const_iterator n_iter = namesBegin;
	 n_iter != namesEnd;
	 ++n_iter) {
      const ALE::Obj<IntSection>& origSection = 
	origSieveMesh->getIntSection(*n_iter);
      assert(!origSection.isNull());
      const ALE::Obj<IntSection>& newSection  = 
	newSieveMesh->getIntSection(*n_iter);
      assert(!newSection.isNull());

      // We assume all integer sections are complete sections
      newSection->setChart(newSieveMesh->getSieve()->getChart());
      DistributionType::distributeSection(origSection, partition, renumbering,
					  sendMeshOverlap, recvMeshOverlap, 
					  newSection);
#if 0 // DEBUGGING
      std::string serialName("Serial ");
      std::string parallelName("Parallel ");
      serialName   += *n_iter;
      parallelName += *n_iter;
      origSection->view(serialName.c_str());
      newSection->view(parallelName.c_str());
#endif
    } // for
  } // if
  if (origSieveMesh->getArrowSections()->size() > 1)
    throw std::logic_error("Need to distribute more arrow sections");
  
  logger.stagePop();
  logger.stagePush("MeshLabels");

  // Distribute labels
  info << journal::at(__HERE__)
       << "Distributing labels." << journal::endl;

  const SieveMesh::labels_type& labels = origSieveMesh->getLabels();
  const SieveMesh::labels_type::const_iterator labelsBegin = labels.begin();
  const SieveMesh::labels_type::const_iterator labelsEnd = labels.end();

  for (SieveMesh::labels_type::const_iterator l_iter = labelsBegin;
       l_iter != labelsEnd;
       ++l_iter) {
    if (newSieveMesh->hasLabel(l_iter->first)) continue;
    const ALE::Obj<SieveMesh::label_type>& origLabel = l_iter->second;
    assert(!origLabel.isNull());
    const ALE::Obj<SieveMesh::label_type>& newLabel  = 
      newSieveMesh->createLabel(l_iter->first);
    assert(!newLabel.isNull());

#ifdef IMESH_NEW_LABELS
    newLabel->setChart(newSieve->getChart());
    // Size the local mesh
    partitioner_type::sizeLocalSieveV(origLabel, partition, renumbering, 
				      newLabel);
    // Create the remote meshes
    DistributionType::completeConesV(origLabel, newLabel, renumbering, 
				     sendMeshOverlap, recvMeshOverlap);
    // Create the local mesh
    partitioner_type::createLocalSieveV(origLabel, partition, renumbering, 
					newLabel);
    newLabel->symmetrize();
#else
	DistributionType::distributeLabelV(newSieveMesh->getSieve(), origLabel, partition, renumbering, sendMeshOverlap, recvMeshOverlap, newLabel);
#if 0 // DEBUGGING
    std::string serialName("Serial ");
    std::string parallelName("Parallel ");
    serialName   += l_iter->first;
    parallelName += l_iter->first;
    origLabel->view(serialName.c_str());
    newLabel->view(parallelName.c_str());
#endif
#endif
  } // for

  logger.stagePop();
  logger.stagePush("MeshOverlap");

  // Create the parallel overlap
  info << journal::at(__HERE__)
       << "Creating the parallel overlap." << journal::endl;

  ALE::Obj<SieveMesh::send_overlap_type> sendParallelMeshOverlap = 
    newSieveMesh->getSendOverlap();
  assert(!sendParallelMeshOverlap.isNull());
  ALE::Obj<SieveMesh::recv_overlap_type> recvParallelMeshOverlap = 
    newSieveMesh->getRecvOverlap();
  assert(!recvParallelMeshOverlap.isNull());

  //   Can I figure this out in a nicer way?
  ALE::SetFromMap<std::map<point_type,point_type> > globalPoints(renumbering);
  ALE::OverlapBuilder<>::constructOverlap(globalPoints, renumbering, 
					  sendParallelMeshOverlap, 
					  recvParallelMeshOverlap);
  newSieveMesh->setCalculatedOverlap(true);

  logger.stagePop();
  logger.stagePop(); // Mesh
  //logger.setDebug(0);
} // distribute


// End of file 
