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
  if (strlen(partitioner) == 0)
    *newMesh = ALE::Distribution<ALE::Mesh>::distributeMesh(origMesh);
  else
    *newMesh = 
      ALE::Distribution<ALE::Mesh>::distributeMesh(origMesh, 0, partitioner);
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
  output->open(mesh, cs);
  output->openTimeStep(t, mesh, cs);
  output->appendCellField(t, "partition", partition, mesh);
  output->closeTimeStep();
  output->close();
} // write


// End of file 
