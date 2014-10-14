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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Distributor.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field<Mesh>
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/DataWriter.hh" // USES DataWriter

#include "journal/info.h" // USES journal::info_t

#include <cstring> // USES strlen()
#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

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
					  const topology::Mesh& origMesh)
{ // distribute
  PYLITH_METHOD_BEGIN;
  
  assert(newMesh);
  newMesh->coordsys(origMesh.coordsys());
  journal::info_t info("distributor");

  PetscDM newDM = NULL;
  PetscErrorCode err = DMPlexDistribute(origMesh.dmMesh(), 0, NULL, &newDM);PYLITH_CHECK_ERROR(err);
  newMesh->dmMesh(newDM);

  PYLITH_METHOD_END;
} // distribute

// ----------------------------------------------------------------------
// Write partitioning info for distributed mesh.
void
pylith::topology::Distributor::write(meshio::DataWriter* const writer,
				     const topology::Mesh& mesh)
{ // write
  PYLITH_METHOD_BEGIN;

  journal::info_t info("distributor");
    
  const int commRank = mesh.commRank();
  if (0 == commRank) {
    info << journal::at(__HERE__)
	 << "Writing partition." << journal::endl;
  } // if

  // Setup and allocate field
  const int fiberDim = 1;
  topology::Field partition(mesh);
  partition.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
  partition.allocate();
  partition.scale(1.0);
  partition.label("partition");
  partition.vectorFieldType(topology::FieldBase::SCALAR);

  PylithScalar rankReal = PylithScalar(commRank);

  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  topology::VecVisitorMesh partitionVisitor(partition);
  PetscScalar* partitionArray = partitionVisitor.localArray();

  for (PetscInt c = cStart; c < cEnd; ++c) {
    const PetscInt off = partitionVisitor.sectionOffset(c);
    assert(fiberDim == partitionVisitor.sectionDof(c));
    partitionArray[off] = rankReal;
  } // for

  //partition->view("PARTITION");
  const PylithScalar t = 0.0;
  const int numTimeSteps = 0;
  writer->open(mesh, numTimeSteps);
  writer->openTimeStep(t, mesh);
  writer->writeCellField(t, partition);
  writer->closeTimeStep();
  writer->close();

  PYLITH_METHOD_END;
} // write

// End of file 
