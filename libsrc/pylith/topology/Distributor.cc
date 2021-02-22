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
// Copyright (c) 2010-2017 University of California, Davis
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
#include "pylith/utils/journals.hh" // pythia::journal

#include <cstring> // USES strlen()
#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::topology::Distributor::Distributor(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::topology::Distributor::~Distributor(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Distribute mesh among processors.
void
pylith::topology::Distributor::distribute(topology::Mesh* const newMesh,
                                          const topology::Mesh& origMesh,
                                          const char* partitionerName) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::info_t info("mesh_distributor");

    assert(newMesh);
    newMesh->setCoordSys(origMesh.getCoordSys());

    const int commRank = origMesh.commRank();
    if (0 == commRank) {
        info << pythia::journal::at(__HERE__)
             << "Partitioning mesh using PETSc '" << partitionerName << "' partitioner." << pythia::journal::endl;
    } // if

    PetscErrorCode err = 0;
    PetscPartitioner partitioner = 0;
    PetscDM dmOrig = origMesh.dmMesh();assert(dmOrig);
    err = DMPlexGetPartitioner(dmOrig, &partitioner);PYLITH_CHECK_ERROR(err);
    err = PetscPartitionerSetType(partitioner, partitionerName);PYLITH_CHECK_ERROR(err);

    if (0 == commRank) {
        info << pythia::journal::at(__HERE__)
             << "Distributing partitioned mesh." << pythia::journal::endl;
    } // if

    PetscDM dmNew = NULL;
    err = DMPlexDistribute(origMesh.dmMesh(), 0, NULL, &dmNew);PYLITH_CHECK_ERROR(err);
    newMesh->dmMesh(dmNew);

    PYLITH_METHOD_END;
} // distribute


// ---------------------------------------------------------------------------------------------------------------------
// Write partitioning info for distributed mesh.
void
pylith::topology::Distributor::write(meshio::DataWriter* const writer,
                                     const topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    const int commRank = mesh.commRank();
    if (0 == commRank) {
        pythia::journal::info_t info("mesh_distributor");
        info << pythia::journal::at(__HERE__)
             << "Writing partition." << pythia::journal::endl;
    } // if

    // Setup and allocate field
    pylith::topology::Field partition(mesh);
    const char* components[1] = {"partition"};
    const int numComponents = 1;
    const int basisOrder = 0;
    const int quadOrder = 0;
    const int dim = -1;
    const double scale = 1.0;
    pylith::topology::FieldBase::CellBasis cellBasis = mesh.isSimplex() ?
                                                       pylith::topology::FieldBase::SIMPLEX_BASIS : pylith::topology::FieldBase::TENSOR_BASIS;
    const bool isBasisContinuous = true;
    partition.subfieldAdd("partition", "partition", pylith::topology::Field::SCALAR, components, numComponents, scale,
                          basisOrder, quadOrder, dim, cellBasis, isBasisContinuous, pylith::topology::Field::POLYNOMIAL_SPACE);
    partition.subfieldsSetup();
    partition.createDiscretization();
    partition.allocate();
    partition.setLabel("partition");

    PylithScalar rankReal = PylithReal(commRank);

    PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    topology::VecVisitorMesh partitionVisitor(partition);
    PetscScalar* partitionArray = partitionVisitor.localArray();

    for (PetscInt c = cStart; c < cEnd; ++c) {
        const PetscInt off = partitionVisitor.sectionOffset(c);
        assert(numComponents == partitionVisitor.sectionDof(c));
        partitionArray[off] = rankReal;
    } // for

    // partition->view("PARTITION");
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
