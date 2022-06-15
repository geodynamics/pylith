// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
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

// This is a copy of DMPlexDistributeOverlap()
static PetscErrorCode DMPlexDistributeOverlap_PyLith(DM dm, DM *dmOverlap)
{
  MPI_Comm     comm;
  PetscMPIInt  size, rank;
  PetscSection rootSection, leafSection;
  IS           rootrank, leafrank;
  DM           dmCoord;
  DMLabel      lblOverlap;
  PetscSF      sfOverlap, sfStratified, sfPoint;
  DMLabel      ovLabels[1], ovExLabel;
  PetscInt     ovValues[1] = {20}, ovExValue = 100;

  PetscCall(PetscObjectGetComm((PetscObject)dm,&comm));
  PetscCallMPI(MPI_Comm_size(comm, &size));
  PetscCallMPI(MPI_Comm_rank(comm, &rank));
  /* Compute point overlap with neighbouring processes on the distributed DM */
  PetscCall(PetscSectionCreate(comm, &rootSection));
  PetscCall(PetscSectionCreate(comm, &leafSection));
  PetscCall(DMPlexDistributeOwnership(dm, rootSection, &rootrank, leafSection, &leafrank));
  PetscCall(DMGetLabel(dm, "fault", &ovLabels[0]));
  PetscCall(DMGetLabel(dm, "material-id", &ovExLabel));
  PetscCall(DMPlexCreateOverlapLabelFromLabels(dm, 1, ovLabels, ovValues, ovExLabel, ovExValue, rootSection, rootrank, leafSection, leafrank, &lblOverlap));
  /* Convert overlap label to stratified migration SF */
  PetscCall(DMPlexPartitionLabelCreateSF(dm, lblOverlap, &sfOverlap));
  PetscCall(DMPlexStratifyMigrationSF(dm, sfOverlap, &sfStratified));
  PetscCall(PetscSFDestroy(&sfOverlap));
  sfOverlap = sfStratified;
  PetscCall(PetscObjectSetName((PetscObject) sfOverlap, "Overlap SF"));
  PetscCall(PetscSFSetFromOptions(sfOverlap));

  PetscCall(PetscSectionDestroy(&rootSection));
  PetscCall(PetscSectionDestroy(&leafSection));
  PetscCall(ISDestroy(&rootrank));
  PetscCall(ISDestroy(&leafrank));

  /* Build the overlapping DM */
  PetscCall(DMPlexCreate(comm, dmOverlap));
  PetscCall(PetscObjectSetName((PetscObject) *dmOverlap, "Parallel Mesh"));
  PetscCall(DMPlexMigrate(dm, sfOverlap, *dmOverlap));
  /* Store the overlap in the new DM */
  PetscCall(DMPlexSetOverlap(*dmOverlap, dm, 1));
  /* Build the new point SF */
  PetscCall(DMPlexCreatePointSF(*dmOverlap, sfOverlap, PETSC_FALSE, &sfPoint));
  PetscCall(DMSetPointSF(*dmOverlap, sfPoint));
  PetscCall(DMGetCoordinateDM(*dmOverlap, &dmCoord));
  if (dmCoord) PetscCall(DMSetPointSF(dmCoord, sfPoint));
  PetscCall(PetscSFDestroy(&sfPoint));
  /* Cleanup overlap partition */
  PetscCall(DMLabelDestroy(&lblOverlap));
  PetscCall(PetscSFDestroy(&sfOverlap));
  return(0);
}

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

    const int commRank = origMesh.getCommRank();
    if (0 == commRank) {
        info << pythia::journal::at(__HERE__)
             << "Partitioning mesh using PETSc '" << partitionerName << "' partitioner." << pythia::journal::endl;
    } // if

    PetscErrorCode err = 0;
    PetscPartitioner partitioner = 0;
    PetscDM dmOrig = origMesh.getDM();assert(dmOrig);
    err = DMPlexGetPartitioner(dmOrig, &partitioner);PYLITH_CHECK_ERROR(err);
    err = PetscPartitionerSetType(partitioner, partitionerName);PYLITH_CHECK_ERROR(err);

    if (0 == commRank) {
        info << pythia::journal::at(__HERE__)
             << "Distributing partitioned mesh." << pythia::journal::endl;
    } // if

    PetscDM dmTmp = NULL, dmNew = NULL;
    const PetscInt overlap = 0;
    err = DMPlexDistribute(origMesh.getDM(), overlap, NULL, &dmTmp);PYLITH_CHECK_ERROR(err);
    err = DMPlexDistributeOverlap_PyLith(dmTmp, &dmNew);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&dmTmp);PYLITH_CHECK_ERROR(err);
    err = DMPlexDistributeSetDefault(dmNew, PETSC_FALSE);PYLITH_CHECK_ERROR(err);
    newMesh->setDM(dmNew);

    PYLITH_METHOD_END;
} // distribute


// ---------------------------------------------------------------------------------------------------------------------
// Write partitioning info for distributed mesh.
void
pylith::topology::Distributor::write(meshio::DataWriter* const writer,
                                     const topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    const int commRank = mesh.getCommRank();
    if (0 == commRank) {
        pythia::journal::info_t info("mesh_distributor");
        info << pythia::journal::at(__HERE__)
             << "Writing partition." << pythia::journal::endl;
    } // if

    // Setup and allocate PETSc vector
    PylithScalar rankReal = PylithReal(commRank);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    PetscVec partitionVec = NULL;
    PylithScalar* partitionArray = NULL;
    PetscErrorCode err;
    err = VecCreate(mesh.getComm(), &partitionVec);PYLITH_CHECK_ERROR(err);
    err = VecSetSizes(partitionVec, cEnd-cStart, PETSC_DECIDE);PYLITH_CHECK_ERROR(err);
    err = VecGetArray(partitionVec, &partitionArray);PYLITH_CHECK_ERROR(err);
    for (PetscInt c = cStart; c < cEnd; ++c) {
        partitionArray[c] = rankReal;
    } // for
    err = VecRestoreArray(partitionVec, &partitionArray);PYLITH_CHECK_ERROR(err);

    const PylithScalar t = 0.0;
    const int numTimeSteps = 0;
    writer->open(mesh, numTimeSteps);
    writer->openTimeStep(t, mesh);
    assert(0); // :TODO: Fix this
    // writer->writeCellField(t, partitionVec);
    writer->closeTimeStep();
    writer->close();

    PYLITH_METHOD_END;
} // write


// End of file
