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
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/faults/FaultCohesive.hh" // USES FaultCohesive
#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/utils/journals.hh" // pythia::journal

#include <cstring> // USES strlen()
#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace topology {
        class _Distributor {
public:

            /** Distribute custom overlap based on PETSc labels.
             *
             * The overlap excludes cohesive cells but includes cells adjacent to faults.
             * This is a custom version of DMPlexDistributeOverlap()
             *
             * @param[out] dmOverlap PETSc DM for the overlap.
             * @param[in] dmMesh PETSc DM for the current mesh.
             * @param[in] faults Array of fault interfaces.
             * @param[in] numFaults Number of fault interfaces.
             *
             * @returns PETSc error code (0==success).
             */
            static
            PetscErrorCode distributeOverlap(PetscDM* dmOverlap,
                                             PetscDM dmMesh,
                                             pylith::faults::FaultCohesive* faults[],
                                             const int numFaults);

        }; // _Distributor
    } // topology
} // pylith

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::topology::Distributor::Distributor(void) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::topology::Distributor::~Distributor(void) {}


// ------------------------------------------------------------------------------------------------
// Distribute mesh among processors.
void
pylith::topology::Distributor::distribute(pylith::topology::Mesh* const newMesh,
                                          const pylith::topology::Mesh& origMesh,
                                          pylith::faults::FaultCohesive* faults[],
                                          const int numFaults,
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
    err = _Distributor::distributeOverlap(&dmNew, dmTmp, faults, numFaults);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&dmTmp);PYLITH_CHECK_ERROR(err);
    err = DMPlexDistributeSetDefault(dmNew, PETSC_FALSE);PYLITH_CHECK_ERROR(err);
    newMesh->setDM(dmNew);

    PYLITH_METHOD_END;
} // distribute


// ------------------------------------------------------------------------------------------------
// Write partitioning info for distributed mesh.
void
pylith::topology::Distributor::write(meshio::DataWriter* const writer,
                                     const topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    if (pylith::utils::MPI::isRoot()) {
        pythia::journal::info_t info("mesh_distributor");
        info << pythia::journal::at(__HERE__)
             << "Writing partition." << pythia::journal::endl;
    } // if

    // Setup and allocate PETSc vector
    const int commRank = mesh.getCommRank();
    PylithScalar rankReal = PylithReal(commRank);

    pylith::topology::Field partitionField(mesh);
    partitionField.setLabel("partition");

    pylith::topology::Field::Description description;
    description.label = "partition";
    description.alias = "partition";
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = "rank";
    description.scale = 1.0;
    description.validator = NULL;

    pylith::topology::Field::Discretization discretization(0, 1);

    partitionField.subfieldAdd(description, discretization);
    partitionField.subfieldsSetup();
    partitionField.createDiscretization();
    partitionField.allocate();
    partitionField.zeroLocal();
    partitionField.createOutputVector();

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    topology::Stratum cellsStratum(dmMesh, pylith::topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    VecVisitorMesh partitionVisitor(partitionField);
    PetscScalar* partitionArray = partitionVisitor.localArray();
    for (PetscInt point = cStart; point < cEnd; ++point) {
        const PetscInt off = partitionVisitor.sectionOffset(point);
        if (partitionVisitor.sectionDof(point) > 0) {
            partitionArray[off] = rankReal;
        } // if
    } // for
    partitionVisitor.clear();
    partitionField.scatterLocalToOutput();

    const int basisOrder = 0;
    pylith::meshio::OutputSubfield* outputField =
        pylith::meshio::OutputSubfield::create(partitionField, mesh, "partition", basisOrder);
    outputField->project(partitionField.getOutputVector());

    const PylithScalar t = 0.0;
    const bool isInfo = true;
    writer->open(mesh, isInfo);
    writer->openTimeStep(t, mesh);
    writer->writeCellField(t, *outputField);
    writer->closeTimeStep();
    writer->close();

    delete outputField;outputField = NULL;

    PYLITH_METHOD_END;
} // write


// ------------------------------------------------------------------------------------------------
// This is a copy of DMPlexDistributeOverlap()
PetscErrorCode
pylith::topology::_Distributor::distributeOverlap(PetscDM* dmOverlap,
                                                  PetscDM dmMesh,
                                                  pylith::faults::FaultCohesive* faults[],
                                                  const int numFaults) {
    PYLITH_METHOD_BEGIN;
    assert(dmOverlap);

    MPI_Comm comm;
    PetscMPIInt size, rank;
    PetscSection rootSection, leafSection;
    PetscIS rootrank, leafrank;
    PetscDM dmCoord;
    PetscDMLabel lblOverlap;
    PetscSF sfOverlap, sfStratified, sfPoint;
    PetscErrorCode err;

    if (0 == numFaults) {
        err = PetscObjectReference((PetscObject)dmMesh);PYLITH_CHECK_ERROR(err);
        *dmOverlap = dmMesh;
        PYLITH_METHOD_RETURN(0);
    } // if

    PetscDMLabel* ovIncludeLabels = (numFaults > 0) ? new PetscDMLabel[numFaults] : NULL;
    PetscInt* ovIncludeLabelValues = (numFaults > 0) ? new PetscInt[numFaults] : NULL;
    PetscDMLabel* ovExcludeLabels = (numFaults > 0) ? new PetscDMLabel[numFaults] : NULL;
    PetscInt* ovExcludeLabelValues = (numFaults > 0) ? new PetscInt[numFaults] : NULL;

    for (int i = 0; i < numFaults; ++i) {
        const char* surfaceLabelName = faults[i]->getSurfaceLabelName();
        err = DMGetLabel(dmMesh, surfaceLabelName, &ovIncludeLabels[i]);PYLITH_CHECK_ERROR(err);
        ovIncludeLabelValues[i] = faults[i]->getSurfaceLabelValue();

        const char* cohesiveLabelName = faults[i]->getCohesiveLabelName();
        err = DMGetLabel(dmMesh, cohesiveLabelName, &ovExcludeLabels[i]);PYLITH_CHECK_ERROR(err);
        ovExcludeLabelValues[i] = faults[i]->getCohesiveLabelValue();
    } // for

    PetscCall(PetscObjectGetComm((PetscObject)dmMesh,&comm));
    PetscCallMPI(MPI_Comm_size(comm, &size));
    PetscCallMPI(MPI_Comm_rank(comm, &rank));
    /* Compute point overlap with neighbouring processes on the distributed DM */
    PetscCall(PetscSectionCreate(comm, &rootSection));
    PetscCall(PetscSectionCreate(comm, &leafSection));
    PetscCall(DMPlexDistributeOwnership(dmMesh, rootSection, &rootrank, leafSection, &leafrank));
    PetscCall(DMPlexCreateOverlapLabelFromLabels(dmMesh, numFaults, ovIncludeLabels, ovIncludeLabelValues,
                                                 numFaults, ovExcludeLabels, ovExcludeLabelValues, rootSection, rootrank, leafSection, leafrank, &lblOverlap));

    delete[] ovIncludeLabels;ovIncludeLabels = NULL;
    delete[] ovIncludeLabelValues;ovIncludeLabelValues = NULL;
    delete[] ovExcludeLabels;ovExcludeLabels = NULL;
    delete[] ovExcludeLabelValues;ovExcludeLabelValues = NULL;

    /* Convert overlap label to stratified migration SF */
    PetscCall(DMPlexPartitionLabelCreateSF(dmMesh, lblOverlap, &sfOverlap));
    PetscCall(DMPlexStratifyMigrationSF(dmMesh, sfOverlap, &sfStratified));
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
    PetscCall(DMPlexMigrate(dmMesh, sfOverlap, *dmOverlap));
    /* Store the overlap in the new DM */
    PetscCall(DMPlexSetOverlap(*dmOverlap, dmMesh, 1));
    /* Build the new point SF */
    PetscCall(DMPlexCreatePointSF(*dmOverlap, sfOverlap, PETSC_FALSE, &sfPoint));
    PetscCall(DMSetPointSF(*dmOverlap, sfPoint));
    PetscCall(DMGetCoordinateDM(*dmOverlap, &dmCoord));
    if (dmCoord) { PetscCall(DMSetPointSF(dmCoord, sfPoint));}
    PetscCall(PetscSFDestroy(&sfPoint));
    /* Cleanup overlap partition */
    PetscCall(DMLabelDestroy(&lblOverlap));
    PetscCall(PetscSFDestroy(&sfOverlap));

    PYLITH_METHOD_RETURN(0);
} // distributeOverlap


// End of file
