// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/topology/Distributor.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/faults/FaultCohesive.hh" // USES FaultCohesive
#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/utils/journals.hh" // pythia::journal

#include "petsc/private/dmpleximpl.h"

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

            static
            void write(pylith::meshio::DataWriter* writer,
                       const pylith::topology::Mesh& mesh);

            // :KLUDGE: Local copy of DMPlexCopy_Internal()
            static
            PetscErrorCode DMPlexCopy(DM dmin,
                                      PetscBool copyPeriodicity,
                                      PetscBool copyOverlap,
                                      DM dmout);


        }; // _Distributor
    } // topology
} // pylith

// ------------------------------------------------------------------------------------------------
static PetscErrorCode
DMPlexGetAdjacency_SupportOnly_Internal(DM dm,
                                        PetscInt p,
                                        PetscInt *adjSize,
                                        PetscInt adj[],
                                        void *ctx) {
    const PetscInt *support = NULL;
    PetscInt maxAdjSize = *adjSize;

    PetscFunctionBeginHot;
    PylithCallPetsc(DMPlexGetSupportSize(dm, p, adjSize));
    PylithCallPetsc(DMPlexGetSupport(dm, p, &support));
    PetscCheck(*adjSize <= maxAdjSize, PETSC_COMM_SELF, PETSC_ERR_PLIB, "Invalid mesh exceeded adjacency allocation (%" PetscInt_FMT ")", maxAdjSize);
    for (PetscInt s = 0; s < *adjSize; ++s) {
        adj[s] = support[s];
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::topology::Distributor::Distributor(void) :
    _writer(nullptr),
    _partitioner("parmetis"),
    _useEdgeWeighting(true) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::topology::Distributor::~Distributor(void) {
    _writer = nullptr; // :KLUDGE: Use shared pointer
}


// ------------------------------------------------------------------------------------------------
// Set data writer.
void
pylith::topology::Distributor::setDataWriter(pylith::meshio::DataWriter* writer) {
    _writer = writer;
} // setDataWriter


// ------------------------------------------------------------------------------------------------
// Set edge weighting.
void
pylith::topology::Distributor::setUseEdgeWeighting(const bool flag) {
    _useEdgeWeighting = flag;
} // setUseEdgeWeighting


// ------------------------------------------------------------------------------------------------
// Set partitioner.
void
pylith::topology::Distributor::setPartitioner(const char* partitioner) {
    if ((0 == strcasecmp(partitioner, "parmetis")) || (0 == strcasecmp(partitioner, "chaco")) || (0 == strcasecmp(partitioner, "simple"))) {
        _partitioner = partitioner;
    } else {
        std::ostringstream msg;
        msg << "Unknown partitioner '" << partitioner << "'. Partitioner must be 'parmetis', 'chaco', or 'simple'.";
        throw std::runtime_error(msg.str());
    } // if/else
} // setPartitioner


// ------------------------------------------------------------------------------------------------
// Distribute mesh.
pylith::topology::Mesh*
pylith::topology::Distributor::distribute(const pylith::topology::Mesh& mesh,
                                          const std::vector<pylith::faults::FaultCohesive*>& faults) const {
    PYLITH_METHOD_BEGIN;


    const int commRank = mesh.getCommRank();
    if (0 == commRank) {
        PYLITH_COMPONENT_INFO("Partitioning mesh using PETSc '" << _partitioner << "' partitioner.");
    } // if

    PetscPartitioner partitioner = nullptr;
    PetscDM dmOrig = mesh.getDM();assert(dmOrig);
    PylithCallPetsc(DMPlexGetPartitioner(dmOrig, &partitioner));
    PylithCallPetsc(PetscPartitionerSetType(partitioner, _partitioner.c_str()));
    if ((_partitioner == std::string("parmetis")) && _useEdgeWeighting) {
        PylithCallPetsc(PetscOptionsSetValue(NULL, "-petscpartitioner_use_vertex_weights", "true"));
        PylithCallPetsc(PetscPartitionerSetFromOptions(partitioner));
    } // if

    if (0 == commRank) {
        PYLITH_COMPONENT_INFO("Distributing partitioner mesh.");
    } // if

    PetscDM dmTmp = NULL, dmNew = NULL;
    const PetscInt overlap = 0;
    PylithCallPetsc(DMPlexDistribute(dmOrig, overlap, NULL, &dmTmp));
    pylith::topology::Mesh* meshNew = nullptr;
    if (dmTmp) {
#if 1
        PylithCallPetsc(Distributor::distributeOverlap(&dmNew, dmTmp, faults));
        PylithCallPetsc(DMDestroy(&dmTmp));
        PylithCallPetsc(DMPlexDistributeSetDefault(dmNew, PETSC_FALSE));
        PylithCallPetsc(DMPlexReorderCohesiveSupports(dmNew));
        PylithCallPetsc(DMViewFromOptions(dmNew, NULL, "-pylith_dist_dm_view"));
        meshNew = new pylith::topology::Mesh(dmNew, mesh);assert(meshNew);
#else
        meshNew = new pylith::topology::Mesh(dmTmp, mesh);assert(meshNew);
#endif
    } else {
        PetscObjectReference(PetscObject(dmOrig));
        meshNew = new pylith::topology::Mesh(dmOrig, mesh);assert(meshNew);
    } // if/else

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        meshNew->view(":mesh_distributed.txt:ascii_info_detail");
        pylith::topology::Mesh* meshExploded = pylith::topology::MeshOps::explode(*meshNew);
        meshExploded->view(":mesh_domain_after_initialize.tex:ascii_latex");
        delete meshExploded;meshExploded = nullptr;
    } // if
    if (_writer) {
        _Distributor::write(_writer, *meshNew);
    } // if

    if (pylith::topology::MeshOps::getNumCells(*meshNew) == 0) {
        std::ostringstream msg;
        msg << "No cells are assigned to process " << commRank << " after distribution. "
            << "Either there are too many processes for the mesh or there is a topology related error.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(meshNew);
} // distribute


// ------------------------------------------------------------------------------------------------
// This is a copy of DMPlexDistributeOverlap()
PetscErrorCode
pylith::topology::Distributor::distributeOverlap(PetscDM* dmOverlap,
                                                 PetscDM dmMesh,
                                                 const std::vector<pylith::faults::FaultCohesive*>& faults) {
    PYLITH_METHOD_BEGIN;
    assert(dmOverlap);

    MPI_Comm comm;
    PetscMPIInt size, rank;
    PetscSection rootSection, leafSection;
    PetscIS rootrank, leafrank;
    PetscDM dmCoord;
    PetscDMLabel lblOverlap;
    PetscSF sfOverlap, sfStratified, sfPoint;
    PetscInt dim;
    const char* meshName;

    const size_t numFaults = faults.size();
    PylithCallPetsc(PetscObjectGetComm((PetscObject)dmMesh,&comm));
    PetscCallMPI(MPI_Comm_size(comm, &size));
    PetscCallMPI(MPI_Comm_rank(comm, &rank));
    if ((0 == numFaults) || (1 == size)) {
        PylithCallPetsc(PetscObjectReference((PetscObject)dmMesh));
        *dmOverlap = dmMesh;
        PYLITH_METHOD_RETURN(0);
    } // if

    PetscDMLabel* ovIncludeLabels = (numFaults > 0) ? new PetscDMLabel[numFaults] : NULL;
    PetscInt* ovIncludeLabelValues = (numFaults > 0) ? new PetscInt[numFaults] : NULL;
    PetscDMLabel* ovExcludeLabels = (numFaults > 0) ? new PetscDMLabel[numFaults] : NULL;
    PetscInt* ovExcludeLabelValues = (numFaults > 0) ? new PetscInt[numFaults] : NULL;

    PylithCallPetsc(DMGetDimension(dmMesh, &dim));
    for (size_t i = 0; i < numFaults; ++i) {
        const char* surfaceLabelName = faults[i]->getSurfaceLabelName();
        PylithCallPetsc(DMGetLabel(dmMesh, surfaceLabelName, &ovIncludeLabels[i]));
        ovIncludeLabelValues[i] = dim - 1;

        const char* cohesiveLabelName = faults[i]->getCohesiveLabelName();
        PylithCallPetsc(DMGetLabel(dmMesh, cohesiveLabelName, &ovExcludeLabels[i]));
        ovExcludeLabelValues[i] = faults[i]->getCohesiveLabelValue();
    } // for

    /* Compute point overlap with neighboring processes on the distributed DM */
    PylithCallPetsc(PetscSectionCreate(comm, &rootSection));
    PylithCallPetsc(PetscSectionCreate(comm, &leafSection));
    PylithCallPetsc(DMPlexSetAdjacencyUser(dmMesh, DMPlexGetAdjacency_SupportOnly_Internal, NULL));
    PylithCallPetsc(DMPlexDistributeOwnership(dmMesh, rootSection, &rootrank, leafSection, &leafrank));
    PylithCallPetsc(DMPlexCreateOverlapLabelFromLabels(dmMesh, numFaults, ovIncludeLabels, ovIncludeLabelValues,
                                                       numFaults, ovExcludeLabels, ovExcludeLabelValues, rootSection, rootrank, leafSection, leafrank, &lblOverlap));
    PylithCallPetsc(DMPlexSetAdjacencyUser(dmMesh, NULL, NULL));

    delete[] ovIncludeLabels;ovIncludeLabels = NULL;
    delete[] ovIncludeLabelValues;ovIncludeLabelValues = NULL;
    delete[] ovExcludeLabels;ovExcludeLabels = NULL;
    delete[] ovExcludeLabelValues;ovExcludeLabelValues = NULL;

    /* Convert overlap label to stratified migration SF */
    PylithCallPetsc(DMPlexPartitionLabelCreateSF(dmMesh, lblOverlap, PETSC_TRUE, &sfOverlap));
    PylithCallPetsc(DMPlexStratifyMigrationSF(dmMesh, sfOverlap, &sfStratified));
    PylithCallPetsc(PetscSFDestroy(&sfOverlap));
    sfOverlap = sfStratified;
    PylithCallPetsc(PetscObjectSetName((PetscObject) sfOverlap, "Overlap SF "));
    PylithCallPetsc(PetscSFSetFromOptions(sfOverlap));

    PylithCallPetsc(PetscSectionDestroy(&rootSection));
    PylithCallPetsc(PetscSectionDestroy(&leafSection));
    PylithCallPetsc(ISDestroy(&rootrank));
    PylithCallPetsc(ISDestroy(&leafrank));

    /* Build the overlapping DM */
    PylithCallPetsc(DMPlexCreate(comm, dmOverlap));
    PylithCallPetsc(PetscObjectGetName((PetscObject) dmMesh, &meshName));
    PylithCallPetsc(PetscObjectSetName((PetscObject) *dmOverlap, meshName));
    PylithCallPetsc(DMPlexMigrate(dmMesh, sfOverlap, *dmOverlap));
    PylithCallPetsc(_Distributor::DMPlexCopy(dmMesh, PETSC_TRUE, PETSC_FALSE, *dmOverlap));

    PetscReal scale = 0.0;
    PylithCallPetsc(DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &scale));
    PylithCallPetsc(DMPlexSetScale(*dmOverlap, PETSC_UNIT_LENGTH, scale));

    /* Store the overlap in the new DM */
    PylithCallPetsc(DMPlexSetOverlap(*dmOverlap, dmMesh, 1));
    /* Build the new point SF */
    PylithCallPetsc(DMPlexCreatePointSF(*dmOverlap, sfOverlap, PETSC_FALSE, &sfPoint));
    PylithCallPetsc(DMSetPointSF(*dmOverlap, sfPoint));
    PylithCallPetsc(DMGetCoordinateDM(*dmOverlap, &dmCoord));
    if (dmCoord) { PylithCallPetsc(DMSetPointSF(dmCoord, sfPoint));}
    PylithCallPetsc(PetscSFDestroy(&sfPoint));
    /* Cleanup overlap partition */
    PylithCallPetsc(DMLabelDestroy(&lblOverlap));
    PylithCallPetsc(PetscSFDestroy(&sfOverlap));

    PYLITH_METHOD_RETURN(0);
} // distributeOverlap


// ------------------------------------------------------------------------------------------------
// Write partitioning info for distributed mesh.
void
pylith::topology::_Distributor::write(meshio::DataWriter* const writer,
                                      const topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    // Setup and allocate PETSc vector
    const int commRank = mesh.getCommRank();
    PylithScalar rankReal = PylithReal(commRank);

    pylith::topology::Field partitionField(mesh);
    partitionField.setLabel("partition ");

    pylith::topology::Field::Description description;
    description.label = "partition ";
    description.alias = "partition ";
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = "rank ";
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
    const int refineLevels = 0;
    pylith::meshio::OutputSubfield* outputField =
        pylith::meshio::OutputSubfield::create(partitionField, mesh, "partition ", basisOrder, refineLevels);
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
PetscErrorCode
pylith::topology::_Distributor::DMPlexCopy(DM dmin,
                                           PetscBool copyPeriodicity,
                                           PetscBool copyOverlap,
                                           DM dmout) {
    const PetscReal     *maxCell, *Lstart, *L;
    VecType vecType;
    MatType matType;
    PetscBool dist, useCeed, balance_partition;
    DMReorderDefaultFlag reorder;
    MatOrderingType otype;


    PetscFunctionBegin;
    if (dmin == dmout) { PetscFunctionReturn(PETSC_SUCCESS);}
    PylithCallPetsc(DMGetVecType(dmin, &vecType));
    PylithCallPetsc(DMSetVecType(dmout, vecType));
    PylithCallPetsc(DMGetMatType(dmin, &matType));
    PylithCallPetsc(DMSetMatType(dmout, matType));
    if (copyPeriodicity) {
        PylithCallPetsc(DMGetPeriodicity(dmin, &maxCell, &Lstart, &L));
        PylithCallPetsc(DMSetPeriodicity(dmout, maxCell, Lstart, L));
        PylithCallPetsc(DMLocalizeCoordinates(dmout));
    }
    PylithCallPetsc(DMPlexDistributeGetDefault(dmin, &dist));
    PylithCallPetsc(DMPlexDistributeSetDefault(dmout, dist));
    PylithCallPetsc(DMPlexReorderGetDefault(dmin, &reorder));
    PylithCallPetsc(DMPlexReorderSetDefault(dmout, reorder));

    PylithCallPetsc(DMReorderSectionGetDefault(dmin, &reorder));
    PylithCallPetsc(DMReorderSectionSetDefault(dmout, reorder));
    PylithCallPetsc(DMReorderSectionGetType(dmin, &otype));
    PylithCallPetsc(DMReorderSectionSetType(dmout, otype));

    PylithCallPetsc(DMPlexGetUseCeed(dmin, &useCeed));
    PylithCallPetsc(DMPlexSetUseCeed(dmout, useCeed));
    PylithCallPetsc(DMPlexGetPartitionBalance(dmin, &balance_partition));
    PylithCallPetsc(DMPlexSetPartitionBalance(dmout, balance_partition));
    PylithCallPetsc(DMPlexCopyFlags(dmin, dmout));
    PetscFunctionReturn(PETSC_SUCCESS);
}


// End of file
