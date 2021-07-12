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

#include "MeshIOPetsc.hh" // implementation of class methods

#include "MeshBuilder.hh" // USES MeshBuilder
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "pylith/utils/array.hh" // USES scalar_array, int_array, string_vector
#include "spatialdata/utils/LineParser.hh" // USES LineParser

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <iomanip> // USES setw(), setiosflags(), resetiosflags()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <fstream> // USES std::ifstream, std::ofstream
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES std::typeid

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    namespace meshio {
    } // meshio
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOPetsc::MeshIOPetsc(void) :
    _useIndexZero(true) { // constructor
    PyreComponent::setName("meshiopetsc");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOPetsc::~MeshIOPetsc(void) { // destructor
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate Petsc and local data structures.
void
pylith::meshio::MeshIOPetsc::deallocate(void) { // deallocate
    PYLITH_METHOD_BEGIN;

    MeshIO::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Read mesh.
void
pylith::meshio::MeshIOPetsc::_read(void) { // _read
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_read()");
    
    assert(_mesh);
    MPI_Comm comm = _mesh->getComm();
    PetscErrorCode err = 0;
    PetscDM dmMesh = NULL;

    err = DMCreate(comm, &dmMesh);PYLITH_CHECK_ERROR(err);
    err = DMSetType(dmMesh, DMPLEX);PYLITH_CHECK_ERROR(err);
    err = DMSetFromOptions(dmMesh);PYLITH_CHECK_ERROR(err);
    err = DMViewFromOptions(dmMesh, NULL, "-dm_view");PYLITH_CHECK_ERROR(err);

    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    const char* const labelName = pylith::topology::Mesh::getCellsLabelName();
    for (PetscInt c = cStart; c < cEnd; ++c) {
        err = DMSetLabelValue(dmMesh, labelName, c, 1);PYLITH_CHECK_ERROR(err);
    } // for

    
    DMLabel faceLabel;
    err = DMGetLabel(dmMesh, "Face Sets", &faceLabel);PYLITH_CHECK_ERROR(err);PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelComplete(dmMesh, faceLabel);PYLITH_CHECK_ERROR(err);

    IS is;
    DMLabel xnegLabel;
    char xnegName[] = "boundary_xneg";
    err = DMCreateLabel(dmMesh, xnegName);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmMesh, xnegName, &xnegLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(faceLabel, 1, &is);PYLITH_CHECK_ERROR(err);
    err = DMLabelSetStratumIS(xnegLabel, 1, is);PYLITH_CHECK_ERROR(err);
    
    DMLabel xposLabel;
    char xposName[] = "boundary_xpos";
    err = DMCreateLabel(dmMesh, xposName);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmMesh, xposName, &xposLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(faceLabel, 3, &is);PYLITH_CHECK_ERROR(err);
    err = DMLabelSetStratumIS(xposLabel, 1, is);PYLITH_CHECK_ERROR(err);
    
    DMLabel domainLabel;
    char domainName[] = "domain_all";
    err = DMCreateLabel(dmMesh, domainName);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmMesh, domainName, &domainLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(faceLabel, 1, &is);PYLITH_CHECK_ERROR(err);
    err = DMLabelSetStratumIS(domainLabel, 1, is);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(faceLabel, 2, &is);PYLITH_CHECK_ERROR(err);
    err = DMLabelSetStratumIS(domainLabel, 1, is);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(faceLabel, 3, &is);PYLITH_CHECK_ERROR(err);
    err = DMLabelSetStratumIS(domainLabel, 1, is);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(faceLabel, 4, &is);PYLITH_CHECK_ERROR(err);
    err = DMLabelSetStratumIS(domainLabel, 1, is);PYLITH_CHECK_ERROR(err);

    _mesh->setDM(dmMesh);

    PYLITH_METHOD_END;
} // read


// ---------------------------------------------------------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIOPetsc::_write(void) const {  
    PYLITH_JOURNAL_LOGICERROR("Writing meshes via MeshIOPetsc not implemented.");
} // write


