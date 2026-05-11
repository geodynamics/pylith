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

#include "pylith/meshio/MeshIO.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps

#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_INFO
#include "pylith/utils/Exceptions.hh" // USES Exception
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept>

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIO::MeshIO(void) :
    _mesh(nullptr),
    _cs(nullptr) {}


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIO::~MeshIO(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::MeshIO::deallocate(void) {
    _mesh = nullptr; // :KLUDGE: Use shared pointer
    _cs = nullptr; // :KLUDGE: Use shared pointer
} // deallocate


// ----------------------------------------------------------------------
// Read mesh from file.
void
pylith::meshio::MeshIO::read(pylith::topology::Mesh* mesh,
                             const bool checkTopology) {
    PYLITH_METHOD_BEGIN;

    assert(mesh);
    assert(!_mesh);

    _mesh = mesh;
    _read();

    // Check for bounding box with positive volume.
    PylithReal cmin[3];
    PylithReal cmax[3];
    PylithCallPetsc(DMGetBoundingBox(_mesh->getDM(), cmin, cmax));
    const PetscInt dim = _mesh->getDimension();
    PylithReal volume = 1.0;
    for (int i = 0; i < dim; ++i) {
        volume *= cmax[i] - cmin[i];
    } // for
    std::ostringstream msg;
    msg << "Domain bounding box:";
    for (int i = 0; i < dim; ++i) {
        msg << "\n    (" << cmin[i] << ", " << cmax[i] << ")";
    } // for
    PYLITH_COMPONENT_INFO_ROOT(pylith::journal::application_flow, msg.str());
    const PetscReal tolerance = 1.0e-8;
    if (volume < tolerance) {
        PYLITH_ERROR(pylith::ValueError, pylith::journal::user_input,
                     "Domain bounding box volume (" << volume << ") is less than minimum tolerance ("
                                                    << tolerance << "). This usually means you are trying to use a 2D mesh in 3D. Check that you are exporting "
                                                    << "your mesh from the mesh generation software correctly and that your have specified the correct "
                                                    << " coordinate system for the problem.");
    } // if

    // Check mesh consistency
    if (checkTopology) {
        pylith::topology::MeshOps::checkTopology(*_mesh);
    } // if

    if (_cs) {
        _mesh->setCoordSys(_cs);
    } // if

    pythia::journal::debug_t debug(pylith::journal::mesh_full_detail);
    if (debug.state()) {
        _mesh->view("::ascii_info_detail");
        _mesh->view(":mesh.tex:ascii_latex");
    } // if

    _mesh = NULL;

    PYLITH_METHOD_END;
} // read


// ----------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIO::write(pylith::topology::Mesh* const mesh) {
    PYLITH_METHOD_BEGIN;

    assert(mesh);
    assert(!_mesh);

    _mesh = mesh;
    _write();
    _mesh = nullptr;

    PYLITH_METHOD_END;
} // write


// ----------------------------------------------------------------------
// Set coordinate system used in mesh.
void
pylith::meshio::MeshIO::setCoordSys(spatialdata::geocoords::CoordSys* const cs) {
    _cs = cs;
} // setCoordSys


// End of file
