// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "Neumann.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/feassemble/IntegratorObserver.hh" // USES IntegratorObserver
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> \
    // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::Neumann::Neumann(void) :
    _scaleName("pressure")
{ // constructor
    _description.label = "unknown";
    _description.vectorFieldType = pylith::topology::FieldBase::OTHER;
    _description.numComponents = 0;
    _description.scale = 1.0;
    _description.validator = NULL;

} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::Neumann::~Neumann(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::Neumann::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    IntegratorBoundary::deallocate();

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Name of scale associated with Neumann boundary condition (e.g., pressure for elasticity).
void
pylith::bc::Neumann::scaleName(const char* value) {
    if (( value == std::string("length")) ||
        ( value == std::string("time")) ||
        ( value == std::string("pressure")) ||
        ( value == std::string("density")) ||
        ( value == std::string("pressure")) ) {
        _scaleName = value;
    } else {
        std::ostringstream msg;
        msg << "Unknown name of scale ("<<value<<") for Neumann boundary condition '" << label() << "'.";
        throw std::runtime_error(msg.str());
    } // if
} // scaleName

// End of file
