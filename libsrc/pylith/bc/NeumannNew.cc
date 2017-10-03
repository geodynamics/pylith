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

#include "NeumannNew.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::NeumannNew::NeumannNew(void) :
    _boundaryMesh(NULL),
    _bcKernel(NULL)
{ // constructor
    _description.label = "unknown";
    _description.vectorFieldType = pylith::topology::Field::OTHER;
    _description.numComponents = 0;
    _description.scale = 1.0;
    _description.validator = NULL;
    const pylith::topology::FieldBase::Discretization defaultInfo = {-1, -1, true, pylith::topology::FieldBase::POLYNOMIAL_SPACE};
    _auxFieldsFEInfo["default"] = defaultInfo;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::NeumannNew::~NeumannNew(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::NeumannNew::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    IntegratorPointwise::deallocate();

    delete _boundaryMesh; _boundaryMesh = NULL;
    _bcKernel = NULL;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::NeumannNew::initialize(const pylith::topology::Field& solution)
{ // initialize
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(solution="<<solution.label()<<")");

    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    _description = info.description;

    _setFEKernelsRHSResidual(solution);

    _boundaryMesh = new topology::Mesh(solution.mesh(), _label.c_str()); assert(_boundaryMesh);
    PetscDM dmBoundary = _boundaryMesh->dmMesh(); assert(dmBoundary);
    pylith::topology::CoordsVisitor::optimizeClosure(dmBoundary);

    delete _auxFields; _auxFields = new pylith::topology::Field(*_boundaryMesh); assert(_auxFields);
    delete _auxFieldsQuery; _auxFieldsQuery = new pylith::topology::FieldQuery(*_auxFields); assert(_auxFieldsQuery);
    _auxFields->label("auxiliary fields");
    _auxFieldsSetup();
    _auxFields->subfieldsSetup();
    _auxFields->allocate();
    _auxFields->zeroLocal();

    if (_auxFieldsDB) {
        assert(_normalizer);
        _auxFieldsQuery->openDB(_auxFieldsDB, _normalizer->lengthScale());
        _auxFieldsQuery->queryDB();
        _auxFieldsQuery->closeDB(_auxFieldsDB);
    } else { // else
        PYLITH_COMPONENT_ERROR("Unknown case for setting up auxiliary fields.");
        throw std::logic_error("Unknown case for setting up auxiliary fields.");
    } // if/else
    _auxFields->view("AUXILIARY FIELDS"); // :DEBUGGING: TEMPORARY

    PYLITH_METHOD_END;
} // initialize


// End of file
