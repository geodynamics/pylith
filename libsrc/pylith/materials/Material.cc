// -*- C++ -*-
//
// ---------------------------------------------------------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ---------------------------------------------------------------------------------------------------------------------
//

#include <portinfo>

#include "Material.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactory.hh" // TEMPORARY

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> \
    // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::Material::Material(const int dimension) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactory),
    _dimension(dimension),
    _materialId(0),
    _descriptiveLabel("") {
    //
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::Material::~Material(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::Material::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::problems::Physics::deallocate();
    delete _auxiliaryFactory;_auxiliaryFactory = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Get spatial dimension of material.
int
pylith::materials::Material::getDimension(void) const {
    return _dimension;
} // getDimension


// ---------------------------------------------------------------------------------------------------------------------
// Set value of label material-id used to identify material cells.
void
pylith::materials::Material::setMaterialId(const int value) {
    PYLITH_COMPONENT_DEBUG("setMmaterialId(value="<<value<<")");

    _materialId = value;
} // setMaterialId


// ---------------------------------------------------------------------------------------------------------------------
// Get value of label material-id used to identify material cells.
int
pylith::materials::Material::getMaterialId(void) const {
    return _materialId;
} // getMaterialId


// ---------------------------------------------------------------------------------------------------------------------
// Set descriptive label of material.
void
pylith::materials::Material::setDescriptiveLabel(const char* value) {
    PYLITH_COMPONENT_DEBUG("setDescriptiveLabel(value="<<value<<")");

    _descriptiveLabel = value;
} // setDescriptiveLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get label of material.
const char*
pylith::materials::Material::getDescriptiveLabel(void) const {
    return _descriptiveLabel.c_str();
} // getDescriptiveLabel


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
pylith::feassemble::Constraint*
createConstraint(const pylith::topology::Field& solution) {
    return NULL;
} // createConstraint


#if 0
// ---------------------------------------------------------------------------------------------------------------------
// Get physical property parameters and initial state (if used) from database.
void
pylith::materials::Material::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("intialize(solution="<<solution.label()<<")");

    // Get cells associated with material
    const pylith::topology::Mesh& mesh = solution.mesh();
    PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
    pylith::topology::CoordsVisitor::optimizeClosure(dmMesh);

    const bool includeOnlyCells = true;
    delete _materialIS;_materialIS = new pylith::topology::StratumIS(dmMesh, "material-id", _id, includeOnlyCells);assert(_materialIS);

    delete _auxField;_auxField = new pylith::topology::Field(mesh);assert(_auxField);
    _auxField->label("auxiliary subfields");
    _auxFieldSetup();
    _auxField->subfieldsSetup();
    pylith::topology::FieldOps::checkDiscretization(solution, *_auxField);
    _auxField->allocate();
    _auxField->zeroLocal();

    assert(_normalizer);
    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory();assert(factory);
    factory->initializeSubfields();

    //_auxField->view("MATERIAL AUXILIARY FIELD"); // :DEBUG: TEMPORARY
    const bool infoOnly = true;
    notifyObservers(0.0, 0, solution, infoOnly);

    PYLITH_METHOD_END;
} // initialize


#endif

// End of file
