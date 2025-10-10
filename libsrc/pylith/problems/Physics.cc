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

#include "pylith/problems/Physics.hh" // implementation of class methods

#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/scales/Scales.hh" // USES Scales

#include "pylith/utils/error.hh" // USES PYLITH_JMETHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> // USES std::runtime_error

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::Physics::Physics(void) :
    _scales(NULL),
    _labelName(pylith::topology::Mesh::cells_label_name),
    _labelValue(1),
    _observers(new pylith::problems::ObserversPhysics) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::Physics::~Physics(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Physics::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    delete _scales;_scales = NULL;
    delete _observers;_observers = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set name of label marking material.
void
pylith::problems::Physics::setLabelName(const char* value) {
    PYLITH_COMPONENT_DEBUG("setLabelName(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for material label.");
    } // if

    _labelName = value;
} // setLabelName


// ------------------------------------------------------------------------------------------------
// Get name of label marking material.
const char*
pylith::problems::Physics::getLabelName(void) const {
    return _labelName.c_str();
} // getLabelName


// ------------------------------------------------------------------------------------------------
// Set value of label marking material.
void
pylith::problems::Physics::setLabelValue(const int value) {
    _labelValue = value;
} // setLabelValue


// ------------------------------------------------------------------------------------------------
// Get value of label marking material.
int
pylith::problems::Physics::getLabelValue(void) const {
    return _labelValue;
} // getLabelValue


// ------------------------------------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::problems::Physics::setScales(const pylith::scales::Scales& dim) {
    PYLITH_COMPONENT_DEBUG("setScales(dim="<<typeid(dim).name()<<")");

    if (!_scales) {
        _scales = new pylith::scales::Scales(dim);
    } else {
        *_scales = dim;
    } // if/else
} // setScales


// ------------------------------------------------------------------------------------------------
/** Get manager of scales used to nondimensionalize problem.
 *
 * @param dim Nondimensionalizer.
 */
const pylith::scales::Scales&
pylith::problems::Physics::getScales(void) const {
    assert(_scales);
    return *_scales;
} // getScales


// ------------------------------------------------------------------------------------------------
// Set formulation for equations.
void
pylith::problems::Physics::setFormulation(const FormulationEnum value) {
    _formulation = value;
} // setFormulation


// ------------------------------------------------------------------------------------------------
// Set database for auxiliary field.
void
pylith::problems::Physics::setAuxiliaryFieldDB(spatialdata::spatialdb::SpatialDB* value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setAuxiliaryFieldDB(value="<<value<<")");

    pylith::feassemble::AuxiliaryFactory* factory = _getAuxiliaryFactory();assert(factory);
    factory->setQueryDB(value);

    PYLITH_METHOD_END;
} // setAuxiliaryFieldDB


// ------------------------------------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::problems::Physics::setAuxiliarySubfieldDiscretization(const char* subfieldName,
                                                              const int basisOrder,
                                                              const int quadOrder,
                                                              const int dimension,
                                                              const pylith::topology::FieldBase::CellBasis cellBasis,
                                                              const pylith::topology::FieldBase::SpaceEnum feSpace,
                                                              const bool isBasisContinuous) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setAuxiliarySubfieldDiscretization(subfieldName="<<subfieldName<<", basisOrder="<<basisOrder<<", quadOrder="<<quadOrder<<", dimension="<<dimension<<", cellBasis="<<cellBasis<<", isBasisContinuous="<<isBasisContinuous<<")");

    pylith::feassemble::AuxiliaryFactory* factory = _getAuxiliaryFactory();assert(factory);
    const bool isFaultOnly = false;
    factory->setSubfieldDiscretization(subfieldName, basisOrder, quadOrder, dimension, isFaultOnly, cellBasis, feSpace,
                                       isBasisContinuous);

    PYLITH_METHOD_END;
} // setAuxSubfieldDiscretization


// ------------------------------------------------------------------------------------------------
// Set discretization information for derived subfield.
void
pylith::problems::Physics::setDerivedSubfieldDiscretization(const char* subfieldName,
                                                            const int basisOrder,
                                                            const int quadOrder,
                                                            const int dimension,
                                                            const pylith::topology::FieldBase::CellBasis cellBasis,
                                                            const pylith::topology::FieldBase::SpaceEnum feSpace,
                                                            const bool isBasisContinuous) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setDerivedSubfieldDiscretization(subfieldName="<<subfieldName<<", basisOrder="<<basisOrder<<", quadOrder="<<quadOrder<<", dimension="<<dimension<<", cellBasis="<<cellBasis<<", isBasisContinuous="<<isBasisContinuous<<")");

    pylith::topology::FieldFactory* factory = _getDerivedFactory();assert(factory);
    const bool isFaultOnly = false;
    factory->setSubfieldDiscretization(subfieldName, basisOrder, quadOrder, dimension, isFaultOnly, cellBasis, feSpace,
                                       isBasisContinuous);

    PYLITH_METHOD_END;
} // setDerivedSubfieldDiscretization


// ----------------------------------------------------------------------
// Register observer to receive notifications.
void
pylith::problems::Physics::registerObserver(pylith::problems::ObserverPhysics* observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("registerObserver(observer="<<typeid(observer).name()<<")");

    assert(_observers);
    _observers->registerObserver(observer);

    PYLITH_METHOD_END;
} // registerObserver


// ----------------------------------------------------------------------
// Remove observer from receiving notifications.
void
pylith::problems::Physics::removeObserver(pylith::problems::ObserverPhysics* observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("removeObserver(observer="<<typeid(observer).name()<<")");

    assert(_observers);
    _observers->removeObserver(observer);

    PYLITH_METHOD_END;
} // removeObserver


// ------------------------------------------------------------------------------------------------
// Get observers receiving notifications of physics updates.
pylith::problems::ObserversPhysics*
pylith::problems::Physics::getObservers(void) {
    return _observers;
} // getObservers


// ------------------------------------------------------------------------------------------------
// Get constants used in kernels (point-wise functions).
const pylith::real_array&
pylith::problems::Physics::getKernelConstants(const PylithReal dt) {
    _updateKernelConstants(dt);

    return _kernelConstants;
} // getKernelConstants


// ------------------------------------------------------------------------------------------------
// Create diagnostic field.
pylith::topology::Field*
pylith::problems::Physics::createDiagnosticField(const pylith::topology::Field& solution,
                                                 const pylith::topology::Mesh& physicsMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDiagnosticField(solution="<<solution.getLabel()<<", physicsMesh=)"<<typeid(physicsMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(NULL);

}


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::problems::Physics::createDerivedField(const pylith::topology::Field& solution,
                                              const pylith::topology::Mesh& physicsMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getLabel()<<", physicsMesh=)"<<typeid(physicsMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(NULL);
} // createDerivedField


// ------------------------------------------------------------------------------------------------
// Update auxiliary field for given time.
void
pylith::problems::Physics::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateAuxiliaryField(auxiliaryField="<<auxiliaryField<<", t="<<t<<") empty method");

    PYLITH_METHOD_END;
} // updateAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Get derived factory associated with physics.
pylith::topology::FieldFactory*
pylith::problems::Physics::_getDerivedFactory(void) {
    return NULL;
} // _getDerivedFactory


// ------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::problems::Physics::_updateKernelConstants(const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_updateKernelConstants(dt="<<dt<<")");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // _updateKernelConstants


// End of file
