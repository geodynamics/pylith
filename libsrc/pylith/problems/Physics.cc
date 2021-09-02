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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Physics.hh" // implementation of class methods

#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_JMETHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::Physics::Physics(void) :
    _normalizer(NULL),
    _observers(new pylith::problems::ObserversPhysics)
{}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::Physics::~Physics(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Physics::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    delete _normalizer;_normalizer = NULL;
    delete _observers;_observers = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::problems::Physics::setNormalizer(const spatialdata::units::Nondimensional& dim) {
    PYLITH_COMPONENT_DEBUG("setNormalizer(dim="<<typeid(dim).name()<<")");

    if (!_normalizer) {
        _normalizer = new spatialdata::units::Nondimensional(dim);
    } else {
        *_normalizer = dim;
    } // if/else
} // setNormalizer


// ---------------------------------------------------------------------------------------------------------------------
/** Get manager of scales used to nondimensionalize problem.
 *
 * @param dim Nondimensionalizer.
 */
const spatialdata::units::Nondimensional&
pylith::problems::Physics::getNormalizer(void) const {
    assert(_normalizer);
    return *_normalizer;
} // getNormalizer


// ---------------------------------------------------------------------------------------------------------------------
// Set formulation for equations.
void
pylith::problems::Physics::setFormulation(const FormulationEnum value) {
    _formulation = value;
} // setFormulation


// ---------------------------------------------------------------------------------------------------------------------
// Set database for auxiliary field.
void
pylith::problems::Physics::setAuxiliaryFieldDB(spatialdata::spatialdb::SpatialDB* value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setAuxiliaryFieldDB(value="<<value<<")");

    pylith::feassemble::AuxiliaryFactory* factory = _getAuxiliaryFactory();assert(factory);
    factory->setQueryDB(value);

    PYLITH_METHOD_END;
} // setAuxiliaryFieldDB


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
// Get observers receiving notifications of physics updates.
pylith::problems::ObserversPhysics*
pylith::problems::Physics::getObservers(void) {
    return _observers;
} // getObservers


// ---------------------------------------------------------------------------------------------------------------------
// Get constants used in kernels (point-wise functions).
const pylith::real_array&
pylith::problems::Physics::getKernelConstants(const PylithReal dt) {
    _updateKernelConstants(dt);

    return _kernelConstants;
} // getKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary field for given time.
void
pylith::problems::Physics::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                const PylithReal t,
                                                const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateAuxiliaryField(auxiliaryField="<<auxiliaryField<<", t="<<t<<", dt="<<dt<<") empty method");

    PYLITH_METHOD_END;
} // updateAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Get derived factory associated with physics.
pylith::topology::FieldFactory*
pylith::problems::Physics::_getDerivedFactory(void) {
    return NULL;
} // _getDerivedFactory


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::problems::Physics::_updateKernelConstants(const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_updateKernelConstants(dt="<<dt<<")");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // _updateKernelConstants


// End of file
