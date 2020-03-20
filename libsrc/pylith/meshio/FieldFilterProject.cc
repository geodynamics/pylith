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

#include "FieldFilterProject.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::layoutsMatch()

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::FieldFilterProject::FieldFilterProject(void) :
    _fieldProj(NULL),
    _passThruFns(NULL),
    _basisOrder(1) { // constructor
    PyreComponent::setName("fieldfilterproject");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::FieldFilterProject::~FieldFilterProject(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::FieldFilterProject::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    FieldFilter::deallocate();

    delete _fieldProj;_fieldProj = NULL;
    delete[] _passThruFns;_passThruFns = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Copy constructor.
pylith::meshio::FieldFilterProject::FieldFilterProject(const FieldFilterProject& f) :
    FieldFilter(f),
    _fieldProj(NULL),
    _passThruFns(NULL),
    _basisOrder(f._basisOrder)
{}


// ---------------------------------------------------------------------------------------------------------------------
// Create copy of filter.
pylith::meshio::FieldFilter*
pylith::meshio::FieldFilterProject::clone(void) const {
    return new FieldFilterProject(*this);
} // clone


// ---------------------------------------------------------------------------------------------------------------------
// Set basis order for projected field.
void
pylith::meshio::FieldFilterProject::basisOrder(const int value) {
    PYLITH_METHOD_BEGIN;

    if (value < 0) {
        std::ostringstream msg;
        msg << "Basis order (" << value << ") for field filter must be nonnegative.";
        throw std::runtime_error(msg.str());
    } // if

    _basisOrder = value;

    PYLITH_METHOD_END;
} // basisOrder


// ---------------------------------------------------------------------------------------------------------------------
// Filter field.
pylith::topology::Field*
pylith::meshio::FieldFilterProject::filter(pylith::topology::Field* fieldIn) {
    PYLITH_METHOD_BEGIN;

    if (!fieldIn) {
        PYLITH_METHOD_RETURN(NULL);
    } // if

    assert(fieldIn);
    if (_fieldProj && !pylith::topology::FieldOps::layoutsMatch(*fieldIn, *_fieldProj)) {
        delete _fieldProj;_fieldProj = NULL;
        delete[] _passThruFns;_passThruFns = NULL;
    } // if

    pylith::topology::Field::Discretization feP1;
    feP1.basisOrder = _basisOrder;
    if (!_fieldProj) {
        _fieldProj = new pylith::topology::Field(fieldIn->mesh());assert(_fieldProj);

        // Set subfield information in projected field, except basis order.
        const pylith::string_vector& subfieldNames = fieldIn->subfieldNames();
        const size_t numSubfields = subfieldNames.size();
        for (size_t i = 0; i < numSubfields; ++i) {
            const pylith::topology::Field::SubfieldInfo& info = fieldIn->subfieldInfo(subfieldNames[i].c_str());
            feP1.isBasisContinuous = info.fe.isBasisContinuous;
            feP1.feSpace = info.fe.feSpace;
            feP1.quadOrder = info.fe.quadOrder;
            feP1.dimension = info.fe.dimension;

            if (info.fe.quadOrder < _basisOrder) {
                PYLITH_COMPONENT_WARNING(
                    "Projecting subfield '"
                        << info.description.label << "' in field '" << fieldIn->getLabel() << " from basis order "
                        << info.fe.basisOrder << " to basis order " << _basisOrder
                        << " with quadrature order " << info.fe.quadOrder << " will result in under integration of the "
                        << "subfield. Accurate projection requires a quadrature order of at least " << _basisOrder << "."
                    );
            } // if

            _fieldProj->subfieldAdd(info.description, feP1);
        } // for
        _fieldProj->subfieldsSetup();
        _fieldProj->createDiscretization();
        _fieldProj->allocate();
    } else {
        // Update subfield information in projected field to match input field, except basis order.
        pylith::topology::Field::Discretization feP1;
        feP1.basisOrder = _basisOrder;

        const pylith::string_vector& subfieldInNames = fieldIn->subfieldNames();
        const pylith::string_vector& subfieldProjNames = _fieldProj->subfieldNames();
        assert(subfieldInNames.size() == subfieldProjNames.size());
        const size_t numSubfields = subfieldInNames.size();
        for (size_t i = 0; i < numSubfields; ++i) {
            const pylith::topology::Field::SubfieldInfo& infoIn = fieldIn->subfieldInfo(subfieldInNames[i].c_str());
            feP1.isBasisContinuous = infoIn.fe.isBasisContinuous;
            feP1.feSpace = infoIn.fe.feSpace;
            feP1.quadOrder = infoIn.fe.quadOrder;
            feP1.dimension = infoIn.fe.dimension;

            if (infoIn.fe.quadOrder < _basisOrder) {
                PYLITH_COMPONENT_WARNING(
                    "Projecting subfield '"
                        << infoIn.description.label << "' in field ''" << fieldIn->getLabel() << "'' from basis order "
                        << infoIn.fe.basisOrder << " to basis order " << _basisOrder
                        << " with quadrature order " << infoIn.fe.quadOrder << " will result in under integration of the "
                        << "subfield. Accurate projection requires a quadrature order of at least " << _basisOrder << "."
                    );
            } // if
            _fieldProj->subfieldUpdate(subfieldProjNames[i].c_str(), subfieldInNames[i].c_str(), infoIn.description, feP1);
        } // for
    } // if/else
    _fieldProj->setLabel(fieldIn->getLabel());
    _fieldProj->dimensionalizeOkay(true);

    if (!_passThruFns) {
        const size_t numSubfields = _fieldProj->subfieldNames().size();
        _passThruFns = (numSubfields > 0) ? new PetscPointFunc[numSubfields] : NULL;
        for (size_t i = 0; i < numSubfields; ++i) {
            _passThruFns[i] = passThruSoln;
        } // for
    } // if

    PetscErrorCode err = 0;
    const PylithReal t = 0.0;
    err = DMProjectFieldLocal(_fieldProj->dmMesh(), t, fieldIn->localVector(), _passThruFns, INSERT_VALUES,
                              _fieldProj->localVector());PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(_fieldProj);
} // filter


// ---------------------------------------------------------------------------------------------------------------------
// Identify function kernel.
void
pylith::meshio::FieldFilterProject::passThruSoln(const PylithInt dim,
                                                 const PylithInt numS,
                                                 const PylithInt numA,
                                                 const PylithInt sOff[],
                                                 const PylithInt sOff_x[],
                                                 const PylithScalar s[],
                                                 const PylithScalar s_t[],
                                                 const PylithScalar s_x[],
                                                 const PylithInt aOff[],
                                                 const PylithInt aOff_x[],
                                                 const PylithScalar a[],
                                                 const PylithScalar a_t[],
                                                 const PylithScalar a_x[],
                                                 const PylithReal t,
                                                 const PylithScalar x[],
                                                 const PylithInt numConstants,
                                                 const PylithScalar constants[],
                                                 PylithScalar field[]) {
    assert(s);
    assert(field);
    assert(1 == numS);

    const PylithInt sEnd = sOff[1];
    for (PylithInt i = sOff[0]; i < sEnd; ++i) {
        field[i] = s[i];
    } // for
} // passThruSoln


// End of file
