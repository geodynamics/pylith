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

// ----------------------------------------------------------------------
const char* pylith::meshio::FieldFilterProject::_pyreComponent = "fieldfilterproject";

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::FieldFilterProject::FieldFilterProject(void) :
    _fieldP1(NULL) {
    PyreComponent::name(_pyreComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::FieldFilterProject::~FieldFilterProject(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::FieldFilterProject::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    FieldFilter::deallocate();

    delete _fieldP1; _fieldP1 = NULL;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::FieldFilterProject::FieldFilterProject(const FieldFilterProject& f) :
    FieldFilter(f),
    _fieldP1(NULL)
{}

// ----------------------------------------------------------------------
// Create copy of filter.
pylith::meshio::FieldFilter*
pylith::meshio::FieldFilterProject::clone(void) const {
    return new FieldFilterProject(*this);
} // clone

// ----------------------------------------------------------------------
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

// ----------------------------------------------------------------------
// Filter field.
pylith::topology::Field*
pylith::meshio::FieldFilterProject::filter(pylith::topology::Field* fieldIn) {
    PYLITH_METHOD_BEGIN;

    if (!fieldIn) {
        PYLITH_METHOD_RETURN(NULL);
    } // if

    assert(fieldIn);
    if (_fieldP1 && !pylith::topology::FieldOps::layoutsMatch(*fieldIn, *_fieldP1)) {
        delete _fieldP1; _fieldP1 = NULL;
        delete _passThruFns; _passThruFns = NULL;
    } // if
    if (!_fieldP1) {
        _fieldP1 = new pylith::topology::Field(fieldIn->mesh()); assert(_fieldP1);

        // Set discretization of all subfields to basis order.
        pylith::topology::Field::Discretization feP1;
        feP1.basisOrder = _basisOrder;

        const pylith::string_vector& subfieldNames = fieldIn->subfieldNames();
        const size_t numSubfields = subfieldNames.size();
        for (size_t i = 0; i < numSubfields; ++i) {
            const pylith::topology::Field::SubfieldInfo& info = fieldIn->subfieldInfo(subfieldNames[i].c_str());
            feP1.isBasisContinuous = info.fe.isBasisContinuous;
            feP1.feSpace = info.fe.feSpace;
            feP1.quadOrder = info.fe.quadOrder;

            if (info.fe.quadOrder < _basisOrder) {
                PYLITH_COMPONENT_WARNING(
                    "Projecting subfield '"
                    << info.description.label << "' in field ''" << fieldIn->label() << "'' from basis order "
                    << info.fe.basisOrder << " to basis order " << _basisOrder
                    << " with quadrature order " << info.fe.quadOrder << " will result in under integration of the "
                    << "subfield. Accurate projection requires a quadrature order of at least " << _basisOrder << "."
                    );
            } // if

            _fieldP1->subfieldAdd(info.description, feP1);
        } // for
        _fieldP1->subfieldsSetup();
        _fieldP1->allocate();
    } // if
    _fieldP1->label(fieldIn->label());
    _fieldP1->dimensionalizeOkay(true);

    if (!_passThruFns) {
        const size_t numSubfields = _fieldP1->subfieldNames().size();
        _passThruFns = (numSubfields > 0) ? new PetscPointFunc[numSubfields] : NULL;
        for (size_t i = 0; i < numSubfields; ++i) {
            _passThruFns[i] = passThruSoln;
        } // for
    } // if

    const PylithReal t = 0.0;
    PetscDM dmFieldP1 = _fieldP1->dmMesh();
    PetscErrorCode err = DMProjectFieldLocal(dmFieldP1, t, _fieldP1->localVector(), _passThruFns, INSERT_VALUES, _fieldP1->localVector());PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(_fieldP1);
} // filter


// ----------------------------------------------------------------------
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

    const PylithInt sEnd = sOff[numS];
    for (PylithInt i = 0; i < sEnd; ++i) {
        field[i] = s[i];
    } // for
} // passThruSoln


// End of file
