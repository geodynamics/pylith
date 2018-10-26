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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "FieldOps.hh" // implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include "petscdm.h" // USES PetscDM
#include "petscds.h" // USES PetscDS

#include <cppunit/extensions/HelperMacros.h>

// ----------------------------------------------------------------------
// Create PetscFE object for discretization.
PetscFE
pylith::topology::FieldOps::createFE(const FieldBase::Discretization& feinfo,
                                     const PetscDM dm,
                                     const bool isSimplex,
                                     const int numComponents) { // createFE
    PYLITH_METHOD_BEGIN;

    const int basisOrder = PetscMax(feinfo.basisOrder, 0);
    const int quadOrder = PetscMax(feinfo.quadOrder > 0 ? feinfo.quadOrder : basisOrder, 0);
    const PetscBool basisContinuity = feinfo.isBasisContinuous ? PETSC_TRUE : PETSC_FALSE;
    const PetscBool useTensor = isSimplex ? PETSC_FALSE : PETSC_TRUE;

    PetscErrorCode err;

    // Get spatial dimension of mesh.
    int dim = 0;
    err = DMGetDimension(dm, &dim);PYLITH_CHECK_ERROR(err);

    // Create space
    PetscSpace space = NULL;
    err = PetscSpaceCreate(PetscObjectComm((PetscObject) dm), &space);PYLITH_CHECK_ERROR(err);assert(space);
    err = PetscSpaceSetType(space, feinfo.feSpace == FieldBase::POLYNOMIAL_SPACE ? PETSCSPACEPOLYNOMIAL : PETSCSPACEPOINT);PYLITH_CHECK_ERROR(err);
    err = PetscSpaceSetNumComponents(space, numComponents);PYLITH_CHECK_ERROR(err);
    err = PetscSpaceSetDegree(space, basisOrder, PETSC_DETERMINE);
    if (feinfo.feSpace == FieldBase::POLYNOMIAL_SPACE) {
        err = PetscSpacePolynomialSetTensor(space, useTensor);PYLITH_CHECK_ERROR(err);
        err = PetscSpaceSetNumVariables(space, dim);PYLITH_CHECK_ERROR(err);
    } // if
    err = PetscSpaceSetUp(space);PYLITH_CHECK_ERROR(err);

    // Create dual space
    PetscDualSpace dualspace = NULL;
    PetscDM dmCell = NULL;
    err = PetscDualSpaceCreate(PetscObjectComm((PetscObject) dm), &dualspace);PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceCreateReferenceCell(dualspace, dim, isSimplex ? PETSC_TRUE : PETSC_FALSE, &dmCell);PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceSetDM(dualspace, dmCell);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&dmCell);PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceSetNumComponents(dualspace, numComponents);PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceSetType(dualspace, PETSCDUALSPACELAGRANGE);PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceLagrangeSetTensor(dualspace, useTensor);PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceSetOrder(dualspace, basisOrder);PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceLagrangeSetContinuity(dualspace, basisContinuity);
    err = PetscDualSpaceSetUp(dualspace);PYLITH_CHECK_ERROR(err);

    // Create element
    PetscFE fe = NULL;
    err = PetscFECreate(PetscObjectComm((PetscObject) dm), &fe);PYLITH_CHECK_ERROR(err);
    err = PetscFESetType(fe, PETSCFEBASIC);PYLITH_CHECK_ERROR(err);
    err = PetscFESetBasisSpace(fe, space);PYLITH_CHECK_ERROR(err);
    err = PetscFESetDualSpace(fe, dualspace);PYLITH_CHECK_ERROR(err);
    err = PetscFESetNumComponents(fe, numComponents);PYLITH_CHECK_ERROR(err);
    err = PetscFESetUp(fe);PYLITH_CHECK_ERROR(err);
    err = PetscSpaceDestroy(&space);PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceDestroy(&dualspace);PYLITH_CHECK_ERROR(err);

    // Create quadrature
    PetscQuadrature quadrature = NULL;
    const int basisNumComponents = 1;
    const int numPoints = quadOrder + 1;
    const PylithReal xRefMin = -1.0;
    const PylithReal xRefMax = +1.0;
    if (isSimplex) {
        err = PetscDTGaussJacobiQuadrature(dim, basisNumComponents, numPoints, xRefMin, xRefMax, &quadrature);PYLITH_CHECK_ERROR(err);
    } else {
        err = PetscDTGaussTensorQuadrature(dim, basisNumComponents, numPoints, xRefMin, xRefMax, &quadrature);PYLITH_CHECK_ERROR(err);
    }
    err = PetscFESetQuadrature(fe, quadrature);PYLITH_CHECK_ERROR(err);
    err = PetscQuadratureDestroy(&quadrature);PYLITH_CHECK_ERROR(err);
    if (feinfo.feSpace == FieldBase::POLYNOMIAL_SPACE) {
        PetscQuadrature faceQuadrature = NULL;
        err = PetscDTGaussJacobiQuadrature(dim-1, basisNumComponents, numPoints, xRefMin, xRefMax, &faceQuadrature);PYLITH_CHECK_ERROR(err);
        err = PetscFESetFaceQuadrature(fe, faceQuadrature);PYLITH_CHECK_ERROR(err);
        err = PetscQuadratureDestroy(&faceQuadrature);PYLITH_CHECK_ERROR(err);
    } else {
        PetscQuadrature faceQuadrature = NULL;
        err = PetscDTGaussJacobiQuadrature(dim-1, basisNumComponents, 0, xRefMin, xRefMax, &faceQuadrature);PYLITH_CHECK_ERROR(err);
        err = PetscFESetFaceQuadrature(fe, faceQuadrature);PYLITH_CHECK_ERROR(err);
        err = PetscQuadratureDestroy(&faceQuadrature);PYLITH_CHECK_ERROR(err);
    }

    PYLITH_METHOD_RETURN(fe);
} // createFE


// ----------------------------------------------------------------------
// Check compatibility of discretization of subfields in the auxiliary field and target field.
void
pylith::topology::FieldOps::checkDiscretization(const pylith::topology::Field& target,
                                                const pylith::topology::Field& auxiliary) {
    PYLITH_METHOD_BEGIN;
    // PYLITH_JOURNAL_DEBUG("checkDiscretization(target="<<target.label()<<", auxiliary="<<auxiliary.label()<<")");

    // Verify that the quadrature order of the target subfields all
    // match and that they match the quadrature order of the auxiliary
    // subfields, because this is assumed by DMPlex integration
    // routines.

    // Get quadrature order in target subfields.
    PetscInt quadOrder = -1;
    { // target subfields
        const pylith::string_vector& subfieldNames = target.subfieldNames();
        const size_t numSubfields = subfieldNames.size();
        for (size_t i = 0; i < numSubfields; ++i) {
            const pylith::topology::Field::SubfieldInfo& sinfo = target.subfieldInfo(subfieldNames[i].c_str());
            if (quadOrder > 0) {
                if (quadOrder != sinfo.fe.quadOrder) {
                    std::ostringstream msg;
                    msg << "Quadrature order of subfields in target field '" << target.label()
                        << "' must all be the same. Expected quadrature order of " << quadOrder << ", but subfield '"
                        << subfieldNames[i] << "' has a quadrature order of " << sinfo.fe.quadOrder << ".";
                    throw std::runtime_error(msg.str());
                } // if
            } else {
                quadOrder = sinfo.fe.quadOrder;
            } // else
        } // for
    } // target subfields

    // Check quadrature order in auxiliary subfields.
    { // auxiliary subfields
        const pylith::string_vector& subfieldNames = auxiliary.subfieldNames();
        const size_t numSubfields = subfieldNames.size();
        for (size_t i = 0; i < numSubfields; ++i) {
            const pylith::topology::Field::SubfieldInfo& sinfo = auxiliary.subfieldInfo(subfieldNames[i].c_str());
            if (quadOrder > 0) {
                if (quadOrder != sinfo.fe.quadOrder) {
                    std::ostringstream msg;
                    msg << "Quadrature order of subfields in auxiliary field '" << auxiliary.label()
                        << "' must all match the quadrature order in the target subfields '" << target.label()
                        << "'. Expected quadrature order of " << quadOrder << ", but subfield '" << subfieldNames[i]
                        << "' has a quadrature order of " << sinfo.fe.quadOrder << ".";
                    throw std::runtime_error(msg.str());
                } // if
            } else {
                quadOrder = sinfo.fe.quadOrder;
            } // else
        } // for
    } // auxiliary subfields

    PYLITH_METHOD_END;
} // checkDiscretization


// ----------------------------------------------------------------------
// Check to see if fields have the same subfields and match in size.
bool
pylith::topology::FieldOps::layoutsMatch(const pylith::topology::Field& fieldA,
                                         const pylith::topology::Field& fieldB) {
    PYLITH_METHOD_BEGIN;
    // PYLITH_JOURNAL_DEBUG("layoutsMatch(fieldA="<<fieldA.label()<<", fieldB="<<fieldB.label()<<")");

    bool isMatch = true;

    // Check to see if fields have same chart and section sizes.
    if (fieldA.chartSize() != fieldB.chartSize()) { isMatch = false; }
    if (fieldA.sectionSize() != fieldB.sectionSize()) { isMatch = false; }

    // Check to see if number of subfields match.
    const pylith::string_vector& subfieldNamesA = fieldA.subfieldNames();
    const pylith::string_vector& subfieldNamesB = fieldB.subfieldNames();
    if (subfieldNamesA.size() != subfieldNamesB.size()) { isMatch = false; }

    // Check to see if subfields have same number of components and discretizations.
    const size_t numSubfields = subfieldNamesA.size();
    for (size_t i = 0; i < numSubfields; ++i) {
        const pylith::topology::Field::SubfieldInfo& infoA = fieldA.subfieldInfo(subfieldNamesA[i].c_str());
        const pylith::topology::Field::SubfieldInfo& infoB = fieldB.subfieldInfo(subfieldNamesB[i].c_str());

        if (infoA.description.numComponents != infoB.description.numComponents) { isMatch = false; }
        if (infoA.fe.basisOrder != infoB.fe.basisOrder) { isMatch = false; }
    } // for

    // Must match across all processors.
    PetscInt matchLocal = isMatch;
    PetscInt matchGlobal = 0;
    PetscErrorCode err = MPI_Allreduce(&matchLocal, &matchGlobal, 1, MPIU_INT, MPI_LOR, fieldA.mesh().comm());PYLITH_CHECK_ERROR(err);
    isMatch = matchGlobal == 1;

    // PYLITH_JOURNAL_DEBUG("layoutsMatch return value="<<isMatch<<".");
    PYLITH_METHOD_RETURN(isMatch);
} // layoutsMatch


// End of file
