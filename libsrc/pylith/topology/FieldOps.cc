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

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include "petscdm.h" // USES PetscDM
#include "petscds.h" // USES PetscDS


// ----------------------------------------------------------------------
// Create PetscFE object for discretization.
PetscFE
pylith::topology::FieldOps::createFE(const FieldBase::Discretization& feinfo,
                                     const PetscDM dm,
                                     const bool isSimplex,
                                     const int numComponents)
{ // createFE
    PYLITH_METHOD_BEGIN;

    const int basisOrder = PetscMax(feinfo.basisOrder, 0);
    const int quadOrder = PetscMax(feinfo.quadOrder > 0 ? feinfo.quadOrder : basisOrder, 1);
    const PetscBool basisContinuity = feinfo.isBasisContinuous ? PETSC_TRUE : PETSC_FALSE;
    const PetscBool useTensor = isSimplex ? PETSC_FALSE : PETSC_TRUE;

    PetscErrorCode err;

    // Get spatial dimension of mesh.
    int dim = 0;
    err = DMGetDimension(dm, &dim); PYLITH_CHECK_ERROR(err);

    // Create space
    PetscSpace space = NULL;
    err = PetscSpaceCreate(PetscObjectComm((PetscObject) dm), &space); PYLITH_CHECK_ERROR(err); assert(space);
    err = PetscSpaceSetType(space, feinfo.feSpace == FieldBase::POLYNOMIAL_SPACE ? PETSCSPACEPOLYNOMIAL : PETSCSPACEPOINT); PYLITH_CHECK_ERROR(err);
    err = PetscSpaceSetNumComponents(space, numComponents); PYLITH_CHECK_ERROR(err);
    err = PetscSpacePolynomialSetTensor(space, useTensor); PYLITH_CHECK_ERROR(err);
    err = PetscSpacePolynomialSetNumVariables(space, dim); PYLITH_CHECK_ERROR(err);
    err = PetscSpaceSetOrder(space, basisOrder);
    err = PetscSpaceSetUp(space); PYLITH_CHECK_ERROR(err);

    // Create dual space
    PetscDualSpace dualspace = NULL;
    PetscDM dmCell = NULL;
    err = PetscDualSpaceCreate(PetscObjectComm((PetscObject) dm), &dualspace); PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceCreateReferenceCell(dualspace, dim, isSimplex ? PETSC_TRUE : PETSC_FALSE, &dmCell); PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceSetDM(dualspace, dmCell); PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&dmCell); PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceSetNumComponents(dualspace, numComponents); PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceSetType(dualspace, PETSCDUALSPACELAGRANGE); PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceLagrangeSetTensor(dualspace, useTensor); PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceSetOrder(dualspace, basisOrder); PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceLagrangeSetContinuity(dualspace, basisContinuity);
    err = PetscDualSpaceSetUp(dualspace); PYLITH_CHECK_ERROR(err);

    // Create element
    PetscFE fe = NULL;
    err = PetscFECreate(PetscObjectComm((PetscObject) dm), &fe); PYLITH_CHECK_ERROR(err);
    err = PetscFESetType(fe, PETSCFEBASIC); PYLITH_CHECK_ERROR(err);
    err = PetscFESetBasisSpace(fe, space); PYLITH_CHECK_ERROR(err);
    err = PetscFESetDualSpace(fe, dualspace); PYLITH_CHECK_ERROR(err);
    err = PetscFESetNumComponents(fe, numComponents); PYLITH_CHECK_ERROR(err);
    err = PetscFESetUp(fe); PYLITH_CHECK_ERROR(err);
    err = PetscSpaceDestroy(&space); PYLITH_CHECK_ERROR(err);
    err = PetscDualSpaceDestroy(&dualspace); PYLITH_CHECK_ERROR(err);

    // Create quadrature
    PetscQuadrature quadrature = NULL;
    const int basisNumComponents = 1;
    const int numPoints = quadOrder + 1;
    const PylithReal xRefMin = -1.0;
    const PylithReal xRefMax = +1.0;
    if (isSimplex) {
        err = PetscDTGaussJacobiQuadrature(dim, basisNumComponents, numPoints, xRefMin, xRefMax, &quadrature); PYLITH_CHECK_ERROR(err);
    } else {
        err = PetscDTGaussTensorQuadrature(dim, basisNumComponents, numPoints, xRefMin, xRefMax, &quadrature); PYLITH_CHECK_ERROR(err);
    }
    err = PetscFESetQuadrature(fe, quadrature); PYLITH_CHECK_ERROR(err);
    err = PetscFESetFaceQuadrature(fe, quadrature); PYLITH_CHECK_ERROR(err);
    err = PetscQuadratureDestroy(&quadrature); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(fe);
} // createFE


// End of file
