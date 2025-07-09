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

#include "pylith/topology/FieldOps.hh" // implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/VisitorSubmesh.hh" // USES VecVisitorSubmesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include "petscdm.h" // USES PetscDM

std::map<pylith::topology::FieldBase::Discretization, pylith::topology::FE> pylith::topology::FieldOps::feStore = std::map<pylith::topology::FieldBase::Discretization, pylith::topology::FE>();

void
pylith::topology::FieldOps::deallocate(void) {
    pylith::topology::FieldOps::feStore.clear();
}


// ------------------------------------------------------------------------------------------------
// Create PetscFE object for discretization.
PetscFE
pylith::topology::FieldOps::createFE(const FieldBase::Discretization& feinfo,
                                     const PetscDM dm,
                                     const int numComponents) {
    PYLITH_METHOD_BEGIN;
    // Get spatial dimension of mesh.
    PetscFE fe;
    PetscInt dim;
    PetscErrorCode err;

    err = DMGetDimension(dm, &dim);PYLITH_CHECK_ERROR(err);
    dim = (feinfo.dimension < 0) ? dim : feinfo.dimension;assert(dim > 0);
    FieldBase::Discretization feKey = FieldBase::Discretization(feinfo.basisOrder, feinfo.quadOrder, dim, numComponents,
                                                                feinfo.isFaultOnly, feinfo.cellBasis, feinfo.feSpace,
                                                                feinfo.isBasisContinuous);
    std::map<FieldBase::Discretization, pylith::topology::FE>::const_iterator hasFE = pylith::topology::FieldOps::feStore.find(feKey);

    if (hasFE == pylith::topology::FieldOps::feStore.end()) {
        const int basisOrder = PetscMax(feKey.basisOrder, 0);
        const int quadOrder = PetscMax(feKey.quadOrder > 0 ? feKey.quadOrder : basisOrder, 0);
        const PetscBool simplexBasis = pylith::topology::FieldBase::SIMPLEX_BASIS == feKey.cellBasis ? PETSC_TRUE : PETSC_FALSE;
        const PetscBool useTensor = pylith::topology::FieldBase::TENSOR_BASIS == feKey.cellBasis ? PETSC_TRUE : PETSC_FALSE;
        const PetscBool basisContinuity = feKey.isBasisContinuous ? PETSC_TRUE : PETSC_FALSE;

        // Create space
        PetscSpace space = NULL;
        err = PetscSpaceCreate(PETSC_COMM_SELF, &space);PYLITH_CHECK_ERROR(err);assert(space);
        err = PetscSpaceSetType(space, feKey.feSpace == FieldBase::POLYNOMIAL_SPACE ?
                                PETSCSPACEPOLYNOMIAL : PETSCSPACEPOINT);PYLITH_CHECK_ERROR(err);
        err = PetscSpaceSetNumComponents(space, numComponents);PYLITH_CHECK_ERROR(err);
        err = PetscSpaceSetDegree(space, basisOrder, PETSC_DETERMINE);
        if (feKey.feSpace == FieldBase::POLYNOMIAL_SPACE) {
            err = PetscSpacePolynomialSetTensor(space, useTensor);PYLITH_CHECK_ERROR(err);
        } // if
        err = PetscSpaceSetNumVariables(space, dim);PYLITH_CHECK_ERROR(err);
        err = PetscSpaceSetUp(space);PYLITH_CHECK_ERROR(err);

        // Create dual space
        PetscDualSpace dualspace = NULL;
        PetscDM dmCell = NULL;
        err = PetscDualSpaceCreate(PETSC_COMM_SELF, &dualspace);PYLITH_CHECK_ERROR(err);
        err = DMPlexCreateReferenceCell(PETSC_COMM_SELF, DMPolytopeTypeSimpleShape(dim, simplexBasis), &dmCell);PYLITH_CHECK_ERROR(err);
        err = PetscDualSpaceSetDM(dualspace, dmCell);PYLITH_CHECK_ERROR(err);
        err = DMDestroy(&dmCell);PYLITH_CHECK_ERROR(err);
        err = PetscDualSpaceSetNumComponents(dualspace, numComponents);PYLITH_CHECK_ERROR(err);
        err = PetscDualSpaceSetType(dualspace, PETSCDUALSPACELAGRANGE);PYLITH_CHECK_ERROR(err);
        err = PetscDualSpaceLagrangeSetTensor(dualspace, useTensor);PYLITH_CHECK_ERROR(err);
        err = PetscDualSpaceSetOrder(dualspace, basisOrder);PYLITH_CHECK_ERROR(err);
        err = PetscDualSpaceLagrangeSetContinuity(dualspace, basisContinuity);
        err = PetscDualSpaceSetUp(dualspace);PYLITH_CHECK_ERROR(err);

        // Create element
        err = PetscFECreate(PETSC_COMM_SELF, &fe);PYLITH_CHECK_ERROR(err);
        err = PetscFESetType(fe, PETSCFEBASIC);PYLITH_CHECK_ERROR(err);
        err = PetscFESetBasisSpace(fe, space);PYLITH_CHECK_ERROR(err);
        err = PetscFESetDualSpace(fe, dualspace);PYLITH_CHECK_ERROR(err);
        err = PetscFESetNumComponents(fe, numComponents);PYLITH_CHECK_ERROR(err);
        err = PetscFESetUp(fe);PYLITH_CHECK_ERROR(err);
        err = PetscSpaceDestroy(&space);PYLITH_CHECK_ERROR(err);
        err = PetscDualSpaceDestroy(&dualspace);PYLITH_CHECK_ERROR(err);

        // Create quadrature
        PetscQuadrature quadrature = NULL;
        PetscQuadrature faceQuadrature = NULL;
        DMPolytopeType ct;
        switch (dim) {
        case 0: ct = DM_POLYTOPE_POINT;break;
        case 1: ct = DM_POLYTOPE_SEGMENT;break;
        case 2: ct = useTensor ? DM_POLYTOPE_QUADRILATERAL : DM_POLYTOPE_TRIANGLE;break;
        case 3: ct = useTensor ? DM_POLYTOPE_HEXAHEDRON : DM_POLYTOPE_TETRAHEDRON;break;
        default: throw std::logic_error("Cannot handle dimension");
        }
        err = PetscDTCreateDefaultQuadrature(ct, quadOrder, &quadrature, &faceQuadrature);PYLITH_CHECK_ERROR(err);
        err = PetscFESetQuadrature(fe, quadrature);PYLITH_CHECK_ERROR(err);
        err = PetscQuadratureDestroy(&quadrature);PYLITH_CHECK_ERROR(err);
        err = PetscFESetFaceQuadrature(fe, faceQuadrature);PYLITH_CHECK_ERROR(err);
        err = PetscQuadratureDestroy(&faceQuadrature);PYLITH_CHECK_ERROR(err);

        assert(feKey.feSpace == FieldBase::POLYNOMIAL_SPACE);
        pylith::topology::FieldOps::feStore.insert(std::pair<FieldBase::Discretization, pylith::topology::FE>(feKey, fe));
    } else {
        throw std::logic_error("FieldOps::createFE() :TODO: Can't reuse PetscFE due to naming of fields, so make a deep copy of fe.");
        fe = hasFE->second._fe;
        err = PetscObjectReference((PetscObject) fe);PYLITH_CHECK_ERROR(err);
    }

    PYLITH_METHOD_RETURN(fe);
} // createFE


// ------------------------------------------------------------------------------------------------
// Check compatibility of discretization of subfields in the auxiliary field and target field.
void
pylith::topology::FieldOps::checkDiscretization(const pylith::topology::Field& target,
                                                const pylith::topology::Field& auxiliary) {
    PYLITH_METHOD_BEGIN;
    // PYLITH_JOURNAL_DEBUG("checkDiscretization(target="<<target.getLabel()<<",
    // auxiliary="<<auxiliary.getLabel()<<")");

    // Verify that the quadrature order of the target subfields all
    // match and that they match the quadrature order of the auxiliary
    // subfields, because this is assumed by DMPlex integration
    // routines.

    // Get quadrature order in target subfields.
    PetscInt quadOrder = -1;
    { // target subfields
        const pylith::string_vector& subfieldNames = target.getSubfieldNames();
        const size_t numSubfields = subfieldNames.size();
        for (size_t i = 0; i < numSubfields; ++i) {
            const pylith::topology::Field::SubfieldInfo& sinfo = target.getSubfieldInfo(subfieldNames[i].c_str());
            if (quadOrder > 0) {
                if (quadOrder != sinfo.fe.quadOrder) {
                    std::ostringstream msg;
                    msg << "Quadrature order of subfields in target field '" << target.getLabel()
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
        const pylith::string_vector& subfieldNames = auxiliary.getSubfieldNames();
        const size_t numSubfields = subfieldNames.size();
        for (size_t i = 0; i < numSubfields; ++i) {
            const pylith::topology::Field::SubfieldInfo& sinfo = auxiliary.getSubfieldInfo(subfieldNames[i].c_str());
            if (quadOrder > 0) {
                if (quadOrder != sinfo.fe.quadOrder) {
                    std::ostringstream msg;
                    msg << "Quadrature order of subfields in auxiliary field '" << auxiliary.getLabel()
                        << "' must all match the quadrature order in the target subfields '" << target.getLabel()
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


// ------------------------------------------------------------------------------------------------
// Check that 'field' contains required fields for 'reason'.
void
pylith::topology::FieldOps::checkSubfieldsExist(const pylith::string_vector& requiredFields,
                                                const std::string& reason,
                                                const pylith::topology::Field& field) {
    pylith::string_vector missingFields;
    size_t numMissing = 0;
    const size_t numRequired = requiredFields.size();
    for (size_t i = 0; i < numRequired; ++i) {
        if (!field.hasSubfield(requiredFields[i].c_str())) {
            missingFields.resize(numMissing+1);
            missingFields[numMissing++] = requiredFields[i];
        } // if
    } // for

    if (numMissing) {
        const std::vector<std::string>& subfieldNames = field.getSubfieldNames();
        std::ostringstream msg;
        msg << "Could not find";
        for (size_t i = 0; i < numMissing; ++i) {
            msg << " '" << missingFields[i] << "'";
            if (i+1 < numMissing) {
                msg << ",";
            }
        }
        msg << " in " << field.getLabel() << " field. Field contains: ";
        for (size_t i = 0; i < subfieldNames.size(); ++i) {
            msg << "'" << subfieldNames[i] << "'";
            if (i+1 < subfieldNames.size()) {
                msg << ",";
            }
        } // for
        msg << "; the missing fields are required for " << reason;
        throw std::runtime_error(msg.str());
    } // if

} // checkSubfieldsExist


// ------------------------------------------------------------------------------------------------
// Get names of subfields extending over entire domain.
pylith::string_vector
pylith::topology::FieldOps::getSubfieldNamesDomain(const pylith::topology::Field& field) {
    PYLITH_METHOD_BEGIN;
    const pylith::string_vector& subfieldNames = field.getSubfieldNames();

    // Restrict fields to those defined over the entire domain
    // (as opposed to those defined over a subset like the fault_lagrange_multiplier).
    PetscDS fieldDS = NULL;
    PetscErrorCode err = DMGetDS(field.getDM(), &fieldDS);PYLITH_CHECK_ERROR(err);
    PylithInt numFields = 0;
    err = PetscDSGetNumFields(fieldDS, &numFields);PYLITH_CHECK_ERROR(err);
    assert(numFields > 0);
    pylith::string_vector subfieldNamesDomain(numFields);
    for (PylithInt iField = 0; iField < numFields; ++iField) {
        PetscObject discretization = NULL;
        err = PetscDSGetDiscretization(fieldDS, iField, &discretization);PYLITH_CHECK_ERROR(err);
        PylithInt fieldIndex = -1;
        err = PetscDSGetFieldIndex(fieldDS, discretization, &fieldIndex);PYLITH_CHECK_ERROR(err);
        assert(fieldIndex >= 0 && size_t(fieldIndex) < subfieldNames.size());
        subfieldNamesDomain[iField] = subfieldNames[fieldIndex];
    } // for

    PYLITH_METHOD_RETURN(subfieldNamesDomain);
} // getSubfieldNamesDomain


// ------------------------------------------------------------------------------------------------
// Check to see if fields have the same subfields and match in size.
bool
pylith::topology::FieldOps::layoutsMatch(const pylith::topology::Field& fieldA,
                                         const pylith::topology::Field& fieldB) {
    PYLITH_METHOD_BEGIN;
    // PYLITH_JOURNAL_DEBUG("layoutsMatch(fieldA="<<fieldA.getLabel()<<", fieldB="<<fieldB.getLabel()<<")");

    bool isMatch = true;

    // Check to see if fields have same chart and section sizes.
    if (fieldA.getChartSize() != fieldB.getChartSize()) { isMatch = false; }
    if (fieldA.getStorageSize() != fieldB.getStorageSize()) { isMatch = false; }

    // Check to see if number of subfields match.
    const pylith::string_vector& subfieldNamesA = fieldA.getSubfieldNames();
    const pylith::string_vector& subfieldNamesB = fieldB.getSubfieldNames();
    if (subfieldNamesA.size() != subfieldNamesB.size()) { isMatch = false; }

    // Check to see if subfields have same number of components and discretizations.
    const size_t numSubfields = subfieldNamesA.size();
    if (isMatch) {
        for (size_t i = 0; i < numSubfields; ++i) {
            const pylith::topology::Field::SubfieldInfo& infoA = fieldA.getSubfieldInfo(subfieldNamesA[i].c_str());
            const pylith::topology::Field::SubfieldInfo& infoB = fieldB.getSubfieldInfo(subfieldNamesB[i].c_str());

            if (infoA.description.numComponents != infoB.description.numComponents) { isMatch = false; }
            if (infoA.fe.cellBasis != infoB.fe.cellBasis) { isMatch = false; }
            if (infoA.fe.basisOrder != infoB.fe.basisOrder) { isMatch = false; }
        } // for
    } // if

    // Must match across all processors.
    PetscInt matchLocal = isMatch;
    PetscInt matchGlobal = 0;
    PetscErrorCode err = MPI_Allreduce(&matchLocal, &matchGlobal, 1, MPIU_INT, MPI_LOR, fieldA.getMesh().getComm());PYLITH_CHECK_ERROR(err);
    isMatch = matchGlobal == 1;

    // PYLITH_JOURNAL_DEBUG("layoutsMatch return value="<<isMatch<<".");
    PYLITH_METHOD_RETURN(isMatch);
} // layoutsMatch


// ------------------------------------------------------------------------------------------------
// Create label for output.
void
pylith::topology::FieldOps::createOutputLabel(const pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;

    const char* outputLabelName = "output";
    PetscDM fieldDM = field->getDM();
    PetscDMLabel outputLabel = NULL;
    PetscErrorCode err = PETSC_SUCCESS;
    err = DMCreateLabel(fieldDM, outputLabelName);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(fieldDM, outputLabelName, &outputLabel);PYLITH_CHECK_ERROR(err);
    pylith::topology::Stratum faultStratum(fieldDM, pylith::topology::Stratum::HEIGHT, 1);
    for (PetscInt point = faultStratum.begin(); point != faultStratum.end(); ++point) {
        err = DMLabelSetValue(outputLabel, point, 1);
    } // for
    err = DMPlexLabelComplete(fieldDM, outputLabel);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
void
pylith::topology::FieldOps::transformVector(PetscVec* outputVector,
                                            const PetscDM& outputDM,
                                            const PetscVec& inputVector,
                                            const PetscDM& inputDM) {
    PYLITH_METHOD_BEGIN;

    assert(outputVector);
    assert(outputDM);
    assert(inputVector);
    assert(inputDM);

    PetscErrorCode err = PETSC_SUCCESS;
    PetscSection inputSection = PETSC_NULLPTR, outputSection = PETSC_NULLPTR;
    err = DMGetGlobalSection(inputDM, &inputSection);PYLITH_CHECK_ERROR(err);
    err = DMGetGlobalSection(outputDM, &outputSection);PYLITH_CHECK_ERROR(err);

    // Verify sizes
    PetscInt outputNumFields = 0, inputNumFields = 0;
    err = PetscSectionGetNumFields(inputSection, &inputNumFields);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetNumFields(outputSection, &outputNumFields);PYLITH_CHECK_ERROR(err);
    assert(inputNumFields == outputNumFields);

    PetscInt inputSize = 0, outputSize = 0;
    err = PetscSectionGetStorageSize(inputSection, &inputSize);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetStorageSize(outputSection, &outputSize);PYLITH_CHECK_ERROR(err);
    assert(inputSize == outputSize);

    PetscInt inputVecSize = 0, outputVecSize = 0;
    err = VecGetLocalSize(inputVector, &inputVecSize);PYLITH_CHECK_ERROR(err);
    err = VecGetLocalSize(*outputVector, &outputVecSize);PYLITH_CHECK_ERROR(err);
    assert(inputVecSize == inputSize);
    assert(outputVecSize == outputSize);

    // Copy values from input vector to output vector
    const PetscScalar* inputArray = PETSC_NULLPTR;
    PetscScalar* outputArray = PETSC_NULLPTR;
    err = VecGetArrayRead(inputVector, &inputArray);PYLITH_CHECK_ERROR(err);
    err = VecGetArray(*outputVector, &outputArray);PYLITH_CHECK_ERROR(err);

    const PetscInt inputSectionLocalOffset = getGlobalSectionLocalOffset(inputSection);
    const PetscInt outputSectionLocalOffset = getGlobalSectionLocalOffset(outputSection);

    // Check whether inputDM and outputDM have the same names for points.
    // We will assume that if inputDM and outputDM have the same range of points,
    // then the names are the same and we ignore the subpoint map.
    PetscDMLabel subpointMap = PETSC_NULLPTR; // Mapping of points in output DM back to input DM
    PetscInt pStartIn, pEndIn, pStartOut, pEndOut;
    PetscBool renamePoints = PETSC_FALSE;
    err = DMPlexGetSubpointMap(outputDM, &subpointMap);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetChart(inputDM, &pStartIn, &pEndIn);
    err = DMPlexGetChart(outputDM, &pStartOut, &pEndOut);
    renamePoints = subpointMap && ((pStartIn != pStartOut) || (pEndIn != pEndOut)) ? PETSC_TRUE : PETSC_FALSE;
    if (renamePoints) {
        PetscIS subpointIS = PETSC_NULLPTR; // Mapping of points in output DM back to input DM
        PetscInt subpointISSize = 0;
        const PetscInt* subpointISPoints = PETSC_NULLPTR;
        err = DMPlexGetSubpointIS(outputDM, &subpointIS);PYLITH_CHECK_ERROR(err);
        err = ISGetLocalSize(subpointIS, &subpointISSize);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(subpointIS, &subpointISPoints);PYLITH_CHECK_ERROR(err);

        PetscInt pStart = 0, pEnd = 0;
        err = PetscSectionGetChart(outputSection, &pStart, &pEnd);
        for (PetscInt iPoint = 0; iPoint < subpointISSize; ++iPoint) {
            PetscInt inputDof = 0, inputOffset = 0;
            PetscInt outputDof = 0, outputOffset = 0;
            const PetscInt inputPoint = subpointISPoints[iPoint];
            const PetscInt outputPoint = pStart + iPoint;
            err = PetscSectionGetDof(inputSection, inputPoint, &inputDof);PYLITH_CHECK_ERROR(err);
            err = PetscSectionGetOffset(inputSection, inputPoint, &inputOffset);PYLITH_CHECK_ERROR(err);
            inputOffset -= inputSectionLocalOffset;
            err = PetscSectionGetDof(outputSection, outputPoint, &outputDof);PYLITH_CHECK_ERROR(err);
            err = PetscSectionGetOffset(outputSection, outputPoint, &outputOffset);PYLITH_CHECK_ERROR(err);
            outputOffset -= outputSectionLocalOffset;

            assert(inputDof == outputDof);
            for (PetscInt iDof = 0; iDof < inputDof; ++iDof) {
                outputArray[outputOffset + iDof] = inputArray[inputOffset + iDof];
            } // for
        } // for
    } else {
        PetscInt pStart = 0, pEnd = 0;
        err = PetscSectionGetChart(outputSection, &pStart, &pEnd);
        for (PetscInt point = pStart; point < pEnd; ++point) {
            PetscInt inputDof = 0, inputOffset = 0;
            PetscInt outputDof = 0, outputOffset = 0;
            const PetscInt inputPoint = point;
            const PetscInt outputPoint = point;
            err = PetscSectionGetDof(inputSection, inputPoint, &inputDof);PYLITH_CHECK_ERROR(err);
            err = PetscSectionGetOffset(inputSection, inputPoint, &inputOffset);PYLITH_CHECK_ERROR(err);
            inputOffset -= inputSectionLocalOffset;
            err = PetscSectionGetDof(outputSection, outputPoint, &outputDof);PYLITH_CHECK_ERROR(err);
            err = PetscSectionGetOffset(outputSection, outputPoint, &outputOffset);PYLITH_CHECK_ERROR(err);
            outputOffset -= outputSectionLocalOffset;

            assert(inputDof == outputDof);
            for (PetscInt iDof = 0; iDof < inputDof; ++iDof) {
                outputArray[outputOffset + iDof] = inputArray[inputOffset + iDof];
            } // for
        } // for
    } // if/else

    err = VecRestoreArrayRead(inputVector, &inputArray);PYLITH_CHECK_ERROR(err);
    err = VecRestoreArray(*outputVector, &outputArray);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
PetscInt
pylith::topology::FieldOps::getGlobalSectionLocalOffset(PetscSection globalSection) {
    PetscInt pStart = 0, pEnd = 0;
    PetscInt localOffset = 0;
    PetscErrorCode err = PETSC_SUCCESS;

    err = PetscSectionGetChart(globalSection, &pStart, &pEnd);
    if (pStart < pEnd) {
        for (PetscInt point = pStart; point < pEnd; ++point) {
            err = PetscSectionGetOffset(globalSection, point, &localOffset);PYLITH_CHECK_ERROR(err);
            if (localOffset >= 0) {
                break;
            } // if
        } // for
    } // if

    return localOffset;
}


// End of file
