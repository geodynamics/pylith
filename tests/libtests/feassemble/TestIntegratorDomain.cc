// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestIntegratorDomain.hh" // Implementation of class methods

#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/testing/FieldTester.hh" // USES FieldTester
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // pythia::journal

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestIntegratorDomain::setUp(void) {
    _integrator = NULL;
    _data = NULL;

    _mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    _solutionFields = NULL;
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate testing data.
void
pylith::feassemble::TestIntegratorDomain::tearDown(void) {
    delete _integrator;_integrator = NULL;
    delete _data;_data = NULL;

    delete _solutionFields;_solutionFields = NULL;
    delete _mesh;_mesh = NULL;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Test getPhysicsDomainMesh(), getAuxiliaryField(), getDerivedField(), getMaterialId(), setMaterialId().
void
pylith::feassemble::TestIntegratorDomain::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    _initializeMin();

    CPPUNIT_ASSERT(_integrator);
    const pylith::topology::Mesh& meshIntegrator = _integrator->getPhysicsDomainMesh();
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT_EQUAL(_mesh->dimension(), meshIntegrator.dimension());
    CPPUNIT_ASSERT_EQUAL(_mesh->numCorners(), meshIntegrator.numCorners());
    CPPUNIT_ASSERT_EQUAL(_mesh->numCells(), meshIntegrator.numCells());
    CPPUNIT_ASSERT_EQUAL(_mesh->numVertices(), meshIntegrator.numVertices());

    const pylith::topology::Field* auxiliaryField = _integrator->getAuxiliaryField();
    CPPUNIT_ASSERT(auxiliaryField);
    CPPUNIT_ASSERT_EQUAL(std::string("auxiliary subfields"), std::string(auxiliaryField->getLabel()));

    const pylith::topology::Field* derivedField = _integrator->getDerivedField();
    CPPUNIT_ASSERT(derivedField);
    CPPUNIT_ASSERT_EQUAL(std::string("derived subfields"), std::string(derivedField->getLabel()));

    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of default label name failed.",
                                 std::string("material-id"), std::string(_integrator->getLabelName()));
    const std::string& labelName = "material-label";
    _integrator->setLabelName(labelName.c_str());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of custom label name failed.",
                                 labelName, std::string(_integrator->getLabelName()));

    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of default label value.", 1, _integrator->getLabelValue());
    _integrator->setLabelValue(_data->materialId);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of custom label value.", _data->materialId, _integrator->getLabelValue());

    PYLITH_METHOD_END;
} // testAccessors


// ----------------------------------------------------------------------
// Test setKernelsRHSResidual(), setKernelsLHSResidual(), setKernelsLHSJacobian().
void
pylith::feassemble::TestIntegratorDomain::testSetKernels(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad Implement testSetKernels().", false);

    PYLITH_METHOD_END;
} // testSetKernels


// ---------------------------------------------------------------------------------------------------------------------
// Test initialize().
void
pylith::feassemble::TestIntegratorDomain::testInitialize(void) {
    PYLITH_METHOD_BEGIN;

    _initializeMin();

    CPPUNIT_ASSERT(_integrator);
    const pylith::topology::Field* auxiliaryField = _integrator->getAuxiliaryField();

    CPPUNIT_ASSERT(_data);

    CPPUNIT_ASSERT_EQUAL(std::string("auxiliary subfields"), std::string(auxiliaryField->getLabel()));
    CPPUNIT_ASSERT_EQUAL(_data->dimension, auxiliaryField->getSpaceDim());

    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT(_data->normalizer);
    const PylithReal lengthScale = _data->normalizer->getLengthScale();
    PylithReal norm = pylith::testing::FieldTester::checkFieldWithDB(*auxiliaryField, _data->auxiliaryDB, lengthScale);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Values in auxiliary field do not match spatial database.", 0.0, norm, tolerance);

    // Verify solution and perturbation fields can be exactly represented by discretization.
    CPPUNIT_ASSERT(_solutionFields);
    pylith::topology::Field& solution = _solutionFields->get("solution");
    norm = pylith::testing::FieldTester::checkFieldWithDB(solution, _data->solutionDB, lengthScale);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Solution field failed representation test.", 0.0, norm, tolerance);

    pylith::topology::Field& perturbation = _solutionFields->get("perturbation");
    norm = pylith::testing::FieldTester::checkFieldWithDB(perturbation, _data->perturbationDB, lengthScale);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Perturbation field failed representation test.", 0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testInitialize


// ---------------------------------------------------------------------------------------------------------------------
// Test poststep().
void
pylith::feassemble::TestIntegratorDomain::testPoststep(void) {
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    // We test updating the state variables in the auxiliary field by
    // passing in the perturbation as the "new" solution and the existing
    // auxiliary field. We test whether the "updated" auxiliary field
    // matches the database with the updated auxiliary field.

    CPPUNIT_ASSERT(_solutionFields);
    pylith::topology::Field& perturbation = _solutionFields->get("perturbation");
    CPPUNIT_ASSERT(_integrator);
    CPPUNIT_ASSERT(_data);
    _integrator->poststep(_data->t, _data->tindex, _data->dt, perturbation);

    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT(_data->normalizer);
    const PylithReal lengthScale = _data->normalizer->getLengthScale();
    PylithReal norm = pylith::testing::FieldTester::checkFieldWithDB(*_integrator->getAuxiliaryField(),
                                                                     _data->auxiliaryUpdateDB, lengthScale);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Updated auxiliary field values do not match spatial database.", 0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testPoststep


// ---------------------------------------------------------------------------------------------------------------------
// Test updateState().
void
pylith::feassemble::TestIntegratorDomain::testUpdateState(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_integrator);
    CPPUNIT_ASSERT(_data);
    _integrator->updateState(_data->t);

    // Nothing to check. updateState() should not do anything.

    PYLITH_METHOD_END;
} // testUpdateState


// ---------------------------------------------------------------------------------------------------------------------
// Test computeRHSResidual(), computeLHSResidual().
void
pylith::feassemble::TestIntegratorDomain::testComputeResidual(void) {
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_solutionFields);
    pylith::topology::Field& solution = _solutionFields->get("solution");
    pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");

    pylith::topology::Field residualRHS(solution);
    residualRHS.setLabel("residual RHS");
    residualRHS.createDiscretization();
    residualRHS.allocate();

    pylith::topology::Field residualLHS(solution);
    residualLHS.setLabel("residual LHS");
    residualLHS.createDiscretization();
    residualLHS.allocate();

#if 0 // :DEBUG:
    PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2"); // :DEBUG:
    DMSetFromOptions(residualRHS.dmMesh()); // :DEBUG:
#endif // :DEBUG:

    CPPUNIT_ASSERT(_data);
    const PylithReal t = _data->t;
    const PylithReal dt = _data->dt;
    CPPUNIT_ASSERT(_integrator);
    _integrator->computeRHSResidual(&residualRHS, t, dt, solution);
    _integrator->computeLHSResidual(&residualLHS, t, dt, solution, solutionDot);

    // We don't use Dirichlet BC, so we must manually zero out the residual values on the boundary.
    _zeroBoundary(&residualRHS);
    _zeroBoundary(&residualLHS);

#if 0 // :DEBUG:
    solution.view("SOLUTION"); // :DEBUG:
    solutionDot.view("SOLUTION_DOT"); // :DEBUG:
    residualRHS.view("RESIDUAL RHS"); // :DEBUG:
    residualLHS.view("RESIDUAL LHS"); // :DEBUG:
#endif // :DEBUG:

    PetscErrorCode err;
    PetscVec residualVec = NULL;
    err = VecDuplicate(residualRHS.localVector(), &residualVec);CPPUNIT_ASSERT(!err);
    err = VecWAXPY(residualVec, -1.0, residualRHS.localVector(), residualLHS.localVector());CPPUNIT_ASSERT(!err);

    PylithReal norm = 0.0;
    PylithReal normRHS = 0.0;
    PylithReal normLHS = 0.0;
    err = VecNorm(residualRHS.localVector(), NORM_2, &normRHS);CPPUNIT_ASSERT(!err);
    err = VecNorm(residualLHS.localVector(), NORM_2, &normLHS);CPPUNIT_ASSERT(!err);
    err = VecNorm(residualVec, NORM_2, &norm);CPPUNIT_ASSERT(!err);
    err = VecDestroy(&residualVec);CPPUNIT_ASSERT(!err);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test of F(s) - G(s) == 0 failed.", 0.0, norm, tolerance);
    CPPUNIT_ASSERT_MESSAGE("RHS and LHS residuals are both exactly zero, which is suspicious.", normRHS > 0.0 || normLHS > 0.0);

    PYLITH_METHOD_END;
} // testComputeResidual


// ---------------------------------------------------------------------------------------------------------------------
// Test computeLHSJacobian().
void
pylith::feassemble::TestIntegratorDomain::testComputeJacobian(void) {
    PYLITH_METHOD_BEGIN;

    // Create linear problem (MMS) with two trial solutions, s,s_dot and p,p_dot.
    //
    // Check that Jf(s,s_dot)*(p - s) = F(p,p_dot) - F(s,s_dot).

    // Call initialize()
    _initializeFull(); // includes setting up auxField

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_solutionFields);
    pylith::topology::Field& solution = _solutionFields->get("solution");
    pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");
    pylith::topology::Field& perturbation = _solutionFields->get("perturbation");
    pylith::topology::Field& perturbationDot = _solutionFields->get("perturbation_dot");

    pylith::topology::Field residual1(solution);
    residual1.setLabel("residual1");
    residual1.createDiscretization();
    residual1.allocate();

    pylith::topology::Field residual2(perturbation);
    residual2.setLabel("residual2");
    residual2.createDiscretization();
    residual2.allocate();

#if 0 // :DEBUG:
    PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2"); // :DEBUG:
    DMSetFromOptions(_solution1->dmMesh()); // :DEBUG:
#endif // :DEBUG:

    CPPUNIT_ASSERT(_data);
    const PylithReal t = _data->t;
    const PylithReal dt = _data->dt;
    const PylithReal s_tshift = _data->s_tshift;
    CPPUNIT_ASSERT(_integrator);
    _integrator->computeLHSResidual(&residual1, t, dt, solution, solutionDot);
    _integrator->computeLHSResidual(&residual2, t, dt, perturbation, perturbationDot);

    // residual1.view("RESIDUAL 1 LHS"); // :DEBUG:
    // residual2.view("RESIDUAL 2 LHS"); // :DEBUG:

    PetscErrorCode err;

    PetscVec residualVec = NULL;
    err = VecDuplicate(residual1.localVector(), &residualVec);CPPUNIT_ASSERT(!err);
    err = VecWAXPY(residualVec, -1.0, residual1.localVector(), residual2.localVector());CPPUNIT_ASSERT(!err);

    PetscVec solutionIncrVec = NULL;
    err = VecDuplicate(solution.localVector(), &solutionIncrVec);CPPUNIT_ASSERT(!err);
    err = VecWAXPY(solutionIncrVec, -1.0, solution.localVector(), perturbation.localVector());CPPUNIT_ASSERT(!err);

    // Compute Jacobian
    PetscMat jacobianMat = NULL;
    err = DMCreateMatrix(solution.dmMesh(), &jacobianMat);CPPUNIT_ASSERT(!err);
    err = MatZeroEntries(jacobianMat);CPPUNIT_ASSERT(!err);
    PetscMat precondMat = jacobianMat; // Use Jacobian == preconditioner

    _integrator->computeLHSJacobian(jacobianMat, precondMat, t, dt, s_tshift, solution, solutionDot);
    CPPUNIT_ASSERT_EQUAL(false, _integrator->needNewLHSJacobian(false));
    err = MatAssemblyBegin(jacobianMat, MAT_FINAL_ASSEMBLY);PYLITH_CHECK_ERROR(err);
    err = MatAssemblyEnd(jacobianMat, MAT_FINAL_ASSEMBLY);PYLITH_CHECK_ERROR(err);

    // result = J*(-solutionIncr) + residual
    PetscVec resultVec = NULL;
    err = VecDuplicate(residualVec, &resultVec);CPPUNIT_ASSERT(!err);
    err = VecZeroEntries(resultVec);CPPUNIT_ASSERT(!err);
    err = VecScale(solutionIncrVec, -1.0);CPPUNIT_ASSERT(!err);
    err = MatMultAdd(jacobianMat, solutionIncrVec, residualVec, resultVec);CPPUNIT_ASSERT(!err);

#if 0 // :DEBUG:
    std::cout << "Solution INCR" << std::endl;
    VecView(solutionIncrVec, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "F2-F1" << std::endl;
    VecView(residualVec, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "RESULT" << std::endl;
    VecView(resultVec, PETSC_VIEWER_STDOUT_SELF);
#endif // :DEBUG:

    PylithReal norm = 0.0, normSolutionIncr = 0.0, normResidual = 0.0;
    err = VecNorm(resultVec, NORM_2, &norm);CPPUNIT_ASSERT(!err);
    err = VecNorm(solutionIncrVec, NORM_2, &normSolutionIncr);CPPUNIT_ASSERT(!err);
    err = VecNorm(residualVec, NORM_2, &normResidual);CPPUNIT_ASSERT(!err);
    err = VecDestroy(&resultVec);CPPUNIT_ASSERT(!err);
    err = VecDestroy(&solutionIncrVec);CPPUNIT_ASSERT(!err);
    err = VecDestroy(&residualVec);CPPUNIT_ASSERT(!err);
    err = MatDestroy(&jacobianMat);CPPUNIT_ASSERT(!err);

    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Check of Jf(s)*(p-s) - (F(p) - F(s)) == 0 failed.", 0.0, norm, tolerance);
    CPPUNIT_ASSERT_MESSAGE("Norm of resulting vector is exactly zero, which is suspicious.", (0 < normResidual && 0 < norm) || (0 == normResidual && 0 == norm));

    PYLITH_METHOD_END;
} // testComputeJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Do minimal initilaization of test data.
void
pylith::feassemble::TestIntegratorDomain::_initializeMin(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_data);

    pylith::meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(_data->meshFilename);
    iohandler.filename(_data->meshFilename);
    iohandler.read(_mesh);CPPUNIT_ASSERT(_mesh);

    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any cells.", _mesh->numCells() > 0);
    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any vertices.", _mesh->numVertices() > 0);

    // Setup coordinates.
    _mesh->setCoordSys(_data->cs);
    CPPUNIT_ASSERT(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    // Setup solution fields.
    delete _solutionFields;_solutionFields = new pylith::topology::Fields(*_mesh);CPPUNIT_ASSERT(_solutionFields);
    _solutionFields->add("solution","solution");
    _solutionFields->add("solution_dot","solution_dot");
    _solutionFields->add("perturbation","perturbation");
    _solutionFields->add("perturbation_dot","perturbation_dot");
    _setupSolutionFields();

    PYLITH_METHOD_END;
} // _initializeMin


// ---------------------------------------------------------------------------------------------------------------------
// Complete initilaization of test data.
void
pylith::feassemble::TestIntegratorDomain::_initializeFull(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);

    CPPUNIT_ASSERT(_solutionFields);
    CPPUNIT_ASSERT(_integrator);
    _integrator->initialize(_solutionFields->get("solution"));

    PYLITH_METHOD_END;
} // _initializeFull


#if 0
// ---------------------------------------------------------------------------------------------------------------------
// Setup and populate solution fields.
void
pylith::feassemble::TestIntegratorDomain::_setupSolutionFields(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_solutionFields);
    CPPUNIT_ASSERT(_data->solutionDiscretizations);
    CPPUNIT_ASSERT(_data->normalizer);

    { // Solution
        pylith::topology::Field& solution = _solutionFields->get("solution");
        pylith::problems::SolutionFactory factory(solution, *_data->normalizer);
        factory.displacement(_data->solutionDiscretizations[0]);
        if (_data->isExplicit) {
            factory.velocity(_data->solutionDiscretizations[1]);
        } // if
        solution.subfieldsSetup();
        solution.createDiscretization();
        solution.allocate();
        factory.setValues(_data->solutionDB);
    } // Solution

    { // Time derivative of solution
        pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");
        pylith::problems::SolutionFactory factory(solutionDot, *_data->normalizer);
        factory.displacementDot(_data->solutionDiscretizations[0]);
        if (_data->isExplicit) {
            factory.velocityDot(_data->solutionDiscretizations[1]);
        } // if
        solutionDot.subfieldsSetup();
        solutionDot.createDiscretization();
        solutionDot.allocate();
        factory.setValues(_data->solutionDB);
    } // Time derivative of solution

    { // Perturbation
        pylith::topology::Field& perturbation = _solutionFields->get("perturbation");
        const pylith::topology::Field& solution = _solutionFields->get("solution");
        perturbation.cloneSection(solution);
        perturbation.createDiscretization();
        perturbation.allocate();
        perturbation.zeroLocal();
        pylith::problems::SolutionFactory factory(perturbation, *_data->normalizer);
        factory.setValues(_data->perturbationDB);
    } // Perturbation

    { // Time derivative perturbation
        pylith::topology::Field& perturbationDot = _solutionFields->get("perturbation_dot");
        const pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");
        perturbationDot.cloneSection(solutionDot);
        perturbationDot.createDiscretization();
        perturbationDot.allocate();
        perturbationDot.zeroLocal();
        pylith::problems::SolutionFactory factory(perturbationDot, *_data->normalizer);
        factory.setValues(_data->perturbationDB);
    } // Time derivative perturbation

    PYLITH_METHOD_END;
} // _setupSolutionFields


#endif

// ---------------------------------------------------------------------------------------------------------------------
// Set field to zero on the boundary.
void
pylith::feassemble::TestIntegratorDomain::_zeroBoundary(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(field);
    CPPUNIT_ASSERT(_data->boundaryLabel);

    PetscDM dmMesh = field->mesh().dmMesh();CPPUNIT_ASSERT(dmMesh);
    PetscDMLabel label = NULL;
    PetscIS pointIS = NULL;
    const PetscInt *points;
    PetscInt numPoints = 0;
    PetscBool hasLabel = PETSC_FALSE;
    PetscErrorCode err;
    err = DMHasLabel(dmMesh, _data->boundaryLabel, &hasLabel);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(hasLabel);
    err = DMGetLabel(dmMesh, _data->boundaryLabel, &label);CPPUNIT_ASSERT(!err);
    err = DMLabelGetStratumIS(label, 1, &pointIS);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(pointIS);
    err = ISGetLocalSize(pointIS, &numPoints);CPPUNIT_ASSERT(!err);
    err = ISGetIndices(pointIS, &points);CPPUNIT_ASSERT(!err);

    pylith::topology::VecVisitorMesh fieldVisitor(*field);
    PylithScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);

    for (PylithInt p = 0; p < numPoints; ++p) {
        const PylithInt p_bc = points[p];

        const PylithInt off = fieldVisitor.sectionOffset(p_bc);
        const PylithInt dof = fieldVisitor.sectionDof(p_bc);
        for (PylithInt i = 0; i < dof; ++i) {
            fieldArray[off+i] = 0.0;
        } // for
    } // for

    err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _zeroBoundary


// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::feassemble::TestIntegratorDomain_Data::TestIntegratorDomain_Data(void) :
    dimension(0),
    meshFilename(0),
    boundaryLabel(NULL),
    materialId(0),
    cs(NULL),

    normalizer(new spatialdata::units::Nondimensional),

    t(0.0),
    dt(0.0),
    tindex(0),
    s_tshift(0.0),

    numSolutionSubfields(0),
    solutionDiscretizations(NULL),
    solutionDB(new spatialdata::spatialdb::UserFunctionDB),
    perturbationDB(new spatialdata::spatialdb::UserFunctionDB),

    numAuxiliarySubfields(0),
    auxiliarySubfields(NULL),
    auxiliaryDiscretizations(NULL),
    auxiliaryDB(new spatialdata::spatialdb::UserFunctionDB),
    auxiliaryUpdateDB(NULL),

    hasLHSJacobianLumpedInv(false) {
    CPPUNIT_ASSERT(normalizer);

    CPPUNIT_ASSERT(solutionDB);
    solutionDB->setLabel("solution");

    CPPUNIT_ASSERT(perturbationDB);
    perturbationDB->setLabel("solution+perturbation");

    CPPUNIT_ASSERT(auxiliaryDB);
    auxiliaryDB->setLabel("auxiliary field");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::feassemble::TestIntegratorDomain_Data::~TestIntegratorDomain_Data(void) {
    delete cs;cs = NULL;
    delete normalizer;normalizer = NULL;
    delete solutionDB;solutionDB = NULL;
    delete auxiliaryDB;auxiliaryDB = NULL;
    delete auxiliaryUpdateDB;auxiliaryUpdateDB = NULL;
} // destructor


// End of file
