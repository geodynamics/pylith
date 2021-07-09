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

#include "TestMaterial.hh" // Implementation of class methods

#include "pylith/materials/Material.hh" // USES Material
#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // pythia::journal

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestMaterial::setUp(void) {
    _mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    _solutionFields = NULL;
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::materials::TestMaterial::tearDown(void) {
    delete _solutionFields;_solutionFields = NULL;
    delete _mesh;_mesh = NULL;
} // tearDown


// ----------------------------------------------------------------------
// Test auxField().
void
pylith::materials::TestMaterial::testAuxField(void) {
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    Material* material = _material();CPPUNIT_ASSERT(material);
    TestMaterial_Data* data = _data();CPPUNIT_ASSERT(data);

    const pylith::topology::Field* auxField = material->auxField();CPPUNIT_ASSERT(auxField);
    for (int i = 0; i < data->numAuxSubfields; ++i) {
        CPPUNIT_ASSERT(auxField->hasSubfield(data->auxSubfields[i]));
    } // for

    CPPUNIT_ASSERT(!auxField->hasSubfield("abc4598245"));

    PYLITH_METHOD_END;
} // testAuxField


// ----------------------------------------------------------------------
// Test auxSubfieldDiscretization().
void
pylith::materials::TestMaterial::testAuxSubfieldDiscretization(void) {
    PYLITH_METHOD_BEGIN;

    const topology::FieldBase::Discretization infoDefault = pylith::topology::Field::Discretization(1, 1, true, pylith::topology::FieldBase::POLYNOMIAL_SPACE);
    const topology::FieldBase::Discretization infoA = pylith::topology::Field::Discretization(1, 2, false, pylith::topology::FieldBase::POLYNOMIAL_SPACE);
    const topology::FieldBase::Discretization infoB = pylith::topology::Field::Discretization(2, 2, true, pylith::topology::FieldBase::POINT_SPACE);

    Material* material = _material();CPPUNIT_ASSERT(material);
    material->auxSubfieldDiscretization("A", infoA.basisOrder, infoA.quadOrder, infoA.isBasisContinuous, infoA.feSpace);
    material->auxSubfieldDiscretization("B", infoB.basisOrder, infoB.quadOrder, infoB.isBasisContinuous, infoB.feSpace);

    CPPUNIT_ASSERT(material->_auxiliaryFactory());
    { // A
        const topology::FieldBase::Discretization& test = material->_auxiliaryFactory()->getSubfieldDiscretization("A");
        CPPUNIT_ASSERT_EQUAL(infoA.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoA.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoA.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoA.feSpace, test.feSpace);
    } // A

    { // B
        const topology::FieldBase::Discretization& test = material->_auxiliaryFactory()->getSubfieldDiscretization("B");
        CPPUNIT_ASSERT_EQUAL(infoB.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoB.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoB.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoB.feSpace, test.feSpace);
    } // B

    { // C (default)
        const topology::FieldBase::Discretization& test = material->_auxiliaryFactory()->getSubfieldDiscretization("C");
        CPPUNIT_ASSERT_EQUAL(infoDefault.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoDefault.feSpace, test.feSpace);
    } // C (default)

    { // default
        const topology::FieldBase::Discretization& test = material->_auxiliaryFactory()->getSubfieldDiscretization("default");
        CPPUNIT_ASSERT_EQUAL(infoDefault.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoDefault.feSpace, test.feSpace);
    } // default

    PYLITH_METHOD_END;
} // testAuxSubfieldDiscretization


// ----------------------------------------------------------------------
// Test auxFieldDB().
void
pylith::materials::TestMaterial::testAuxFieldDB(void) {
    PYLITH_METHOD_BEGIN;

    const std::string label = "test db";
    spatialdata::spatialdb::UserFunctionDB db;
    db.setLabel(label.c_str());

    Material* material = _material();CPPUNIT_ASSERT(material);
    material->auxFieldDB(&db);

    CPPUNIT_ASSERT(material->_auxiliaryFactory());
    CPPUNIT_ASSERT(material->_auxiliaryFactory()->queryDB());
    CPPUNIT_ASSERT_EQUAL(label, std::string(material->_auxiliaryFactory()->queryDB()->getLabel()));

    PYLITH_METHOD_END;
} // testAuxFieldDB


// ----------------------------------------------------------------------
// Test normalizer().
void
pylith::materials::TestMaterial::testNormalizer(void) {
    PYLITH_METHOD_BEGIN;

    spatialdata::units::Nondimensional normalizer;
    const double scale = 5.0;
    normalizer.setLengthScale(scale);

    Material* material = _material();CPPUNIT_ASSERT(material);
    material->normalizer(normalizer);
    CPPUNIT_ASSERT_EQUAL(scale, material->_normalizer->getLengthScale());

    PYLITH_METHOD_END;
} // testNormalizer


// ----------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::materials::TestMaterial::testVerifyConfiguration(void) {
    PYLITH_METHOD_BEGIN;

    // Call verifyConfiguration()
    Material* material = _material();CPPUNIT_ASSERT(material);
    CPPUNIT_ASSERT(_solutionFields);
    material->verifyConfiguration(_solutionFields->get("solution"));

    // Nothing to test.

    PYLITH_METHOD_END;
} // testVerifyConfiguration


// ----------------------------------------------------------------------
// Test dimension(), id(), and getLabel().
void
pylith::materials::TestMaterial::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    Material* material = _material();CPPUNIT_ASSERT(material);
    TestMaterial_Data* data = _data();CPPUNIT_ASSERT(data);

    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of Material::dimension() failed.", data->dimension, material->dimension());

    const int matId = 1234;
    material->id(matId);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of Material::id() failed.", matId, material->id());

    const std::string& matLabel = "xyz";
    material->setLabel(matLabel.c_str());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of Material::getLabel() failed.", matLabel, std::string(material->getLabel()));

    PYLITH_METHOD_END;
} // testAccessors


// ----------------------------------------------------------------------
// Test initialize().
void
pylith::materials::TestMaterial::testInitialize(void) {
    PYLITH_METHOD_BEGIN;

    // Call initialize()
    _initializeFull(); // includes setting up auxField

    Material* material = _material();CPPUNIT_ASSERT(material);
    const pylith::topology::Field* auxField = material->auxField();CPPUNIT_ASSERT(auxField);

    // material->_auxiliaryField->view("AUX FIELDS"); // :DEBUGGING:

    // Check result
    TestMaterial_Data* data = _data();CPPUNIT_ASSERT(data);
    CPPUNIT_ASSERT_EQUAL(std::string("auxiliary subfields"), std::string(auxField->getLabel()));
    CPPUNIT_ASSERT_EQUAL(data->dimension, auxField->getSpaceDim());

    PylithReal norm = 0.0;
    PylithReal t = 0.0;
    const PetscDM dm = auxField->dmMesh();CPPUNIT_ASSERT(dm);
    pylith::topology::FieldQuery query(*auxField);
    query.initializeWithDefaultQueryFns();
    CPPUNIT_ASSERT(data->normalizer);
    query.openDB(data->auxDB, data->normalizer->getLengthScale());
    PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, query.functions(), (void**)query.contextPtrs(), auxField->localVector(), &norm);CPPUNIT_ASSERT(!err);
    query.closeDB(data->auxDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test of auxiliary field values failed.", 0.0, norm, tolerance);

#if 1
    // Verify solution and perturbation fields can be exactly represented by discretization.
    norm = 0.0;
    t = 0.0;

    pylith::topology::Field& solution = _solutionFields->get("solution");
    // solution.view("SOLUTION"); // :DEBUG:
    const PetscDM dmSoln = solution.dmMesh();CPPUNIT_ASSERT(dmSoln);
    pylith::topology::FieldQuery solnQuery(solution);
    solnQuery.initializeWithDefaultQueryFns();
    CPPUNIT_ASSERT(data->normalizer);
    solnQuery.openDB(data->solnDB, data->normalizer->getLengthScale());
    err = DMPlexComputeL2DiffLocal(dmSoln, t, solnQuery.functions(), (void**)solnQuery.contextPtrs(), solution.localVector(), &norm);CPPUNIT_ASSERT(!err);
    solnQuery.closeDB(data->solnDB);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Discretized solution field failed representation test.", 0.0, norm, tolerance);

    pylith::topology::Field& perturbation = _solutionFields->get("perturbation");
    // perturbation.view("PERTURBATION"); // :DEBUG:
    const PetscDM dmPerturb = perturbation.dmMesh();CPPUNIT_ASSERT(dmPerturb);
    pylith::topology::FieldQuery perturbQuery(perturbation);
    perturbQuery.initializeWithDefaultQueryFns();
    CPPUNIT_ASSERT(data->normalizer);
    perturbQuery.openDB(data->perturbDB, data->normalizer->getLengthScale());
    err = DMPlexComputeL2DiffLocal(dmPerturb, t, perturbQuery.functions(), (void**)perturbQuery.contextPtrs(), perturbation.localVector(), &norm);CPPUNIT_ASSERT(!err);
    perturbQuery.closeDB(data->perturbDB);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Discretized perturbation field failed representation test.", 0.0, norm, tolerance);
#endif

    PYLITH_METHOD_END;
} // testInitialize


// ----------------------------------------------------------------------
// Test computeResidual().
void
pylith::materials::TestMaterial::testComputeResidual(void) {
    PYLITH_METHOD_BEGIN;

    // Call initialize()
    _initializeFull(); // includes setting up auxField

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_solutionFields);
    pylith::topology::Field& solution = _solutionFields->get("solution");
    pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");

    pylith::topology::Field residualRHS(*_mesh);
    residualRHS.cloneSection(solution);
    residualRHS.setLabel("residual RHS");
    residualRHS.createDiscretization();
    residualRHS.allocate();

    pylith::topology::Field residualLHS(*_mesh);
    residualLHS.cloneSection(solution);
    residualLHS.setLabel("residual LHS");
    residualLHS.createDiscretization();
    residualLHS.allocate();

    Material* material = _material();CPPUNIT_ASSERT(material);
    TestMaterial_Data* data = _data();CPPUNIT_ASSERT(data);

#if 0 // :DEBUG:
    PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2"); // :DEBUG:
    DMSetFromOptions(residualRHS.dmMesh()); // :DEBUG:
#endif // :DEBUG:

    const PylithReal t = data->t;
    const PylithReal dt = data->dt;
    material->computeRHSResidual(&residualRHS, t, dt, solution);
    material->computeLHSResidual(&residualLHS, t, dt, solution, solutionDot);

    // We don't use Dirichlet BC, so we must manually zero out the residual values for constrained DOF.
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
    // Avoid trivial satisfaction of norm with zero values.
    CPPUNIT_ASSERT_MESSAGE("RHS and LHS residuals are both exactly zero, which is suspicious.", normRHS > 0.0 || normLHS > 0.0);

    PYLITH_METHOD_END;
} // testComputeResidual


// ----------------------------------------------------------------------
// Test computeJacobian().
void
pylith::materials::TestMaterial::testComputeJacobian(void) {
    PYLITH_METHOD_BEGIN;

    // Create linear problem (MMS) with two trial solutions, s and p.
    //
    // Check that Jg(s)*(p - s) = G(p) - G(s).

    // Call initialize()
    _initializeFull();

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_solutionFields);
    pylith::topology::Field& solution = _solutionFields->get("solution");
    pylith::topology::Field& perturbation = _solutionFields->get("perturbation");

    Material* material = _material();CPPUNIT_ASSERT(material);
    TestMaterial_Data* data = _data();CPPUNIT_ASSERT(data);

    pylith::topology::Field residual1(*_mesh);
    residual1.cloneSection(solution);
    residual1.setLabel("residual1");
    residual1.createDiscretization();
    residual1.allocate();

    pylith::topology::Field residual2(*_mesh);
    residual2.cloneSection(perturbation);
    residual2.setLabel("residual2");
    residual2.createDiscretization();
    residual2.allocate();

#if 0 // :DEBUG:
    PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "3"); // :DEBUG:
    DMSetFromOptions(_solution1->dmMesh()); // :DEBUG:
#endif // :DEBUG:

    const PylithReal t = data->t;
    const PylithReal dt = data->dt;
    material->computeLHSResidual(&residual1, t, dt, solution);
    material->computeLHSResidual(&residual2, t, dt, perturbation);

    // residual1.view("RESIDUAL 1 RHS"); // :DEBUG:
    // residual2.view("RESIDUAL 2 RHS"); // :DEBUG:

    // Compute Jacobian
    PetscErrorCode err;
    PetscMat jacobianMat = NULL;
    err = DMCreateMatrix(solution.dmMesh(), &jacobianMat);CPPUNIT_ASSERT(!err);
    err = MatZeroEntries(jacobianMat);CPPUNIT_ASSERT(!err);
    PetscMat precondMat = jacobianMat; // Use Jacobian == preconditioner

    material->computeLHSJacobian(jacobianMat, precondMat, t, dt, solution);
    CPPUNIT_ASSERT_EQUAL(false, material->needNewLHSJacobian());
    // _zeroBoundary(&residual1);
    // _zeroBoundary(&residual2, jacobianMat);
    err = MatAssemblyBegin(jacobianMat, MAT_FINAL_ASSEMBLY);PYLITH_CHECK_ERROR(err);
    err = MatAssemblyEnd(jacobianMat, MAT_FINAL_ASSEMBLY);PYLITH_CHECK_ERROR(err);

    // Check that J(s)*(p - s) = G(p) - G(s).

    PetscVec residualVec = NULL;
    err = VecDuplicate(residual1.localVector(), &residualVec);CPPUNIT_ASSERT(!err);
    err = VecWAXPY(residualVec, -1.0, residual1.localVector(), residual2.localVector());CPPUNIT_ASSERT(!err);

    PetscVec solnIncrVec = NULL;
    err = VecDuplicate(solution.localVector(), &solnIncrVec);CPPUNIT_ASSERT(!err);
    err = VecWAXPY(solnIncrVec, -1.0, solution.localVector(), perturbation.localVector());CPPUNIT_ASSERT(!err);

    // result = Jg*(-solnIncr) + residual
    PetscVec resultVec = NULL;
    err = VecDuplicate(residualVec, &resultVec);CPPUNIT_ASSERT(!err);
    err = VecZeroEntries(resultVec);CPPUNIT_ASSERT(!err);
    err = VecScale(solnIncrVec, -1.0);CPPUNIT_ASSERT(!err);
    err = MatMultAdd(jacobianMat, solnIncrVec, residualVec, resultVec);CPPUNIT_ASSERT(!err);

#if 0 // :DEBUG:
    std::cout << "SOLN INCR" << std::endl;
    VecView(solnIncrVec, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "G2-G1" << std::endl;
    VecView(residualVec, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "RESULT" << std::endl;
    VecView(resultVec, PETSC_VIEWER_STDOUT_SELF);
#endif // :DEBUG:

    PylithReal norm = 0.0;
    err = VecNorm(resultVec, NORM_2, &norm);CPPUNIT_ASSERT(!err);
    err = VecDestroy(&resultVec);CPPUNIT_ASSERT(!err);
    err = VecDestroy(&solnIncrVec);CPPUNIT_ASSERT(!err);
    err = VecDestroy(&residualVec);CPPUNIT_ASSERT(!err);
    err = MatDestroy(&jacobianMat);CPPUNIT_ASSERT(!err);

    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Check of Jg(s)*(p-s) - (G(p) - G(s)) == 0 failed.", 0.0, norm, tolerance);
    CPPUNIT_ASSERT_MESSAGE("Norm of resulting vector is exactly zero, which is suspicious.", norm > 0.0);

    PYLITH_METHOD_END;
} // testComputeJacobian


// ----------------------------------------------------------------------
// Test computeLHSJacobianImplicit().
void
pylith::materials::TestMaterial::testComputeLHSJacobianImplicit(void) {
    PYLITH_METHOD_BEGIN;

    Material* material = _material();CPPUNIT_ASSERT(material);
    const TestMaterial_Data* data = _data();CPPUNIT_ASSERT(data);
    if (data->isExplicit) {
        PYLITH_METHOD_END;
    } // if

    // Create linear problem (MMS) with two trial solutions, s,s_dor and p,p_dot.
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

    pylith::topology::Field residual1(*_mesh);
    residual1.cloneSection(solution);
    residual1.setLabel("residual1");
    residual1.createDiscretization();
    residual1.allocate();

    pylith::topology::Field residual2(*_mesh);
    residual2.cloneSection(perturbation);
    residual2.setLabel("residual2");
    residual2.createDiscretization();
    residual2.allocate();

#if 0 // :DEBUG:
    PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2"); // :DEBUG:
    DMSetFromOptions(_solution1->dmMesh()); // :DEBUG:
#endif // :DEBUG:

    const PylithReal t = data->t;
    const PylithReal dt = data->dt;
    const PylithReal s_tshift = data->s_tshift;
    material->computeLHSResidual(&residual1, t, dt, solution, solutionDot);
    material->computeLHSResidual(&residual2, t, dt, perturbation, perturbationDot);

    // residual1.view("RESIDUAL 1 LHS"); // :DEBUG:
    // residual2.view("RESIDUAL 2 LHS"); // :DEBUG:

    PetscErrorCode err;

    PetscVec residualVec = NULL;
    err = VecDuplicate(residual1.localVector(), &residualVec);CPPUNIT_ASSERT(!err);
    err = VecWAXPY(residualVec, -1.0, residual1.localVector(), residual2.localVector());CPPUNIT_ASSERT(!err);

    PetscVec solnIncrVec = NULL;
    err = VecDuplicate(solution.localVector(), &solnIncrVec);CPPUNIT_ASSERT(!err);
    err = VecWAXPY(solnIncrVec, -1.0, solution.localVector(), perturbation.localVector());CPPUNIT_ASSERT(!err);

    // Compute Jacobian
    PetscMat jacobianMat = NULL;
    err = DMCreateMatrix(solution.dmMesh(), &jacobianMat);CPPUNIT_ASSERT(!err);
    err = MatZeroEntries(jacobianMat);CPPUNIT_ASSERT(!err);
    PetscMat precondMat = jacobianMat; // Use Jacobian == preconditioner

    material->computeLHSJacobianImplicit(jacobianMat, precondMat, t, dt, s_tshift, solution, solutionDot);
    CPPUNIT_ASSERT_EQUAL(false, material->needNewLHSJacobian());
    err = MatAssemblyBegin(jacobianMat, MAT_FINAL_ASSEMBLY);PYLITH_CHECK_ERROR(err);
    err = MatAssemblyEnd(jacobianMat, MAT_FINAL_ASSEMBLY);PYLITH_CHECK_ERROR(err);

    // result = J*(-solnIncr) + residual
    PetscVec resultVec = NULL;
    err = VecDuplicate(residualVec, &resultVec);CPPUNIT_ASSERT(!err);
    err = VecZeroEntries(resultVec);CPPUNIT_ASSERT(!err);
    err = VecScale(solnIncrVec, -1.0);CPPUNIT_ASSERT(!err);
    err = MatMultAdd(jacobianMat, solnIncrVec, residualVec, resultVec);CPPUNIT_ASSERT(!err);

#if 0 // :DEBUG:
    std::cout << "SOLN INCR" << std::endl;
    VecView(solnIncrVec, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "F2-F1" << std::endl;
    VecView(residualVec, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "RESULT" << std::endl;
    VecView(resultVec, PETSC_VIEWER_STDOUT_SELF);
#endif // :DEBUG:

    PylithReal norm = 0.0, normSolnIncr = 0.0, normResidual = 0.0;
    err = VecNorm(resultVec, NORM_2, &norm);CPPUNIT_ASSERT(!err);
    err = VecNorm(solnIncrVec, NORM_2, &normSolnIncr);CPPUNIT_ASSERT(!err);
    err = VecNorm(residualVec, NORM_2, &normResidual);CPPUNIT_ASSERT(!err);
    err = VecDestroy(&resultVec);CPPUNIT_ASSERT(!err);
    err = VecDestroy(&solnIncrVec);CPPUNIT_ASSERT(!err);
    err = VecDestroy(&residualVec);CPPUNIT_ASSERT(!err);
    err = MatDestroy(&jacobianMat);CPPUNIT_ASSERT(!err);

    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Check of Jf(s)*(p-s) - (F(p) - F(s)) == 0 failed.", 0.0, norm, tolerance);
    CPPUNIT_ASSERT_MESSAGE("Norm of resulting vector is exactly zero, which is suspicious.", (0 < normResidual && 0 < norm) || (0 == normResidual && 0 == norm));

    PYLITH_METHOD_END;
} // testComputeLHSJacobianImplicit


// ----------------------------------------------------------------------
// Test computeLHSJacobianExplicit().
void
pylith::materials::TestMaterial::testComputeLHSJacobianInverseExplicit(void) {
    PYLITH_METHOD_BEGIN;

    Material* material = _material();CPPUNIT_ASSERT(material);
    TestMaterial_Data* data = _data();CPPUNIT_ASSERT(data);
    if (!data->isExplicit) {
        PYLITH_METHOD_END;
    } // if

    CPPUNIT_ASSERT_MESSAGE("Test not implemented.", false); // :TODO: ADD MORE HERE

    PYLITH_METHOD_END;
} // testComputeLHSJacobianInverseExplicit


// ----------------------------------------------------------------------
// Test updateStateVars().
void
pylith::materials::TestMaterial::testUpdateStateVars(void) {
    PYLITH_METHOD_BEGIN;

    TestMaterial_Data* data = _data();CPPUNIT_ASSERT(data);
    if (!data->auxUpdateDB) {
        PYLITH_METHOD_END;
    } // if

    // Call initialize()
    _initializeFull(); // includes setting up auxField

    // We test updating the state variables in the auxiliary field by
    // passing the perturbation as the "new" solution and the existing
    // auxiliary field. We test whether the "updated" auxiliary field
    // matches the database with the updated auxiliary field.

    Material* material = _material();CPPUNIT_ASSERT(material);
    CPPUNIT_ASSERT(_solutionFields);
    pylith::topology::Field& perturbation = _solutionFields->get("perturbation");
#if 0
    material->_auxiliaryField->view("INITIAL_AUX FIELDS"); // :DEBUGGING:
#endif
    material->_updateStateVars(data->t, data->dt, perturbation);

    const pylith::topology::Field* auxField = material->auxField();CPPUNIT_ASSERT(auxField);
    material->_auxiliaryField->view("UPDATED_AUX FIELDS"); // :DEBUGGING:

    // Check updated auxiliary field.
    PylithReal norm = 0.0;
    PylithReal t = 0.0;
    const PetscDM dm = auxField->dmMesh();CPPUNIT_ASSERT(dm);
    pylith::topology::FieldQuery query(*auxField);
    query.initializeWithDefaultQueryFns();
    CPPUNIT_ASSERT(data->normalizer);
    query.openDB(data->auxUpdateDB, data->normalizer->getLengthScale());
#if 0 // :DEBUG:
    PetscOptionsSetValue(NULL, "-dm_plex_print_l2", "1"); // :DEBUG:
    DMSetFromOptions(dm); // :DEBUG:
#endif
    PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, query.functions(), (void**)query.contextPtrs(), auxField->localVector(), &norm);CPPUNIT_ASSERT(!err);
    query.closeDB(data->auxUpdateDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Check of updated auxiliary field values failed.", 0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testUpdateStateVars


// ----------------------------------------------------------------------
// Do minimal initilaization of test data.
void
pylith::materials::TestMaterial::_initializeMin(void) {
    PYLITH_METHOD_BEGIN;

    Material* material = _material();CPPUNIT_ASSERT(material);
    TestMaterial_Data* data = _data();CPPUNIT_ASSERT(data);

    pylith::meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(data->meshFilename);
    iohandler.filename(data->meshFilename);
    iohandler.read(_mesh);CPPUNIT_ASSERT(_mesh);

    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any cells.", _mesh->numCells() > 0);
    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any vertices.", _mesh->numVertices() > 0);

    // Setup coordinates.
    _mesh->setCoordSys(data->cs);
    CPPUNIT_ASSERT(data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *data->normalizer);

    // id and label initialized in derived class
    material->normalizer(*data->normalizer);
    material->gravityField(data->gravityField);

    // Setup solution fields.
    delete _solutionFields;_solutionFields = new pylith::topology::Fields(*_mesh);CPPUNIT_ASSERT(_solutionFields);
    _solutionFields->add("solution","solution");
    _solutionFields->add("solution_dot","solution_dot");
    _solutionFields->add("perturbation","perturbation");
    _solutionFields->add("perturbation_dot","perturbation_dot");
    this->_setupSolutionFields();

    PYLITH_METHOD_END;
} // _initializeMin


// ----------------------------------------------------------------------
// Complete initilaization of test data.
void
pylith::materials::TestMaterial::_initializeFull(void) {
    PYLITH_METHOD_BEGIN;

    Material* material = _material();CPPUNIT_ASSERT(material);
    TestMaterial_Data* data = _data();CPPUNIT_ASSERT(data);
    CPPUNIT_ASSERT(_mesh);

    // Set auxiliary fields spatial database.
    material->auxFieldDB(data->auxDB);

    for (int i = 0; i < data->numAuxSubfields; ++i) {
        const pylith::topology::FieldBase::Discretization& info = data->auxDiscretizations[i];
        material->auxSubfieldDiscretization(data->auxSubfields[i], info.basisOrder, info.quadOrder, info.isBasisContinuous, info.feSpace);
    } // for

    CPPUNIT_ASSERT(_solutionFields);
    material->initialize(_solutionFields->get("solution"));

    PYLITH_METHOD_END;
} // _initializeFull


// ----------------------------------------------------------------------
// Set field to zero on the boundary.
void
pylith::materials::TestMaterial::_zeroBoundary(pylith::topology::Field* field) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(field);
    TestMaterial_Data* data = _data();CPPUNIT_ASSERT(data);
    CPPUNIT_ASSERT(data->boundaryLabel);

    PetscDM dmMesh = field->mesh().dmMesh();CPPUNIT_ASSERT(dmMesh);
    PetscDMLabel label = NULL;
    PetscIS pointIS = NULL;
    const PetscInt *points;
    PetscInt numPoints = 0;
    PetscBool hasLabel = PETSC_FALSE;
    PetscErrorCode err;
    err = DMHasLabel(dmMesh, data->boundaryLabel, &hasLabel);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(hasLabel);
    err = DMGetLabel(dmMesh, data->boundaryLabel, &label);CPPUNIT_ASSERT(!err);
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


// ----------------------------------------------------------------------
// Constructor
pylith::materials::TestMaterial_Data::TestMaterial_Data(void) :
    dimension(0),
    meshFilename(0),
    boundaryLabel(NULL),
    cs(NULL),
    gravityField(NULL),

    normalizer(new spatialdata::units::Nondimensional),

    t(0.0),
    dt(0.0),
    s_tshift(0.0),
    perturbation(1.0e-4),

    numSolnSubfields(0),
    solnDiscretizations(NULL),
    solnDB(new spatialdata::spatialdb::UserFunctionDB),
    perturbDB(new spatialdata::spatialdb::UserFunctionDB),

    numAuxSubfields(0),
    auxSubfields(NULL),
    auxDiscretizations(NULL),
    auxDB(new spatialdata::spatialdb::UserFunctionDB),
    auxUpdateDB(NULL),

    isExplicit(false) { // constructor
    CPPUNIT_ASSERT(normalizer);

    CPPUNIT_ASSERT(solnDB);
    solnDB->setLabel("solution");

    CPPUNIT_ASSERT(perturbDB);
    perturbDB->setLabel("solution+perturbation");

    CPPUNIT_ASSERT(auxDB);
    auxDB->setLabel("auxiliary field");
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::materials::TestMaterial_Data::~TestMaterial_Data(void) {
    delete cs;cs = NULL;
    delete gravityField;gravityField = NULL;
    delete normalizer;normalizer = NULL;
    delete solnDB;solnDB = NULL;
    delete auxDB;auxDB = NULL;
    delete auxUpdateDB;auxUpdateDB = NULL;
} // destructor


// End of file
