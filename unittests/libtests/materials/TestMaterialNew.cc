// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestMaterialNew.hh" // Implementation of class methods

#include "pylith/materials/MaterialNew.hh" // USES MaterialNew
#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/spatialdb/SimpleGridDB.hh" // USES SimpleDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "journal/debug.h" // USES journal::debug_t

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestMaterialNew::setUp(void)
{ // setUp
    _mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    _solutionFields = NULL;
    _auxDB = NULL;
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::materials::TestMaterialNew::tearDown(void)
{ // tearDown
    delete _solutionFields; _solutionFields = NULL;
    delete _mesh; _mesh = NULL;
    delete _auxDB; _auxDB = NULL;
} // tearDown


// ----------------------------------------------------------------------
// Test hasAuxField().
void
pylith::materials::TestMaterialNew::testHasAuxField(void)
{ // testHasAuxField
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

    for (int i = 0; i < data->numAuxFields; ++i) {
        CPPUNIT_ASSERT(material->hasAuxField(data->auxFields[i]));
    } // for

    CPPUNIT_ASSERT(!material->hasAuxField("abc4598245"));

    PYLITH_METHOD_END;
} // testHaxAuxField


// ----------------------------------------------------------------------
// Test auxFieldDiscretization().
void
pylith::materials::TestMaterialNew::testAuxFieldsDiscretization(void)
{ // testAuxFieldsDiscretization
    PYLITH_METHOD_BEGIN;

    const topology::FieldBase::Discretization infoDefault = {-1, -1, true, pylith::topology::FieldBase::POLYNOMIAL_SPACE};
    const topology::FieldBase::Discretization infoA = {1, 2, false, pylith::topology::FieldBase::POLYNOMIAL_SPACE};
    const topology::FieldBase::Discretization infoB = {2, 2, true, pylith::topology::FieldBase::POINT_SPACE};

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    material->auxFieldDiscretization("A", infoA.basisOrder, infoA.quadOrder, infoA.isBasisContinuous, infoA.feSpace);
    material->auxFieldDiscretization("B", infoB.basisOrder, infoB.quadOrder, infoB.isBasisContinuous, infoB.feSpace);

    { // A
        const topology::FieldBase::Discretization& test = material->auxFieldDiscretization("A");
        CPPUNIT_ASSERT_EQUAL(infoA.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoA.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoA.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoA.feSpace, test.feSpace);
    } // A

    { // B
        const topology::FieldBase::Discretization& test = material->auxFieldDiscretization("B");
        CPPUNIT_ASSERT_EQUAL(infoB.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoB.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoB.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoB.feSpace, test.feSpace);
    } // B

    { // C (default)
        const topology::FieldBase::Discretization& test = material->auxFieldDiscretization("C");
        CPPUNIT_ASSERT_EQUAL(infoDefault.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoDefault.feSpace, test.feSpace);
    } // C (default)

    { // default
        const topology::FieldBase::Discretization& test = material->auxFieldDiscretization("default");
        CPPUNIT_ASSERT_EQUAL(infoDefault.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoDefault.feSpace, test.feSpace);
    } // default

    PYLITH_METHOD_END;
} // testAuxFieldsDiscretization


// ----------------------------------------------------------------------
// Test auxFieldsDB().
void
pylith::materials::TestMaterialNew::testAuxFieldsDB(void)
{ // testAuxFieldsDB
    PYLITH_METHOD_BEGIN;

    const std::string label = "test db";
    spatialdata::spatialdb::SimpleGridDB db;
    db.label(label.c_str());

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    material->auxFieldsDB(&db);

    CPPUNIT_ASSERT(material->_auxFieldsDB);
    CPPUNIT_ASSERT_EQUAL(label, std::string(material->_auxFieldsDB->label()));

    PYLITH_METHOD_END;
} // testAuxFieldsDB


// ----------------------------------------------------------------------
// Test normalizer().
void
pylith::materials::TestMaterialNew::testNormalizer(void)
{ // testNormalizer
    PYLITH_METHOD_BEGIN;

    spatialdata::units::Nondimensional normalizer;
    const double scale = 5.0;
    normalizer.lengthScale(scale);

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    material->normalizer(normalizer);
    CPPUNIT_ASSERT_EQUAL(scale, material->_normalizer->lengthScale());

    PYLITH_METHOD_END;
} // testNormalizer


// ----------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::materials::TestMaterialNew::testVerifyConfiguration(void)
{ // testVerifyConfiguration
    PYLITH_METHOD_BEGIN;

    // Call verifyConfiguration()
    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    CPPUNIT_ASSERT(_solutionFields);
    material->verifyConfiguration(_solutionFields->get("solution"));

    // Nothing to test.

    PYLITH_METHOD_END;
} // testVerifyConfiguration


// ----------------------------------------------------------------------
// Test dimension().
void
pylith::materials::TestMaterialNew::testDimension(void)
{ // testDimension
    PYLITH_METHOD_BEGIN;

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);
    CPPUNIT_ASSERT_EQUAL(data->dimension, material->dimension());

    PYLITH_METHOD_END;
} // testDimension


// ----------------------------------------------------------------------
// Test id().
void
pylith::materials::TestMaterialNew::testId(void)
{ // testId
    PYLITH_METHOD_BEGIN;

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);
    CPPUNIT_ASSERT_EQUAL(data->materialId, material->id());

    PYLITH_METHOD_END;
} // testId


// ----------------------------------------------------------------------
// Test label().
void
pylith::materials::TestMaterialNew::testLabel(void)
{ // testLabel
    PYLITH_METHOD_BEGIN;

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);
    CPPUNIT_ASSERT_EQUAL(std::string(data->materialLabel), std::string(material->label()));

    PYLITH_METHOD_END;
} // testLabel


// ----------------------------------------------------------------------
// Test initialize().
void
pylith::materials::TestMaterialNew::testInitialize(void)
{ // testInitialize
    PYLITH_METHOD_BEGIN;

    // Call initialize()
    _initializeFull(); // includes setting up auxFields

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    const pylith::topology::Field& auxFields = material->auxFields();

    //material->_auxFields->view("AUX FIELDS"); // :DEBUGGING:

    // Check result
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);
    CPPUNIT_ASSERT_EQUAL(std::string("auxiliary fields"), std::string(auxFields.label()));
    CPPUNIT_ASSERT_EQUAL(data->dimension, auxFields.spaceDim());

    PylithReal norm = 0.0;
    PylithReal t = 0.0;
    const PetscDM dm = auxFields.dmMesh(); CPPUNIT_ASSERT(dm);
    pylith::topology::FieldQuery* query = material->_auxFieldsQuery;
    CPPUNIT_ASSERT(data->normalizer);
    query->openDB(_auxDB, data->normalizer->lengthScale());

    PetscErrorCode err = DMComputeL2Diff(dm, t, query->functions(), (void**)query->contextPtrs(), auxFields.localVector(), &norm); CPPUNIT_ASSERT(!err);
    query->closeDB(_auxDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testInitialize


// ----------------------------------------------------------------------
// Test computeResidual().
void
pylith::materials::TestMaterialNew::testComputeResidual(void)
{ // testComputeResidual
    PYLITH_METHOD_BEGIN;

    // Call initialize()
    _initializeFull(); // includes setting up auxFields

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_solutionFields);
    pylith::topology::Field& solution = _solutionFields->get("solution");
    pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");

    pylith::topology::Field residualRHS(*_mesh);
    residualRHS.cloneSection(solution);
    residualRHS.label("residual RHS");
    residualRHS.allocate();

    pylith::topology::Field residualLHS(*_mesh);
    residualLHS.cloneSection(solution);
    residualLHS.label("residual LHS");
    residualLHS.allocate();

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

#if 0 // DEBUGGING
    PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2"); // DEBUGGING
    DMSetFromOptions(residualRHS.dmMesh()); // DEBUGGING
#endif

    const PylithReal t = data->t;
    const PylithReal dt = data->dt;
    material->computeRHSResidual(&residualRHS, t, dt, solution);
    material->computeLHSResidual(&residualLHS, t, dt, solution, solutionDot);

    // We don't use Dirichlet BC, so we must manually zero out the residual values for constrained DOF.
    _zeroBoundary(&residualRHS);
    _zeroBoundary(&residualLHS);

#if 0 // DEBUGGING
    solution.view("SOLUTION"); // DEBUGGING
    solutionDot.view("SOLUTION_DOT"); // DEBUGGING
    residualRHS.view("RESIDUAL RHS"); // DEBUGGING
    residualLHS.view("RESIDUAL LHS"); // DEBUGGING
#endif

    PetscErrorCode err;
    PetscVec residualVec = NULL;
    err = VecDuplicate(residualRHS.localVector(), &residualVec); CPPUNIT_ASSERT(!err);
    err = VecWAXPY(residualVec, -1.0, residualRHS.localVector(), residualLHS.localVector()); CPPUNIT_ASSERT(!err);

    PylithReal norm = 0.0, normRHS = 0.0, normLHS = 0.0;
    err = VecNorm(residualRHS.localVector(), NORM_2, &normRHS); CPPUNIT_ASSERT(!err);
    err = VecNorm(residualLHS.localVector(), NORM_2, &normLHS); CPPUNIT_ASSERT(!err);
    err = VecNorm(residualVec, NORM_2, &norm); CPPUNIT_ASSERT(!err);
    err = VecDestroy(&residualVec); CPPUNIT_ASSERT(!err);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT(normRHS > 0.0 || normLHS > 0.0); // Avoid trivial satisfaction of norm with zero values.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testComputeResidual


// ----------------------------------------------------------------------
// Test computeRHSJacobian().
void
pylith::materials::TestMaterialNew::testComputeRHSJacobian(void)
{ // testComputeRHSJacobian
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

    pylith::topology::Field residual1(*_mesh);
    residual1.cloneSection(solution);
    residual1.label("residual1");
    residual1.allocate();

    pylith::topology::Field residual2(*_mesh);
    residual2.cloneSection(perturbation);
    residual2.label("residual2");
    residual2.allocate();

#if 0 // DEBUGGING
    PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2"); // DEBUGGING
    DMSetFromOptions(_solution1->dmMesh()); // DEBUGGING
#endif

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

    const PylithReal t = data->t;
    const PylithReal dt = data->dt;
    material->computeRHSResidual(&residual1, t, dt, solution);
    material->computeRHSResidual(&residual2, t, dt, perturbation);

    //residual1.view("RESIDUAL 1 RHS"); // DEBUGGING
    //residual2.view("RESIDUAL 1 RHS"); // DEBUGGING

    // Check that J(s)*(p - s) = G(p) - G(s).

    PetscErrorCode err;

    PetscVec residualVec = NULL;
    err = VecDuplicate(residual1.localVector(), &residualVec); CPPUNIT_ASSERT(!err);
    err = VecWAXPY(residualVec, -1.0, residual1.localVector(), residual2.localVector()); CPPUNIT_ASSERT(!err);

    PetscVec solnIncrVec = NULL;
    err = VecDuplicate(solution.localVector(), &solnIncrVec); CPPUNIT_ASSERT(!err);
    err = VecWAXPY(solnIncrVec, -1.0, solution.localVector(), perturbation.localVector()); CPPUNIT_ASSERT(!err);

    // Compute Jacobian
    PetscMat jacobianMat = NULL;
    err = DMCreateMatrix(solution.dmMesh(), &jacobianMat); CPPUNIT_ASSERT(!err);
    err = MatZeroEntries(jacobianMat); CPPUNIT_ASSERT(!err);
    PetscMat precondMat = jacobianMat; // Use Jacobian == preconditioner

    material->computeRHSJacobian(jacobianMat, precondMat, t, dt, solution);
    CPPUNIT_ASSERT_EQUAL(false, material->needNewRHSJacobian());
    err = MatAssemblyBegin(jacobianMat, MAT_FINAL_ASSEMBLY); PYLITH_CHECK_ERROR(err);
    err = MatAssemblyEnd(jacobianMat, MAT_FINAL_ASSEMBLY); PYLITH_CHECK_ERROR(err);

    // result = Jg*(-solnIncr) + residual
    PetscVec resultVec = NULL;
    err = VecDuplicate(residualVec, &resultVec); CPPUNIT_ASSERT(!err);
    err = VecZeroEntries(resultVec); CPPUNIT_ASSERT(!err);
    err = VecScale(solnIncrVec, -1.0); CPPUNIT_ASSERT(!err);
    err = MatMultAdd(jacobianMat, solnIncrVec, residualVec, resultVec); CPPUNIT_ASSERT(!err);

#if 0 // DEBUGGING
    std::cout << "SOLN INCR" << std::endl;
    VecView(solnIncrVec, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "G2-G1" << std::endl;
    VecView(residualVec, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "RESULT" << std::endl;
    VecView(resultVec, PETSC_VIEWER_STDOUT_SELF);
#endif

    PylithReal norm = 0.0;
    err = VecNorm(resultVec, NORM_2, &norm); CPPUNIT_ASSERT(!err);
    err = VecDestroy(&resultVec); CPPUNIT_ASSERT(!err);
    err = VecDestroy(&solnIncrVec); CPPUNIT_ASSERT(!err);
    err = VecDestroy(&residualVec); CPPUNIT_ASSERT(!err);
    err = MatDestroy(&jacobianMat); CPPUNIT_ASSERT(!err);

    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
    CPPUNIT_ASSERT(norm > 0.0); // Norm exactly equal to zero almost certainly means test is satisfied trivially.

    PYLITH_METHOD_END;
} // testComputeRHSJacobian


// ----------------------------------------------------------------------
// Test computeLHSJacobianImplicit().
void
pylith::materials::TestMaterialNew::testComputeLHSJacobianImplicit(void)
{ // testComputeLHSJacobianImplicit
    PYLITH_METHOD_BEGIN;

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    const TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);
    if (data->isExplicit) {
        PYLITH_METHOD_END;
    } // if

    // Create linear problem (MMS) with two solutions, s_1 and s_2.
    //
    // Check that Jf(s_1)*(s_2 - s_1) = F(s_2) - F(s_1).

    // Call initialize()
    _initializeFull(); // includes setting up auxFields

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_solutionFields);
    pylith::topology::Field& solution = _solutionFields->get("solution");
    pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");
    pylith::topology::Field& perturbation = _solutionFields->get("perturbation");
    pylith::topology::Field& perturbationDot = _solutionFields->get("perturbation_dot");

    pylith::topology::Field residual1(*_mesh);
    residual1.cloneSection(solution);
    residual1.label("residual1");
    residual1.allocate();

    pylith::topology::Field residual2(*_mesh);
    residual2.cloneSection(perturbation);
    residual2.label("residual2");
    residual2.allocate();

#if 0 // DEBUGGING
    PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2"); // DEBUGGING
    DMSetFromOptions(_solution1->dmMesh()); // DEBUGGING
#endif


    const PylithReal t = data->t;
    const PylithReal dt = data->dt;
    const PylithReal tshift = data->tshift;
    material->computeLHSResidual(&residual1, t, dt, solution, solutionDot);
    material->computeLHSResidual(&residual2, t, dt, perturbation, perturbationDot);

    // Check that Jf(s)*(p - s) = F(p) - F(s).

    //residual1.view("RESIDUAL 1 LHS"); // DEBUGGING
    //residual2.view("RESIDUAL 2 LHS"); // DEBUGGING

    PetscErrorCode err;

    PetscVec residualVec = NULL;
    err = VecDuplicate(residual1.localVector(), &residualVec); CPPUNIT_ASSERT(!err);
    err = VecWAXPY(residualVec, -1.0, residual1.localVector(), residual2.localVector()); CPPUNIT_ASSERT(!err);

    PetscVec solnIncrVec = NULL;
    err = VecDuplicate(solution.localVector(), &solnIncrVec); CPPUNIT_ASSERT(!err);
    err = VecWAXPY(solnIncrVec, -1.0, solution.localVector(), perturbation.localVector()); CPPUNIT_ASSERT(!err);

    // Compute Jacobian
    PetscMat jacobianMat = NULL;
    err = DMCreateMatrix(solution.dmMesh(), &jacobianMat); CPPUNIT_ASSERT(!err);
    err = MatZeroEntries(jacobianMat); CPPUNIT_ASSERT(!err);
    PetscMat precondMat = jacobianMat; // Use Jacobian == preconditioner

    material->computeLHSJacobianImplicit(jacobianMat, precondMat, t, dt, tshift, solution, solutionDot);
    CPPUNIT_ASSERT_EQUAL(false, material->needNewLHSJacobian());
    err = MatAssemblyBegin(jacobianMat, MAT_FINAL_ASSEMBLY); PYLITH_CHECK_ERROR(err);
    err = MatAssemblyEnd(jacobianMat, MAT_FINAL_ASSEMBLY); PYLITH_CHECK_ERROR(err);

    // result = J*(-solnIncr) + residual
    PetscVec resultVec = NULL;
    err = VecDuplicate(residualVec, &resultVec); CPPUNIT_ASSERT(!err);
    err = VecZeroEntries(resultVec); CPPUNIT_ASSERT(!err);
    err = VecScale(solnIncrVec, -1.0); CPPUNIT_ASSERT(!err);
    err = MatMultAdd(jacobianMat, solnIncrVec, residualVec, resultVec); CPPUNIT_ASSERT(!err);

#if 0 // DEBUGGING
    std::cout << "SOLN INCR" << std::endl;
    VecView(solnIncrVec, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "F2-F1" << std::endl;
    VecView(residualVec, PETSC_VIEWER_STDOUT_SELF);
    std::cout << "RESULT" << std::endl;
    VecView(resultVec, PETSC_VIEWER_STDOUT_SELF);
#endif

    PylithReal norm = 0.0, normSolnIncr = 0.0, normResidual = 0.0;
    err = VecNorm(resultVec, NORM_2, &norm); CPPUNIT_ASSERT(!err);
    err = VecNorm(solnIncrVec, NORM_2, &normSolnIncr); CPPUNIT_ASSERT(!err);
    err = VecNorm(residualVec, NORM_2, &normResidual); CPPUNIT_ASSERT(!err);
    err = VecDestroy(&resultVec); CPPUNIT_ASSERT(!err);
    err = VecDestroy(&solnIncrVec); CPPUNIT_ASSERT(!err);
    err = VecDestroy(&residualVec); CPPUNIT_ASSERT(!err);
    err = MatDestroy(&jacobianMat); CPPUNIT_ASSERT(!err);

    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
    CPPUNIT_ASSERT((0 < normResidual && 0 < norm) || (0 == normResidual && 0 == norm));

    PYLITH_METHOD_END;
} // testComputeLHSJacobianImplicit


// ----------------------------------------------------------------------
// Test computeLHSJacobianExplicit().
void
pylith::materials::TestMaterialNew::testComputeLHSJacobianInverseExplicit(void)
{ // testComputeLHSJacobianInverseExplicit
    PYLITH_METHOD_BEGIN;

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);
    if (!data->isExplicit) {
        PYLITH_METHOD_END;
    } // if

    CPPUNIT_ASSERT_MESSAGE("Test not implemented.", false); // :TODO: ADD MORE HERE

    PYLITH_METHOD_END;
} // testComputeLHSJacobianInverseExplicit


// ----------------------------------------------------------------------
// Test updateStateVars().
void
pylith::materials::TestMaterialNew::testUpdateStateVars(void)
{ // testUpdateStateVars
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT_MESSAGE("Test not implemented.", false); // :TODO: ADD MORE HERE

    PYLITH_METHOD_END;
} // testUpdateStateVars


// ----------------------------------------------------------------------
// Do minimal initilaization of test data.
void
pylith::materials::TestMaterialNew::_initializeMin(void)
{ // _initializeMin
    PYLITH_METHOD_BEGIN;

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

    pylith::meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(data->meshFilename);
    iohandler.filename(data->meshFilename);
    iohandler.read(_mesh); CPPUNIT_ASSERT(_mesh);

    // Setup coordinates.
    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_mesh->dimension());
    cs.initialize();
    _mesh->coordsys(&cs);
    CPPUNIT_ASSERT(data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *data->normalizer);

    material->id(data->materialId);
    material->label(data->materialLabel);
    material->normalizer(*data->normalizer);

    // Setup solution fields.
    delete _solutionFields; _solutionFields = new pylith::topology::Fields(*_mesh);CPPUNIT_ASSERT(_solutionFields);
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
pylith::materials::TestMaterialNew::_initializeFull(void)
{ // _initializeFull
    PYLITH_METHOD_BEGIN;

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);
    CPPUNIT_ASSERT(_mesh);

    // Set auxiliary fields spatial database.
    delete _auxDB; _auxDB = new spatialdata::spatialdb::SimpleGridDB; CPPUNIT_ASSERT(_auxDB);
    CPPUNIT_ASSERT(data->auxDBFilename);
    _auxDB->filename(data->auxDBFilename);
    _auxDB->label("IsotropicLinearElasciticityPlaneStrain auxiliary fields database");
    _auxDB->queryType(spatialdata::spatialdb::SimpleGridDB::LINEAR);
    material->auxFieldsDB(_auxDB);

    for (int i = 0; i < data->numAuxFields; ++i) {
        const pylith::topology::FieldBase::Discretization& info = data->auxDiscretizations[i];
        material->auxFieldDiscretization(data->auxFields[i], info.basisOrder, info.quadOrder, info.isBasisContinuous, info.feSpace);
    } // for

    CPPUNIT_ASSERT(_solutionFields);
    material->initialize(_solutionFields->get("solution"));

    PYLITH_METHOD_END;
} // _initializeFull


// ----------------------------------------------------------------------
// Set field to zero on the boundary.
void
pylith::materials::TestMaterialNew::_zeroBoundary(pylith::topology::Field* field)
{ // _zeroBoundary
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(field);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);
    CPPUNIT_ASSERT(data->boundaryLabel);

    PetscDM dmMesh = field->mesh().dmMesh(); CPPUNIT_ASSERT(dmMesh);
    PetscDMLabel label = NULL;
    PetscIS pointIS = NULL;
    const PetscInt *points;
    PetscInt numPoints = 0;
    PetscBool hasLabel = PETSC_FALSE;
    PetscErrorCode err;
    err = DMHasLabel(dmMesh, data->boundaryLabel, &hasLabel); CPPUNIT_ASSERT(!err); CPPUNIT_ASSERT(hasLabel);
    err = DMGetLabel(dmMesh, data->boundaryLabel, &label); CPPUNIT_ASSERT(!err);
    err = DMLabelGetStratumIS(label, 1, &pointIS); CPPUNIT_ASSERT(!err); CPPUNIT_ASSERT(pointIS);
    err = ISGetLocalSize(pointIS, &numPoints); CPPUNIT_ASSERT(!err);
    err = ISGetIndices(pointIS, &points); CPPUNIT_ASSERT(!err);

    pylith::topology::VecVisitorMesh fieldVisitor(*field);
    PylithScalar* fieldArray = fieldVisitor.localArray(); CPPUNIT_ASSERT(fieldArray);

    for (PetscInt p = 0; p < numPoints; ++p) {
        const PetscInt p_bc = points[p];

        const PylithInt off = fieldVisitor.sectionOffset(p_bc);
        const PylithInt dof = fieldVisitor.sectionDof(p_bc);
        for (PylithInt i = 0; i < dof; ++i) {
            fieldArray[off+i] = 0.0;
        } // for
    } // for

    err = ISRestoreIndices(pointIS, &points); PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&pointIS); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _zeroBoundary


// ----------------------------------------------------------------------
// Constructor
pylith::materials::TestMaterialNew_Data::TestMaterialNew_Data(void) :
    dimension(0),
    meshFilename(0),
    materialLabel(NULL),
    materialId(0),
    boundaryLabel(NULL),

    normalizer(new spatialdata::units::Nondimensional),

    t(0.0),
    dt(0.0),
    tshift(0.0),

    numSolnFields(0),
    solnDiscretizations(NULL),
    solnDBFilename(NULL),
    pertDBFilename(NULL),

    numAuxFields(0),
    auxFields(NULL),
    auxDiscretizations(NULL),
    auxDBFilename(NULL),

    isExplicit(false)
{ // constructor
    CPPUNIT_ASSERT(normalizer);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::materials::TestMaterialNew_Data::~TestMaterialNew_Data(void)
{ // destructor
    delete normalizer; normalizer = NULL;
} // destructor


// End of file
