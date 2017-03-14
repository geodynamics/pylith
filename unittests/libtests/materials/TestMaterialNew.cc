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

#include "TestMaterialNew_Data.hh" // USES IsotropicLinearElasticityPlaneStrainData

#include "pylith/materials/MaterialNew.hh" // USES MaterialNew
#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "journal/debug.h" // USES journal::debug_t

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestMaterialNew::setUp(void)
{ // setUp
    _mesh = new topology::Mesh();
    _solution1 = NULL;
    _solution2 = NULL;
    _solution1Dot = NULL;
    _solution2Dot = NULL;
    _auxDB = NULL;
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::materials::TestMaterialNew::tearDown(void)
{ // tearDown
    delete _solution1; _solution1 = NULL;
    delete _solution2; _solution2 = NULL;
    delete _solution1Dot; _solution1Dot = NULL;
    delete _solution2Dot; _solution2Dot = NULL;
    delete _mesh; _mesh = NULL;
    delete _auxDB; _auxDB = NULL;
} // tearDown


// ----------------------------------------------------------------------
// Test _setFEKernelsRHSResidual().
void
pylith::materials::TestMaterialNew::test_setFEKernelsRHSResidual(void)
{ // test_setFEKernelsRHSResidual
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    CPPUNIT_ASSERT(_solution1);
    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

    material->_setFEKernelsRHSResidual(*_solution1);

    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(_solution1->dmMesh(), &prob); CPPUNIT_ASSERT(!err); CPPUNIT_ASSERT(prob);

    const int numSolnFields = data->numSolnFields;
    const int numFieldKernels = data->numKernelsResidual;
    for (int iField=0; iField < numSolnFields; ++iField) {
        const int indexK = iField*numFieldKernels;

        const PetscPointFunc* kernelsE = data->kernelsRHSResidual; CPPUNIT_ASSERT(kernelsE);
        PetscPointFunc g0 = NULL;
        PetscPointFunc g1 = NULL;
        err = PetscDSGetResidual(prob, iField, &g0, &g1); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT(kernelsE[indexK+0] == g0);
        CPPUNIT_ASSERT(kernelsE[indexK+1] == g1);
    } // for

    PYLITH_METHOD_END;
} // test_setFEKernelsRHSResidual


// ----------------------------------------------------------------------
// Test _setFEKernelsRHSJacobian().
void
pylith::materials::TestMaterialNew::test_setFEKernelsRHSJacobian(void)
{ // test_setFEKernelsRHSJacobian
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    CPPUNIT_ASSERT(_solution1);
    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

    material->_setFEKernelsRHSJacobian(*_solution1);

    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(_solution1->dmMesh(), &prob); CPPUNIT_ASSERT(!err); CPPUNIT_ASSERT(prob);

    const int numSolnFields = data->numSolnFields;
    const int numFieldKernels = data->numKernelsJacobian;
    for (int iField=0; iField < numSolnFields; ++iField) {
        for (int jField=0; jField < numSolnFields; ++jField) {
            const int indexK = iField*numSolnFields*numFieldKernels + jField*numFieldKernels;

            const PetscPointJac* kernelsE = data->kernelsRHSJacobian; CPPUNIT_ASSERT(kernelsE);
            PetscPointJac Jg0 = NULL;
            PetscPointJac Jg1 = NULL;
            PetscPointJac Jg2 = NULL;
            PetscPointJac Jg3 = NULL;
            err = PetscDSGetJacobian(prob, iField, jField, &Jg0, &Jg1, &Jg2, &Jg3); CPPUNIT_ASSERT(!err);
            CPPUNIT_ASSERT(kernelsE[indexK+0] == Jg0);
            CPPUNIT_ASSERT(kernelsE[indexK+1] == Jg1);
            CPPUNIT_ASSERT(kernelsE[indexK+2] == Jg2);
            CPPUNIT_ASSERT(kernelsE[indexK+3] == Jg3);
        } // for
    } // for

    PYLITH_METHOD_END;
} // test_setFEKernelsRHSJacobian


// ----------------------------------------------------------------------
// Test _setFEKernelsLHSResidual().
void
pylith::materials::TestMaterialNew::test_setFEKernelsLHSResidual(void)
{ // test_setFEKernelsLHSResidual
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    CPPUNIT_ASSERT(_solution1);
    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

    material->_setFEKernelsLHSResidual(*_solution1);

    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(_solution1->dmMesh(), &prob); CPPUNIT_ASSERT(!err); CPPUNIT_ASSERT(prob);

    const int numSolnFields = data->numSolnFields;
    const int numFieldKernels = data->numKernelsResidual;
    for (int iField=0; iField < numSolnFields; ++iField) {
        const int indexK = iField*numFieldKernels;

        const PetscPointFunc* kernelsE = data->kernelsLHSResidual; CPPUNIT_ASSERT(kernelsE);
        PetscPointFunc f0 = NULL;
        PetscPointFunc f1 = NULL;
        err = PetscDSGetResidual(prob, iField, &f0, &f1); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT(kernelsE[indexK+0] == f0);
        CPPUNIT_ASSERT(kernelsE[indexK+1] == f1);
    } // for

    PYLITH_METHOD_END;
} // test_setFEKernelsLHSResidual


// ----------------------------------------------------------------------
// Test _setFEKernelsLHSJacobianImplicit().
void
pylith::materials::TestMaterialNew::test_setFEKernelsLHSJacobianImplicit(void)
{ // test_setFEKernelsLHSJacobianImplicit
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    CPPUNIT_ASSERT(_solution1);
    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

    material->_setFEKernelsLHSJacobianImplicit(*_solution1);

    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(_solution1->dmMesh(), &prob); CPPUNIT_ASSERT(!err); CPPUNIT_ASSERT(prob);

    const int numSolnFields = data->numSolnFields;
    const int numFieldKernels = data->numKernelsJacobian;
    for (int iField=0; iField < numSolnFields; ++iField) {
        for (int jField=0; jField < numSolnFields; ++jField) {
            const int indexK = iField*numSolnFields*numFieldKernels + jField*numFieldKernels;

            const PetscPointJac* kernelsE = data->kernelsLHSJacobianImplicit; CPPUNIT_ASSERT(kernelsE);
            PetscPointJac Jf0 = NULL;
            PetscPointJac Jf1 = NULL;
            PetscPointJac Jf2 = NULL;
            PetscPointJac Jf3 = NULL;
            err = PetscDSGetJacobian(prob, iField, jField, &Jf0, &Jf1, &Jf2, &Jf3); CPPUNIT_ASSERT(!err);
            CPPUNIT_ASSERT(kernelsE[indexK+0] == Jf0);
            CPPUNIT_ASSERT(kernelsE[indexK+1] == Jf1);
            CPPUNIT_ASSERT(kernelsE[indexK+2] == Jf2);
            CPPUNIT_ASSERT(kernelsE[indexK+3] == Jf3);
        } // for
    } // for

    PYLITH_METHOD_END;
} // test_setFEKernelsLHSJacobianImplicit


// ----------------------------------------------------------------------
// Test _setFEKernelsLHSJacobianExplicit().
void
pylith::materials::TestMaterialNew::test_setFEKernelsLHSJacobianExplicit(void)
{ // test_setFEKernelsLHSJacobianExplicit
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    CPPUNIT_ASSERT(_solution1);
    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

    material->_setFEKernelsLHSJacobianExplicit(*_solution1);

    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(_solution1->dmMesh(), &prob); CPPUNIT_ASSERT(!err); CPPUNIT_ASSERT(prob);

    const int numSolnFields = data->numSolnFields;
    const int numFieldKernels = data->numKernelsJacobian;
    for (int iField=0; iField < numSolnFields; ++iField) {
        for (int jField=0; jField < numSolnFields; ++jField) {
            const int indexK = iField*numSolnFields*numFieldKernels + jField*numFieldKernels;

            const PetscPointJac* kernelsE = data->kernelsLHSJacobianExplicit; CPPUNIT_ASSERT(kernelsE);
            PetscPointJac Jf0 = NULL;
            PetscPointJac Jf1 = NULL;
            PetscPointJac Jf2 = NULL;
            PetscPointJac Jf3 = NULL;
            err = PetscDSGetJacobian(prob, iField, jField, &Jf0, &Jf1, &Jf2, &Jf3); CPPUNIT_ASSERT(!err);
            CPPUNIT_ASSERT(kernelsE[indexK+0] == Jf0);
            CPPUNIT_ASSERT(kernelsE[indexK+1] == Jf1);
            CPPUNIT_ASSERT(kernelsE[indexK+2] == Jf2);
            CPPUNIT_ASSERT(kernelsE[indexK+3] == Jf3);
        } // for
    } // for

    PYLITH_METHOD_END;
} // test_setFEKernelsLHSJacobianExplicit


// ----------------------------------------------------------------------
// Test hasAuxField().
void
pylith::materials::TestMaterialNew::testHasAuxField(void)
{ // testHasAuxField
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

    for (int i=0; i < data->numAuxFields; ++i) {
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

    const topology::FieldBase::DiscretizeInfo infoDefault = {-1, -1, true};
    const topology::FieldBase::DiscretizeInfo infoA = {1, 2, false};
    const topology::FieldBase::DiscretizeInfo infoB = {2, 2, true};

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    material->auxFieldDiscretization("A", infoA.basisOrder, infoA.quadOrder, infoA.isBasisContinuous);
    material->auxFieldDiscretization("B", infoB.basisOrder, infoB.quadOrder, infoB.isBasisContinuous);

    { // A
        const topology::FieldBase::DiscretizeInfo& test = material->auxFieldDiscretization("A");
        CPPUNIT_ASSERT_EQUAL(test.basisOrder, infoA.basisOrder);
        CPPUNIT_ASSERT_EQUAL(test.quadOrder, infoA.quadOrder);
        CPPUNIT_ASSERT_EQUAL(test.isBasisContinuous, infoA.isBasisContinuous);
    } // A

    { // B
        const topology::FieldBase::DiscretizeInfo& test = material->auxFieldDiscretization("B");
        CPPUNIT_ASSERT_EQUAL(test.basisOrder, infoB.basisOrder);
        CPPUNIT_ASSERT_EQUAL(test.quadOrder, infoB.quadOrder);
        CPPUNIT_ASSERT_EQUAL(test.isBasisContinuous, infoB.isBasisContinuous);
    } // B

    { // C (default)
        const topology::FieldBase::DiscretizeInfo& test = material->auxFieldDiscretization("C");
        CPPUNIT_ASSERT_EQUAL(test.basisOrder, infoDefault.basisOrder);
        CPPUNIT_ASSERT_EQUAL(test.quadOrder, infoDefault.quadOrder);
        CPPUNIT_ASSERT_EQUAL(test.isBasisContinuous, infoDefault.isBasisContinuous);
    } // C (default)

    { // default
        const topology::FieldBase::DiscretizeInfo& test = material->auxFieldDiscretization("default");
        CPPUNIT_ASSERT_EQUAL(test.basisOrder, infoDefault.basisOrder);
        CPPUNIT_ASSERT_EQUAL(test.quadOrder, infoDefault.quadOrder);
        CPPUNIT_ASSERT_EQUAL(test.isBasisContinuous, infoDefault.isBasisContinuous);
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
    spatialdata::spatialdb::SimpleDB db;
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

    _initializeMin();

    // Call verifyConfiguration()
    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    CPPUNIT_ASSERT(_mesh);
    material->verifyConfiguration(*_mesh);

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
    query->openDB(_auxDB, data->lengthScale);

    PetscErrorCode err = DMComputeL2Diff(dm, t, query->functions(), (void**)query->contextPtrs(), auxFields.globalVector(), &norm); CPPUNIT_ASSERT(!err);
    query->closeDB(_auxDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testInitialize


// ----------------------------------------------------------------------
// Test computeRHSResidual().
void
pylith::materials::TestMaterialNew::testComputeResidual(void)
{ // testComputeResidual
    PYLITH_METHOD_BEGIN;

    // Call initialize()
    _initializeFull(); // includes setting up auxFields

    CPPUNIT_ASSERT(_mesh);
    pylith::topology::Field residualRHS(*_mesh);
    residualRHS.cloneSection(*_solution1);
    residualRHS.label("residual RHS");
    residualRHS.allocate();
    residualRHS.zeroLocal();

    pylith::topology::Field residualLHS(*_mesh);
    residualLHS.cloneSection(*_solution1);
    residualLHS.label("residual LHS");
    residualLHS.allocate();
    residualLHS.zeroLocal();

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

    const PylithReal t = data->t;
    const PylithReal dt = data->dt;
    material->computeRHSResidual(residualRHS.localVector(), t, dt, *_solution1);
    material->computeLHSResidual(residualLHS.localVector(), t, dt, *_solution1, _solution1Dot->localVector());

    _zeroBoundary(&residualRHS);
    _zeroBoundary(&residualLHS);

    // Scatter local to global.
    residualRHS.complete();
    residualLHS.complete();

    PetscVec residualVec = NULL;
    PetscErrorCode err;
    err = VecDuplicate(residualRHS.globalVector(), &residualVec); CPPUNIT_ASSERT(!err);
    err = VecZeroEntries(residualVec); CPPUNIT_ASSERT(!err);
    err = VecWAXPY(residualVec, -1.0, residualLHS.globalVector(), residualRHS.globalVector()); CPPUNIT_ASSERT(!err);

    //residualRHS.view("RESIDUAL RHS"); // DEBUGGING
    //residualLHS.view("RESIDUAL LHS"); // DEBUGGING

    PylithReal norm = 0.0;
    err = VecNorm(residualVec, NORM_2, &norm); CPPUNIT_ASSERT(!err);
    err = VecDestroy(&residualVec); CPPUNIT_ASSERT(!err);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
    CPPUNIT_ASSERT(norm > 0.0); // Norm exactly equal to zero almost certainly means test is satisfied trivially.

    PYLITH_METHOD_END;
} // testComputeResidual


// ----------------------------------------------------------------------
// Test computeRHSJacobian().
void
pylith::materials::TestMaterialNew::testComputeRHSJacobian(void)
{ // testComputeRHSJacobian
    PYLITH_METHOD_BEGIN;

    // Create linear problem (MMS) with two solutions, s_1 and s_2.
    //
    // Check that Jg(s_1)*(s_2 - s_1) = G(s_2) - G(s_1).

    // Call initialize()
    _initializeFull(); // includes setting up auxFields

    CPPUNIT_ASSERT(_mesh);
    pylith::topology::Field residual1(*_mesh);
    residual1.cloneSection(*_solution1);
    residual1.label("residual1");
    residual1.allocate();
    residual1.zeroLocal();

    pylith::topology::Field residual2(*_mesh);
    residual2.cloneSection(*_solution2);
    residual2.label("residual2");
    residual2.allocate();
    residual2.zeroLocal();

#if 0 // DEBUGGING
    PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2");
    DMSetFromOptions(_solution1->dmMesh());
#endif

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

    const PylithReal t = data->t;
    const PylithReal dt = data->dt;
    material->computeRHSResidual(residual1.localVector(), t, dt, *_solution1);
    material->computeRHSResidual(residual2.localVector(), t, dt, *_solution2);

    // Scatter local to global.
    _solution1->createScatter(_solution1->mesh());
    _solution2->createScatter(_solution2->mesh());
    _solution1->scatterLocalToContext();
    _solution2->scatterLocalToContext();
    residual1.complete();
    residual2.complete();

    //residual1.view("RESIDUAL 1 RHS"); // DEBUGGING
    //residual2.view("RESIDUAL 1 RHS"); // DEBUGGING

    // Check that J(s_1)*(s_2 - s_1) = G(s_2) - G(s_1).

    PetscVec residualVec = NULL;
    PetscErrorCode err;
    err = VecDuplicate(residual1.globalVector(), &residualVec); CPPUNIT_ASSERT(!err);
    err = VecZeroEntries(residualVec); CPPUNIT_ASSERT(!err);
    err = VecWAXPY(residualVec, -1.0, residual1.globalVector(), residual2.globalVector()); CPPUNIT_ASSERT(!err);

    PetscVec solnIncrVec = NULL;
    err = VecDuplicate(_solution1->globalVector(), &solnIncrVec); CPPUNIT_ASSERT(!err);
    err = VecZeroEntries(solnIncrVec); CPPUNIT_ASSERT(!err);
    err = VecWAXPY(solnIncrVec, -1.0, _solution1->globalVector(), _solution2->globalVector()); CPPUNIT_ASSERT(!err);

    // Compute Jacobian
    PetscMat jacobianMat = NULL;
    err = DMCreateMatrix(_solution1->dmMesh(), &jacobianMat); CPPUNIT_ASSERT(!err);
    err = MatZeroEntries(jacobianMat); CPPUNIT_ASSERT(!err);
    PetscMat precondMat = jacobianMat; // Use Jacobian == preconditioner

    material->computeRHSJacobian(jacobianMat, precondMat, t, dt, *_solution1);
    CPPUNIT_ASSERT_EQUAL(false, material->needNewRHSJacobian());
    err = MatAssemblyBegin(jacobianMat, MAT_FINAL_ASSEMBLY); PYLITH_CHECK_ERROR(err);
    err = MatAssemblyEnd(jacobianMat, MAT_FINAL_ASSEMBLY); PYLITH_CHECK_ERROR(err);

    // result = Jg*(-solnIncr) + residual
    PetscVec resultVec = NULL;
    err = VecDuplicate(_solution1->globalVector(), &resultVec); CPPUNIT_ASSERT(!err);
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

    // Create linear problem (MMS) with two solutions, s_1 and s_2.
    //
    // Check that Jf(s_1)*(s_2 - s_1) = F(s_2) - F(s_1).

    // Call initialize()
    _initializeFull(); // includes setting up auxFields

    CPPUNIT_ASSERT(_mesh);
    pylith::topology::Field residual1(*_mesh);
    residual1.cloneSection(*_solution1);
    residual1.label("residual1");
    residual1.allocate();
    residual1.zeroLocal();

    pylith::topology::Field residual2(*_mesh);
    residual2.cloneSection(*_solution2);
    residual2.label("residual2");
    residual2.allocate();
    residual2.zeroLocal();

#if 0 // DEBUGGING
    PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2");
    DMSetFromOptions(_solution1->dmMesh());
#endif

    MaterialNew* material = _material(); CPPUNIT_ASSERT(material);
    TestMaterialNew_Data* data = _data(); CPPUNIT_ASSERT(data);

    const PylithReal t = data->t;
    const PylithReal dt = data->dt;
    const PylithReal tshift = data->tshift;
    material->computeLHSResidual(residual1.localVector(), t, dt, *_solution1, _solution1Dot->localVector());
    material->computeLHSResidual(residual2.localVector(), t, dt, *_solution2, _solution2Dot->localVector());

    // Scatter local to global.
    _solution1->createScatter(_solution1->mesh());
    _solution2->createScatter(_solution2->mesh());
    _solution1->scatterLocalToContext();
    _solution2->scatterLocalToContext();
    residual1.complete();
    residual2.complete();

    // Check that Jf(s_1)*(s_2 - s_1) = F(s_2) - F(s_1).

    //residual1.view("RESIDUAL 1 LHS"); // DEBUGGING
    //residual2.view("RESIDUAL 2 LHS"); // DEBUGGING

    PetscVec residualVec = NULL;
    PetscErrorCode err;
    err = VecDuplicate(residual1.globalVector(), &residualVec); CPPUNIT_ASSERT(!err);
    err = VecZeroEntries(residualVec); CPPUNIT_ASSERT(!err);
    err = VecWAXPY(residualVec, -1.0, residual1.globalVector(), residual2.globalVector()); CPPUNIT_ASSERT(!err);

    PetscVec solnIncrVec = NULL;
    err = VecDuplicate(_solution1->globalVector(), &solnIncrVec); CPPUNIT_ASSERT(!err);
    err = VecZeroEntries(solnIncrVec); CPPUNIT_ASSERT(!err);
    err = VecWAXPY(solnIncrVec, -1.0, _solution1->globalVector(), _solution2->globalVector()); CPPUNIT_ASSERT(!err);

    // Compute Jacobian
    PetscMat jacobianMat = NULL;
    err = DMCreateMatrix(_solution1->dmMesh(), &jacobianMat); CPPUNIT_ASSERT(!err);
    err = MatZeroEntries(jacobianMat); CPPUNIT_ASSERT(!err);
    PetscMat precondMat = jacobianMat; // Use Jacobian == preconditioner

    material->computeLHSJacobianImplicit(jacobianMat, precondMat, t, dt, tshift, *_solution1, _solution1Dot->localVector());
    CPPUNIT_ASSERT_EQUAL(false, material->needNewLHSJacobian());
    err = MatAssemblyBegin(jacobianMat, MAT_FINAL_ASSEMBLY); PYLITH_CHECK_ERROR(err);
    err = MatAssemblyEnd(jacobianMat, MAT_FINAL_ASSEMBLY); PYLITH_CHECK_ERROR(err);

    // result = J*(-solnIncr) + residual
    PetscVec resultVec = NULL;
    err = VecDuplicate(_solution1->globalVector(), &resultVec); CPPUNIT_ASSERT(!err);
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
} // testComputeLHSJacobianImplicit


// ----------------------------------------------------------------------
// Get material.
pylith::materials::MaterialNew*
pylith::materials::TestMaterialNew::_material(void)
{ // material
    throw std::logic_error("Implement TestMaterialNew::material in derived class.");
} // material


// ----------------------------------------------------------------------
// Get test data.
pylith::materials::TestMaterialNew_Data*
pylith::materials::TestMaterialNew::_data(void)
{ // data
    throw std::logic_error("Implement TestMaterialNew::data in derived class.");
} // data


// ----------------------------------------------------------------------
// Test computeLHSJacobianExplicit().
void
pylith::materials::TestMaterialNew::testComputeLHSJacobianInverseExplicit(void)
{ // testComputeLHSJacobianInverseExplicit
    PYLITH_METHOD_BEGIN;

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

    meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(data->meshFilename);
    iohandler.filename(data->meshFilename);
    iohandler.read(_mesh); CPPUNIT_ASSERT(_mesh);

    // Setup coordinates.
    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_mesh->dimension());
    cs.initialize();
    _mesh->coordsys(&cs);

    // Setup scales.
    spatialdata::units::Nondimensional normalizer;
    normalizer.lengthScale(data->lengthScale);
    normalizer.pressureScale(data->pressureScale);
    normalizer.timeScale(data->timeScale);
    normalizer.densityScale(data->densityScale);
    topology::MeshOps::nondimensionalize(_mesh, normalizer);

    material->id(data->materialId);
    material->label(data->materialLabel);
    material->normalizer(normalizer);

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

    _initializeMin();

    // Set auxiliary fields spatial database.
    delete _auxDB; _auxDB = new spatialdata::spatialdb::SimpleDB; CPPUNIT_ASSERT(_auxDB);
    spatialdata::spatialdb::SimpleIOAscii dbIO;
    CPPUNIT_ASSERT(data->auxDBFilename);
    dbIO.filename(data->auxDBFilename);
    _auxDB->ioHandler(&dbIO);
    _auxDB->label("IsotropicLinearElasciticityPlaneStrain auxiliary fields database");
    _auxDB->queryType(spatialdata::spatialdb::SimpleDB::LINEAR);
    material->auxFieldsDB(_auxDB);

    for (int i=0; i < data->numAuxFields; ++i) {
        const pylith::topology::FieldBase::DiscretizeInfo& info = data->auxDiscretizations[i];
        material->auxFieldDiscretization(data->auxFields[i], info.basisOrder, info.quadOrder, info.isBasisContinuous);
    } // for


    // Create solution field 1.
    CPPUNIT_ASSERT(data->solnDBFilename);
    bool isClone = false;

    delete _solution1; _solution1 = new pylith::topology::Field(*_mesh);
    _solution1->label("solution1");
    _setupSolutionField(_solution1, data->solnDBFilename, isClone);

    // Create time derivative of solution field 1.
    delete _solution1Dot; _solution1Dot = new pylith::topology::Field(*_mesh);
    _solution1Dot->label("solution1_dot");
    _setupSolutionField(_solution1Dot, data->solnDBFilename);

    // Create solution field 2; solution 2 = solution 1 + perturbation
    CPPUNIT_ASSERT(data->pertDBFilename);
    isClone = true;

    delete _solution2; _solution2 = new pylith::topology::Field(*_mesh);
    _solution2->cloneSection(*_solution1);
    _solution2->label("solution2");
    _setupSolutionField(_solution2, data->pertDBFilename, isClone);

    // Create time derivative of solution field 2.
    delete _solution2Dot; _solution2Dot = new pylith::topology::Field(*_mesh);
    _solution2Dot->cloneSection(*_solution1Dot);
    _solution2Dot->label("solution2_dot");
    _setupSolutionField(_solution2Dot, data->pertDBFilename, isClone);

    CPPUNIT_ASSERT(_solution1);
    material->initialize(*_solution1);

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
        for (PylithInt i=0; i < dof; ++i) {
            fieldArray[off+i] = 0.0;
        } // for
    } // for

    err = ISRestoreIndices(pointIS, &points); PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&pointIS); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _zeroBoundary


// ----------------------------------------------------------------------
// Setup and populate solution field.
void
pylith::materials::TestMaterialNew::_setupSolutionField(pylith::topology::Field* field,
                                                        const char* dbFilename,
                                                        const bool isClone)
{ // _setupSolutionField
    throw std::logic_error("Implement TestMaterialNew::_setupSolutionField() in derived class.");
} // _setupSolutionField


// End of file
