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

#include "TestElasticityExplicitTri3.hh" // Implementation of class methods

#include "pylith/feassemble/ElasticityExplicitTri3.hh" // USES ElasticityExplicitTri3
#include "data/ElasticityExplicitData2DLinear.hh"
#include "pylith/feassemble/GeometryTri2D.hh" // USES GeometryTri2D

#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <math.h> // USES fabs()

#include <stdexcept> // USES std::exception

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitTri3::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _quadrature = new Quadrature();CPPUNIT_ASSERT(_quadrature);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  _data = new ElasticityExplicitData2DLinear();
  CPPUNIT_ASSERT(_data);
  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"), std::string(_data->matType));
  _gravityField = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::feassemble::TestElasticityExplicitTri3::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
  delete _material; _material = 0;
  delete _gravityField; _gravityField = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::feassemble::TestElasticityExplicitTri3::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  ElasticityExplicitTri3 integrator;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test timeStep().
void
pylith::feassemble::TestElasticityExplicitTri3::testTimeStep(void)
{ // testTimeStep
  PYLITH_METHOD_BEGIN;

  ElasticityExplicitTri3 integrator;

  const PylithScalar dt1 = 2.0;
  integrator.timeStep(dt1);
  CPPUNIT_ASSERT_EQUAL(dt1, integrator._dt);
  integrator.timeStep(dt1);
  CPPUNIT_ASSERT_EQUAL(dt1, integrator._dtm1);
  CPPUNIT_ASSERT_EQUAL(dt1, integrator._dt);

  PYLITH_METHOD_END;
} // testTimeStep

// ----------------------------------------------------------------------
// Test material().
void
pylith::feassemble::TestElasticityExplicitTri3::testMaterial(void)
{ // testMaterial
  PYLITH_METHOD_BEGIN;

  ElasticityExplicitTri3 integrator;

  materials::ElasticPlaneStrain material;
  const int id = 3;
  const std::string label("my material");
  material.id(id);
  material.label(label.c_str());
  integrator.material(&material);
  CPPUNIT_ASSERT_EQUAL(id, integrator._material->id());
  CPPUNIT_ASSERT_EQUAL(label, std::string(integrator._material->label()));
  CPPUNIT_ASSERT_EQUAL(integrator._dt, integrator._material->timeStep());
  const PylithScalar dt = 2.0;
  integrator.timeStep(dt);
  CPPUNIT_ASSERT_EQUAL(dt, integrator._material->timeStep());

  PYLITH_METHOD_END;
} // testMaterial

// ----------------------------------------------------------------------
// Test needNewJacobian().
void
pylith::feassemble::TestElasticityExplicitTri3::testNeedNewJacobian(void)
{ // testNeedNewJacobian
  PYLITH_METHOD_BEGIN;

  ElasticityExplicitTri3 integrator;

  materials::ElasticPlaneStrain material;
  integrator.material(&material);
  CPPUNIT_ASSERT_EQUAL(true, integrator.needNewJacobian());
  integrator._needNewJacobian = false;
  CPPUNIT_ASSERT_EQUAL(false, integrator.needNewJacobian());  

  PYLITH_METHOD_END;
} // testNeedNewJacobian

// ----------------------------------------------------------------------
// Test initialize().
void 
pylith::feassemble::TestElasticityExplicitTri3::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  ElasticityExplicitTri3 integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::feassemble::TestElasticityExplicitTri3::testIntegrateResidual(void)
{ // testIntegrateResidual
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  ElasticityExplicitTri3 integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  topology::Field& residual = fields.get("residual");
  const PylithScalar t = 1.0;
  integrator.integrateResidual(residual, t, &fields);

  const PylithScalar* valsE = _data->valsResidual;

#if 0 // DEBUGGING
  residual.view("RESIDUAL");
  std::cout << "EXPECTED RESIDUAL" << std::endl;
  const int size = _data->spaceDim * _data->numVertices;
  for (int i=0; i < size; ++i)
    std::cout << "  " << valsE[i] << std::endl;
#endif // DEBUGGING

  const PetscDM dmMesh = mesh.dmMesh();
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  CPPUNIT_ASSERT_EQUAL(_data->numVertices, verticesStratum.size());

  topology::VecVisitorMesh residualVisitor(residual);
  const PetscScalar* residualArray = residualVisitor.localArray();CPPUNIT_ASSERT(residualArray);

  const PylithScalar accScale = _data->lengthScale / pow(_data->timeScale, 2);
  const PylithScalar residualScale = _data->densityScale * accScale*pow(_data->lengthScale, _data->spaceDim);

  const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
  for (PetscInt v = vStart, index = 0; v < vEnd; ++v) {
    const PetscInt off = residualVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(_data->spaceDim, residualVisitor.sectionDof(v));

    for (int d=0; d < _data->spaceDim; ++d, ++index) {
      if (fabs(valsE[index]) > 1.0)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, residualArray[off+d]/valsE[index]*residualScale, tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[index], residualArray[off+d]*residualScale, tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::feassemble::TestElasticityExplicitTri3::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  ElasticityExplicitTri3 integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);
  integrator._needNewJacobian = true;

  const int spaceDim = _data->spaceDim;
  const PylithScalar lengthScale = _data->lengthScale;

  topology::Field jacobian(mesh);
  jacobian.label("Jacobian");
  jacobian.vectorFieldType(topology::FieldBase::VECTOR);
  jacobian.subfieldAdd("displacement", spaceDim, topology::Field::VECTOR, lengthScale);
  jacobian.subfieldAdd("lagrange_multiplier", spaceDim, topology::Field::VECTOR);

  jacobian.subfieldsSetup();
  jacobian.setupSolnChart();
  jacobian.setupSolnDof(spaceDim);
  jacobian.allocate();
  jacobian.zeroAll();

  const PylithScalar t = 1.0;
  integrator.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, integrator.needNewJacobian());
  jacobian.complete();

  const PylithScalar* valsE = _data->valsJacobian;
  const int numBasis = _data->numVertices;

#if 0 // DEBUGGING
  jacobian.view("JACOBIAN");
  std::cout << "\n\nJACOBIAN FULL" << std::endl;
  const int n = numBasis*spaceDim;
  for (int i=0; i < n; ++i)
    std::cout << "  " << valsE[i] << "\n";
#endif // DEBUGGING

  const PetscDM dmMesh = mesh.dmMesh();
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  CPPUNIT_ASSERT_EQUAL(_data->numVertices, verticesStratum.size());

  topology::VecVisitorMesh jacobianVisitor(jacobian);
  const PetscScalar* jacobianArray = jacobianVisitor.localArray();CPPUNIT_ASSERT(jacobianArray);

  const PylithScalar jacobianScale = _data->densityScale / pow(_data->timeScale, 2) * pow(_data->lengthScale, _data->spaceDim);

  const PylithScalar tolerance = (sizeof(double) == sizeof(PylithScalar)) ? 1.0e-06 : 1.0e-05;
  for (PetscInt v = vStart, index = 0; v < vEnd; ++v) {
    const PetscInt off = jacobianVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(_data->spaceDim, jacobianVisitor.sectionDof(v));

    for (int d=0; d < _data->spaceDim; ++d, ++index) {
      if (fabs(valsE[index]) > 1.0)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, jacobianArray[off+d]/valsE[index]*jacobianScale, tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[index], jacobianArray[off+d]*jacobianScale, tolerance);
    } // for
  } // for

  PYLITH_METHOD_END;
} // testIntegrateJacobian

// ----------------------------------------------------------------------
// Test updateStateVars().
void 
pylith::feassemble::TestElasticityExplicitTri3::testUpdateStateVars(void)
{ // testUpdateStateVars
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  ElasticityExplicitTri3 integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  const PylithScalar t = 1.0;
  integrator.updateStateVars(t, &fields);

  PYLITH_METHOD_END;
} // testUpdateStateVars

// ----------------------------------------------------------------------
// Test StableTimeStep().
void
pylith::feassemble::TestElasticityExplicitTri3::testStableTimeStep(void)
{ // testStableTimeStep
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  topology::Mesh mesh;
  ElasticityExplicitTri3 integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  const PylithScalar dtStable = integrator.stableTimeStep(mesh);
  const PylithScalar tolerance = 1.0e-6;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dtStable/_data->dtStableExplicit*_data->timeScale, tolerance);

  PYLITH_METHOD_END;
} // testStableTimeStep

// Initialize elasticity integrator.
void
pylith::feassemble::TestElasticityExplicitTri3::_initialize(topology::Mesh* mesh,
							    ElasticityExplicitTri3* const integrator,
							    topology::SolutionFields* fields)
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(integrator);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_quadrature);
  CPPUNIT_ASSERT(_material);

  const int spaceDim = _data->spaceDim;
  const PylithScalar dt = _data->dt;

  // Setup mesh
  PetscDM dmMesh;

  // Cells and vertices
  const PetscBool interpolate = PETSC_TRUE;
  PetscErrorCode err;
  err = DMPlexCreateFromCellList(PETSC_COMM_WORLD, _data->cellDim, _data->numCells, _data->numVertices, _data->numBasis, interpolate, _data->cells, _data->spaceDim, _data->vertices, &dmMesh);PYLITH_CHECK_ERROR(err);
  mesh->dmMesh(dmMesh, "domain");

  // Material ids
  PetscInt cStart, cEnd, c;
  err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    err = DMSetLabelValue(dmMesh, "material-id", c, _data->matId);PYLITH_CHECK_ERROR(err);
  } // for
 
  // Setup quadrature
  _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			  _data->basisDerivRef, _data->numQuadPts,
			  _data->numBasis, _data->cellDim,
			  _data->quadPts, _data->numQuadPts, _data->cellDim,
			  _data->quadWts, _data->numQuadPts,
			  spaceDim);

  // Setup coordinate system.
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);

  // Setup scales.
  const PylithScalar timeScale = _data->timeScale;
  const PylithScalar lengthScale = _data->lengthScale;
  const PylithScalar velScale = lengthScale / timeScale;
  const PylithScalar accScale = lengthScale / (timeScale*timeScale);
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.densityScale(_data->densityScale);
  normalizer.timeScale(_data->timeScale);
  topology::MeshOps::nondimensionalize(mesh, normalizer);

  // Setup material
  spatialdata::spatialdb::SimpleIOAscii iohandler;
  iohandler.filename(_data->matDBFilename);
  spatialdata::spatialdb::SimpleDB dbProperties;
  dbProperties.ioHandler(&iohandler);
  
  _material->id(_data->matId);
  _material->label(_data->matLabel);
  _material->dbProperties(&dbProperties);
  _material->normalizer(normalizer);

  integrator->quadrature(_quadrature);
  integrator->gravityField(_gravityField);
  integrator->timeStep(_data->dt / _data->timeScale);
  integrator->material(_material);
  integrator->initialize(*mesh);

  // Setup fields
  CPPUNIT_ASSERT(fields);
  fields->add("residual", "residual");
  fields->add("dispIncr(t->t+dt)", "displacement_increment");
  fields->add("disp(t)", "displacement");
  fields->add("disp(t-dt)", "displacement");
  fields->add("velocity(t)", "velocity");
  fields->add("acceleration(t)", "acceleration");
  fields->solutionName("dispIncr(t->t+dt)");
  
  topology::Field& residual = fields->get("residual");
  residual.subfieldAdd("displacement", spaceDim, topology::Field::VECTOR, lengthScale);
  residual.subfieldAdd("lagrange_multiplier", spaceDim, topology::Field::VECTOR);

  residual.subfieldsSetup();
  residual.setupSolnChart();
  residual.setupSolnDof(spaceDim);
  residual.allocate();
  residual.zeroAll();
  fields->copyLayout("residual");

  topology::VecVisitorMesh dispTVisitor(fields->get("disp(t)"));
  PetscScalar* dispTArray = dispTVisitor.localArray();CPPUNIT_ASSERT(dispTArray);

  topology::VecVisitorMesh dispTmdtVisitor(fields->get("disp(t-dt)"));
  PetscScalar* dispTmdtArray = dispTmdtVisitor.localArray();CPPUNIT_ASSERT(dispTmdtArray);

  topology::VecVisitorMesh dispTIncrVisitor(fields->get("dispIncr(t->t+dt)"));
  PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();CPPUNIT_ASSERT(dispTIncrArray);

  topology::VecVisitorMesh velVisitor(fields->get("velocity(t)"));
  PetscScalar* velArray = velVisitor.localArray();CPPUNIT_ASSERT(velArray);

  topology::VecVisitorMesh accVisitor(fields->get("acceleration(t)"));
  PetscScalar* accArray = accVisitor.localArray();CPPUNIT_ASSERT(accArray);

  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  
  for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
    const PetscInt dtoff = dispTVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dispTVisitor.sectionDof(v));

    const PetscInt dmoff = dispTmdtVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dispTmdtVisitor.sectionDof(v));

    const PetscInt dioff = dispTIncrVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, dispTIncrVisitor.sectionDof(v));

    const PetscInt voff = velVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, velVisitor.sectionDof(v));

    const PetscInt aoff = accVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(spaceDim, accVisitor.sectionDof(v));

    for(int iDim=0; iDim < spaceDim; ++iDim) {
      dispTArray[dtoff+iDim] = _data->fieldT[iVertex*spaceDim+iDim] / lengthScale;
      dispTmdtArray[dmoff+iDim] = _data->fieldTmdt[iVertex*spaceDim+iDim] / lengthScale;
      dispTIncrArray[dioff+iDim] = _data->fieldTIncr[iVertex*spaceDim+iDim] / lengthScale;

      velArray[voff+iDim] = (_data->fieldTIncr[iVertex*spaceDim+iDim] +
			     _data->fieldT[iVertex*spaceDim+iDim] -
			     _data->fieldTmdt[iVertex*spaceDim+iDim]) / (2.0*dt) / velScale;
      accArray[aoff+iDim] = (_data->fieldTIncr[iVertex*spaceDim+iDim] -
			     _data->fieldT[iVertex*spaceDim+iDim] +
			     _data->fieldTmdt[iVertex*spaceDim+iDim]) / (dt*dt) / accScale;
    } // for
  } // for

  PYLITH_METHOD_END;
} // _initialize


// End of file 
