// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestElasticityExplicit.hh" // Implementation of class methods

#include "pylith/feassemble/ElasticityExplicit.hh" // USES ElasticityExplicit
#include "data/ElasticityExplicitData.hh" // USES IntegratorData

#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES Mesh
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
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicit );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityExplicit::setUp(void)
{ // setUp
  _quadrature = new Quadrature<topology::Mesh>();
  _data = 0;
  _material = 0;
  _gravityField = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::feassemble::TestElasticityExplicit::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
  delete _material; _material = 0;
  delete _gravityField; _gravityField = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::feassemble::TestElasticityExplicit::testConstructor(void)
{ // testConstructor
  ElasticityExplicit integrator;
} // testConstructor

// ----------------------------------------------------------------------
// Test timeStep().
void
pylith::feassemble::TestElasticityExplicit::testTimeStep(void)
{ // testTimeStep
  ElasticityExplicit integrator;

  const double dt1 = 2.0;
  integrator.timeStep(dt1);
  CPPUNIT_ASSERT_EQUAL(dt1, integrator._dt);
  integrator.timeStep(dt1);
  CPPUNIT_ASSERT_EQUAL(dt1, integrator._dtm1);
  CPPUNIT_ASSERT_EQUAL(dt1, integrator._dt);
} // testTimeStep

// ----------------------------------------------------------------------
// Test material().
void
pylith::feassemble::TestElasticityExplicit::testMaterial(void)
{ // testMaterial
  ElasticityExplicit integrator;

  materials::ElasticIsotropic3D material;
  const int id = 3;
  const std::string label("my material");
  material.id(id);
  material.label(label.c_str());
  integrator.material(&material);
  CPPUNIT_ASSERT_EQUAL(id, integrator._material->id());
  CPPUNIT_ASSERT_EQUAL(label, std::string(integrator._material->label()));
  CPPUNIT_ASSERT_EQUAL(integrator._dt, integrator._material->timeStep());
  const double dt = 2.0;
  integrator.timeStep(dt);
  CPPUNIT_ASSERT_EQUAL(dt, integrator._material->timeStep());
} // testMaterial

// ----------------------------------------------------------------------
// Test needNewJacobian().
void
pylith::feassemble::TestElasticityExplicit::testNeedNewJacobian(void)
{ // testNeedNewJacobian
  ElasticityExplicit integrator;

  materials::ElasticIsotropic3D material;
  integrator.material(&material);
  CPPUNIT_ASSERT_EQUAL(true, integrator.needNewJacobian());
  integrator._needNewJacobian = false;
  CPPUNIT_ASSERT_EQUAL(false, integrator.needNewJacobian());  
} // testNeedNewJacobian

// ----------------------------------------------------------------------
// Test useSolnIncr().
void
pylith::feassemble::TestElasticityExplicit::testUseSolnIncr(void)
{ // testUseSolnIncr
  ElasticityExplicit integrator;

  materials::ElasticIsotropic3D material;
  integrator.material(&material);
  CPPUNIT_ASSERT_EQUAL(false, integrator._useSolnIncr);
  try {
    integrator.useSolnIncr(false);

    // Should have thrown exception, so don't make it here.
    CPPUNIT_ASSERT(false);
  } catch (const std::logic_error& err) {
    // Expect logic error so don't do anything.
  } catch (...) {
    CPPUNIT_ASSERT(false);
  } // try/catch
} // testUseSolnIncr

// ----------------------------------------------------------------------
// Test initialize().
void 
pylith::feassemble::TestElasticityExplicit::testInitialize(void)
{ // testInitialize
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  ElasticityExplicit integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);

} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::feassemble::TestElasticityExplicit::testIntegrateResidual(void)
{ // testIntegrateResidual
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  ElasticityExplicit integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  topology::Field<topology::Mesh>& residual = fields.get("residual");
  const double t = 1.0;
  integrator.integrateResidual(residual, t, &fields);

  const double* valsE = _data->valsResidual;
  const int sizeE = _data->spaceDim * _data->numVertices;

  const ALE::Obj<RealSection>& residualSection = residual.section();
  CPPUNIT_ASSERT(!residualSection.isNull());
  const double* vals = residualSection->restrictSpace();
  const int size = residualSection->sizeWithBC();
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

#if 0
  residual.view("RESIDUAL");
  std::cout << "EXPECTED RESIDUAL" << std::endl;
  for (int i=0; i < size; ++i)
    std::cout << "  " << valsE[i] << std::endl;
#endif

  const double tolerance = 1.0e-06;
  for (int i=0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i], vals[i], tolerance);
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::feassemble::TestElasticityExplicit::testIntegrateResidualLumped(void)
{ // testIntegrateResidualLumped
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  ElasticityExplicit integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  topology::Field<topology::Mesh>& residual = fields.get("residual");
  const double t = 1.0;
  integrator.integrateResidualLumped(residual, t, &fields);

  const double* valsE = _data->valsResidualLumped;
  const int sizeE = _data->spaceDim * _data->numVertices;

  const ALE::Obj<RealSection>& residualSection = residual.section();
  CPPUNIT_ASSERT(!residualSection.isNull());
  const double* vals = residualSection->restrictSpace();
  const int size = residualSection->sizeWithBC();
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

#if 0
  residual.view("RESIDUAL");
  std::cout << "EXPECTED RESIDUAL" << std::endl;
  for (int i=0; i < size; ++i)
    std::cout << "  " << valsE[i] << std::endl;
#endif

  const double tolerance = 1.0e-06;
  for (int i=0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i], vals[i], tolerance);
} // testIntegrateResidualLumped

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::feassemble::TestElasticityExplicit::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  ElasticityExplicit integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);
  integrator._needNewJacobian = true;

  topology::Jacobian jacobian(fields.solution());

  const double t = 1.0;
  integrator.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, integrator.needNewJacobian());
  jacobian.assemble("final_assembly");

  const double* valsE = _data->valsJacobian;
  const int nrowsE = _data->numVertices * _data->spaceDim;
  const int ncolsE = _data->numVertices * _data->spaceDim;

  const PetscMat jacobianMat = jacobian.matrix();

  int nrows = 0;
  int ncols = 0;
  MatGetSize(jacobianMat, &nrows, &ncols);
  CPPUNIT_ASSERT_EQUAL(nrowsE, nrows);
  CPPUNIT_ASSERT_EQUAL(ncolsE, ncols);

  PetscMat jDense;
  PetscMat jSparseAIJ;
  MatConvert(jacobianMat, MATSEQAIJ, MAT_INITIAL_MATRIX, &jSparseAIJ);
  MatConvert(jSparseAIJ, MATSEQDENSE, MAT_INITIAL_MATRIX, &jDense);

  double_array vals(nrows*ncols);
  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  MatGetValues(jDense, nrows, &rows[0], ncols, &cols[0], &vals[0]);
  const double tolerance = 1.0e-06;
  for (int iRow=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol) {
      const int index = ncols*iRow+iCol;
      if (fabs(valsE[index]) > 1.0)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[index]/valsE[index], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[index], vals[index], tolerance);
    } // for
  MatDestroy(jDense);
  MatDestroy(jSparseAIJ);
} // testIntegrateJacobian

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::feassemble::TestElasticityExplicit::testIntegrateJacobianLumped(void)
{ // testIntegrateJacobianLumped
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  ElasticityExplicit integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);
  integrator._needNewJacobian = true;

  topology::Field<topology::Mesh> jacobian(mesh);
  jacobian.label("Jacobian");
  jacobian.vectorFieldType(topology::FieldBase::VECTOR);
  jacobian.newSection(topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
  jacobian.allocate();

  const double t = 1.0;
  integrator.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, integrator.needNewJacobian());
  jacobian.complete();

  const double* valsE = _data->valsJacobianLumped;

#if 0 // DEBUGGING
  // TEMPORARY
  jacobian.view("JACOBIAN");
  std::cout << "\n\nJACOBIAN FULL" << std::endl;
  const int n = numBasis*spaceDim;
  for (int r=0; r < n; ++r) {
    for (int c=0; c < n; ++c) 
      std::cout << "  " << valsMatrixE[r*n+c];
    std::cout << "\n";
  } // for
#endif // DEBUGGING

  const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
  CPPUNIT_ASSERT(!jacobianSection.isNull());
  const double* vals = jacobianSection->restrictSpace();
  const int size = jacobianSection->sizeWithBC();
  const int sizeE = _data->numVertices * _data->spaceDim;
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  const double tolerance = 1.0e-06;
  for (int i=0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i], vals[i], tolerance);
} // testIntegrateJacobianLumped

// ----------------------------------------------------------------------
// Test updateStateVars().
void 
pylith::feassemble::TestElasticityExplicit::testUpdateStateVars(void)
{ // testUpdateStateVars
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  ElasticityExplicit integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  const double t = 1.0;
  integrator.updateStateVars(t, &fields);
} // testUpdateStateVars

// ----------------------------------------------------------------------
// Test StableTimeStep().
void
pylith::feassemble::TestElasticityExplicit::testStableTimeStep(void)
{ // testStableTimeStep
  topology::Mesh mesh;
  ElasticityExplicit integrator;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  const double stableTimeStep = integrator.stableTimeStep(mesh);
  CPPUNIT_ASSERT_EQUAL(1.0e+30, stableTimeStep);
} // testStableTimeStep

// ----------------------------------------------------------------------
// Initialize elasticity integrator.
void
pylith::feassemble::TestElasticityExplicit::_initialize(
					 topology::Mesh* mesh,
					 ElasticityExplicit* const integrator,
					 topology::SolutionFields* fields)
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != integrator);
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != _quadrature);
  CPPUNIT_ASSERT(0 != _material);

  const int spaceDim = _data->spaceDim;
  const double dt = _data->dt;

  // Setup mesh
  mesh->createSieveMesh(_data->cellDim);
  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  ALE::Obj<SieveMesh::sieve_type> sieve = 
    new SieveMesh::sieve_type(mesh->comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  // Cells and vertices
  const bool interpolate = false;
  ALE::Obj<ALE::Mesh::sieve_type> s = 
    new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());
  
  ALE::SieveBuilder<ALE::Mesh>::buildTopology(s, 
					      _data->cellDim, _data->numCells,
                                              const_cast<int*>(_data->cells), 
					      _data->numVertices,
                                              interpolate, _data->numBasis);
  std::map<ALE::Mesh::point_type,ALE::Mesh::point_type> renumbering;
  ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
  sieveMesh->setSieve(sieve);
  sieveMesh->stratify();
  ALE::SieveBuilder<SieveMesh>::buildCoordinates(sieveMesh, spaceDim, 
						 _data->vertices);

  // Material ids
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->heightStratum(0);
  CPPUNIT_ASSERT(!cells.isNull());
  const ALE::Obj<SieveMesh::label_type>& labelMaterials = 
    sieveMesh->createLabel("material-id");
  CPPUNIT_ASSERT(!labelMaterials.isNull());
  int i = 0;
  for(SieveMesh::label_sequence::iterator e_iter=cells->begin(); 
      e_iter != cells->end();
      ++e_iter)
    sieveMesh->setValue(labelMaterials, *e_iter, _data->matId);

  // Setup quadrature
  _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			  _data->basisDerivRef, _data->numQuadPts,
			  _data->numBasis, _data->cellDim,
			  _data->quadPts, _data->numQuadPts, _data->cellDim,
			  _data->quadWts, _data->numQuadPts,
			  spaceDim);

  spatialdata::units::Nondimensional normalizer;
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(spaceDim);
  cs.initialize();
  mesh->coordsys(&cs);
  mesh->nondimensionalize(normalizer);

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
  integrator->timeStep(_data->dt);
  integrator->material(_material);
  integrator->initialize(*mesh);

  // Setup fields
  CPPUNIT_ASSERT(0 != fields);
  fields->add("residual", "residual");
  fields->add("dispIncr(t->t+dt)", "displacement_increment");
  fields->add("disp(t)", "displacement");
  fields->add("disp(t-dt)", "displacement");
  fields->add("acceleration(t)", "acceleration");
  fields->solutionName("dispIncr(t->t+dt)");
  
  topology::Field<topology::Mesh>& residual = fields->get("residual");
  residual.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  residual.allocate();
  residual.zero();
  fields->copyLayout("residual");

  const int fieldSize = spaceDim * _data->numVertices;
  const ALE::Obj<RealSection>& dispIncrSection = 
    fields->get("dispIncr(t->t+dt)").section();
  const ALE::Obj<RealSection>& dispTSection = 
    fields->get("disp(t)").section();
  const ALE::Obj<RealSection>& dispTmdtSection = 
    fields->get("disp(t-dt)").section();
  const ALE::Obj<RealSection>& accSection = 
    fields->get("acceleration(t)").section();
  CPPUNIT_ASSERT(!dispIncrSection.isNull());
  CPPUNIT_ASSERT(!dispTSection.isNull());
  CPPUNIT_ASSERT(!dispTmdtSection.isNull());
  CPPUNIT_ASSERT(!accSection.isNull());

  double_array accVertex(spaceDim);


  const int offset = _data->numCells;
  for (int iVertex=0; iVertex < _data->numVertices; ++iVertex) {
    dispIncrSection->updatePoint(iVertex+offset, 
				 &_data->fieldTIncr[iVertex*spaceDim]);
    dispTSection->updatePoint(iVertex+offset, 
			      &_data->fieldT[iVertex*spaceDim]);
    dispTmdtSection->updatePoint(iVertex+offset, 
				 &_data->fieldTmdt[iVertex*spaceDim]);

    for (int iDim=0; iDim < spaceDim; ++iDim)
      accVertex[iDim] = (_data->fieldTIncr[iVertex*spaceDim+iDim] -
			 _data->fieldT[iVertex*spaceDim+iDim] +
			 _data->fieldTmdt[iVertex*spaceDim+iDim]) / (dt*dt);
    accSection->updatePoint(iVertex+offset, &accVertex[0]);

  } // for
} // _initialize


// End of file 
