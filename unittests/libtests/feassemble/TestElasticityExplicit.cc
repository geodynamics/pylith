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
#include "data/IntegratorData.hh" // USES IntegratorData

#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

#include <math.h> // USES fabs()

#include <stdexcept> // TEMPORARY

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicit );

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
// Test StableTimeStep().
void
pylith::feassemble::TestElasticityExplicit::testStableTimeStep(void)
{ // testStableTimeStep
  ElasticityExplicit integrator;

  const double dt1 = 2.0;
  integrator.timeStep(dt1);
  const double stableTimeStep = integrator.stableTimeStep();
  CPPUNIT_ASSERT_EQUAL(dt1, stableTimeStep);
} // testStableTimeStep

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
} // testMaterial

// ----------------------------------------------------------------------
// Test updateState().
void 
pylith::feassemble::TestElasticityExplicit::testUpdateState(void)
{ // testUpdateState
  CPPUNIT_ASSERT(0 != _data);

  ALE::Obj<ALE::Mesh> mesh;
  ElasticityExplicit integrator;
  topology::FieldsManager fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  const ALE::Obj<real_section_type>& dispT = fields.getReal("dispT");
  CPPUNIT_ASSERT(!dispT.isNull());
  integrator.updateState(dispT, mesh);
} // testUpdateState

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::feassemble::TestElasticityExplicit::testIntegrateResidual(void)
{ // testIntegrateResidual
  CPPUNIT_ASSERT(0 != _data);

  ALE::Obj<ALE::Mesh> mesh;
  ElasticityExplicit integrator;
  topology::FieldsManager fields(mesh);
  _initialize(&mesh, &integrator, &fields);

  const ALE::Obj<real_section_type>& residual = fields.getReal("residual");
  CPPUNIT_ASSERT(!residual.isNull());
  integrator.integrateResidual(residual, &fields, mesh);

  const double* valsE = _data->valsResidual;
  const int sizeE = _data->spaceDim * _data->numVertices;

  const double* vals = residual->restrict();
  const int size = residual->sizeWithBC();
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  const double tolerance = 1.0e-06;
  for (int i=0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i], vals[i], tolerance);
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::feassemble::TestElasticityExplicit::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  CPPUNIT_ASSERT(0 != _data);

  ALE::Obj<ALE::Mesh> mesh;
  ElasticityExplicit integrator;
  topology::FieldsManager fields(mesh);
  _initialize(&mesh, &integrator, &fields);
  mesh->getFactory()->clear();

  const ALE::Obj<pylith::real_section_type>& dispTpdt = 
    fields.getReal("dispTpdt");
  CPPUNIT_ASSERT(!dispTpdt.isNull());

  PetscMat jacobian;
  PetscErrorCode err = MeshCreateMatrix(mesh, dispTpdt, MATMPIBAIJ, &jacobian);
  CPPUNIT_ASSERT(0 == err);

  integrator.integrateJacobian(&jacobian, &fields, mesh);
  err = MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY);
  CPPUNIT_ASSERT(0 == err);
  err = MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY);
  CPPUNIT_ASSERT(0 == err);

  const double* valsE = _data->valsJacobian;
  const int nrowsE = _data->numVertices * _data->spaceDim;
  const int ncolsE = _data->numVertices * _data->spaceDim;

  int nrows = 0;
  int ncols = 0;
  MatGetSize(jacobian, &nrows, &ncols);
  CPPUNIT_ASSERT_EQUAL(nrowsE, nrows);
  CPPUNIT_ASSERT_EQUAL(ncolsE, ncols);

  PetscMat jDense;
  PetscMat jSparseAIJ;
  MatConvert(jacobian, MATSEQAIJ, MAT_INITIAL_MATRIX, &jSparseAIJ);
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
// Initialize elasticity integrator.
void
pylith::feassemble::TestElasticityExplicit::_initialize(
					 ALE::Obj<ALE::Mesh>* mesh,
					 ElasticityExplicit* const integrator,
					 topology::FieldsManager* fields)
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != integrator);
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != _quadrature);
  CPPUNIT_ASSERT(0 != _material);

  // Setup mesh
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(_data->spaceDim);
  cs.initialize();
  *mesh = new ALE::Mesh(PETSC_COMM_WORLD, _data->cellDim);
  CPPUNIT_ASSERT(!mesh->isNull());
  ALE::Obj<sieve_type> sieve = new sieve_type((*mesh)->comm());
  CPPUNIT_ASSERT(!sieve.isNull());
  const bool interpolate = false;
  ALE::SieveBuilder<Mesh>::buildTopology(sieve, _data->cellDim, 
	       _data->numCells, const_cast<int*>(_data->cells), 
	       _data->numVertices, interpolate, _data->numBasis);
  (*mesh)->setSieve(sieve);
  (*mesh)->stratify();
  ALE::SieveBuilder<Mesh>::buildCoordinates((*mesh), _data->spaceDim,
					    _data->vertices);
  const ALE::Obj<Mesh::label_type>& labelMaterials = 
    (*mesh)->createLabel("material-id");  
  int i = 0;
  const ALE::Obj<Mesh::label_sequence>& cells = (*mesh)->heightStratum(0);
  CPPUNIT_ASSERT(!cells.isNull());
  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter)
    (*mesh)->setValue(labelMaterials, *c_iter, _data->matId);

  // Setup quadrature
  _quadrature->initialize(_data->basis, _data->basisDeriv, _data->quadPts,
			  _data->quadWts, _data->cellDim, _data->numBasis,
			  _data->numQuadPts, _data->spaceDim);

  // Setup material
  spatialdata::spatialdb::SimpleIOAscii iohandler;
  iohandler.filename(_data->matDBFilename);
  spatialdata::spatialdb::SimpleDB db;
  db.ioHandler(&iohandler);
  
  _material->id(_data->matId);
  _material->label(_data->matLabel);
  _material->timestep(_data->dt);
  _material->db(&db);
  _material->initialize(*mesh, &cs, _quadrature);

  integrator->quadrature(_quadrature);
  integrator->timeStep(_data->dt);
  integrator->material(_material);

  // Setup fields
  CPPUNIT_ASSERT(0 != fields);
  fields->addReal("residual");
  fields->addReal("dispTpdt");
  fields->addReal("dispT");
  fields->addReal("dispTmdt");
  const char* history[] = { "dispTpdt", "dispT", "dispTmdt" };
  const int historySize = 3;
  fields->createHistory(history, historySize);
  
  const ALE::Obj<real_section_type>& residual = fields->getReal("residual");
  CPPUNIT_ASSERT(!residual.isNull());
  residual->setFiberDimension((*mesh)->depthStratum(0), _data->spaceDim);
  (*mesh)->allocate(residual);
  residual->zero();
  fields->copyLayout("residual");

  const int fieldSize = _data->spaceDim * _data->numVertices;
  const ALE::Obj<real_section_type>& dispTpdt = fields->getReal("dispTpdt");
  const ALE::Obj<real_section_type>& dispT = fields->getReal("dispT");
  const ALE::Obj<real_section_type>& dispTmdt = fields->getReal("dispTmdt");
  CPPUNIT_ASSERT(!dispTpdt.isNull());
  CPPUNIT_ASSERT(!dispT.isNull());
  CPPUNIT_ASSERT(!dispTmdt.isNull());
  const int offset = _data->numCells;
  for (int iVertex=0; iVertex < _data->numVertices; ++iVertex) {
    dispTpdt->updatePoint(iVertex+offset, 
			  &_data->fieldTpdt[iVertex*_data->spaceDim]);
    dispT->updatePoint(iVertex+offset, 
		       &_data->fieldT[iVertex*_data->spaceDim]);
    dispTmdt->updatePoint(iVertex+offset, 
			  &_data->fieldTmdt[iVertex*_data->spaceDim]);
  } // for
} // _initialize


// End of file 
