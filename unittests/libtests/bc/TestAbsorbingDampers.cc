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

#include "TestAbsorbingDampers.hh" // Implementation of class methods

#include "pylith/bc/AbsorbingDampers.hh" // USES AbsorbingDampers

#include "data/AbsorbingDampersData.hh" // USES AbsorbingDampersData

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestAbsorbingDampers );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestAbsorbingDampers::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestAbsorbingDampers::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestAbsorbingDampers::testConstructor(void)
{ // testConstructor
  AbsorbingDampers bc;
} // testConstructor

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestAbsorbingDampers::testInitialize(void)
{ // testInitialize
  ALE::Obj<Mesh> mesh;
  AbsorbingDampers bc;
  topology::FieldsManager fields(mesh);
  _initialize(&mesh, &bc, &fields);

  CPPUNIT_ASSERT(0 != _data);
  
  const ALE::Obj<Mesh>& boundaryMesh = bc._boundaryMesh;

  // Check boundary mesh
  CPPUNIT_ASSERT(!boundaryMesh.isNull());

  const int cellDim = boundaryMesh->getDimension();
  const ALE::Obj<sieve_type>& sieve = boundaryMesh->getSieve();
  const ALE::Obj<Mesh::label_sequence>& cells = boundaryMesh->heightStratum(1);
  const int numVertices = boundaryMesh->depthStratum(0)->size();
  const int numCells = cells->size();

  CPPUNIT_ASSERT_EQUAL(_data->cellDim, cellDim);
  CPPUNIT_ASSERT_EQUAL(_data->numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(_data->numCells, numCells);

  //boundaryMesh->view("BOUNDARY MESH");

  const int boundaryDepth = boundaryMesh->depth()-1; // depth of bndry cells  
  int iCell = 0;
  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    const int numCorners = (boundaryMesh->getDimension() > 0) ?
      sieve->nCone(*c_iter, boundaryDepth)->size() : 1;
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numCorners);

    const ALE::Obj<sieve_type::traits::coneSequence>& cone = 
      sieve->cone(*c_iter);
    for(sieve_type::traits::coneSequence::iterator v_iter = cone->begin();
	v_iter != cone->end();
	++v_iter)
      CPPUNIT_ASSERT_EQUAL(_data->cells[iCell++], *v_iter);
  } // for

  // Check damping constants
  const int sizeE = _data->numCells * _data->numQuadPts * _data->spaceDim;
  const double* valsE = _data->dampingConsts;

  const int size = bc._dampingConsts->sizeWithBC();
  const double* vals = bc._dampingConsts->restrict();

  //bc._dampingConsts->view("DAMPING CONSTS");

  CPPUNIT_ASSERT_EQUAL(sizeE, size);
  const double tolerance = 1.0e-06;
  for (int i=0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i], vals[i], tolerance);
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestAbsorbingDampers::testIntegrateResidual(void)
{ // testIntegrateResidual
  CPPUNIT_ASSERT(0 != _data);

  ALE::Obj<ALE::Mesh> mesh;
  AbsorbingDampers bc;
  topology::FieldsManager fields(mesh);
  _initialize(&mesh, &bc, &fields);

  const ALE::Obj<real_section_type>& residual = fields.getReal("residual");
  CPPUNIT_ASSERT(!residual.isNull());
  const double t = 1.0;
  bc.integrateResidual(residual, t, &fields, mesh);

  const double* valsE = _data->valsResidual;
  const int totalNumVertices = mesh->depthStratum(0)->size();
  const int sizeE = _data->spaceDim * totalNumVertices;

  const double* vals = residual->restrict();
  const int size = residual->sizeWithBC();
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  //residual->view("RESIDUAL");

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
pylith::bc::TestAbsorbingDampers::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  CPPUNIT_ASSERT(0 != _data);

  ALE::Obj<ALE::Mesh> mesh;
  AbsorbingDampers bc;
  topology::FieldsManager fields(mesh);
  _initialize(&mesh, &bc, &fields);
  bc._needNewJacobian = true;

  const ALE::Obj<pylith::real_section_type>& dispTpdt = 
    fields.getReal("dispTpdt");
  CPPUNIT_ASSERT(!dispTpdt.isNull());

  PetscMat jacobian;
  PetscErrorCode err = MeshCreateMatrix(mesh, dispTpdt, MATMPIBAIJ, &jacobian);
  CPPUNIT_ASSERT(0 == err);

  const double t = 1.0;
  bc.integrateJacobian(&jacobian, t, &fields, mesh);
  CPPUNIT_ASSERT_EQUAL(false, bc.needNewJacobian());

  err = MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY);
  CPPUNIT_ASSERT(0 == err);
  err = MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY);
  CPPUNIT_ASSERT(0 == err);

  const double* valsE = _data->valsJacobian;
  const int totalNumVertices = mesh->depthStratum(0)->size();
  const int nrowsE = totalNumVertices * _data->spaceDim;
  const int ncolsE = totalNumVertices * _data->spaceDim;

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

#if 0
  std::cout << "JACOBIAN\n";
  for (int iRow=0, i=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol, ++i)
      std::cout << "  iRow: " << iRow << ", iCol: " << iCol << ", value: " << vals[i] << std::endl;
#endif

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
void
pylith::bc::TestAbsorbingDampers::_initialize(
					ALE::Obj<Mesh>* mesh,
					AbsorbingDampers* const bc,
					topology::FieldsManager* fields) const
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != bc);
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != _quadrature);

  try {
    // Setup mesh
    meshio::MeshIOAscii iohandler;
    iohandler.filename(_data->meshFilename);
    //iohandler.debug(true);
    iohandler.read(mesh);
    CPPUNIT_ASSERT(!mesh->isNull());

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim((*mesh)->getDimension());
    cs.initialize();

    // Setup quadrature
    _quadrature->initialize(_data->basis, _data->basisDerivRef, _data->quadPts,
			    _data->quadWts, _data->cellDim, _data->numBasis,
			    _data->numQuadPts, _data->spaceDim);

    spatialdata::spatialdb::SimpleDB db("TestAbsorbingDampers");
    spatialdata::spatialdb::SimpleIOAscii dbIO;
    dbIO.filename(_data->spatialDBFilename);
    db.ioHandler(&dbIO);
    db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

    const double upDirVals[] = { 0.0, 0.0, 1.0 };
    double_array upDir(upDirVals, 3);

    bc->quadrature(_quadrature);
    bc->timeStep(_data->dt);
    bc->label(_data->label);
    bc->db(&db);
    bc->initialize(*mesh, &cs, upDir);

    //bc->_boundaryMesh->view("BOUNDARY MESH");

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
    
    const int totalNumVertices = (*mesh)->depthStratum(0)->size();
    const int numMeshCells = (*mesh)->heightStratum(0)->size();
    const int fieldSize = _data->spaceDim * totalNumVertices;
    const ALE::Obj<real_section_type>& dispTpdt = fields->getReal("dispTpdt");
    const ALE::Obj<real_section_type>& dispT = fields->getReal("dispT");
    const ALE::Obj<real_section_type>& dispTmdt = fields->getReal("dispTmdt");
    CPPUNIT_ASSERT(!dispTpdt.isNull());
    CPPUNIT_ASSERT(!dispT.isNull());
    CPPUNIT_ASSERT(!dispTmdt.isNull());
    const int offset = numMeshCells;
    for (int iVertex=0; iVertex < totalNumVertices; ++iVertex) {
      dispTpdt->updatePoint(iVertex+offset, 
			    &_data->fieldTpdt[iVertex*_data->spaceDim]);
      dispT->updatePoint(iVertex+offset, 
			 &_data->fieldT[iVertex*_data->spaceDim]);
      dispTmdt->updatePoint(iVertex+offset, 
			    &_data->fieldTmdt[iVertex*_data->spaceDim]);
    } // for
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _initialize


// End of file 
