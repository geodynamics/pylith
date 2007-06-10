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

#include "TestFaultCohesiveKin.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "data/CohesiveKinData.hh" // USES CohesiveKinData

#include "pylith/faults/EqKinSrc.hh" // USES EqKinSrc
#include "pylith/faults/BruneSlipFn.hh" // USES BruneSlipFn
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

#include <stdexcept> // USES runtime_error

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKin );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKin::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = 0;
  _eqsrc = new EqKinSrc();
  _slipfn = new BruneSlipFn();
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::faults::TestFaultCohesiveKin::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
  delete _eqsrc; _eqsrc = 0;
  delete _slipfn; _slipfn = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestFaultCohesiveKin::testConstructor(void)
{ // testConstructor
  FaultCohesiveKin fault;
} // testConstructor

// ----------------------------------------------------------------------
// Test eqsrc().
void
pylith::faults::TestFaultCohesiveKin::testEqsrc(void)
{ // testEqsrc
  FaultCohesiveKin fault;

  EqKinSrc eqsrc;
  fault.eqsrc(&eqsrc);
  CPPUNIT_ASSERT(&eqsrc == fault._eqsrc);
} // testEqsrc

// ----------------------------------------------------------------------
// Test useLagrangeConstraints().
void
pylith::faults::TestFaultCohesiveKin::testUseLagrangeConstraints(void)
{ // testUseLagrangeConstraints
  FaultCohesiveKin fault;
  CPPUNIT_ASSERT_EQUAL(true, fault._useLagrangeConstraints());
} // testUseLagrangeConstraints

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestFaultCohesiveKin::testInitialize(void)
{ // testInitialize
  ALE::Obj<Mesh> mesh;
  FaultCohesiveKin fault;
  _initialize(&mesh, &fault);
  
  // Check set of constraint vertices
  CPPUNIT_ASSERT_EQUAL(_data->numConstraintVert,
		       int(fault._constraintVert.size()));
  const std::set<Mesh::point_type>::const_iterator vertConstraintBegin = 
    fault._constraintVert.begin();
  const std::set<Mesh::point_type>::const_iterator vertConstraintEnd = 
    fault._constraintVert.end();
  int iVertex = 0;
  for (std::set<Mesh::point_type>::const_iterator v_iter=vertConstraintBegin;
       v_iter != vertConstraintEnd;
       ++v_iter, ++iVertex)
    CPPUNIT_ASSERT_EQUAL(_data->constraintVertices[iVertex],
			 *v_iter);

  // Check orientation
  iVertex = 0;
  const int cellDim = _data->cellDim;
  const int spaceDim = _data->spaceDim;
  const int orientationSize = spaceDim*spaceDim;
  for (std::set<Mesh::point_type>::const_iterator v_iter=vertConstraintBegin;
       v_iter != vertConstraintEnd;
       ++v_iter) {
    const int fiberDim = fault._orientation->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(orientationSize, fiberDim);
    const real_section_type::value_type* vertexOrient = 
      fault._orientation->restrictPoint(*v_iter);
    CPPUNIT_ASSERT(0 != vertexOrient);
    const double tolerance = 1.0e-06;
    for (int i=0; i < orientationSize; ++i) {
      const int index = iVertex*orientationSize+i;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->orientation[index],
				   vertexOrient[i], tolerance);
    } // for
  } // for
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::faults::TestFaultCohesiveKin::testIntegrateResidual(void)
{ // testIntegrateResidual
  ALE::Obj<Mesh> mesh;
  FaultCohesiveKin fault;
  _initialize(&mesh, &fault);

  // Setup fields
  topology::FieldsManager fields(mesh);
  fields.addReal("residual");
  fields.addReal("dispTpdt");
  fields.addReal("dispT");
  const char* history[] = { "dispTpdt", "dispT" };
  const int historySize = 2;
  fields.createHistory(history, historySize);
  
  const ALE::Obj<real_section_type>& residual = fields.getReal("residual");
  CPPUNIT_ASSERT(!residual.isNull());
  const int spaceDim = _data->spaceDim;
  residual->setFiberDimension(mesh->depthStratum(0), spaceDim);
  mesh->allocate(residual);
  residual->zero();
  fields.copyLayout("residual");

  const ALE::Obj<real_section_type>& dispTpdt = fields.getReal("dispTpdt");
  const ALE::Obj<real_section_type>& dispT = fields.getReal("dispT");
  CPPUNIT_ASSERT(!dispTpdt.isNull());
  CPPUNIT_ASSERT(!dispT.isNull());

  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const Mesh::label_sequence::iterator vBegin = vertices->begin();
  const Mesh::label_sequence::iterator vEnd = vertices->end();
  int iVertex = 0;
  for (Mesh::label_sequence::iterator v_iter=vBegin;
       v_iter != vEnd;
       ++v_iter, ++iVertex) {
    dispTpdt->updatePoint(*v_iter, &_data->fieldTpdt[iVertex*spaceDim]);
    dispT->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);
  } // for
  
  //dispT->view("DISP T");

  // Call integrateResidual()
  fault.integrateResidual(residual, &fields, mesh);

  //residual->view("RESIDUAL");

  // Check values
  const double* valsE = _data->valsResidual;
  iVertex = 0;
  const int fiberDimE = spaceDim;
  const double tolerance = 1.0e-06;
  for (Mesh::label_sequence::iterator v_iter=vBegin;
       v_iter != vEnd;
       ++v_iter, ++iVertex) {
    const int fiberDim = residual->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
    const real_section_type::value_type* vals = 
      residual->restrictPoint(*v_iter);
    for (int i=0; i < fiberDimE; ++i) {
      const int index = iVertex*spaceDim+i;
      if (valsE[index] > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[index], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[index], vals[i], tolerance);
    } // for
  } // for
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::faults::TestFaultCohesiveKin::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  ALE::Obj<Mesh> mesh;
  FaultCohesiveKin fault;
  _initialize(&mesh, &fault);

  // Setup fields
  topology::FieldsManager fields(mesh);
  fields.addReal("residual");
  fields.addReal("dispTpdt");
  fields.addReal("dispT");
  const char* history[] = { "dispTpdt", "dispT" };
  const int historySize = 2;
  fields.createHistory(history, historySize);
  
  const ALE::Obj<real_section_type>& residual = fields.getReal("residual");
  CPPUNIT_ASSERT(!residual.isNull());
  const int spaceDim = _data->spaceDim;
  residual->setFiberDimension(mesh->depthStratum(0), spaceDim);
  mesh->allocate(residual);
  residual->zero();
  fields.copyLayout("residual");

  const ALE::Obj<real_section_type>& dispTpdt = fields.getReal("dispTpdt");
  const ALE::Obj<real_section_type>& dispT = fields.getReal("dispT");
  CPPUNIT_ASSERT(!dispTpdt.isNull());
  CPPUNIT_ASSERT(!dispT.isNull());

  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const Mesh::label_sequence::iterator vBegin = vertices->begin();
  const Mesh::label_sequence::iterator vEnd = vertices->end();
  int iVertex = 0;
  for (Mesh::label_sequence::iterator v_iter=vBegin;
       v_iter != vEnd;
       ++v_iter, ++iVertex) {
    dispTpdt->updatePoint(*v_iter, &_data->fieldTpdt[iVertex*spaceDim]);
    dispT->updatePoint(*v_iter, &_data->fieldT[iVertex*spaceDim]);
  } // for
  
  PetscMat jacobian;
  PetscErrorCode err = MeshCreateMatrix(mesh, dispTpdt, MATMPIBAIJ, &jacobian);
  CPPUNIT_ASSERT(0 == err);

  fault.integrateJacobian(&jacobian, &fields, mesh);
  CPPUNIT_ASSERT_EQUAL(false, fault.needNewJacobian());

  err = MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY);
  CPPUNIT_ASSERT(0 == err);
  err = MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY);
  CPPUNIT_ASSERT(0 == err);

  //MatView(jacobian, PETSC_VIEWER_STDOUT_WORLD);

  const double* valsE = _data->valsJacobian;
  const int nrowsE = dispT->sizeWithBC();
  const int ncolsE = nrowsE;

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
// Test setConstraintSizes().
void
pylith::faults::TestFaultCohesiveKin::testSetConstraintSizes(void)
{ // testSetConstraintSizes
  // Make sure the fiber dimension at each point is equal to the
  // spatial dimension of the mesh, because the Lagrange multiplier
  // formation does not eliminate any DOF from the system of
  // equations.

  ALE::Obj<Mesh> mesh;
  FaultCohesiveKin fault;
  _initialize(&mesh, &fault);
  
  const int fiberDim = 3;
  const ALE::Obj<real_section_type>& field = 
    new real_section_type(mesh->comm(), mesh->debug());
  CPPUNIT_ASSERT(!field.isNull());

  field->setFiberDimension(mesh->depthStratum(0), fiberDim);
  fault.setConstraintSizes(field, mesh);
  mesh->allocate(field);

  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const Mesh::label_sequence::iterator vBegin = vertices->begin();
  const Mesh::label_sequence::iterator vEnd = vertices->end();
  for (Mesh::label_sequence::iterator v_iter=vBegin;
       v_iter != vEnd;
       ++v_iter) {
    CPPUNIT_ASSERT_EQUAL(fiberDim, field->getFiberDimension(*v_iter));
    CPPUNIT_ASSERT_EQUAL(0, field->getConstraintDimension(*v_iter));
  } // for
} // testSetConstraintSizes

// ----------------------------------------------------------------------
// Test setField().
void
pylith::faults::TestFaultCohesiveKin::testSetField(void)
{ // testSetField
  ALE::Obj<Mesh> mesh;
  FaultCohesiveKin fault;
  _initialize(&mesh, &fault);

  // Setup fields
  topology::FieldsManager fields(mesh);
  fields.addReal("disp");
  
  const ALE::Obj<real_section_type>& disp = fields.getReal("disp");
  CPPUNIT_ASSERT(!disp.isNull());
  const int spaceDim = _data->spaceDim;
  disp->setFiberDimension(mesh->depthStratum(0), spaceDim);
  mesh->allocate(disp);
  disp->zero();

  const double t = 2.134;
  fault.setField(t, disp, mesh);

  // Check values
  const double* valsE = _data->valsSlip;
  const int fiberDimE = spaceDim;
  const double tolerance = 1.0e-06;

  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  CPPUNIT_ASSERT(!vertices.isNull());
  const Mesh::label_sequence::iterator vBegin = vertices->begin();
  const Mesh::label_sequence::iterator vEnd = vertices->end();
  int iVertex = 0;
  for (Mesh::label_sequence::iterator v_iter=vBegin;
       v_iter != vEnd;
       ++v_iter, ++iVertex) {
    const int fiberDim = disp->getFiberDimension(*v_iter);
    CPPUNIT_ASSERT_EQUAL(fiberDimE, fiberDim);
    const real_section_type::value_type* vals = 
      disp->restrictPoint(*v_iter);
    for (int i=0; i < fiberDimE; ++i) {
      const int index = iVertex*spaceDim+i;
      if (valsE[index] > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[index], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[index], vals[i], tolerance);
    } // for
  } // for
} // testSetField

// ----------------------------------------------------------------------
// Initialize FaultCohesiveKin interface condition.
void
pylith::faults::TestFaultCohesiveKin::_initialize(ALE::Obj<ALE::Mesh>* mesh,
					FaultCohesiveKin* const fault) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != fault);
  CPPUNIT_ASSERT(0 != _quadrature);
  CPPUNIT_ASSERT(0 != _eqsrc);
  CPPUNIT_ASSERT(0 != _slipfn);

  try {
    meshio::MeshIOAscii iohandler;
    iohandler.filename(_data->meshFilename);
    iohandler.read(mesh);
    CPPUNIT_ASSERT(!mesh->isNull());
    (*mesh)->getFactory()->clear();

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim((*mesh)->getDimension());
    cs.initialize();

    _quadrature->initialize(_data->verticesRef, 
			    _data->basis, _data->basisDeriv, _data->quadPts,
			    _data->quadWts, _data->cellDim, _data->numBasis,
			    _data->numQuadPts, _data->spaceDim);

    // Setup earthquake source
    spatialdata::spatialdb::SimpleDB dbFinalSlip("final slip");
    spatialdata::spatialdb::SimpleIOAscii ioFinalSlip;
    ioFinalSlip.filename(_data->finalSlipFilename);
    dbFinalSlip.ioHandler(&ioFinalSlip);
  
    spatialdata::spatialdb::SimpleDB dbSlipTime("slip time");
    spatialdata::spatialdb::SimpleIOAscii ioSlipTime;
    ioSlipTime.filename(_data->slipTimeFilename);
    dbSlipTime.ioHandler(&ioSlipTime);
  
    spatialdata::spatialdb::SimpleDB dbPeakRate("peak rate");
    spatialdata::spatialdb::SimpleIOAscii ioPeakRate;
    ioPeakRate.filename(_data->peakRateFilename);
    dbPeakRate.ioHandler(&ioPeakRate);

    _slipfn->dbFinalSlip(&dbFinalSlip);
    _slipfn->dbSlipTime(&dbSlipTime);
    _slipfn->dbPeakRate(&dbPeakRate);
  
    _eqsrc->slipfn(_slipfn);
  
    fault->id(_data->id);
    fault->label(_data->label);
    fault->quadrature(_quadrature);
    fault->eqsrc(_eqsrc);
    fault->adjustTopology(*mesh);

    const double upDirVals[] = { 0.0, 0.0, 1.0 };
    double_array upDir(upDirVals, 3);

    fault->initialize(*mesh, &cs, upDir); 
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _initialize


// End of file 
