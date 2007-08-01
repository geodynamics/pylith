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

#include "TestNeumann.hh" // Implementation of class methods

#include "pylith/bc/Neumann.hh" // USES Neumann

#include "data/NeumannData.hh" // USES NeumannData
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumann );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestNeumann::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestNeumann::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::bc::TestNeumann::testConstructor(void)
{ // testConstructor
  Neumann bc;
} // testConstructor

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestNeumann::testInitialize(void)
{ // testInitialize
#if 0
  ALE::Obj<Mesh> mesh;
  Neumann bc;
  topology::FieldsManager fields(mesh);
  _initialize(&mesh, &bc, &fields);

  CPPUNIT_ASSERT(0 != _data);

  // Check submesh
  CPPUNIT_ASSERT_EQUAL(_data.cellDim, _boundaryMesh->getDimension());
  CPPUNIT_ASSERT_EQUAL(_data.numBoundaryCells, cells->size());
  int iCell = 0;
  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    const int numCorners = sieve->nCone(*c_iter, mesh->depth())->size();
    CPPUNIT_ASSERT_EQUAL(_data.numCorners[iCell++], numCorners);
    const ALE::Obj<sieve_type::traits::coneSequence>& cone =
      sieve->cone(*c_iter);
    for(sieve_type::traits::coneSequence::iterator v_iter = cone->begin();
        v_iter != cone->end();
        ++v_iter)
      CPPUNIT_ASSERT_EQUAL(_data.cells[i++], *v_iter);
  } // for

  // Check traction values
  numQuadPts = _data->numQuadPts;
  spaceDim = _data->spaceDim;
  fiberDim = numQuadPts * spaceDim;
  double_array tractionCell(fiberDim);
  int index = 0;
  const double tolerance = 1.0e-06;

  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {

    _boundaryMesh->restrict(_tractionGlobal, *c_iter,
			    &tractionCell[0], tractionCell.size());

    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      for (int iDim =0; iDim < spaceDim; ++iDim) {
	CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->tractionCell[index],
				     tractionCell[iQuad*spaceDim+iDim],
				     tolerance);
	++index;
      } // for
    } // for
  } // for
#endif

} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestNeumann::testIntegrateResidual(void)
{ // testIntegrateResidual
#if 0
  CPPUNIT_ASSERT(0 != _data);

  ALE::Obj<Mesh> mesh;
  Neumann bc;
  Neumann integrator;
  topology::FieldsManager fields(mesh);
  _initialize(&mesh, &bc, &fields);

  const ALE::Obj<real_section_type>& residual = fields->getReal("residual");
  CPPUNIT_ASSERT(!residual.isNull());
  const double t = 0.0;
  integrator.integrateResidual(residual, t, &fields, mesh);

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

#endif
} // testIntegrateResidual

// ----------------------------------------------------------------------
void
pylith::bc::TestNeumann::_initialize(ALE::Obj<Mesh>* mesh,
				     Neumann* const bc) const

// Maybe need to deal with fields inside integrateResidual
				     // Neumann* const bc,
				     // topology::FieldsManager* fields) const
{ // _initialize
#if 0
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != bc);
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _quadrature);

  // Set up mesh
  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);
  CPPUNIT_ASSERT(!mesh->isNull());
  (*mesh)->getFactory()->clear();

  // Set up coordinates
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim((*mesh)->getDimension());
  cs.initialize();

  // Set up fields
  fields->addReal("residual");
  fields->addReal("solution");
  const ALE::Obj<real_section_type>& residual = fields->getReal("residual");
 CPPUNIT_ASSERT(!residual.isNull());
 residual->setFiberDimension((*mesh)->depthStratum(0), _data->spaceDim);
  (*mesh)->allocate(residual);
  residual->zero();
  fields->copyLayout("residual");
  const ALE::Obj<real_section_type>& solution = fields.getReal("solution");
  CPPUNIT_ASSERT(!solution.isNull());

  // Set up quadrature
  _quadrature->initialize(_data->basis, _data->basisDeriv, _data->quadPts,
			  _data->quadWts, _data->cellDim, _data->numBasis,
			  _data->numQuadPts, _data->spaceDim);

  // Set up database
  spatialdata::spatialdb::SimpleDB db("TestNeumann");
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename(_data->dbFilename);
  db.ioHandler(&dbIO);

  const double upDirVals[] = { 0.0, 0.0, 1.0 };
  double_array upDir(upDirVals, 3);

  bc->id(_data->id);
  bc->label(_data->label);
  bc->db(&db);
  bc->initialize(*mesh, &cs, upDir);
#endif
} // _initialize


// End of file 
