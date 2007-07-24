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
  ALE::Obj<Mesh> mesh;
  Neumann bc;
  _initialize(&mesh, &bc);

  CPPUNIT_ASSERT(0 != _data);

  // Check submesh
  CPPUNIT_ASSERT_EQUAL(_data.boundaryCellDim, _boundaryMesh->getDimension());
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
  const ALE::Obj<ALE::Mesh::label_sequence>& vertices =
    _boundaryMesh->heightStratum(0);
  const int numVertices = vertices->size();
  CPPUNIT_ASSERT_EQUAL(_data.numVertices, numVertices);
  int i = 0;
  const int spaceDim = _data.spaceDim;
  for(Mesh::label_sequance::iterator v_iter =
	vertices->begin();
      v_iter != vertices->end();
      ++v_iter) {
    const Mesh::real_section_type::value_type *tractionVals =
      _tractionGlobal->restrictPoint(*v_iter);
    const double tolerance = 1.0e-06;
    for(int iDim=0; iDim < spaceDim; ++iDim)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->tractionVals[i], tractionVals[iDim],
				   tolerance);
  } // for

} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestNeumann::testIntegrateResidual(void)
{ // testIntegrateResidual
  ALE::Obj<Mesh> mesh;
  Neumann bc;

  _initialize(&mesh, &bc);
  const ALE::Obj<real_section_type>& residual = fields->getReal("residual");
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  fields->addReal("dispTBctpdt");
  fields->copyLayout("residual");
  const double t = 0.0;

  CPPUNIT_ASSERT(0 != _data);

  const int numCells = mesh->heightStratum(0)->size();
  const int offset = numCells;
  int iConstraint = 0;
  for (Mesh::label_sequence::iterator v_iter = vertices->begin();
       v_iter != vertices->end();
       ++v_iter) {
    if (*v_iter != _data->constrainedPoints[iConstraint] + offset) {
      CPPUNIT_ASSERT_EQUAL(_data->numDOF, field->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(0, field->getConstraintDimension(*v_iter));
    } else {
      CPPUNIT_ASSERT_EQUAL(_data->numDOF, field->getFiberDimension(*v_iter));
      CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, 
			   field->getConstraintDimension(*v_iter));
      ++iConstraint;
    } // if/else
  } // for
} // testIntegrateResidual

// ----------------------------------------------------------------------
void
pylith::bc::TestNeumann::_initialize(ALE::Obj<Mesh>* mesh,
				       Neumann* const bc) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != bc);
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
} // _initialize


// End of file 
