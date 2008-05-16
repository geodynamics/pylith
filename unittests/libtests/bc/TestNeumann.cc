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

#include <stdexcept> // USES std::runtime_erro

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
  topology::FieldsManager fields(mesh);
  _initialize(&mesh, &bc, &fields);

  CPPUNIT_ASSERT(0 != _data);

  const ALE::Obj<Mesh>& boundaryMesh = bc._boundaryMesh;

  // Check boundary mesh
  CPPUNIT_ASSERT(!boundaryMesh.isNull());

  const int cellDim = boundaryMesh->getDimension();
  const ALE::Obj<Mesh::label_sequence>& cells = boundaryMesh->heightStratum(1);
  const int numBoundaryVertices = boundaryMesh->depthStratum(0)->size();
  const int numBoundaryCells = cells->size();

  CPPUNIT_ASSERT_EQUAL(_data->cellDim, cellDim);
  CPPUNIT_ASSERT_EQUAL(_data->numBoundaryVertices, numBoundaryVertices);
  CPPUNIT_ASSERT_EQUAL(_data->numBoundaryCells, numBoundaryCells);

  // boundaryMesh->view("BOUNDARY MESH");

  const int boundaryDepth = boundaryMesh->depth()-1; // depth of boundary cells
  const ALE::Obj<real_section_type>& coordinates =
    mesh->getRealSection("coordinates");
  // coordinates->view("Mesh coordinates from TestNeumann::testInitialize");

  const int spaceDim = _data->spaceDim;
  const int numBasis = bc._quadrature->numBasis();
  const int cellVertSize = _data->numCorners * spaceDim;
  double_array cellVertices(cellVertSize);

  const double tolerance = 1.0e-06;

  // check cell vertices
  int iCell = 0;
  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    const int numCorners = boundaryMesh->getNumCellCorners(*c_iter, boundaryDepth);
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numCorners);

    boundaryMesh->restrict(coordinates, *c_iter, &cellVertices[0],
			   cellVertices.size());
    double vert =0.0;
    double vertE =0.0;
    // std::cout << "c_iter " << *c_iter << " vertex coords:" << std::endl;
    for(int iVert = 0; iVert < numCorners; ++iVert) {
      for(int iDim = 0; iDim < spaceDim; ++iDim) {
	vertE = _data->cellVertices[iDim+spaceDim*iVert+iCell*cellVertSize];
	vert = cellVertices[iDim+spaceDim*iVert];
        // std::cout << "  " << cellVertices[iDim+spaceDim*iVert];
	if (fabs(vertE) > 1.0)
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vert/vertE, tolerance);
	else
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(vert, vertE, tolerance);
      } // for
      // std::cout << std::endl;
    } // for
    iCell++;
  } // for

  // Check traction values
  const int numQuadPts = _data->numQuadPts;
  const int fiberDim = numQuadPts * spaceDim;
  double_array tractionsCell(fiberDim);
  int index = 0;

  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {

    bc._boundaryMesh->restrict(bc._tractions, *c_iter,
			       &tractionsCell[0], tractionsCell.size());

    // std::cout << "Tractions at quadrature points: " << std::endl;
    // std::cout << "Computed    Expected" << std::endl;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      for (int iDim =0; iDim < spaceDim; ++iDim) {
	const double tractionsCellData = _data->tractionsCell[index];
        // std::cout << "  " << tractionsCell[iQuad*spaceDim+iDim] << "   " << tractionsCellData << std::endl;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(tractionsCellData,
				     tractionsCell[iQuad*spaceDim+iDim],
				     tolerance);
	++index;
      } // for
      // std::cout << std::endl;
    } // for
  } // for

} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestNeumann::testIntegrateResidual(void)
{ // testIntegrateResidual
  CPPUNIT_ASSERT(0 != _data);

  ALE::Obj<Mesh> mesh;
  Neumann bc;
  topology::FieldsManager fields(mesh);
  _initialize(&mesh, &bc, &fields);

  const ALE::Obj<real_section_type>& residual = fields.getReal("residual");
  CPPUNIT_ASSERT(!residual.isNull());

  const int spaceDim = _data->spaceDim;

  const double t = 0.0;
  bc.integrateResidual(residual, t, &fields, mesh);

  const double* valsE = _data->valsResidual;
  const int totalNumVertices = mesh->depthStratum(0)->size();
  const int sizeE = _data->spaceDim * totalNumVertices;

  const double* vals = residual->restrict();
  const int size = residual->sizeWithBC();
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  //residual->view("RESIDUAL");

  const double tolerance = 1.0e-06;
  // std::cout << "computed residuals: " << std::endl;
  for (int i=0; i < size; ++i)
    // std::cout << "  " << vals[i] << std::endl;
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i], vals[i], tolerance);
  
} // testIntegrateResidual

// ----------------------------------------------------------------------
void
pylith::bc::TestNeumann::_initialize(ALE::Obj<Mesh>* mesh,
				     Neumann* const bc,
				     topology::FieldsManager* fields) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != bc);
  CPPUNIT_ASSERT(0 != _quadrature);

  try {
    // Set up mesh
    meshio::MeshIOAscii iohandler;
    iohandler.filename(_data->meshFilename);
    iohandler.read(mesh);
    CPPUNIT_ASSERT(!mesh->isNull());

    // Set up coordinates
    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim((*mesh)->getDimension());
    cs.initialize();

    // Set up quadrature
    _quadrature->initialize(_data->basis, _data->basisDerivRef, _data->quadPts,
			    _data->quadWts, _data->cellDim, _data->numBasis,
			    _data->numQuadPts, _data->spaceDim);

    // Set up database
    spatialdata::spatialdb::SimpleDB db("TestNeumann");
    spatialdata::spatialdb::SimpleIOAscii dbIO;
    dbIO.filename(_data->spatialDBFilename);
    db.ioHandler(&dbIO);
    db.queryType(spatialdata::spatialdb::SimpleDB::LINEAR);

    const double upDirVals[] = { 0.0, 0.0, 1.0 };
    double_array upDir(upDirVals, 3);

    bc->quadrature(_quadrature);
    bc->label(_data->label);
    bc->db(&db);
    bc->initialize(*mesh, &cs, upDir);

    // Set up fields
    CPPUNIT_ASSERT(0 != fields);
    fields->addReal("residual");
    fields->addReal("dispTBctpdt");

    const ALE::Obj<real_section_type>& residual = fields->getReal("residual");
    CPPUNIT_ASSERT(!residual.isNull());
    residual->setChart((*mesh)->getSieve()->getChart());
    residual->setFiberDimension((*mesh)->depthStratum(0), _data->spaceDim);
    (*mesh)->allocate(residual);
    residual->zero();
    fields->copyLayout("residual");
    const ALE::Obj<real_section_type>& dispTBctpdt = 
      fields->getReal("dispTBctpdt");
    CPPUNIT_ASSERT(!dispTBctpdt.isNull());
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _initialize


// End of file 
