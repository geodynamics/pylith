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

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

#include <stdexcept> // USES std::runtime_erro

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestNeumann );

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveMesh;
typedef pylith::topology::SubMesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealSection SubRealSection;
typedef pylith::topology::SubMesh::RestrictVisitor RestrictVisitor;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestNeumann::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = new feassemble::Quadrature<topology::SubMesh>();
  CPPUNIT_ASSERT(0 != _quadrature);
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
  topology::Mesh mesh;
  Neumann bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  CPPUNIT_ASSERT(0 != _data);

  const topology::SubMesh& boundaryMesh = *bc._boundaryMesh;
  const ALE::Obj<SieveSubMesh>& submesh = boundaryMesh.sieveMesh();

  // Check boundary mesh
  CPPUNIT_ASSERT(!submesh.isNull());

  const int cellDim = boundaryMesh.dimension();
  const int numCorners = _data->numCorners;
  const int spaceDim = _data->spaceDim;
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = submesh->heightStratum(1);
  const int numBoundaryVertices = submesh->depthStratum(0)->size();
  const int numBoundaryCells = cells->size();

  CPPUNIT_ASSERT_EQUAL(_data->cellDim, cellDim);
  CPPUNIT_ASSERT_EQUAL(_data->numBoundaryVertices, numBoundaryVertices);
  CPPUNIT_ASSERT_EQUAL(_data->numBoundaryCells, numBoundaryCells);

  const int boundaryDepth = submesh->depth()-1; // depth of boundary cells
  const ALE::Obj<SubRealSection>& coordinates =
    submesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates, numCorners*spaceDim);
  // coordinates->view("Mesh coordinates from TestNeumann::testInitialize");

  const int numBasis = bc._quadrature->numBasis();
  const int cellVertSize = _data->numCorners * spaceDim;
  double_array cellVertices(cellVertSize);

  const double tolerance = 1.0e-06;

  // check cell vertices
  int iCell = 0;
  for(SieveSubMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    const int numCorners = submesh->getNumCellCorners(*c_iter, boundaryDepth);
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numCorners);

    coordsVisitor.clear();
    submesh->restrictClosure(*c_iter, coordsVisitor);
    double vert =0.0;
    double vertE =0.0;
    const double* cellVertices = coordsVisitor.getValues();
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
  const ALE::Obj<SubRealSection>& tractionSection = bc._tractions->section();

  for(SieveSubMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    tractionSection->restrictPoint(*c_iter,
				   &tractionsCell[0], tractionsCell.size());
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
      for (int iDim =0; iDim < spaceDim; ++iDim) {
	const double tractionsCellData = _data->tractionsCell[index];
	CPPUNIT_ASSERT_DOUBLES_EQUAL(tractionsCellData,
				     tractionsCell[iQuad*spaceDim+iDim],
				     tolerance);
	++index;
      } // for
  } // for

} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestNeumann::testIntegrateResidual(void)
{ // testIntegrateResidual
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  Neumann bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  topology::Field<topology::Mesh>& residual = fields.get("residual");
  const double t = 0.0;
  bc.integrateResidual(residual, t, &fields);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  CPPUNIT_ASSERT(!sieveMesh->depthStratum(0).isNull());

  const double* valsE = _data->valsResidual;
  const int totalNumVertices = sieveMesh->depthStratum(0)->size();
  const int sizeE = _data->spaceDim * totalNumVertices;

  const ALE::Obj<RealSection>& residualSection = residual.section();
  CPPUNIT_ASSERT(!residualSection.isNull());

  const double* vals = residualSection->restrictSpace();
  const int size = residualSection->sizeWithBC();
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  //residual.view("RESIDUAL");

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
pylith::bc::TestNeumann::_initialize(topology::Mesh* mesh,
				     Neumann* const bc,
				     topology::SolutionFields* fields) const
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != mesh);
  CPPUNIT_ASSERT(0 != bc);
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _quadrature);

  try {
    // Set up mesh
    meshio::MeshIOAscii iohandler;
    iohandler.filename(_data->meshFilename);
    iohandler.read(mesh);

    // Set up coordinates
    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(mesh->dimension());
    cs.initialize();
    mesh->coordsys(&cs);

    // Set up quadrature
    _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			    _data->basisDerivRef, _data->numQuadPts, 
			    _data->numBasis, _data->cellDim,
			    _data->quadPts, _data->numQuadPts, _data->cellDim,
			    _data->quadWts, _data->numQuadPts,
			    _data->spaceDim);

    // Set up database
    spatialdata::spatialdb::SimpleDB db("TestNeumann");
    spatialdata::spatialdb::SimpleIOAscii dbIO;
    dbIO.filename(_data->spatialDBFilename);
    db.ioHandler(&dbIO);
    db.queryType(spatialdata::spatialdb::SimpleDB::LINEAR);

    const double upDir[] = { 0.0, 0.0, 1.0 };

    bc->quadrature(_quadrature);
    bc->label(_data->label);
    bc->db(&db);
    bc->initialize(*mesh, upDir);

    // Set up fields
    CPPUNIT_ASSERT(0 != fields);
    fields->add("residual");
    fields->add("dispTBctpdt");
    fields->solutionName("dispTBctpdt");

    topology::Field<topology::Mesh>& residual = fields->get("residual");
    const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
    CPPUNIT_ASSERT(!sieveMesh.isNull());
    const ALE::Obj<SieveMesh::label_sequence>& vertices = 
      sieveMesh->depthStratum(0);
    CPPUNIT_ASSERT(!vertices.isNull());
    residual.newSection(vertices, _data->spaceDim);
    residual.allocate();
    residual.zero();

    fields->copyLayout("residual");
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _initialize


// End of file 
