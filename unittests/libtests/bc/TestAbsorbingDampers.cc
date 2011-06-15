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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestAbsorbingDampers.hh" // Implementation of class methods

#include "pylith/bc/AbsorbingDampers.hh" // USES AbsorbingDampers

#include "data/AbsorbingDampersData.hh" // USES AbsorbingDampersData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii

#include <stdexcept> // USES std::runtime_erro

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestAbsorbingDampers );

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveMesh;
typedef pylith::topology::SubMesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealUniformSection SubRealUniformSection;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestAbsorbingDampers::setUp(void)
{ // setUp
  _data = 0;
  _quadrature = new feassemble::Quadrature<topology::SubMesh>();
  CPPUNIT_ASSERT(0 != _quadrature);
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
// Test db().
void
pylith::bc::TestAbsorbingDampers::testDB(void)
{ // testDB
  const std::string& label = "my db";
  spatialdata::spatialdb::SimpleDB db(label.c_str());
  AbsorbingDampers bc;
  bc.db(&db);
  
  CPPUNIT_ASSERT(0 != bc._db);
  CPPUNIT_ASSERT_EQUAL(label, std::string(bc._db->label()));
} // testDB
    
// ----------------------------------------------------------------------
// Test initialize().
void
pylith::bc::TestAbsorbingDampers::testInitialize(void)
{ // testInitialize
  topology::Mesh mesh;
  AbsorbingDampers bc;
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
  const int numVertices = submesh->depthStratum(0)->size();
  const int numCells = cells->size();
  const int boundaryDepth = submesh->depth()-1; // depth of boundary cells

  CPPUNIT_ASSERT_EQUAL(_data->cellDim, cellDim);
  CPPUNIT_ASSERT_EQUAL(_data->numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(_data->numCells, numCells);

  const ALE::Obj<SieveMesh::sieve_type>& sieve = submesh->getSieve();
  ALE::ISieveVisitor::PointRetriever<SieveSubMesh::sieve_type> pV(sieve->getMaxConeSize());
  int dp = 0;
  for(SieveSubMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    const int numCorners = submesh->getNumCellCorners(*c_iter, boundaryDepth);
    CPPUNIT_ASSERT_EQUAL(_data->numCorners, numCorners);

    sieve->cone(*c_iter, pV);
    const SieveSubMesh::point_type *cone = pV.getPoints();
    for(int p = 0; p < pV.getSize(); ++p, ++dp) {
      CPPUNIT_ASSERT_EQUAL(_data->cells[dp], cone[p]);
    }
    pV.clear();
  } // for

  // Check damping constants
  const int numQuadPts = _data->numQuadPts;
  const int fiberDim = numQuadPts * spaceDim;
  double_array dampersCell(fiberDim);
  int index = 0;
  CPPUNIT_ASSERT(0 != bc._parameters);
  const ALE::Obj<SubRealUniformSection>& dampersSection = 
    bc._parameters->section();

  const double tolerance = 1.0e-06;
  for(SieveSubMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    dampersSection->restrictPoint(*c_iter,
				  &dampersCell[0], dampersCell.size());
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
      for (int iDim =0; iDim < spaceDim; ++iDim) {
	const double dampersCellData = _data->dampingConsts[index];
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, 
				     dampersCell[iQuad*spaceDim+iDim]/dampersCellData,
				     tolerance);
	++index;
      } // for
  } // for
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::bc::TestAbsorbingDampers::testIntegrateResidual(void)
{ // testIntegrateResidual
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  const topology::SubMesh& boundaryMesh = *bc._boundaryMesh;
  const ALE::Obj<SieveSubMesh>& submesh = boundaryMesh.sieveMesh();

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

  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  const topology::SubMesh& boundaryMesh = *bc._boundaryMesh;
  const ALE::Obj<SieveSubMesh>& submesh = boundaryMesh.sieveMesh();

  topology::Field<topology::Mesh>& solution = fields.solution();
  const ALE::Obj<RealSection>& solutionSection = solution.section();
  CPPUNIT_ASSERT(!solutionSection.isNull());

  topology::Jacobian jacobian(solution);

  const double t = 1.0;
  bc.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, bc.needNewJacobian());
  jacobian.assemble("final_assembly");

  CPPUNIT_ASSERT(!sieveMesh->depthStratum(0).isNull());

  const double* valsE = _data->valsJacobian;
  const int totalNumVertices = sieveMesh->depthStratum(0)->size();
  const int nrowsE = totalNumVertices * _data->spaceDim;
  const int ncolsE = totalNumVertices * _data->spaceDim;

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

#if 0
  std::cout << "JACOBIAN\n";
  for (int iRow=0, i=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol, ++i)
      std::cout << "  iRow: " << iRow << ", iCol: " << iCol << ", value: " << vals[i] << ", valueE: " << valsE[i] << std::endl;
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
  MatDestroy(&jDense);
  MatDestroy(&jSparseAIJ);
} // testIntegrateJacobian

// ----------------------------------------------------------------------
// Test integrateJacobianLumped().
void
pylith::bc::TestAbsorbingDampers::testIntegrateJacobianLumped(void)
{ // testIntegrateJacobianLumped
  CPPUNIT_ASSERT(0 != _data);

  topology::Mesh mesh;
  AbsorbingDampers bc;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &bc, &fields);

  topology::Field<topology::Mesh> jacobian(mesh);
  jacobian.label("Jacobian");
  jacobian.vectorFieldType(topology::FieldBase::VECTOR);
  jacobian.newSection(topology::FieldBase::VERTICES_FIELD, _data->spaceDim);
  jacobian.allocate();

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());

  const topology::SubMesh& boundaryMesh = *bc._boundaryMesh;
  const ALE::Obj<SieveSubMesh>& submesh = boundaryMesh.sieveMesh();
  CPPUNIT_ASSERT(!submesh.isNull());

  topology::Field<topology::Mesh>& solution = fields.solution();
  const ALE::Obj<RealSection>& solutionSection = solution.section();
  CPPUNIT_ASSERT(!solutionSection.isNull());

  const double t = 1.0;
  bc.integrateJacobian(&jacobian, t, &fields);
  CPPUNIT_ASSERT_EQUAL(false, bc.needNewJacobian());
  jacobian.complete();

  const double* valsMatrixE = _data->valsJacobian;
  const int totalNumVertices = sieveMesh->depthStratum(0)->size();
  const int sizeE = totalNumVertices * _data->spaceDim;
  double_array valsE(sizeE);
  const int spaceDim = _data->spaceDim;
  for (int iVertex=0; iVertex < totalNumVertices; ++iVertex)
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      const int indexRow = (iVertex*spaceDim+iDim)*totalNumVertices*spaceDim;
      double value = 0.0;
      for (int jVertex=0; jVertex < totalNumVertices; ++jVertex)
	value += valsMatrixE[indexRow + jVertex*spaceDim+iDim];
      valsE[iVertex*spaceDim+iDim] = value;
    } // for

#if 0 // DEBUGGING
  jacobian.view("JACOBIAN");
  std::cout << "\n\nJACOBIAN FULL" << std::endl;
  const int n = totalNumVertices*spaceDim;
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
  CPPUNIT_ASSERT_EQUAL(sizeE, size);

  const double tolerance = 1.0e-06;
  for (int i=0; i < size; ++i)
    if (fabs(valsE[i]) > 1.0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vals[i]/valsE[i], tolerance);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsE[i], vals[i], tolerance);

} // testIntegrateJacobianLumped

// ----------------------------------------------------------------------
void
pylith::bc::TestAbsorbingDampers::_initialize(topology::Mesh* mesh,
					      AbsorbingDampers* const bc,
					      topology::SolutionFields* fields) const
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
    iohandler.read(mesh);

    // Set coordinate system
    spatialdata::geocoords::CSCart cs;
    spatialdata::units::Nondimensional normalizer;
    cs.setSpaceDim(mesh->dimension());
    cs.initialize();
    mesh->coordsys(&cs);
    mesh->nondimensionalize(normalizer);

    // Set up quadrature
    _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			    _data->basisDerivRef, _data->numQuadPts, 
			    _data->numBasis, _data->cellDim,
			    _data->quadPts, _data->numQuadPts, _data->cellDim,
			    _data->quadWts, _data->numQuadPts,
			    _data->spaceDim);

    // Set up database
    spatialdata::spatialdb::SimpleDB db("TestAbsorbingDampers");
    spatialdata::spatialdb::SimpleIOAscii dbIO;
    dbIO.filename(_data->spatialDBFilename);
    db.ioHandler(&dbIO);
    db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

    const double upDir[] = { 0.0, 0.0, 1.0 };

    bc->quadrature(_quadrature);
    bc->timeStep(_data->dt);
    bc->label(_data->label);
    bc->db(&db);
    bc->createSubMesh(*mesh);
    bc->initialize(*mesh, upDir);

    //bc->_boundaryMesh->view("BOUNDARY MESH");

    // Setup fields
    CPPUNIT_ASSERT(0 != fields);
    fields->add("residual", "residual");
    fields->add("dispIncr(t->t+dt)", "displacement_increment");
    fields->add("disp(t)", "displacement");
    fields->add("disp(t-dt)", "displacement");
    fields->add("velocity(t)", "velocity");
    fields->solutionName("dispIncr(t->t+dt)");
  
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

    const int totalNumVertices = sieveMesh->depthStratum(0)->size();
    const int numMeshCells = sieveMesh->heightStratum(0)->size();
 int fieldSize = _data->spaceDim * totalNumVertices;
    const ALE::Obj<RealSection>& dispTIncrSection = 
      fields->get("dispIncr(t->t+dt)").section();
    const ALE::Obj<RealSection>& dispTSection = 
      fields->get("disp(t)").section();
    const ALE::Obj<RealSection>& dispTmdtSection = 
      fields->get("disp(t-dt)").section();
    const ALE::Obj<RealSection>& velSection = 
      fields->get("velocity(t)").section();
    CPPUNIT_ASSERT(!dispTIncrSection.isNull());
    CPPUNIT_ASSERT(!dispTSection.isNull());
    CPPUNIT_ASSERT(!dispTmdtSection.isNull());
    CPPUNIT_ASSERT(!velSection.isNull());
    const int offset = numMeshCells;
    const int spaceDim = _data->spaceDim;
    const double dt = _data->dt;
    double_array velVertex(spaceDim);
    for (int iVertex=0; iVertex < totalNumVertices; ++iVertex) {
      dispTIncrSection->updatePoint(iVertex+offset, 
				    &_data->fieldTIncr[iVertex*spaceDim]);
      dispTSection->updatePoint(iVertex+offset, 
				&_data->fieldT[iVertex*spaceDim]);
      dispTmdtSection->updatePoint(iVertex+offset, 
				   &_data->fieldTmdt[iVertex*spaceDim]);
      
      for (int iDim=0; iDim < spaceDim; ++iDim)
	velVertex[iDim] = (_data->fieldTIncr[iVertex*spaceDim+iDim] +
			   _data->fieldT[iVertex*spaceDim+iDim] -
			   _data->fieldTmdt[iVertex*spaceDim+iDim]) / (2*dt);

      velSection->updatePoint(iVertex+offset, &velVertex[0]);
    } // for
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _initialize


// End of file 
