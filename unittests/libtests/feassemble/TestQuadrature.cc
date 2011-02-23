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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestQuadrature.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

#include "pylith/feassemble/GeometryLine1D.hh" // USES GeometryLine1D
#include "pylith/feassemble/GeometryTri2D.hh" // USES GeometryTri2D

#include "data/QuadratureData2DLinear.hh" // USES QuadratureData2DLinear

#include <string.h> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature );

// ----------------------------------------------------------------------
// Test copy constuctor.
void
pylith::feassemble::TestQuadrature::testCopyConstructor(void)
{ // testClone
  // Semi-random values manually set to check cloning
  const double minJacobianE = 1.0;
  const bool checkConditioning = true;
  const int cellDimE = 1;
  const int numBasisE = 2;
  const int numQuadPtsE = 1;
  const int spaceDimE = 1;
  const double basisE[] = { 0.2, 0.4 };
  const double basisDerivE[] = { 0.8, 1.6 };
  const double quadPtsRefE[] = { 3.2 };
  const double quadWtsE[] = { 6.4 };
  const double quadPtsE[] = { 12.8 };
  const double jacobianE[] = { 2.56 };
  const double jacobianInvE[] = { 5.12 };
  const double jacobianDetE[] = { 10.24 };
  GeometryLine1D geometry;

  // Set values
  Quadrature<topology::Mesh> qOrig;
  qOrig.refGeometry(&geometry);
  qOrig.minJacobian(minJacobianE);
  qOrig.checkConditioning(checkConditioning);
  qOrig.initialize(basisE, numQuadPtsE, numBasisE,
		   basisDerivE, numQuadPtsE, numBasisE, cellDimE,
		   quadPtsRefE, numQuadPtsE, cellDimE,
		   quadWtsE, numQuadPtsE,
		   spaceDimE);

  // Copy
  Quadrature<topology::Mesh> qCopy(qOrig);
  
  // Check copy
  CPPUNIT_ASSERT(0 == qCopy._engine);
  CPPUNIT_ASSERT_EQUAL(minJacobianE, qCopy._minJacobian);
  CPPUNIT_ASSERT_EQUAL(checkConditioning, qCopy._checkConditioning);
  CPPUNIT_ASSERT_EQUAL(cellDimE, qCopy.cellDim());
  CPPUNIT_ASSERT_EQUAL(numBasisE, qCopy.numBasis());
  CPPUNIT_ASSERT_EQUAL(numQuadPtsE, qCopy.numQuadPts());
  CPPUNIT_ASSERT_EQUAL(spaceDimE, qCopy.spaceDim());

  const double_array& basis = qCopy.basis();
  size_t size = numBasisE * numQuadPtsE;
  CPPUNIT_ASSERT_EQUAL(size, basis.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisE[i], basis[i]);

  const double_array& basisDerivRef = qCopy._basisDerivRef;
  size = numBasisE * numQuadPtsE * spaceDimE;
  CPPUNIT_ASSERT_EQUAL(size, basisDerivRef.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDerivE[i], basisDerivRef[i]);

  const double_array& quadPtsRef = qCopy._quadPtsRef;
  size = numQuadPtsE * cellDimE;
  CPPUNIT_ASSERT_EQUAL(size, quadPtsRef.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsRefE[i], quadPtsRef[i]);

  const double_array& quadWts = qCopy.quadWts();
  size = numQuadPtsE;
  CPPUNIT_ASSERT_EQUAL(size, quadWts.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadWtsE[i], quadWts[i]);

  CPPUNIT_ASSERT_EQUAL(geometry.cellDim(), qCopy.refGeometry().cellDim());
  CPPUNIT_ASSERT_EQUAL(geometry.spaceDim(), qCopy.refGeometry().spaceDim());
  CPPUNIT_ASSERT_EQUAL(geometry.numCorners(), qCopy.refGeometry().numCorners());
} // testCopyConstructor

// ----------------------------------------------------------------------
// Test checkConditioning()
void
pylith::feassemble::TestQuadrature::testCheckConditioning(void)
{ // testCheckConditioning
  Quadrature<topology::Mesh> q;

  CPPUNIT_ASSERT_EQUAL(false, q.checkConditioning());
  q.checkConditioning(true);
  CPPUNIT_ASSERT_EQUAL(true, q.checkConditioning());
  q.checkConditioning(false);
  CPPUNIT_ASSERT_EQUAL(false, q.checkConditioning());
} // testCheckConditioning

// ----------------------------------------------------------------------
// Test quadPts(), basisDeriv(), jacobian(), and jacobianDet().
void
pylith::feassemble::TestQuadrature::testEngineAccessors(void)
{ // testEngineAccessors
  const int cellDim = 2;
  const int numBasis = 5;
  const int numQuadPts = 1;
  const int spaceDim = 3;
  const double basis[] = { 
    1.1, 1.2, 1.3, 1.4, 1.5
  };
  const double basisDerivRef[] = {
    2.1, 2.2, 2.3,
    2.4, 2.5, 2.6,
    2.7, 2.8, 2.9,
    2.10, 2.11, 2.12,
    2.13, 2.14, 2.15,
  };
  const double quadPtsRef[] = { 3.1, 3.2, 3.3 };
  const double quadWts[] = { 4.0 };

  QuadratureRefCell refCell;
  refCell.initialize(basis, numQuadPts, numBasis,
		     basisDerivRef, numQuadPts, numBasis, cellDim,
		     quadPtsRef, numQuadPts, cellDim,
		     quadWts, numQuadPts,
		     spaceDim);

  Quadrature1D engine(refCell);
  engine.initialize();

  Quadrature<topology::Mesh> q;
  q._engine = engine.clone();

  size_t size = 0;

  size = numQuadPts * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, q.quadPts().size());

  size = numQuadPts * cellDim * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, q.jacobian().size());

  size = numQuadPts;
  CPPUNIT_ASSERT_EQUAL(size, q.jacobianDet().size());

  size = numQuadPts * numBasis * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, q.basisDeriv().size());
} // testEngineAccessors

// ----------------------------------------------------------------------
// Test computeGeometry() and retrieveGeometry()
void
pylith::feassemble::TestQuadrature::testComputeGeometry(void)
{ // testComputeGeometry
  typedef pylith::topology::Mesh::SieveMesh SieveMesh;

  QuadratureData2DLinear data;
  const int cellDim = data.cellDim;
  const int numBasis = data.numBasis;
  const int numQuadPts = data.numQuadPts;
  const int spaceDim = data.spaceDim;

  const int numCells = data.numCells;
  const double* vertCoords = data.vertices;
  const double* quadPtsE = data.quadPts;
  const double* jacobianE = data.jacobian;
  const double* jacobianDetE = data.jacobianDet;
  const double* basisDerivE = data.basisDeriv;

  const double minJacobian = 1.0e-06;

  // Create mesh with test cell
  topology::Mesh mesh(data.cellDim);
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  ALE::Obj<SieveMesh::sieve_type> sieve = 
    new SieveMesh::sieve_type(mesh.comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  // Cells and vertices
  const bool interpolate = false;
  ALE::Obj<SieveFlexMesh::sieve_type> s = 
    new SieveFlexMesh::sieve_type(sieve->comm(), sieve->debug());
  
  ALE::SieveBuilder<SieveFlexMesh>::buildTopology(s, cellDim, numCells,
                                              const_cast<int*>(data.cells), 
					      data.numVertices,
                                              interpolate, numBasis);
  std::map<SieveFlexMesh::point_type,SieveFlexMesh::point_type> renumbering;
  ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
  sieveMesh->setSieve(sieve);
  sieveMesh->stratify();
  ALE::SieveBuilder<SieveMesh>::buildCoordinates(sieveMesh, spaceDim, 
						 data.vertices);

  // Setup quadrature and compute geometry
  GeometryTri2D geometry;
  Quadrature<topology::Mesh> quadrature;
  quadrature.refGeometry(&geometry);
  quadrature.minJacobian(minJacobian);
  quadrature.initialize(data.basis, numQuadPts, numBasis,
			data.basisDerivRef, numQuadPts, numBasis, cellDim,
			data.quadPtsRef, numQuadPts, cellDim,
			data.quadWts, numQuadPts,
			spaceDim);

  const ALE::Obj<SieveMesh::label_sequence>& cells = sieveMesh->heightStratum(0);
  CPPUNIT_ASSERT(!cells.isNull());
  quadrature.initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  quadrature.computeGeometry(mesh, cells);
#else
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<topology::Mesh::RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates, 
						coordinatesCell.size(),
						&coordinatesCell[0]);
#endif

  size_t size = 0;

  // Check values from computeGeometry()
  const double tolerance = 1.0e-06;
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end(); 
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
#if defined(PRECOMPUTE_GEOMETRY)
    quadrature.retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
    quadrature.computeGeometry(coordinatesCell, *c_iter);    
#endif

    const double_array& quadPts = quadrature.quadPts();
    size = numQuadPts * spaceDim;
    CPPUNIT_ASSERT_EQUAL(size, quadPts.size());
    for (size_t i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(quadPtsE[i], quadPts[i], tolerance);

    const double_array& jacobian = quadrature.jacobian();
    size = numQuadPts * cellDim * spaceDim;
    CPPUNIT_ASSERT_EQUAL(size, jacobian.size());
    for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianE[i], jacobian[i], tolerance);

    const double_array& jacobianDet = quadrature.jacobianDet();
    size = numQuadPts;
    CPPUNIT_ASSERT_EQUAL(size, jacobianDet.size());
    for (size_t i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianDetE[i], jacobianDet[i], tolerance);
    
    const double_array& basisDeriv = quadrature.basisDeriv();
    size = numQuadPts * numBasis * spaceDim;
    CPPUNIT_ASSERT_EQUAL(size, basisDeriv.size());
    for (size_t i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(basisDerivE[i], basisDeriv[i], tolerance);
  } // for

  // Check clear()
  quadrature.clear();
  
  CPPUNIT_ASSERT(0 == quadrature._geometryFields);
  CPPUNIT_ASSERT(0 == quadrature._engine);

  // Make sure caling clear without data doesn't generate errors 
  quadrature.clear();
} // testComputeGeometry

// ----------------------------------------------------------------------
// Test computeGeometry() will coordinates and cell.
void
pylith::feassemble::TestQuadrature::testComputeGeometryCell(void)
{ // testComputeGeometryCell
  typedef pylith::topology::Mesh::SieveMesh SieveMesh;

  QuadratureData2DLinear data;
  const int cellDim = data.cellDim;
  const int numBasis = data.numBasis;
  const int numQuadPts = data.numQuadPts;
  const int spaceDim = data.spaceDim;

  const int numCells = data.numCells;
  double_array vertCoords(data.vertices, numBasis*spaceDim);
  const double* quadPtsE = data.quadPts;
  const double* jacobianE = data.jacobian;
  const double* jacobianDetE = data.jacobianDet;
  const double* basisDerivE = data.basisDeriv;

  const double minJacobian = 1.0e-06;

#if 0
  // Create mesh with test cell
  topology::Mesh mesh(data.cellDim);
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  CPPUNIT_ASSERT(!sieveMesh.isNull());
  ALE::Obj<SieveMesh::sieve_type> sieve = 
    new SieveMesh::sieve_type(mesh.comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  // Cells and vertices
  const bool interpolate = false;
  ALE::Obj<SieveFlexMesh::sieve_type> s = 
    new SieveFlexMesh::sieve_type(sieve->comm(), sieve->debug());
  
  ALE::SieveBuilder<SieveFlexMesh>::buildTopology(s, cellDim, numCells,
                                              const_cast<int*>(data.cells), 
					      data.numVertices,
                                              interpolate, numBasis);
  std::map<SieveFlexMesh::point_type,SieveFlexMesh::point_type> renumbering;
  ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
  sieveMesh->setSieve(sieve);
  sieveMesh->stratify();
  ALE::SieveBuilder<SieveMesh>::buildCoordinates(sieveMesh, spaceDim, 
						 data.vertices);
#endif

  // Setup quadrature and compute geometry
  GeometryTri2D geometry;
  Quadrature<topology::Mesh> quadrature;
  quadrature.refGeometry(&geometry);
  quadrature.minJacobian(minJacobian);
  quadrature.initialize(data.basis, numQuadPts, numBasis,
			data.basisDerivRef, numQuadPts, numBasis, cellDim,
			data.quadPtsRef, numQuadPts, cellDim,
			data.quadWts, numQuadPts,
			spaceDim);

#if 0
  const ALE::Obj<SieveMesh::label_sequence>& cells = sieveMesh->heightStratum(0);
  CPPUNIT_ASSERT(!cells.isNull());
#endif

  quadrature.initializeGeometry();
  quadrature.computeGeometry(vertCoords, 0);
  
  size_t size = 0;

  // Check values from computeGeometry()
  const double tolerance = 1.0e-06;

  const double_array& quadPts = quadrature.quadPts();
  size = numQuadPts * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, quadPts.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(quadPtsE[i], quadPts[i], tolerance);
  
  const double_array& jacobian = quadrature.jacobian();
  size = numQuadPts * cellDim * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, jacobian.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianE[i], jacobian[i], tolerance);
  
  const double_array& jacobianDet = quadrature.jacobianDet();
  size = numQuadPts;
  CPPUNIT_ASSERT_EQUAL(size, jacobianDet.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianDetE[i], jacobianDet[i], tolerance);
  
  const double_array& basisDeriv = quadrature.basisDeriv();
  size = numQuadPts * numBasis * spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, basisDeriv.size());
  for (size_t i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(basisDerivE[i], basisDeriv[i], tolerance);

  quadrature.clear();

  CPPUNIT_ASSERT(0 == quadrature._geometryFields);
  CPPUNIT_ASSERT(0 == quadrature._engine);
} // testComputeGeometryCell


// End of file 
