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

#include "TestQuadrature.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature1D.hh" // USES Quadrature1D
#include "pylith/feassemble/GeometryLine1D.hh" // USES GeometryLine1D

#include "data/QuadratureData.hh" // USES QuadratureData

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature );

// ----------------------------------------------------------------------
// Test clone
void
pylith::feassemble::TestQuadrature::testClone(void)
{ // testClone
  // Semi-random values manually set to check cloning
  const double minJacobianE = 1.0;
  const int cellDimE = 1;
  const int numBasisE = 2;
  const int numQuadPtsE = 1;
  const int spaceDimE = 1;
  const double verticesE[] = { 0.12, 0.42 };
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
  Quadrature1D qOrig;
  qOrig._minJacobian = minJacobianE;
  qOrig._cellDim = cellDimE;
  qOrig._numBasis = numBasisE;
  qOrig._numQuadPts = numQuadPtsE;
  qOrig._spaceDim = spaceDimE;

  size_t size = 2;
  qOrig._vertices.resize(size);
  memcpy(&qOrig._vertices[0], verticesE, size*sizeof(double));
  
  size = 2;
  qOrig._basis.resize(size);
  memcpy(&qOrig._basis[0], basisE, size*sizeof(double));
  
  size = 2;
  qOrig._basisDerivRef.resize(size);
  memcpy(&qOrig._basisDerivRef[0], basisDerivE, size*sizeof(double));

  size = 1;
  qOrig._quadPtsRef.resize(size);
  memcpy(&qOrig._quadPtsRef[0], quadPtsRefE, size*sizeof(double));

  size = 1;
  qOrig._quadWts.resize(size);
  memcpy(&qOrig._quadWts[0], quadWtsE, size*sizeof(double));

  size = 1;
  qOrig._quadPts.resize(size);
  memcpy(&qOrig._quadPts[0], quadPtsE, size*sizeof(double));

  size = 1;
  qOrig._jacobian.resize(size);
  memcpy(&qOrig._jacobian[0], jacobianE, size*sizeof(double));

  size = 1;
  qOrig._jacobianInv.resize(size);
  memcpy(&qOrig._jacobianInv[0], jacobianInvE, size*sizeof(double));

  size = 1;
  qOrig._jacobianDet.resize(size);
  memcpy(&qOrig._jacobianDet[0], jacobianDetE, size*sizeof(double));

  qOrig._geometry = geometry.clone();

  // Clone
  const Quadrature* qCopy = qOrig.clone();

  // Check clone
  CPPUNIT_ASSERT(0 != qCopy);

  CPPUNIT_ASSERT_EQUAL(minJacobianE, qCopy->_minJacobian);
  CPPUNIT_ASSERT_EQUAL(cellDimE, qCopy->cellDim());
  CPPUNIT_ASSERT_EQUAL(numBasisE, qCopy->numBasis());
  CPPUNIT_ASSERT_EQUAL(numQuadPtsE, qCopy->numQuadPts());
  CPPUNIT_ASSERT_EQUAL(spaceDimE, qCopy->spaceDim());

  const double_array& vertices = qCopy->vertices();
  size = numBasisE * cellDimE;
  CPPUNIT_ASSERT_EQUAL(size, vertices.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(verticesE[i], vertices[i]);

  const double_array& basis = qCopy->basis();
  size = numBasisE * numQuadPtsE;
  CPPUNIT_ASSERT_EQUAL(size, basis.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisE[i], basis[i]);

  const double_array& basisDeriv = qCopy->basisDeriv();
  size = numBasisE * numQuadPtsE * spaceDimE;
  CPPUNIT_ASSERT_EQUAL(size, basisDeriv.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDerivE[i], basisDeriv[i]);

  const double_array& quadPtsRef = qCopy->_quadPtsRef;
  size = numQuadPtsE * cellDimE;
  CPPUNIT_ASSERT_EQUAL(size, quadPtsRef.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsRefE[i], quadPtsRef[i]);

  const double_array& quadWts = qCopy->quadWts();
  size = numQuadPtsE;
  CPPUNIT_ASSERT_EQUAL(size, quadWts.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadWtsE[i], quadWts[i]);

  const double_array& quadPts = qCopy->quadPts();
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, quadPts.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsE[i], quadPts[i]);

  const double_array& jacobian = qCopy->_jacobian;
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobian.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianE[i], jacobian[i]);

  const double_array& jacobianInv = qCopy->jacobianInv();
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobianInv.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianInvE[i], jacobianInv[i]);

  const double_array& jacobianDet = qCopy->jacobianDet();
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobianDet.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianDetE[i], jacobianDet[i]);

  CPPUNIT_ASSERT(0 != qCopy->_geometry);
  CPPUNIT_ASSERT_EQUAL(geometry.cellDim(), qCopy->_geometry->cellDim());
  CPPUNIT_ASSERT_EQUAL(geometry.spaceDim(), qCopy->_geometry->spaceDim());
  CPPUNIT_ASSERT_EQUAL(geometry.numCorners(), qCopy->_geometry->numCorners());

  delete qCopy; qCopy = 0;
} // testCopy

// ----------------------------------------------------------------------
// Test minJacobian()
void
pylith::feassemble::TestQuadrature::testMinJacobian(void)
{ // testMinJacobian
  Quadrature1D q;
  const double min = 1.0;
  q.minJacobian(min);
  CPPUNIT_ASSERT_EQUAL(min, q._minJacobian);
} // testMinJacobian

// ----------------------------------------------------------------------
// Test refGeometry()
void
pylith::feassemble::TestQuadrature::testRefGeometry(void)
{ // testRefGeometry
  GeometryLine1D geometry;
  Quadrature1D quadrature;
  quadrature.refGeometry(&geometry);
  const CellGeometry& test = quadrature.refGeometry();

  CPPUNIT_ASSERT_EQUAL(geometry.cellDim(), test.cellDim());
  CPPUNIT_ASSERT_EQUAL(geometry.spaceDim(), test.spaceDim());
  CPPUNIT_ASSERT_EQUAL(geometry.numCorners(), test.numCorners());
} // testRefGeometry

// ----------------------------------------------------------------------
// Test initialize()
void
pylith::feassemble::TestQuadrature::testInitialize(void)
{ // initialize
  
  const int cellDim = 1;
  const int numBasis = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double vertices[] = { -1.0, 1.0 };
  const double basis[] = { 0.5, 0.5 };
  const double basisDeriv[] = { -0.5, 0.5 };
  const double quadPtsRef[] = { 0.0 };
  const double quadWts[] = { 2.0 };
  const double minJacobian = 1.0;

  Quadrature1D q;
  q.initialize(vertices, basis, basisDeriv, quadPtsRef, quadWts,
	       cellDim, numBasis, numQuadPts, spaceDim);
  
  CPPUNIT_ASSERT_EQUAL(cellDim, q._cellDim);
  CPPUNIT_ASSERT_EQUAL(numBasis, q._numBasis);
  CPPUNIT_ASSERT_EQUAL(numQuadPts, q._numQuadPts);
  CPPUNIT_ASSERT_EQUAL(spaceDim, q._spaceDim);

  size_t size = numBasis * cellDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(vertices[i], q._vertices[i]);

  size = numBasis * numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basis[i], q._basis[i]);

  size = numBasis * numQuadPts * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDeriv[i], q._basisDerivRef[i]);

  size = numQuadPts * cellDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsRef[i], q._quadPtsRef[i]);

  size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadWts[i], q._quadWts[i]);

  // Make sure Jacobian stuff has been allocated
  size = numQuadPts*cellDim*spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, q._jacobian.size());
  
  size = numQuadPts*spaceDim*cellDim;
  CPPUNIT_ASSERT_EQUAL(size, q._jacobianInv.size());
  
  size = numQuadPts;
  CPPUNIT_ASSERT_EQUAL(size, q._jacobianDet.size());
  
  size = numQuadPts*spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, q._quadPts.size());
} // initialize

// ----------------------------------------------------------------------
// Test initialize() & computeGeometry()
void
pylith::feassemble::TestQuadrature::_testComputeGeometry(Quadrature* pQuad,
					      const QuadratureData& data) const
{ // testComputeGeometry
  const int cellDim = data.cellDim;
  const int numBasis = data.numBasis;
  const int numQuadPts = data.numQuadPts;
  const int spaceDim = data.spaceDim;
  const double* verticesRef = data.verticesRef;
  const double* basis = data.basis;
  const double* basisDeriv = data.basisDeriv;
  const double* quadPtsRef = data.quadPtsRef;
  const double* quadWts = data.quadWts;

  const int numVertices = data.numVertices;
  const int numCells = data.numCells;
  const double* vertCoords = data.vertices;
  const int* cells = data.cells;
  const double* quadPts = data.quadPts;
  const double* jacobian = data.jacobian;
  const double* jacobianInv = data.jacobianInv;
  const double* jacobianDet = data.jacobianDet;

  const double minJacobian = 1.0e-06;

  pQuad->minJacobian(minJacobian);
  pQuad->initialize(verticesRef, basis, basisDeriv, quadPtsRef, quadWts,
		    cellDim, numBasis, numQuadPts, spaceDim);

  // Create mesh with test cell
  ALE::Obj<Mesh> mesh = new Mesh(PETSC_COMM_WORLD, cellDim);
  CPPUNIT_ASSERT(!mesh.isNull());
  ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  const bool interpolate = false;
  ALE::SieveBuilder<Mesh>::buildTopology(sieve, cellDim, numCells,
		     (int*) cells, numVertices, interpolate, numBasis);
  mesh->setSieve(sieve);
  mesh->stratify();
  ALE::SieveBuilder<Mesh>::buildCoordinates(mesh, spaceDim, vertCoords);
  
  // Check values from computeGeometry()
  const ALE::Obj<Mesh::label_sequence>& cellsMesh = mesh->heightStratum(0);
  CPPUNIT_ASSERT(!cellsMesh.isNull());
  const Mesh::label_sequence::iterator e_iter = cellsMesh->begin(); 
  const ALE::Obj<Mesh::real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  CPPUNIT_ASSERT(!coordinates.isNull());
  pQuad->computeGeometry(mesh, coordinates, *e_iter);

  CPPUNIT_ASSERT(1 == numCells);

  const double tolerance = 1.0e-06;
  int size = numQuadPts * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(quadPts[i], pQuad->_quadPts[i], tolerance);

  size = numQuadPts * cellDim * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobian[i], pQuad->_jacobian[i], 
				 tolerance);

  size = numQuadPts * spaceDim * cellDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianInv[i], pQuad->_jacobianInv[i], 
				 tolerance);

  size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianDet[i], pQuad->_jacobianDet[i], 
				 tolerance);
} // testComputeGeometry


// End of file 
