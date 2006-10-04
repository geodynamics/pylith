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

// ----------------------------------------------------------------------
// Test clone
void
pylith::feassemble::TestQuadrature::testClone(void)
{ // testClone
  // Semi-random values manually set to check cloning
  const double minJacobian = 1.0;
  const int cellDim = 1;
  const int numCorners = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double basis[] = { 0.2, 0.4 };
  const double basisDeriv[] = { 0.8, 1.6 };
  const double quadPtsRef[] = { 3.2 };
  const double quadWts[] = { 6.4 };
  const double quadPts[] = { 12.8 };
  const double jacobian[] = { 2.56 };
  const double jacobianInv[] = { 5.12 };
  const double jacobianDet[] = { 10.24 };

  // Set values
  Quadrature1D qOrig;
  qOrig._minJacobian = minJacobian;
  qOrig._cellDim = cellDim;
  qOrig._numCorners = numCorners;
  qOrig._numQuadPts = numQuadPts;
  qOrig._spaceDim = spaceDim;

  int size = 2;
  qOrig._basis = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._basis, basis, size*sizeof(double));
  
  size = 2;
  qOrig._basisDeriv = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._basisDeriv, basisDeriv, size*sizeof(double));

  size = 1;
  qOrig._quadPtsRef = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._quadPtsRef, quadPtsRef, size*sizeof(double));

  size = 1;
  qOrig._quadWts = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._quadWts, quadWts, size*sizeof(double));

  size = 1;
  qOrig._quadPts = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._quadPts, quadPts, size*sizeof(double));

  size = 1;
  qOrig._jacobian = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._jacobian, jacobian, size*sizeof(double));

  size = 1;
  qOrig._jacobianInv = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._jacobianInv, jacobianInv, size*sizeof(double));

  size = 1;
  qOrig._jacobianDet = (size > 0) ? new double[size] : 0;
  memcpy(qOrig._jacobianDet, jacobianDet, size*sizeof(double));

  // Clone
  const Quadrature* qCopy = qOrig.clone();

  // Check clone
  CPPUNIT_ASSERT(0 != qCopy);

  CPPUNIT_ASSERT_EQUAL(minJacobian, qCopy->_minJacobian);
  CPPUNIT_ASSERT_EQUAL(cellDim, qCopy->_cellDim);
  CPPUNIT_ASSERT_EQUAL(numCorners, qCopy->_numCorners);
  CPPUNIT_ASSERT_EQUAL(numQuadPts, qCopy->_numQuadPts);
  CPPUNIT_ASSERT_EQUAL(spaceDim, qCopy->_spaceDim);

  CPPUNIT_ASSERT(0 != qCopy->_basis);
  size = numCorners * numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basis[i], qCopy->_basis[i]);

  CPPUNIT_ASSERT(0 != qCopy->_basisDeriv);
  size = numCorners * numQuadPts * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDeriv[i], qCopy->_basisDeriv[i]);

  CPPUNIT_ASSERT(0 != qCopy->_quadPtsRef);
  size = numQuadPts * cellDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsRef[i], qCopy->_quadPtsRef[i]);

  CPPUNIT_ASSERT(0 != qCopy->_quadWts);
  size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadWts[i], qCopy->_quadWts[i]);

  size = 1;

  CPPUNIT_ASSERT(0 != qCopy->_quadPts);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPts[i], qCopy->_quadPts[i]);
  CPPUNIT_ASSERT(0 != qCopy->_jacobian);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobian[i], qCopy->_jacobian[i]);
  CPPUNIT_ASSERT(0 != qCopy->_jacobianInv);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianInv[i], qCopy->_jacobianInv[i]);
  CPPUNIT_ASSERT(0 != qCopy->_jacobianDet);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianDet[i], qCopy->_jacobianDet[i]);
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
// Test initialize()
void
pylith::feassemble::TestQuadrature::testInitialize(void)
{ // initialize
  
  const int cellDim = 1;
  const int numCorners = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double basis[] = { 0.5, 0.5 };
  const double basisDeriv[] = { -0.5, 0.5 };
  const double quadPtsRef[] = { 0.0 };
  const double quadWts[] = { 2.0 };
  const double minJacobian = 1.0;

  Quadrature1D q;
  q.initialize(basis, basisDeriv, quadPtsRef, quadWts,
	       cellDim, numCorners, numQuadPts, spaceDim);
  
  CPPUNIT_ASSERT_EQUAL(cellDim, q._cellDim);
  CPPUNIT_ASSERT_EQUAL(numCorners, q._numCorners);
  CPPUNIT_ASSERT_EQUAL(numQuadPts, q._numQuadPts);
  CPPUNIT_ASSERT_EQUAL(spaceDim, q._spaceDim);

  int size = numCorners * numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basis[i], q._basis[i]);

  size = numCorners * numQuadPts * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDeriv[i], q._basisDeriv[i]);

  size = numQuadPts * cellDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsRef[i], q._quadPtsRef[i]);

  size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadWts[i], q._quadWts[i]);

  // Make sure Jacobian stuff has been allocated
  CPPUNIT_ASSERT(0 != q._jacobian);
  CPPUNIT_ASSERT(0 != q._jacobianInv);
  CPPUNIT_ASSERT(0 != q._jacobianDet);
  CPPUNIT_ASSERT(0 != q._quadPts);
} // initialize

// ----------------------------------------------------------------------
// Test initialize() & computeGeometry()
void
pylith::feassemble::TestQuadrature::_testComputeGeometry(
					    Quadrature* pQuad,
					    const QuadratureData& data) const
{ // testComputeGeometry
  const int cellDim = data.cellDim;
  const int numCorners = data.numCorners;
  const int numQuadPts = data.numQuadPts;
  const int spaceDim = data.spaceDim;
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

  CPPUNIT_ASSERT(0 != basis);
  CPPUNIT_ASSERT(0 != basisDeriv);
  CPPUNIT_ASSERT(0 != quadPtsRef);
  CPPUNIT_ASSERT(0 != quadWts);
  CPPUNIT_ASSERT(0 != vertCoords);
  CPPUNIT_ASSERT(0 != cells);
  CPPUNIT_ASSERT(0 != quadPts);
  CPPUNIT_ASSERT(0 != jacobian);
  CPPUNIT_ASSERT(0 != jacobianInv);
  CPPUNIT_ASSERT(0 != jacobianDet);

  pQuad->minJacobian(minJacobian);
  pQuad->initialize(basis, basisDeriv, quadPtsRef, quadWts,
		    cellDim, numCorners, numQuadPts, spaceDim);

  // Create mesh with test cell
  typedef ALE::Mesh::topology_type topology_type;
  typedef topology_type::sieve_type sieve_type;
  ALE::Obj<ALE::Mesh> mesh = new ALE::Mesh(PETSC_COMM_WORLD, cellDim);
  ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());
  ALE::Obj<topology_type> topology = new topology_type(mesh->comm());

  const bool interpolate = false;
  ALE::New::SieveBuilder<sieve_type>::buildTopology(sieve, cellDim, numCells,
		     (int*) cells, numVertices, interpolate, numCorners);
  sieve->stratify();
  topology->setPatch(0, sieve);
  topology->stratify();
  mesh->setTopology(topology);
  ALE::New::SieveBuilder<sieve_type>::buildCoordinates(
		    mesh->getRealSection("coordinates"), spaceDim, vertCoords);
  
  // Check values from _computeGeometry()
  const ALE::Mesh::topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type::label_sequence>& elements = 
    topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator e_iter = elements->begin(); 
  const ALE::Obj<ALE::Mesh::real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  pQuad->_computeGeometry(coordinates, *e_iter);

  CPPUNIT_ASSERT(1 == numCells);

  const double tolerance = 1.0e-06;
  int size = numQuadPts * spaceDim;
  CPPUNIT_ASSERT(0 != pQuad->_quadPts);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(quadPts[i], pQuad->_quadPts[i], tolerance);

  size = numQuadPts * cellDim * spaceDim;
  CPPUNIT_ASSERT(0 != pQuad->_jacobian);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobian[i], pQuad->_jacobian[i], tolerance);

  size = numQuadPts * spaceDim * cellDim;
  CPPUNIT_ASSERT(0 != pQuad->_jacobianInv);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianInv[i], pQuad->_jacobianInv[i], 
				 tolerance);

  size = numQuadPts;
  CPPUNIT_ASSERT(0 != pQuad->_jacobianDet);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianDet[i], pQuad->_jacobianDet[i], 
				 tolerance);
} // testQuadratic

// End of file 
