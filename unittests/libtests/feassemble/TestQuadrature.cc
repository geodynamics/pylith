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
  const double basisQuadE[] = { 0.2, 0.4 };
  const double basisDerivQuadE[] = { 0.8, 1.6 };
  const double quadPtsRefE[] = { 3.2 };
  const double quadWtsE[] = { 6.4 };
  const double quadPtsE[] = { 12.8 };
  const double jacobianQuadE[] = { 2.56 };
  const double jacobianInvQuadE[] = { 5.12 };
  const double jacobianDetQuadE[] = { 10.24 };

  // Set values
  Quadrature1D qOrig;
  qOrig._minJacobian = minJacobianE;
  qOrig._cellDim = cellDimE;
  qOrig._numBasis = numBasisE;
  qOrig._numQuadPts = numQuadPtsE;
  qOrig._spaceDim = spaceDimE;

  size_t size = 2;
  qOrig._basisQuad.resize(size);
  memcpy(&qOrig._basisQuad[0], basisQuadE, size*sizeof(double));
  
  size = 2;
  qOrig._basisDerivQuad.resize(size);
  memcpy(&qOrig._basisDerivQuad[0], basisDerivQuadE, size*sizeof(double));

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
  qOrig._jacobianQuad.resize(size);
  memcpy(&qOrig._jacobianQuad[0], jacobianQuadE, size*sizeof(double));

  size = 1;
  qOrig._jacobianInvQuad.resize(size);
  memcpy(&qOrig._jacobianInvQuad[0], jacobianInvQuadE, size*sizeof(double));

  size = 1;
  qOrig._jacobianDetQuad.resize(size);
  memcpy(&qOrig._jacobianDetQuad[0], jacobianDetQuadE, size*sizeof(double));

  // Clone
  const Quadrature* qCopy = qOrig.clone();

  // Check clone
  CPPUNIT_ASSERT(0 != qCopy);

  CPPUNIT_ASSERT_EQUAL(minJacobianE, qCopy->_minJacobian);
  CPPUNIT_ASSERT_EQUAL(cellDimE, qCopy->cellDim());
  CPPUNIT_ASSERT_EQUAL(numBasisE, qCopy->numBasis());
  CPPUNIT_ASSERT_EQUAL(numQuadPtsE, qCopy->numQuadPts());
  CPPUNIT_ASSERT_EQUAL(spaceDimE, qCopy->spaceDim());

  const double_array& basisQuad = qCopy->basisQuad();
  size = numBasisE * numQuadPtsE;
  CPPUNIT_ASSERT_EQUAL(size, basisQuad.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisQuadE[i], basisQuad[i]);

  const double_array& basisDerivQuad = qCopy->basisDerivQuad();
  size = numBasisE * numQuadPtsE * spaceDimE;
  CPPUNIT_ASSERT_EQUAL(size, basisDerivQuad.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDerivQuadE[i], basisDerivQuad[i]);

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

  const double_array& jacobianQuad = qCopy->_jacobianQuad;
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobianQuad.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianQuadE[i], jacobianQuad[i]);

  const double_array& jacobianInvQuad = qCopy->jacobianInvQuad();
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobianInvQuad.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianInvQuadE[i], jacobianInvQuad[i]);

  const double_array& jacobianDetQuad = qCopy->jacobianDetQuad();
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobianDetQuad.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianDetQuadE[i], jacobianDetQuad[i]);

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
// Test initialize()
void
pylith::feassemble::TestQuadrature::testInitialize(void)
{ // initialize
  
  const int cellDim = 1;
  const int numBasis = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double basisVert[] = { 0.5, 0.5, 0.3, 0.7 };
  const double basisDerivVert[] = { -0.5, 0.5, -0.4, 0.7 };
  const double basisQuad[] = { 0.5, 0.5 };
  const double basisDerivQuad[] = { -0.5, 0.5 };
  const double quadPtsRef[] = { 0.0 };
  const double quadWts[] = { 2.0 };
  const double minJacobian = 1.0;

  Quadrature1D q;
  q.initialize(basisVert, basisDerivVert, 
	       basisQuad, basisDerivQuad, quadPtsRef, quadWts,
	       cellDim, numBasis, numQuadPts, spaceDim);
  
  CPPUNIT_ASSERT_EQUAL(cellDim, q._cellDim);
  CPPUNIT_ASSERT_EQUAL(numBasis, q._numBasis);
  CPPUNIT_ASSERT_EQUAL(numQuadPts, q._numQuadPts);
  CPPUNIT_ASSERT_EQUAL(spaceDim, q._spaceDim);

  const int numVertices = numBasis;

  size_t size = numBasis * numVertices;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisVert[i], q._basisVert[i]);

  size = numBasis * numVertices * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDerivVert[i], q._basisDerivVert[i]);

  size = numBasis * numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisQuad[i], q._basisQuad[i]);

  size = numBasis * numQuadPts * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDerivQuad[i], q._basisDerivQuad[i]);

  size = numQuadPts * cellDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsRef[i], q._quadPtsRef[i]);

  size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadWts[i], q._quadWts[i]);

  // Make sure Jacobian stuff has been allocated
  size = numVertices*cellDim*spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, q._jacobianVert.size());
  
  size = numVertices;
  CPPUNIT_ASSERT_EQUAL(size, q._jacobianDetVert.size());
  
  // Make sure Jacobian stuff has been allocated
  size = numQuadPts*cellDim*spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, q._jacobianQuad.size());
  
  size = numQuadPts*spaceDim*cellDim;
  CPPUNIT_ASSERT_EQUAL(size, q._jacobianInvQuad.size());
  
  size = numQuadPts;
  CPPUNIT_ASSERT_EQUAL(size, q._jacobianDetQuad.size());
  
  size = numQuadPts*spaceDim;
  CPPUNIT_ASSERT_EQUAL(size, q._quadPts.size());
} // initialize

// ----------------------------------------------------------------------
// Test initialize() & computeGeometryVert()
void
pylith::feassemble::TestQuadrature::_testComputeGeometryVert(
					    Quadrature* pQuad,
					    const QuadratureData& data) const
{ // testComputeGeometry
  const int cellDim = data.cellDim;
  const int numBasis = data.numBasis;
  const int numQuadPts = data.numQuadPts;
  const int spaceDim = data.spaceDim;
  const double* basisVert = data.basisVert;
  const double* basisDerivVert = data.basisDerivVert;
  const double* basisQuad = data.basisQuad;
  const double* basisDerivQuad = data.basisDerivQuad;
  const double* quadPtsRef = data.quadPtsRef;
  const double* quadWts = data.quadWts;

  const int numVertices = data.numVertices;
  const int numCells = data.numCells;
  const double* vertCoords = data.vertices;
  const int* cells = data.cells;
  const double* quadPts = data.quadPts;
  const double* jacobianVert = data.jacobianVert;
  const double* jacobianDetVert = data.jacobianDetVert;

  const double minJacobian = 1.0e-06;

  pQuad->minJacobian(minJacobian);
  pQuad->initialize(basisVert, basisDerivVert, 
		    basisQuad, basisDerivQuad, quadPtsRef, quadWts,
		    cellDim, numBasis, numQuadPts, spaceDim);

  // Create mesh with test cell
  typedef ALE::Mesh Mesh;
  typedef ALE::Mesh::sieve_type sieve_type;
  ALE::Obj<Mesh> mesh = new Mesh(PETSC_COMM_WORLD, cellDim);
  ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());

  const bool interpolate = false;
  ALE::SieveBuilder<Mesh>::buildTopology(sieve, cellDim, numCells,
		     (int*) cells, numVertices, interpolate, numBasis);
  mesh->setSieve(sieve);
  mesh->stratify();
  ALE::SieveBuilder<Mesh>::buildCoordinates(mesh, spaceDim, vertCoords);
  
  // Check values from computeGeometryVert()
  const ALE::Obj<Mesh::label_sequence>& elements = mesh->heightStratum(0);
  const Mesh::label_sequence::iterator e_iter = elements->begin(); 
  const ALE::Obj<Mesh::real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  pQuad->computeGeometryVert(mesh, coordinates, *e_iter);

  CPPUNIT_ASSERT(1 == numCells);

  const double tolerance = 1.0e-06;
  int size = numVertices * cellDim * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianVert[i], pQuad->_jacobianVert[i], 
				 tolerance);

  size = numVertices;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianDetVert[i], 
				 pQuad->_jacobianDetVert[i], 
				 tolerance);
} // testComputeGeometryVert

// ----------------------------------------------------------------------
// Test initialize() & computeGeometry()
void
pylith::feassemble::TestQuadrature::_testComputeGeometryQuad(
					    Quadrature* pQuad,
					    const QuadratureData& data) const
{ // testComputeGeometry
  const int cellDim = data.cellDim;
  const int numBasis = data.numBasis;
  const int numQuadPts = data.numQuadPts;
  const int spaceDim = data.spaceDim;
  const double* basisVert = data.basisVert;
  const double* basisDerivVert = data.basisDerivVert;
  const double* basisQuad = data.basisQuad;
  const double* basisDerivQuad = data.basisDerivQuad;
  const double* quadPtsRef = data.quadPtsRef;
  const double* quadWts = data.quadWts;

  const int numVertices = data.numVertices;
  const int numCells = data.numCells;
  const double* vertCoords = data.vertices;
  const int* cells = data.cells;
  const double* quadPts = data.quadPts;
  const double* jacobianQuad = data.jacobianQuad;
  const double* jacobianInvQuad = data.jacobianInvQuad;
  const double* jacobianDetQuad = data.jacobianDetQuad;

  const double minJacobian = 1.0e-06;

  pQuad->minJacobian(minJacobian);
  pQuad->initialize(basisVert, basisDerivVert, 
		    basisQuad, basisDerivQuad, quadPtsRef, quadWts,
		    cellDim, numBasis, numQuadPts, spaceDim);

  // Create mesh with test cell
  typedef ALE::Mesh Mesh;
  typedef ALE::Mesh::sieve_type sieve_type;
  ALE::Obj<Mesh> mesh = new Mesh(PETSC_COMM_WORLD, cellDim);
  ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());

  const bool interpolate = false;
  ALE::SieveBuilder<Mesh>::buildTopology(sieve, cellDim, numCells,
		     (int*) cells, numVertices, interpolate, numBasis);
  mesh->setSieve(sieve);
  mesh->stratify();
  ALE::SieveBuilder<Mesh>::buildCoordinates(mesh, spaceDim, vertCoords);
  
  // Check values from computeGeometry()
  const ALE::Obj<Mesh::label_sequence>& elements = mesh->heightStratum(0);
  const Mesh::label_sequence::iterator e_iter = elements->begin(); 
  const ALE::Obj<Mesh::real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  pQuad->computeGeometryQuad(mesh, coordinates, *e_iter);

  CPPUNIT_ASSERT(1 == numCells);

  const double tolerance = 1.0e-06;
  int size = numQuadPts * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(quadPts[i], pQuad->_quadPts[i], tolerance);

  size = numQuadPts * cellDim * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianQuad[i], pQuad->_jacobianQuad[i], 
				 tolerance);

  size = numQuadPts * spaceDim * cellDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianInvQuad[i], 
				 pQuad->_jacobianInvQuad[i], 
				 tolerance);

  size = numQuadPts;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianDetQuad[i], 
				 pQuad->_jacobianDetQuad[i], 
				 tolerance);
} // testComputeGeometryQuad


// End of file 
