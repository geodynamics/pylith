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

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "TestQuadrature.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature1D.hh" // USES Quadrature1D
#include "pylith/feassemble/GeometryLine1D.hh" // USES GeometryLine1D

#include "data/QuadratureData.hh" // USES QuadratureData

#include <string.h> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature );

// ----------------------------------------------------------------------
// Test clone
void
pylith::feassemble::TestQuadrature::testClone(void)
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
  Quadrature1D<topology::Mesh> qOrig;
  qOrig._minJacobian = minJacobianE;
  qOrig._checkConditioning = checkConditioning;
  qOrig._cellDim = cellDimE;
  qOrig._numBasis = numBasisE;
  qOrig._numQuadPts = numQuadPtsE;
  qOrig._spaceDim = spaceDimE;

  size_t size = 2;
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
  const Quadrature<topology::Mesh>* qCopy = qOrig.clone();

  // Check clone
  CPPUNIT_ASSERT(0 != qCopy);

  CPPUNIT_ASSERT_EQUAL(minJacobianE, qCopy->_minJacobian);
  CPPUNIT_ASSERT_EQUAL(checkConditioning, qCopy->_checkConditioning);
  CPPUNIT_ASSERT_EQUAL(cellDimE, qCopy->cellDim());
  CPPUNIT_ASSERT_EQUAL(numBasisE, qCopy->numBasis());
  CPPUNIT_ASSERT_EQUAL(numQuadPtsE, qCopy->numQuadPts());
  CPPUNIT_ASSERT_EQUAL(spaceDimE, qCopy->spaceDim());

  const double_array& basis = qCopy->basis();
  size = numBasisE * numQuadPtsE;
  CPPUNIT_ASSERT_EQUAL(size, basis.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisE[i], basis[i]);

  const double_array& basisDerivRef = qCopy->_basisDerivRef;
  size = numBasisE * numQuadPtsE * spaceDimE;
  CPPUNIT_ASSERT_EQUAL(size, basisDerivRef.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(basisDerivE[i], basisDerivRef[i]);

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
} // testClone

// ----------------------------------------------------------------------
// Test checkConditioning()
void
pylith::feassemble::TestQuadrature::testCheckConditioning(void)
{ // testCheckConditioning
  Quadrature1D<topology::Mesh> q;
  CPPUNIT_ASSERT_EQUAL(false, q.checkConditioning());
  q.checkConditioning(true);
  CPPUNIT_ASSERT_EQUAL(true, q.checkConditioning());
  q.checkConditioning(false);
  CPPUNIT_ASSERT_EQUAL(false, q.checkConditioning());
} // testCheckConditioning

#if 0
// ----------------------------------------------------------------------
// Test computeGeometry() and retrieveGeometry() for meshes.
void
pylith::feassemble::TestQuadrature::_testComputeGeometry(Quadrature<topology::Mesh>* pQuad,
					      const QuadratureData& data) const
{ // testComputeGeometry
  const int cellDim = data.cellDim;
  const int numBasis = data.numBasis;
  const int numQuadPts = data.numQuadPts;
  const int spaceDim = data.spaceDim;
  const double* basis = data.basis;
  const double* basisDerivRef = data.basisDerivRef;
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
  pQuad->initialize(basis, basisDerivRef, quadPtsRef, quadWts,
		    cellDim, numBasis, numQuadPts, spaceDim);

  // Create mesh with test cell
  ALE::Obj<Mesh> mesh = new Mesh(PETSC_COMM_WORLD, cellDim);
  CPPUNIT_ASSERT(!mesh.isNull());
  ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  const bool interpolate = false;
  ALE::Obj<ALE::Mesh::sieve_type> s = new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());

  ALE::SieveBuilder<ALE::Mesh>::buildTopology(s, cellDim, numCells,
		     (int*) cells, numVertices, interpolate, numBasis);
  std::map<Mesh::point_type,Mesh::point_type> renumbering;
  ALE::ISieveConverter::convertSieve(*s, *sieve, renumbering);
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

  size = numQuadPts * numBasis * spaceDim;
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(basisDeriv[i], pQuad->_basisDeriv[i], 
				 tolerance);

} // testComputeGeometry
#endif

// End of file 
