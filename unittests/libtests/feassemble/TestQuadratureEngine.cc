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

#include "TestQuadratureEngine.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature1D.hh" // USES Quadrature1D

#include "data/QuadratureData.hh" // USES QuadratureData

#include <string.h> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadratureEngine );

// ----------------------------------------------------------------------
// Test copy constructor.
void
pylith::feassemble::TestQuadratureEngine::testCopyConstructor(void)
{ // testClone
  // Semi-random values manually set to check cloning
  const double quadPtsE[] = { 12.8 };
  const double jacobianE[] = { 2.56 };
  const double jacobianInvE[] = { 5.12 };
  const double jacobianDetE[] = { 10.24 };
  const double basisDerivE[] = { 0.8, 1.6 };

  Quadrature1D engineOrig;
  int size = 0;

  // Set values
  size = 1;
  engineOrig._quadPts.resize(size);
  memcpy(&engineOrig._quadPts[0], quadPtsE, size*sizeof(double));

  size = 1;
  engineOrig._jacobian.resize(size);
  memcpy(&engineOrig._jacobian[0], jacobianE, size*sizeof(double));

  size = 1;
  engineOrig._jacobianInv.resize(size);
  memcpy(&engineOrig._jacobianInv[0], jacobianInvE, size*sizeof(double));

  size = 1;
  engineOrig._jacobianDet.resize(size);
  memcpy(&engineOrig._jacobianDet[0], jacobianDetE, size*sizeof(double));

  // Copy
  Quadrature1D engine(engineOrig);

  const double_array& quadPts = engine.quadPts();
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, quadPts.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPtsE[i], quadPts[i]);

  const double_array& jacobian = engine._jacobian;
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobian.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianE[i], jacobian[i]);

  const double_array& jacobianInv = engine.jacobianInv();
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobianInv.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianInvE[i], jacobianInv[i]);

  const double_array& jacobianDet = engine.jacobianDet();
  size = 1;
  CPPUNIT_ASSERT_EQUAL(size, jacobianDet.size());
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianDetE[i], jacobianDet[i]);
} // testCopyConstructor

// ----------------------------------------------------------------------
// Test computeGeometry().
void
pylith::feassemble::TestQuadratureEngine::_testComputeGeometry(Quadrature<topology::Mesh>* pQuad,
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

  QuadratureRefCell quadRefCell;
  quadRefCell.minJacobian(minJacobian);
  quadRefCell.initialize(basis, basisDerivRef, quadPtsRef, quadWts,
			 cellDim, numBasis, numQuadPts, spaceDim);

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


// End of file 
